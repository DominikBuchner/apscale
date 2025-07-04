import streamlit as st
import dask.dataframe as dd
import pandas as pd
from pathlib import Path
import os


def filter_sequence_metadata(read_data_to_modify: object) -> object:
    """Function to filter the sequence metadata table with any filter applied in the GUI.

    Args:
        read_data_to_modify (object): Path to the read storage.

    Returns:
        object: Dask dataframe with the applied filters
    """
    # load the sequence metadata, merge with gbif taxonomy and gbif validation
    sequence_metadata = dd.read_hdf(read_data_to_modify, key="sequence_metadata")

    # add gbif taxonomy
    try:
        gbif_taxonomy = dd.read_hdf(read_data_to_modify, key="gbif_taxonomy")
        gbif_taxonomy = gbif_taxonomy.dropna().reset_index()
        sequence_metadata = sequence_metadata.merge(
            gbif_taxonomy, left_on="hash_idx", right_on="hash_idx", how="left"
        )
    except KeyError:
        pass

    try:
        # add gbif validation for sequences and sequence_groups
        gbif_validation_species = dd.read_hdf(
            read_data_to_modify, key="gbif_validation_species"
        )
        # merge both to the sequence metadata
        sequence_metadata = sequence_metadata.merge(
            gbif_validation_species, left_on="hash_idx", right_on="hash_idx", how="left"
        ).rename(columns={"gbif_validation": "gbif_validation_species"})
        # for saving in parquet cannot have mixed bool columns
        sequence_metadata["gbif_validation_species"] = sequence_metadata[
            "gbif_validation_species"
        ].astype(str)
    except KeyError:
        pass

    try:
        gbif_validation_groups = dd.read_hdf(
            read_data_to_modify, key="gbif_validation_species_groups"
        )

        # merge both to the sequence metadata
        sequence_metadata = sequence_metadata.merge(
            gbif_validation_groups, left_on="hash_idx", right_on="hash_idx", how="left"
        ).rename(columns={"gbif_validation": "gbif_validation_species_groups"})
        sequence_metadata["gbif_validation_species_groups"] = sequence_metadata[
            "gbif_validation_species_groups"
        ].astype(str)
    except KeyError:
        pass

    # let the user select columns to filter
    columns_to_filter = st.multiselect("Apply filter to:", sequence_metadata.columns)

    # show a dropdown for each selected filter
    filters = {}

    for col in columns_to_filter:
        unique_values = sequence_metadata[col].dropna().unique()
        dtype = sequence_metadata[col].dtype
        if dtype == "string" or dtype == "bool":
            selected_values = st.multiselect(
                f"Keep records that where {col} matches", options=unique_values
            )
            filters[col] = selected_values
        if dtype == "int" or dtype == "float":

            min_value = st.number_input(
                f"Keep records of {col} which are greater or equal than",
                min_value=sequence_metadata[col].min().compute(),
                max_value=sequence_metadata[col].max().compute(),
                value=sequence_metadata[col].min().compute(),
            )

            max_value = st.number_input(
                f"Keep records of {col} which are smaller or equal than",
                min_value=min_value,
                max_value=sequence_metadata[col].max().compute(),
                value=sequence_metadata[col].max().compute(),
            )

            filters[col] = (min_value, max_value)

    # only perform if any filter is selected
    if not filters:
        return sequence_metadata
    # filter the sequence metadata, generate a mask first
    mask = None

    for column, condition in filters.items():

        # skip empty filters
        if not condition:
            continue
        if isinstance(condition, list) and not condition:
            continue
        if isinstance(condition, tuple) and len(condition) != 2:
            continue

        # this handles numeric values
        if isinstance(condition, tuple) and len(condition) == 2:
            col_mask = (sequence_metadata[column] >= condition[0]) & (
                sequence_metadata[column] <= condition[1]
            )
        # this handles all categories
        elif isinstance(condition, list) and condition:
            col_mask = sequence_metadata[column].isin(condition)
        else:
            # unsupported types go here for now
            continue
        # add everything to the mask
        mask = col_mask if mask is None else mask & col_mask

    if mask is not None:
        return sequence_metadata[mask]
    else:
        return sequence_metadata


def create_data_split(read_data_to_modify: object) -> list:
    sample_metadata = dd.read_hdf(read_data_to_modify, key="sample_metadata")

    columns_to_select = sample_metadata.select_dtypes(
        include=["object", "string"]
    ).columns

    help = """This option will only display columns that do not hold numeric data, as splitting on those will result
    in one read table per sample. Leave empty if a full read table is needed."""

    levels_to_split = st.multiselect(
        label="Columns to split read table on:", options=columns_to_select, help=help
    )

    # collect the sample filters here
    if levels_to_split:
        unique_rows = sample_metadata[levels_to_split].drop_duplicates().compute()
        sample_splits = unique_rows.to_dict(orient="records")
        return sample_splits
    else:
        return []


def export_read_tables(
    read_data_to_modify: str,
    sequence_metadata: object,
    sample_split: list,
    metadata_to_export: list,
    sequence_type: str,
    output_format: str,
):
    """Function to export read tables with the selected filters split by the selected sample metadata.

    Args:
        read_data_to_modify (str): Path to the read storage.
        sequence_metadata (object): Dask dataframe that holds the (filtered) sequence metadata
        sample_split (list): List of dictionarys for filter to apply when splitting output frames by sample metadata.
        metadata_to_export (list): Sequence metadata columns to export.
        sequence_type (str): Wether to export sequences or sequence groups.
        output_format (str): String that defines the output format
    """
    # load the respective count data
    if sequence_type == "sequences":
        read_count_data = dd.read_hdf(read_data_to_modify, key="read_count_data")
    else:
        read_count_data = dd.read_hdf(
            read_data_to_modify, key="sequence_group_read_count_data"
        )

    # add the hash_idx manually
    metadata_to_export = ["hash_idx"] + metadata_to_export

    # only keep the selected columns
    sequence_metadata = sequence_metadata[metadata_to_export]

    # merge with the filtered sequence metadata
    read_count_data = read_count_data.merge(
        sequence_metadata, left_on="hash_idx", right_on="hash_idx", how="right"
    ).dropna(subset=["sample_idx"])

    # merge with the sample metadata
    sample_metadata = dd.read_hdf(read_data_to_modify, key="sample_metadata")
    read_count_data = read_count_data.merge(
        sample_metadata, left_on="sample_idx", right_on="sample_idx", how="left"
    )

    # compute data needed for sorting
    hash_idx_order = read_count_data[["hash_idx", "read_count"]]
    hash_idx_order = (
        hash_idx_order.groupby(["hash_idx"])
        .sum()
        .reset_index()
        .sort_values("read_count", ascending=False)
    )
    hash_idx_order = hash_idx_order.assign(order=1)
    hash_idx_order["order"] = hash_idx_order["order"].cumsum() - 1
    hash_idx_order = hash_idx_order.drop("read_count", axis=1)
    hash_idx_order = hash_idx_order.compute()

    # create a filtering for all subframes if splitting is enabled
    if sample_split:
        masks = []
        for sub_table in sample_split:
            mask = None
            for column, condition in sub_table.items():
                col_mask = read_count_data[column] == condition
                mask = col_mask if mask is None else mask & col_mask
            masks.append(mask)
        frames = [read_count_data[mask] for mask in masks]

        # generate folder names here, as they correspond directly to the maskes
        folder_names = []
        for sub_table in sample_split:
            name = []
            for column, condition in sub_table.items():
                name.append("_".join([column, condition]))
            folder_names.append(name)
        folder_names = ["-".join(name) for name in folder_names]
    else:
        frames = [read_count_data]
        folder_names = ["read_table"]

    # the read tables will be saved here
    save_folder = Path(read_data_to_modify).parents[1]

    # start computing the actual read tables
    for frame, folder in zip(frames, folder_names):
        # find out how many samples are in the current frame
        n_samples = frame[["sample_idx", "sample_name"]].drop_duplicates()

        # compute the batch size to repartition frames if needed
        if output_format == "Excel":
            n_chunks = 10_000_000 // len(hash_idx_order)
            n_chunks = len(n_samples) // n_chunks + 1
        if output_format == "Parquet":
            n_chunks = 100_000_000 // len(hash_idx_order)
            n_chunks = len(n_samples) // n_chunks + 1

        frame = frame.repartition(npartitions=n_chunks)

        # go over the individual partitions and create read tables
        for idx, frame_part in enumerate(frame.to_delayed()):
            # load the samples for that chunk
            frame_df = (
                dd.from_delayed(frame_part)[["sample_idx", "sample_name"]]
                .drop_duplicates()
                .compute()
            )

            # compute all sequence data to generate the current sub table
            all_sequence_data = sequence_metadata.compute()

            sample_sequence_matrix = pd.merge(frame_df, all_sequence_data, how="cross")

            # fetch the read counts for the current part
            read_counts_part = read_count_data[
                read_count_data["sample_idx"].isin(frame_part["sample_idx"])
            ].compute()

            # merge the readcounts into the sample sequence matrix
            sample_sequence_matrix = sample_sequence_matrix.merge(
                read_counts_part[["sample_idx", "hash_idx", "read_count"]],
                on=["sample_idx", "hash_idx"],
                how="left",
            )

            # fill missing values with 0
            sample_sequence_matrix["read_count"] = (
                sample_sequence_matrix["read_count"].fillna(0).astype(int)
            )

            # find out which columns to use for pivot
            pivot_index = [
                col
                for col in sample_sequence_matrix.columns.to_list()
                if col not in ["sample_idx", "sample_name", "read_count"]
            ]

            sample_sequence_matrix[pivot_index] = sample_sequence_matrix[
                pivot_index
            ].fillna("")

            # pivot the table to create a classic samples x sequences table
            wide_read_table = sample_sequence_matrix.pivot_table(
                index=pivot_index,
                columns="sample_name",
                values="read_count",
                fill_value=0,
            )

            # flatten the column index, reset index for saving
            wide_read_table.columns.name = None
            wide_read_table = wide_read_table.reset_index()

            # add the sorting table to sort the read table correctly
            wide_read_table = wide_read_table.merge(
                hash_idx_order,
                left_on="hash_idx",
                right_on="hash_idx",
                how="left",
            )

            # sort the table, drop the empty seq
            wide_read_table = (
                wide_read_table.sort_values(by="order", axis=0)
                .drop(labels=["hash_idx", "order"], axis=1)
                .head(-1)
            )

            # create an output folder for the export
            output_folder = save_folder.joinpath(folder)
            output_folder.mkdir(parents=False, exist_ok=True)

            # write the output file
            if output_format == "Excel":
                output_file = output_folder.joinpath(f"{folder}_read_table.xlsx")
                wide_read_table.to_excel(output_file, index=False)
            else:
                output_file = output_folder.joinpath(
                    f"{folder}_read_table.parquet.snappy"
                )
                wide_read_table.to_parquet(output_file)

            # update the status
            st.toast(f"{output_file.name} has been saved.")


def main():
    # prevent page from scroling up on click
    st.markdown(
        """
    <style>
        * {
        overflow-anchor: none !important;
        }
    </style>""",
        unsafe_allow_html=True,
    )

    # add helptext
    help = """This module can be used to export read tables in either parquet or excel format.
    The read tables can be seperated by metadata fields and will contain the GBIF taxonomy and GBIF
    validation column if those steps where performed. The metadata to be display / exported can
    also be selected."""

    # define the title
    st.title("Export read tables", help=help)
    st.header("Select filters for sequence metadata")

    # potentially filter the metadata. If no filter is selected, the full table will be returned.
    sequence_metadata = filter_sequence_metadata(
        st.session_state["read_data_to_modify"]
    )
    st.divider()

    st.header("Split output tables by sample metadata")
    sample_splits = create_data_split(st.session_state["read_data_to_modify"])

    st.divider()
    # let the user decide the output format
    st.header("Choose output format")
    output_file_format = st.radio(
        label="Output format",
        options=["Excel", "Parquet"],
        index=1,
        horizontal=True,
        key="file_format",
    )

    st.divider()

    st.header("Export read tables")
    output_sequence_type = st.radio(
        label="Choose sequence type",
        options=["sequences", "sequence groups"],
        index=0,
        horizontal=True,
        key="sequence_type",
    )

    # remove hash idx from the output column, as it must be present anyways
    column_choices = [col for col in sequence_metadata.columns if col != "hash_idx"]

    # sequence metadata to include in output table
    metadata_to_export = st.multiselect(
        label="Select sequence metadata to include in the output table",
        options=column_choices,
        default=column_choices,
    )

    export = st.button(label="Export read tables", type="primary")

    if export:
        # export read tables in the apscale analyze folder
        with st.spinner("Computing read tables. Hang on."):
            export_read_tables(
                st.session_state["read_data_to_modify"],
                sequence_metadata,
                sample_splits,
                metadata_to_export,
                output_sequence_type,
                output_file_format,
            )


main()

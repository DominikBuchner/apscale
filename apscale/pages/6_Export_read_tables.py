import streamlit as st
import dask.dataframe as dd


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
    gbif_taxonomy = dd.read_hdf(read_data_to_modify, key="gbif_taxonomy")
    gbif_taxonomy = gbif_taxonomy.dropna().reset_index()
    sequence_metadata = sequence_metadata.merge(
        gbif_taxonomy, left_on="hash_idx", right_on="hash_idx", how="left"
    )

    # add gbif validation for sequences and sequence_groups
    gbif_validation_species = dd.read_hdf(
        read_data_to_modify, key="gbif_validation_species"
    )
    gbif_validation_groups = dd.read_hdf(
        read_data_to_modify, key="gbif_validation_species_groups"
    )

    # merge both to the sequence metadata
    sequence_metadata = sequence_metadata.merge(
        gbif_validation_species, left_on="hash_idx", right_on="hash_idx", how="left"
    ).rename(columns={"gbif_validation": "gbif_validation_species"})

    # merge both to the sequence metadata
    sequence_metadata = sequence_metadata.merge(
        gbif_validation_groups, left_on="hash_idx", right_on="hash_idx", how="left"
    ).rename(columns={"gbif_validation": "gbif_validation_species_groups"})

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
    print(read_data_to_modify)


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
        index=1,
        horizontal=True,
        key="sequence_type",
    )

    # sequence metadata to include in output table
    metadata_to_export = st.multiselect(
        label="Select sequence metadata to include in the output table",
        options=sequence_metadata.columns,
        default=sequence_metadata.columns,
    )

    export = st.button(label="Export read tables", type="primary")

    if export:
        # export read tables in the apscale analyze folder
        export_read_tables(
            st.session_state["read_data_to_modify"],
            sequence_metadata,
            sample_splits,
            metadata_to_export,
            output_sequence_type,
            output_file_format,
        )


main()

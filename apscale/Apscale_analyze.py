import streamlit as st
import apscale, subprocess, sys, shutil
from pathlib import Path
import pandas as pd


def initialize_read_storage(read_data_storage_path: str, project: str) -> None:
    """Functin to initialize the read storage for the analyze module.
    If it already has been initialized, nothing happens, otherwise the read storage
    from step 11 will be copied here.

    Args:
        read_data_storage_path (str): Path to the read storage
        project (str): Path to the project we are working in
    """
    # check if the file exists already
    read_data_storage_name = read_data_storage_path.name
    read_data_to_modify = project.joinpath("12_analyze", "data", read_data_storage_name)

    # add the read_data_to_modify to the session state
    if "read_data_to_modify" not in st.session_state:
        st.session_state["read_data_to_modify"] = read_data_to_modify

    if read_data_to_modify.is_file():
        pass
    else:
        shutil.copyfile(read_data_storage_path, read_data_to_modify)


def reset_metadata(read_data_storage_path, read_data_to_modify, project):
    """Function to reset the metadata file.

    Args:
        read_data_storage_path (_type_): Path to the original read store.
        read_data_to_modify (_type_): Path to the new read store to generate.
    """
    shutil.copyfile(read_data_storage_path, read_data_to_modify)

    # remove all old parquet files from the data folder
    data_folder = Path(project).joinpath("12_analyze", "data")
    files = data_folder.glob("*.parquet.snappy")

    for file in files:
        file.unlink()


def collect_metadata_columns(read_data_to_modify: str) -> tuple:
    """Function to collect the number of columns from the metadata tables.

    Args:
        read_data_to_modify (str): Path to the read data storage with added metadata

    Returns:
        tuple: Number of sample metadata columns, number of sequence metadata columns
    """
    sample_cols, sequence_cols = 0, 0

    with pd.HDFStore(read_data_to_modify, mode="r") as store:
        if "/sample_metadata" in store.keys():
            sample_cols = len(store.get_storer("sample_metadata").table.colnames) - 2
        if "/sequence_metadata" in store.keys():
            sequence_cols = (
                len(store.get_storer("sequence_metadata").table.colnames) - 3
            )

    return sample_cols, sequence_cols


def main(project=Path.cwd()):
    # configure the sidebar
    st.title("Welcome to the Apscale analysis module")
    try:
        project = Path(sys.argv[1])
    except IndexError:
        project = project

    # add the project to the session state so it can be accessed on other pages
    if "project" not in st.session_state:
        st.session_state["project"] = project

    # show the current project
    st.write(f"Current project: **{project.name}**")

    # calculate the number of samples and sequences in the project
    project_name = "_".join(project.name.split("_")[:-1])
    read_data_storage_path = project.joinpath(
        "11_read_table", "data", f"read_data_storage_{project_name}.h5.lz"
    )

    # add the read data storage path to the session state
    if "read_data_storage_path" not in st.session_state:
        st.session_state["read_data_storage_path"] = read_data_storage_path
    # collect number of sequences and number of samples
    with pd.HDFStore(read_data_storage_path, mode="r") as store:
        number_of_sequences = store.get_storer("sequence_data").nrows
        number_of_samples = store.get_storer("sample_data").nrows
        # also collect the number of sequence groups if there are any
        try:
            number_of_sequence_groups = store.get_storer("sequence_group_data").nrows
        except KeyError:
            number_of_sequence_groups = 0

    # if the data storage is not already in the 12_analyze folder, put a copy there
    initialize_read_storage(read_data_storage_path, project)

    st.write(
        f"This project contains **{number_of_sequences - 1}** sequences and **{number_of_samples}** samples."
    )
    st.write(
        f"This project contains **{number_of_sequence_groups - 1}** sequence groups."
    )

    st.write(
        "This module can be used to **add metadata** to your read storage and **export read tables**."
    )

    st.write(
        "This module will grow over time with additional functions we consider useful for metabarcoding analysis."
    )

    # add a divier to display metadata stats
    st.divider()

    # try to collect the number of sample metadata columns and sequence metadata columns
    sample_metadata_cols, sequence_metadata_cols = collect_metadata_columns(
        st.session_state["read_data_to_modify"]
    )

    st.write(
        f"This project currently contains **{sample_metadata_cols} sample metadata columns** and **{sequence_metadata_cols} sequence metadata columns**."
    )

    # add optional data (GBIF species names, GBIF validation)
    with pd.HDFStore(st.session_state["read_data_to_modify"], mode="r") as store:
        keys = {
            "/gbif_taxonomy": "not performed",
            "/gbif_validation_species": "not performed",
            "/gbif_validation_species_groups": "not performed",
        }

        for key in keys.keys():
            if key in store.keys():
                keys[key] = "performed"

        message = f"""GBIF taxonomy module **was {keys["/gbif_taxonomy"]}**.  
        GBIF validation for species **was {keys["/gbif_validation_species"]}**.  
        GBIF validation for species groups **was {keys["/gbif_validation_species_groups"]}**."""

        st.write(message)

    # add a divier to reset all metadata
    st.divider()

    # add the button to reset the metadata
    data_reset = st.button(
        label="Reset all metadata",
        help="Resets all metadata by creating a clean copy of read storage. Old metadata will be lost!",
        type="primary",
    )

    if data_reset:
        reset_metadata(
            st.session_state["read_data_storage_path"],
            st.session_state["read_data_to_modify"],
            project,
        )
        st.info("Metadata has been reset", icon="ℹ️")
        st.rerun()


if __name__ == "__main__":
    main()

import streamlit as st
import duckdb, sys, shutil, time
from pathlib import Path
import pandas as pd


def initialize_read_store(read_data_store_path, project):
    """Functin to initialize the read store for the analyze module.
    If it already has been initialized, nothing happens, otherwise the read store
    from step 11 will be copied here.

    Args:
        read_data_store_path (str): Path to the read store
        project (str): Path to the project we are working in
    """
    read_data_store_name = read_data_store_path.name
    read_data_to_modify = project.joinpath("12_analyze", "data", read_data_store_name)

    # add the read_data_to_modify to the session state
    if "read_data_to_modify" not in st.session_state:
        st.session_state["read_data_to_modify"] = read_data_to_modify
    if read_data_to_modify.is_file():
        pass
    else:
        shutil.copyfile(read_data_store_path, read_data_to_modify)


def get_basic_stats(read_data_to_modify: str) -> tuple:
    """Function to collect basic stats from an apscale project.

    Args:
        read_data_to_modify (str): The read store to look into.

    Returns:
        tuple: Returns a tuple of nr_of_samples, nr_of_sequences, nr_of_sequence_groups
    """
    # establish the connection
    read_data_store = duckdb.connect(read_data_to_modify)

    # collect the number of samples and sequences
    nr_of_samples = read_data_store.execute(
        "SELECT COUNT(*) FROM sample_data"
    ).fetchone()[0]

    nr_of_sequences = (
        read_data_store.execute("SELECT COUNT(*) FROM sequence_data").fetchone()[0] - 1
    )
    nr_of_groups = read_data_store.execute(
        "SELECT COUNT(*) FROM group_data"
    ).fetchone()[0]

    read_data_store.close()

    return nr_of_samples, nr_of_sequences, nr_of_groups


def get_metadata_info(read_data_to_modify: str) -> tuple:
    """Function to get metadata info from the read store.

    Args:
        read_data_to_modify (str): Read store to work with.

    Returns:
        tuple: sample metadata columns, sequence metadata columns
    """
    # connect to the duckdb database
    read_data_store = duckdb.connect(read_data_to_modify)

    (
        sample_metadata_cols,
        sequence_metadata_cols,
        perf_gbif_taxonomy,
        perf_gbif_validation,
    ) = (0, 0, "not performed", "not performed")

    # try to collect the sample metadata
    try:
        sample_metadata = read_data_store.execute(
            "SELECT * FROM sample_metadata LIMIT 1"
        ).df()
        sample_metadata_cols = len(sample_metadata.columns) - 2
    except duckdb.CatalogException:
        pass

    # try to collect the sequence metadata
    try:
        sequence_metadata = read_data_store.execute(
            "SELECT * FROM sequence_metadata LIMIT 1"
        ).df()
        sequence_metadata_cols = len(sequence_metadata.columns) - 5

        # check if one of the gbif validation modules has been run already
        if "gbif_taxonomy" in sequence_metadata.columns:
            perf_gbif_taxonomy = "performed"
        if "gbif_validation" in sequence_metadata.columns:
            perf_gbif_validation = "performed"
    except duckdb.CatalogException:
        pass

    read_data_store.close()

    return (
        sample_metadata_cols,
        sequence_metadata_cols,
        perf_gbif_taxonomy,
        perf_gbif_validation,
    )


def reset_metadata(read_data_storage_path, read_data_to_modify, project):
    """Function to reset the metadata file.

    Args:
        read_data_storage_path (_type_): Path to the original read store.
        read_data_to_modify (_type_): Path to the new read store to generate.
    """
    shutil.copyfile(read_data_storage_path, read_data_to_modify)


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

    # get the duckdb database path
    project_name = "_".join(project.name.split("_")[:-1])
    read_data_store_path = project.joinpath(
        "11_read_table", "data", f"{project_name}_read_data_storage.duckdb"
    )

    # add the read data storage path to the session state
    if "read_data_storage_path" not in st.session_state:
        st.session_state["read_data_store_path"] = read_data_store_path

    # initialize (copy) the read store to the analysis module, to have a file to play with.
    initialize_read_store(read_data_store_path, project)

    # initialize the duckdb connection to gather basic stats
    nr_of_samples, nr_of_sequences, nr_of_groups = get_basic_stats(
        st.session_state["read_data_to_modify"]
    )

    st.write(
        f"This project contains **{nr_of_samples}** samples, **{nr_of_sequences}** sequences and **{nr_of_groups}** sequence groups."
    )

    st.write(
        "This module can be used to **add metadata** to your read store and **export read tables**."
    )

    st.write(
        "This module will grow over time with additional functions we consider useful for metabarcoding analysis."
    )

    # add a divier to display metadata stats
    st.divider()

    (
        sample_metadata_cols,
        sequence_metadata_cols,
        perf_gbif_taxonomy,
        perf_gbif_validation,
    ) = get_metadata_info(st.session_state["read_data_to_modify"])

    st.write(
        f"This project currently contains **{sample_metadata_cols} sample metadata fields** and **{sequence_metadata_cols} sequence metadata fields**."
    )

    message = f"""
        GBIF taxonomy module was **{perf_gbif_taxonomy}**.  
        GBIF validation was **{perf_gbif_validation}**.  
        """

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
            st.session_state["read_data_store_path"],
            st.session_state["read_data_to_modify"],
            project,
        )
        st.info("Metadata has been reset", icon="ℹ️")
        time.sleep(1)
        st.rerun()


if __name__ == "__main__":
    main()

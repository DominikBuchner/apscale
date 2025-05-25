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

    # if the data storage is not already in the 12_analyze folder, put a copy there
    initialize_read_storage(read_data_storage_path, project)

    st.write(
        f"This project contains **{number_of_sequences}** sequences and **{number_of_samples}** samples."
    )

    st.write(
        "This module can be used to **add metadata** to your read storage and **export read tables**."
    )

    st.write(
        "This module will grow over time with additional functions we consider useful for metabarcoding analysis."
    )


if __name__ == "__main__":
    main()

import streamlit as st
import duckdb, sys, shutil
from pathlib import Path


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


if __name__ == "__main__":
    main()

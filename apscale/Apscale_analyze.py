import streamlit as st
import apscale, subprocess, sys
from pathlib import Path
import pandas as pd


def main(project=Path.cwd()):
    # configure the sidebar
    st.title("Welcome to the Apscale analysis module")
    try:
        project = Path(sys.argv[1])
    except IndexError:
        project = project

    # show the current project
    st.write(f"Current project: **{project.name}**")

    # calculate the number of samples and sequences in the project
    project_name = "_".join(project.name.split("_")[:-1])
    read_data_storage_path = project.joinpath(
        "11_read_table", "data", f"read_data_storage_{project_name}.h5.lz"
    )
    # collect number of sequences and number of samples
    with pd.HDFStore(read_data_storage_path, mode="r") as store:
        number_of_sequences = store.get_storer("sequence_data").nrows
        number_of_samples = store.get_storer("sample_data").nrows

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

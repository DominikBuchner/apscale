import streamlit as st
from pathlib import Path
import pandas as pd


def get_file_type_and_read(file_name: str):
    """Function to guess the file type from the file name and return a dataframe with the data

    Args:
        file_name (str): _description_
    """
    extension = "".join(Path(file_name.name).suffixes)

    # define the allowed extensions
    if extension in [".parquet.snappy", ".csv", ".xlsx"]:
        if extension == ".parquet.snappy":
            return pd.read_parquet(file_name)
        elif extension == ".csv":
            return pd.read_csv(file_name)
        elif extension == ".xlsx":
            return pd.read_excel(file_name)
    else:
        return False


def main():
    # define the title
    st.title("Add sample metadata")

    # define the help text
    help_text = """Metadata tables should have one column with sample IDs that are similar to
                   the sample names in the read table. Any number of additional columns can be
                   added as well (e.g. sampling time, longitude, lattitude, habitat, treatments...)"""

    # add some text for context
    st.text(
        "Metadata can be added by uploading a table that holds all metadata information.",
        help=help_text,
    )
    sample_metadata = st.file_uploader(
        "Metadata table upload:", type=[".parquet.snappy", ".csv", ".xlsx"]
    )

    if sample_metadata:
        # read the sample metadata as dataframe
        sample_metadata = get_file_type_and_read(sample_metadata)

        # ask for the sample identifier
        sample_idx = st.selectbox(
            "Please select the column with the sample identifier:",
            options=sample_metadata.columns,
        )


main()

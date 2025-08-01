import streamlit as st
import pandas as pd
import duckdb
from pathlib import Path


def sample_metadata_present(read_data_to_modify: str) -> bool:
    """Function to check if metadata is already present in the dataset.

    Args:
        read_data_to_modify (str): Read data store to work with.

    Returns:
        bool: True if present
    """
    # establish the duckdb connection
    read_data_store = duckdb.connect(read_data_to_modify)

    # look for the sample metadata table
    try:
        read_data_store.execute("SELECT * FROM sample_metadata LIMIT 1")
    except duckdb.CatalogException:
        return False


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


def collect_column_names(metadata_parquet: str) -> list:
    # initialize a duckdb connection to memory
    memory = duckdb.connect(":memory:")
    column_names = memory.execute(
        f"SELECT * FROM read_parquet('{metadata_parquet}') LIMIT 1"
    ).df()
    memory.close()
    return column_names.columns.to_list()


def get_missing_values(
    read_data_to_modify: str, metadata_parquet: str, sample_identifier: str
):
    """Function that return missing values for the selected sample identifier

    Args:
        read_data_to_modify (str): The read store that holds the data.
        metadata_parquet (str): Uploaded parquet that holds teh metadata.
        sample_identifier (str): string with the sample identifier.
    """
    # connect to the database
    read_data_store = duckdb.connect(read_data_to_modify)

    try:
        missing_values = read_data_store.execute(
            f"""
            SELECT sd.sample
            FROM sample_data AS sd
            JOIN read_parquet('{metadata_parquet}') AS mp
                ON sd.sample = mp.{sample_identifier}
            WHERE sd.sample IS DISTINCT FROM mp.{sample_identifier}
            """
        ).df()
    except duckdb.ConversionException:
        st.write("Please select a column that contains strings.")

    return missing_values


def main():
    # define the title
    st.title("Add sample metadata")

    # define the help text
    help_text = """Metadata tables should have one column with sample IDs that are similar to
                   the sample names in the read table. Any number of additional columns can be
                   added as well (e.g. sampling time, longitude, lattitude, habitat, treatments...)"""

    # try to find old metadata
    project = Path(st.session_state["project"])

    metadata_parquet = project.joinpath(
        "12_analyze", "data", "sample_metadata.parquet.snappy"
    )

    # check if sample metadata is present
    metadata_present = sample_metadata_present(st.session_state["read_data_to_modify"])
    parquet_present = metadata_parquet.is_file()

    # if no metadata is present and no parquet is present display the upload button
    if not metadata_present and not parquet_present:
        # display the upload button / handle the upload
        st.text(
            "Metadata can be added by uploading a table that holds all metadata information.",
            help=help_text,
        )
        # collect the sample metadata from file upload
        sample_metadata = st.file_uploader(
            "Metadata table upload:", type=[".parquet.snappy", ".csv", ".xlsx"]
        )

        if sample_metadata:
            with st.spinner("Processing input data:"):
                sample_metadata = get_file_type_and_read(sample_metadata)
                sample_metadata.to_parquet(metadata_parquet)
                st.rerun()

    if parquet_present:
        # notify that the metadata has been found
        st.success("Sample metadata successfully uploaded!", icon="✅")
        # collect the column names
        metadata_fields = collect_column_names(metadata_parquet)
        # add a button to reupload sample metadata
        if st.button("Reupload metadata", type="primary"):
            # delete the old metadata, rerun to go back to upload
            metadata_parquet.unlink()
            st.rerun()
            # give structure to the page
        st.divider()
        # add a selectbox for the sample identifier
        st.subheader("Select and modify metadata")
        # select the sample identifier
        sample_identifier = st.selectbox(
            "Select the sample identifier column:", options=metadata_fields
        )
        # when a sample identifier is selected, check if all values from sample data can be found
        if sample_identifier:
            missing_values = get_missing_values(
                st.session_state["read_data_to_modify"],
                metadata_parquet,
                sample_identifier,
            )


main()

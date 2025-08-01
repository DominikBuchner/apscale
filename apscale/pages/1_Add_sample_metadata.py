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
            WHERE sd.sample NOT IN (
                SELECT {sample_identifier}
                FROM read_parquet('{metadata_parquet}')
            )
            """
        ).df()
    except (duckdb.ConversionException, duckdb.BinderException):
        st.write("Please select a column that contains strings.")
        # return a dummy df
        return pd.DataFrame([[1]], columns=["A"])
    if not missing_values.empty:
        st.text("Samples IDs missing in the provided metadata column:")
        st.write(missing_values)
    else:
        st.success(
            "The provided identifier matches the samples in the read store!",
            icon="✅",
        )

    return missing_values


def display_column_datatypes(metadata_parquet: str):
    # read the first row of the parquet via duckdb to collect the datatypes and names
    memory = duckdb.connect(":memory:")
    memory.execute(
        f"CREATE VIEW temp_view AS SELECT * FROM read_parquet('{metadata_parquet}')"
    )
    info = memory.execute("PRAGMA table_info('temp_view')").fetchdf()
    col_to_types = dict(zip(info["name"], info["type"]))
    human_read_types = {
        "VARCHAR": "string",
        "DOUBLE": "floating point number",
        "BIGINT": "integer",
        "BOOLEAN": "boolean",
        "TIMESTAMP_NS": "timestamp",
    }

    # invert for user selected back conversion
    inverted_types = {v: k for k, v in human_read_types.items()}

    # store the fields to include here
    include_columns = {}
    selected_dtypes = {}

    for field_name in info["name"]:
        col1, col2, col3 = st.columns([2, 2, 2])
        with col1:
            st.write(f"**{field_name}**")
        with col2:
            include_columns[field_name] = st.checkbox(
                "Include", value=True, key=field_name
            )
        with col3:
            selected_dtypes[field_name] = st.multiselect(
                "Data type",
                options=inverted_types.keys(),
                key=f"{field_name}_type",
                max_selections=1,
                default=human_read_types[col_to_types[field_name]],
                disabled=not include_columns[field_name],
                placeholder="Choose datatype",
            )

    # generate a preview
    preview_columns = info["name"].to_list()
    preview_columns = ", ".join(
        f'"{name}"' for name in preview_columns if include_columns[name]
    )
    preview = memory.execute(
        f"""
        SELECT {preview_columns} FROM read_parquet('{metadata_parquet}') LIMIT 5
        """
    ).df()
    st.write("**Current selection:**")
    st.write(preview)

    print(selected_dtypes)
    # return if no field is empty
    if all(selected_dtypes.values()):
        return True
    else:
        return False


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

            if missing_values.empty:
                st.write("Please select the correct datatype for each column.")
                valid_input = display_column_datatypes(metadata_parquet)


main()

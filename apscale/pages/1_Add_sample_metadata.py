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
        return True
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
    column_names = [
        col for col in column_names.columns.to_list() if not col.startswith("__")
    ]
    return column_names


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
                WHERE "{sample_identifier}" IS NOT NULL
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
        f"CREATE TEMPORARY VIEW temp_view AS SELECT * FROM read_parquet('{metadata_parquet}')"
    )
    info = memory.execute("PRAGMA table_info('temp_view')").fetchdf()
    info = info.loc[~info["name"].str.startswith("__")]
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
    st.write("**Current selection preview:**")
    st.write(preview)

    # return if no field is empty
    if all(selected_dtypes.values()):
        # transform single selection list values from multiselect back to regular values
        selected_dtypes = {k: v[0] for k, v in selected_dtypes.items()}
        selected_dtypes = {k: inverted_types[v] for k, v in selected_dtypes.items()}
        selected_dtypes = {
            k: selected_dtypes[k]
            for k in selected_dtypes.keys()
            if k in preview_columns
        }
        return True, selected_dtypes
    else:
        return False, {}


def save_to_read_store(
    read_data_to_modify: str,
    metadata_parquet: str,
    columns_dict: dict,
    sample_identifier: str,
):
    """Function to save the uploaded metadata to the read store.

    Args:
        read_data_to_modify (str): Path to the read data to modify.
        metadata_parquet(str): Path to the uploaded parquet metadata.
        columns_dict (dict): Dict that holds column names and data types
        sample_identifier (str): Column that holds the sample names matching the sample data
    """
    # establish the connection
    read_data_store = duckdb.connect(read_data_to_modify)

    # create the view into the parquet file
    select_clause = ", \n".join(
        [f'CAST("{col}" AS {dtype}) AS "{col}"' for col, dtype in columns_dict.items()]
    )

    read_data_store.execute(
        f"""
    CREATE OR REPLACE TEMPORARY VIEW parquet_data AS
    SELECT {select_clause}
    FROM read_parquet('{metadata_parquet}')
    """
    )

    # build a selector for the columns included
    columns_selector = ", ".join(
        ["sample_idx", "sample"]
        + [f'"{key}"' for key in columns_dict.keys() if key != sample_identifier]
    )

    # join with the sample idx
    read_data_store.execute(
        f"""
    CREATE OR REPLACE TABLE sample_metadata AS 
    (
        SELECT {columns_selector} FROM sample_data AS sd
        LEFT JOIN parquet_data AS pd
            ON sd.sample = pd.{sample_identifier}
    )
    """
    )

    # remove the parquet
    metadata_parquet.unlink()

    # close the connection in the end
    read_data_store.close()


def display_metadata_preview(read_data_to_modifiy, n_rows) -> pd.DataFrame:
    """Function to display a preview of the metadata

    Args:
        read_data_to_modifiy (_type_): Path to the read store.
        n_rows (_type_): Number of rows to return.

    Returns:
        pd.DataFrame: Returns the preview as a dataframe to directly display in streamlit
    """
    read_data_to_modify = duckdb.connect(read_data_to_modifiy)
    preview = read_data_to_modify.execute(
        f"SELECT * FROM sample_metadata LIMIT {n_rows}"
    ).df()
    read_data_to_modify.close()

    return preview


def reset_metadata_table(read_data_to_modify):
    """Function to reset the metadata from the duckdb table.

    Args:
        read_data_to_modify (_type_): Path to the read data store.
    """
    read_data_to_modify = duckdb.connect(read_data_to_modify)
    # remove the sample metadata
    read_data_to_modify.execute("DROP TABLE sample_metadata")
    read_data_to_modify.close()


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
                # here we get the columns and dtypes for storing the data in duckdb
                valid_input, columns_dict = display_column_datatypes(metadata_parquet)

                if valid_input:
                    save_to_store = st.button("Save to read data store", type="primary")
                else:
                    st.button("Save to read data store", type="primary", disabled=True)
                    save_to_store = False

                if save_to_store:
                    save_to_read_store(
                        st.session_state["read_data_to_modify"],
                        metadata_parquet,
                        columns_dict,
                        sample_identifier,
                    )
                    st.rerun()

    if metadata_present:
        # display metadata
        st.write("Sample metadata is already present in the read data store.")
        preview_rows = st.number_input("Rows to display in preview", min_value=5)

        if preview_rows:
            preview = display_metadata_preview(
                st.session_state["read_data_to_modify"], preview_rows
            )
            st.write(preview)

        # button to reset metadata
        reset_metadata = st.button(
            "Reset sample metadata",
            help="This will remove the current metadata table",
            type="primary",
        )

        if reset_metadata:
            reset_metadata_table(st.session_state["read_data_to_modify"])
            st.rerun()


main()

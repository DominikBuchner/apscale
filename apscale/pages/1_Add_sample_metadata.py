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


def check_valid_sample_identifier(sample_names: list) -> list:
    """Function if all sample names from the read storage can be found in the selected sample identifier.

    Args:
        sample_names (list): list of the sample names from the sample metadata

    Returns:
        list: Returns an empty list if all sample names are found, return the missing names otherwise
    """
    # collect all sample names from the read storage
    storage_sample_names = pd.read_hdf(
        st.session_state["read_data_storage_path"],
        key="sample_data",
        columns=["sample_name"],
    )
    difference = set(storage_sample_names["sample_name"].unique()).difference(
        sample_names
    )

    return difference


def add_data_to_read_storage(
    read_data_storage_path: str,
    read_data_to_modify: str,
    data_to_add: str,
    sample_identifier: str,
):
    """Function to add the metadata to the read storage.

    Args:
        read_data_storage_path (str): Path to the original read_store
        read_data_to_modify (str): Path to the copy of the read store that can be modified
        data_to_add (str): the data to add to the read store
        sample_identifier (str): The sample identifier that identifies the sample in the data to add
    """
    # read the data from the storage
    sample_data = pd.read_hdf(read_data_storage_path, key="sample_data")

    # preserve the index
    sample_data = sample_data.reset_index()

    # create the new data by a merge
    updated_data = sample_data.merge(
        data_to_add,
        left_on="sample_name",
        right_on=sample_identifier,
        how="left",
    ).drop(columns=[sample_identifier])

    # reset the index
    updated_data = updated_data.set_index("sample_idx")

    # remove the sample data from the hdf, add the updated data
    with pd.HDFStore(read_data_to_modify) as store:
        if "sample_metadata" in store.keys():
            del store["sample_metadata"]

    # add the updated data
    updated_data.to_hdf(
        read_data_to_modify, key="sample_metadata", format="table", data_columns=True
    )

    # return true on success
    return True


def main():
    # define the title
    st.title("Add sample metadata")

    # define the help text
    help_text = """Metadata tables should have one column with sample IDs that are similar to
                   the sample names in the read table. Any number of additional columns can be
                   added as well (e.g. sampling time, longitude, lattitude, habitat, treatments...)"""

    # try to find old metadata
    project = Path(st.session_state["project"])
    sample_metadata_path = project.joinpath(
        "12_analyze", "data", "sample_metadata.parquet.snappy"
    )
    try:
        sample_metadata = pd.read_parquet(sample_metadata_path)
        # notify that the metadata has been found
        st.success("Sample metadata found!", icon="✅")
        # add a button to reupload sample metadata
        if st.button("Reupload metadata", type="primary"):
            # delete the old metadata, rerun to go back to upload
            sample_metadata_path.unlink()
            st.rerun()
        # give structure to the page
        st.divider()
        # add a selectbox for the sample identifier
        st.subheader("Select and modify metadata")
        # selec the sample identifier
        sample_identifier = st.selectbox(
            "Select the sample identifier column:", options=sample_metadata.columns
        )
        # check if all sample names present in the hdf are present in the metadata folder
        if sample_identifier:
            missing_samples = check_valid_sample_identifier(
                sample_metadata[sample_identifier]
            )
            if missing_samples:
                st.warning(
                    "The provided sample IDs do not match the read storage.", icon="⚠️"
                )
                st.text("Missing samples ids:")
                # create missing samples as markdown
                missing_samples_text = ""
                for sample_id in missing_samples:
                    missing_samples_text += f"- {sample_id}\n"
                st.write(missing_samples_text)
            else:
                st.success("All sample IDs match the read storage.", icon="✅")
                # select the metadata to be added
                metadata_fields = [
                    col for col in sample_metadata.columns if col != sample_identifier
                ]
                metadata_fields = st.multiselect(
                    "Select all metadata fields to add to the read storage:",
                    options=metadata_fields,
                    default=metadata_fields,
                )
                st.divider()
                # tab the input fields
                st.subheader("File preview")
                # define the columns to preview
                preview_columns = [sample_identifier] + metadata_fields
                sample_metadata = sample_metadata[preview_columns]
                st.write(sample_metadata)

                # more structure to the page
                st.divider()

                # finally add the read data to the read storage
                if st.button("Add sample metadata to read storage", type="primary"):
                    # create a copy of the read store in folder 12 and add the data
                    success = add_data_to_read_storage(
                        st.session_state["read_data_storage_path"],
                        st.session_state["read_data_to_modify"],
                        sample_metadata,
                        sample_identifier,
                    )
                    if success:
                        st.success("Data saved successfully.", icon="✅")

    except FileNotFoundError:
        # add some text for context
        st.text(
            "Metadata can be added by uploading a table that holds all metadata information.",
            help=help_text,
        )
        # collect the sample metadata from file upload
        sample_metadata = st.file_uploader(
            "Metadata table upload:", type=[".parquet.snappy", ".csv", ".xlsx"]
        )

        # once sample metadata has been dropped, remove the file uploader
        if sample_metadata:
            # save it to parquet to keep it between runs
            sample_metadata = get_file_type_and_read(sample_metadata)
            sample_metadata.to_parquet(sample_metadata_path)
            st.rerun()


main()

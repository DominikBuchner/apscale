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


def check_valid_sequence_identifier(sequence_ids: list) -> list:
    """Function if all sequence ids from the read storage can be found in the selected sequence identifier.

    Args:
        sequence_ids (list): list of the sequence ids from the sequence metadata

    Returns:
        list: Returns an empty list if all sequence ids are found, return the missing ids otherwise
    """
    # collect all sequence ids from the read storage
    storage_sequence_ids = pd.read_hdf(
        st.session_state["read_data_storage_path"],
        key="sequence_data",
        columns=["hash"],
    )

    difference = set(storage_sequence_ids["hash"].unique()).difference(sequence_ids)

    if difference == {"empty_seq"}:
        return set()
    else:
        return difference


def add_data_to_read_storage(
    read_data_storage_path: str,
    read_data_to_modify: str,
    data_to_add: str,
    sequence_identifier: str,
):
    """Function to add the metadata to the read storage.

    Args:
        read_data_storage_path (str): Path to the original read_store
        read_data_to_modify (str): Path to the copy of the read store that can be modified
        data_to_add (str): the data to add to the read store
        sample_identifier (str): The sequence identifier that identifies the sample in the data to add
    """
    # read the data from the storage
    sequence_data = pd.read_hdf(read_data_storage_path, key="sequence_data")

    # preserve the index
    sequence_data = sequence_data.reset_index()

    # create the new data by a merge
    updated_data = sequence_data.merge(
        data_to_add,
        left_on="hash",
        right_on=sequence_identifier,
        how="left",
    ).drop(columns=[sequence_identifier])

    # reset the index
    updated_data = updated_data.set_index("hash_idx")

    # remove the sample data from the hdf, add the updated data
    with pd.HDFStore(read_data_to_modify) as store:
        del store["sequence_data"]

    # add the updated data
    updated_data.to_hdf(
        read_data_to_modify, key="sequence_data", format="table", data_columns=True
    )

    # return true on success
    return True


def main():
    # define the title
    st.title("Add sequence metadata")

    # define the help text
    help_text = """Metadata tables should have one column with sequence IDs that are similar to
                   the sequence IDs in the read table. Any number of additional columns can be
                   added as well (e.g. taxonomy, trait data, ...)"""

    # try to find old metadata
    project = Path(st.session_state["project"])
    sequence_metadata_path = project.joinpath(
        "12_analyze", "data", "sequence_metadata.parquet.snappy"
    )
    try:
        sequence_metadata = pd.read_parquet(sequence_metadata_path)
        # notify that the metadata has been found
        st.success("Sequence metadata found!", icon="✅")
        # add a button to reupload sample metadata
        if st.button("Reupload metadata", type="primary"):
            # delete the old metadata, rerun to go back to upload
            sequence_metadata_path.unlink()
            st.rerun()
        # give structure to the page
        st.divider()
        # add a selectbox for the sample identifier
        st.subheader("Select and modify metadata")
        # selec the sequence identifier
        sequence_identifier = st.selectbox(
            "Select the sequence identifier column", options=sequence_metadata.columns
        )
        # check if all sample names present in the hdf are present in the metadata folder
        if sequence_identifier:
            missing_sequences = check_valid_sequence_identifier(
                sequence_metadata[sequence_identifier]
            )
            if missing_sequences:
                st.warning(
                    "The provided sequence IDs do not match the read storage.", icon="⚠️"
                )
                st.text("Missing sequence ids:")
                # create missing sequences as markdown
                missing_sequences_text = ""
                for sequence_id in missing_sequences:
                    missing_sequences_text += f"- {sequence_id}\n"
                st.write(missing_sequences_text)
            else:
                st.success("All sequence IDs match the read storage", icon="✅")
                # select the metadata to be added
                metadata_fields = [
                    col
                    for col in sequence_metadata.columns
                    if col != sequence_identifier
                ]
                metadata_fields = st.multiselect(
                    "Select all metadata fields to add to the read storage",
                    options=metadata_fields,
                    default=metadata_fields,
                )
                st.divider()
                # tab the input fields
                st.subheader("File preview")
                # define the columns to preview
                preview_columns = [sequence_identifier] + metadata_fields
                sequence_metadata = sequence_metadata[preview_columns]
                st.write(sequence_metadata)

                # more structure to the page
                st.divider()

                # finally add the read data to the read storage
                if st.button("Add sequence metadata to read storage", type="primary"):
                    # create a copy of the read store in folder 12 and add the data
                    success = add_data_to_read_storage(
                        st.session_state["read_data_storage_path"],
                        st.session_state["read_data_to_modify"],
                        sequence_metadata,
                        sequence_identifier,
                    )
                    if success:
                        st.success("Data saved successfully.", icon="✅")

    except FileNotFoundError:
        # add some text for context
        st.text(
            "Metadata can be added by uploading a table that holds all metadata information.",
            help=help_text,
        )
        # collect the sequence metadata from file upload
        sequence_metadata = st.file_uploader(
            "Metadata table upload:", type=[".parquet.snappy", ".csv", ".xlsx"]
        )

        # once sequence metadata has been dropped, remove the file uploader
        if sequence_metadata:
            # save it to parquet to keep it between runs
            sequence_metadata = get_file_type_and_read(sequence_metadata)
            sequence_metadata.to_parquet(sequence_metadata_path)
            st.rerun()


main()

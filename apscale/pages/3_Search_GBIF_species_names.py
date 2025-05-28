import streamlit as st
import pandas as pd


def seq_metadata_available(read_data_to_modify: str) -> bool:
    """Function to quickly check if the sequence metadata is available.

    Returns:
        bool: True if available, else False
    """
    with pd.HDFStore(read_data_to_modify) as store:
        if "/sequence_metadata" in store.keys():
            return store.get_storer("sequence_metadata").table.colnames
        else:
            return []


def file_preview(read_data_to_modify: str, species_col: str) -> object:
    """Function to generate a filtered df

    Args:
        read_data_to_modify (str): Path to the sequence metadata

    Returns:
        object: Dataframe with the selected column
    """
    sequence_metadata = pd.read_hdf(
        read_data_to_modify, key="sequence_metadata", columns=[species_col]
    )

    return sequence_metadata


def main():
    # header
    st.title("Query the GBIF backbone to harmonize species names")

    # check if sequene metadata is available
    col_names = seq_metadata_available(st.session_state["read_data_to_modify"])
    # remove index from col names
    col_names = [name for name in col_names if name != "index"]

    if col_names:
        st.write("Please select the column that holds the species names.")
        species_column = st.selectbox("Species column label", options=col_names)
        if species_column:
            # generate a file preview
            preview = file_preview(
                st.session_state["read_data_to_modify"], species_column
            )
            st.header("File preview")
            st.write(preview.iloc[1:, :].head(50))
            # get some basic stats
            st.write(
                f"The selected column contains {preview[species_column].nunique()} unique values."
            )
            start_query = st.button("Query GBIF API", type="primary")
            if start_query:
                # extract a list of values for to query the API for [name1, name2, name3]
                pass
                # loop over the list, display progress, return dict with new names

                # map the name to the original dataframe and save it under new key
    else:
        st.write("Please add sequence metadata first.")


main()

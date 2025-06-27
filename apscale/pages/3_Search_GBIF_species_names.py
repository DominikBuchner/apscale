import streamlit as st
import pandas as pd
from pygbif import species
import numpy as np


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


def query_gbif(species_name: str) -> str:
    """Function to query the GBIF species backbone to correct names against it.

    Args:
        species_name (str): Species name as a string.

    Returns:
        str: Species name from GBIF backbone
    """
    api_resp = species.name_backbone(species_name, timeout=60)

    # try to fetch the new name if there is any, else return an empty string
    try:
        new_name = api_resp["species"]
    except KeyError:
        new_name = np.nan

    return new_name


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
                f"The selected column contains **{preview[species_column].nunique()}** unique values."
            )
            start_query = st.button("Query GBIF API", type="primary")
            if start_query:
                # extract a list of values for to query the API for [name1, name2, name3]
                species_names = preview[species_column].dropna().unique()

                # loop over the list, display progress, return dict with new names
                pbar = st.progress(0, text="Looking up GBIF names.")

                # gather results here
                results = {}
                for i, name in enumerate(species_names):
                    results[name] = query_gbif(name)
                    progress = int((i + 1) / len(species_names) * 100)
                    pbar.progress(progress, text=f"Processing {name}")

                # map the name to the original dataframe and save it under new key
                output = preview[[species_column]].copy()
                output["gbif_taxonomy"] = output[species_column].map(results)
                output = output.drop(columns=species_column)

                # save to read storage under a new key
                output.to_hdf(
                    st.session_state["read_data_to_modify"],
                    key="gbif_taxonomy",
                    mode="a",
                    format="table",
                )

                # give user output
                st.success("GBIF taxonomy successfully saved to read storage.")
    else:
        st.write("Please add sequence metadata first.")


main()

import streamlit as st
import dask.dataframe as dd


def filter_sequence_metadata(read_data_to_modify: object) -> dict:
    # load the sequence metadata, merge with gbif taxonomy and gbif validation
    sequence_metadata = dd.read_hdf(read_data_to_modify, key="sequence_metadata")

    # let the user select columns to filter
    columns_to_filter = st.multiselect("Apply filter to:", sequence_metadata.columns)

    # show a dropdown for each selected filter
    filters = {}

    for col in columns_to_filter:
        unique_values = sequence_metadata[col].dropna().unique()
        dtype = sequence_metadata[col].dtype
        if dtype == "string":
            selected_values = st.multiselect(
                f"Keep records that where {col} matches", options=unique_values
            )
        if dtype == "int" or dtype == "float":

            min_value = st.number_input(
                f"Keep records of {col} which are greater or equal than",
                min_value=sequence_metadata[col].min().compute(),
                max_value=sequence_metadata[col].max().compute(),
            )

            max_value = st.number_input(
                f"Keep records of {col} which are smaller or equal than",
                min_value=min_value,
                max_value=sequence_metadata[col].max().compute(),
            )

    st.write(sequence_metadata.head(30))


def main():
    # prevent page from scroling up on click
    st.markdown(
        """
    <style>
        * {
        overflow-anchor: none !important;
        }
    </style>""",
        unsafe_allow_html=True,
    )

    # add helptext
    help = """This module can be used to export read tables in either parquet or excel format.
    The read tables can be seperated by metadata fields and will contain the GBIF taxonomy and GBIF
    validation column if those steps where performed. The metadata to be display / exported can
    also be selected."""

    # define the title
    st.title("Export read tables", help=help)
    st.header("Select filters for sequence metadata")
    sequence_metadata_filter = filter_sequence_metadata(
        st.session_state["read_data_to_modify"]
    )

    st.divider()

    st.header("Split output tables by sample metadata")

    st.divider()
    # let the user decide the output format
    st.header("Choose output format")
    output_file_format = st.radio(
        label="Output format", options=["Excel", "Parquet"], index=1, horizontal=True
    )

    st.header("Export the read tables as:")
    st.divider()

    st.header("Export read tables")


main()

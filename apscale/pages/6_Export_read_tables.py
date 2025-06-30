import streamlit as st


def main():
    # add helptext
    help = """This module can be used to export read tables in either parquet or excel format.
    The read tables can be seperated by metadata fields and will contain the GBIF taxonomy and GBIF
    validation column if those steps where performed. The metadata to be display / exported can
    also be selected."""

    # define the title
    st.title("Export read tables", help=help)


main()

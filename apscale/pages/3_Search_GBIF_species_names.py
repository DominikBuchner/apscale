import streamlit as st
import duckdb
import pandas as pd
from pygbif import species
from joblib import Parallel, delayed


def sequence_metadata_available(read_data_to_modify: str) -> bool:
    """Function to check if sequence metadata is available

    Args:
        read_data_to_modify (str): Read data store to work with

    Returns:
        bool: True if available, else False
    """
    read_data_to_modify = duckdb.connect(read_data_to_modify)

    # look for the sample metadata tabke
    try:
        info = read_data_to_modify.execute(
            "PRAGMA table_info('sequence_metadata')"
        ).df()
        info = info.loc[info["type"] == "VARCHAR"]
        string_columns = info["name"].to_list()
        read_data_to_modify.close()
        return string_columns
    except duckdb.CatalogException:
        read_data_to_modify.close()
        return []


def generate_file_preview(
    read_data_to_modify: str, species_column: str
) -> pd.DataFrame:
    # create the duckdb connection
    read_data_to_modify = duckdb.connect(read_data_to_modify)
    # select all columns that are not null
    preview_table = read_data_to_modify.execute(
        f"""
        SELECT 
            DISTINCT smd.{species_column}
        FROM sequence_metadata as smd
        WHERE 
            smd."{species_column}" IS NOT NULL
            AND smd.hash != 'dummy_hash'
        ORDER BY smd.{species_column}
        LIMIT 50
        """
    ).df()

    species_count = read_data_to_modify.execute(
        f"""
        SELECT 
            COUNT(DISTINCT smd.{species_column})
        FROM sequence_metadata as smd
        WHERE 
            smd."{species_column}" IS NOT NULL
            AND smd.hash != 'dummy_hash'
        """
    ).fetchone()[0]

    st.write(f"There are **{species_count} distinct values** in the selected column.")

    return preview_table


def api_request(species_name: str) -> str:
    api_response = species.name_backbone(species_name, timeout=60)

    # try to fetch the new name, else return a NULL value
    try:
        species_name = api_response["species"]
    except KeyError:
        species_name = pd.NA

    return species_name


def query_gbif(read_data_to_modify: str, species_column: str):
    # create a temp path to save intermediate results
    temp_folder = read_data_to_modify.parent.parent.joinpath("temp")
    temp_folder.mkdir(exist_ok=True)

    # connect to database
    read_data_to_modify = duckdb.connect(read_data_to_modify)

    preview_table = read_data_to_modify.execute(
        f"""
        SELECT 
            DISTINCT smd.{species_column}
        FROM sequence_metadata as smd
        WHERE 
            smd."{species_column}" IS NOT NULL
            AND smd.hash != 'dummy_hash'
        """
    )

    chunk_count = 1

    # loop over the data in chunks, unknown how large this may grow
    while True:
        # fetch a chunk
        chunk = read_data_to_modify.fetch_df_chunk(100)
        if chunk.empty:
            break
        else:
            # create a filename for the chunk
            output_name = temp_folder.joinpath(
                f"gbif_chunk_{chunk_count}.parquet.snappy"
            )
            # execute the API here and save the result to parquet
            species_names = chunk[species_column].to_list()
            corrected_names = Parallel(n_jobs=-2)(
                delayed(api_request)(name) for name in species_names
            )
            # create an output dataframe
            chunk_df = pd.DataFrame(
                data=zip(species_names, corrected_names),
                columns=[species_column, "gbif_taxonomy"],
            )
            # save to intermediate parquet untill all chunks are finished
            chunk_df.to_parquet(output_name)
            chunk_count += 1

    # process parquet files and include them into the database
    parquet_files = temp_folder.joinpath("gbif_chunk_*.parquet.snappy")

    # create a new view across the parquet files
    read_data_to_modify.execute(
        f"""
        CREATE OR REPLACE VIEW gbif_tax_parquet AS
        SELECT * 
        FROM read_parquet('{parquet_files}')
        """
    )

    # drop a potential existing column
    try:
        read_data_to_modify.execute(
            "ALTER TABLE sequence_metadata DROP COLUMN gbif_taxonomy"
        )
    except duckdb.BinderException:
        # add a fresh column
        read_data_to_modify.execute(
            "ALTER TABLE sequence_metadata ADD COLUMN gbif_taxonomy TEXT"
        )

    # add to the sequence metadata
    read_data_to_modify.execute(
        f"""
    UPDATE sequence_metadata AS smd
    SET gbif_taxonomy = gtp.gbif_taxonomy
    FROM gbif_tax_parquet AS gtp
    WHERE smd."{species_column}" = gtp."{species_column}"
    """
    )

    # close the read store
    read_data_to_modify.close()

    # remove parquet files and temp folder
    for file in temp_folder.glob("gbif_chunk_*"):
        if file.is_file():
            file.unlink()

    temp_folder.rmdir()


def check_gbif_taxonomy(read_data_to_modify):
    # connect to database
    read_data_to_modify = duckdb.connect(read_data_to_modify)
    info = read_data_to_modify.execute("PRAGMA table_info(sequence_metadata)").df()
    read_data_to_modify.close()

    if "gbif_taxonomy" in info["name"].to_list():
        return True
    else:
        return False


def display_gbif_taxonomy(read_data_to_modify, n_rows):
    read_data_to_modify = duckdb.connect(read_data_to_modify)
    preview = read_data_to_modify.execute(
        f"""SELECT * FROM sequence_metadata
        WHERE hash != 'dummy_hash' AND gbif_taxonomy IS NOT NULL
        LIMIT {n_rows}"""
    ).df()
    read_data_to_modify.close()

    return preview


def reset_gbif_taxonomy(read_data_to_modify):
    read_data_to_modify = duckdb.connect(read_data_to_modify)
    read_data_to_modify.execute(
        "ALTER TABLE sequence_metadata DROP COLUMN gbif_taxonomy"
    )
    read_data_to_modify.close()


def main():
    st.title("Query the GBIF backbone to harmonize species names")

    # check if sequence metadata is available
    seq_data_col_names = sequence_metadata_available(
        st.session_state["read_data_to_modify"]
    )

    gbif_taxonomy_available = check_gbif_taxonomy(
        st.session_state["read_data_to_modify"]
    )

    if gbif_taxonomy_available:
        # display data, option to reset
        st.write("GBIF taxonomy has already been added.")
        preview_rows = st.number_input("Rows to display in preview", min_value=5)
        st.write(
            display_gbif_taxonomy(
                st.session_state["read_data_to_modify"], n_rows=preview_rows
            )
        )
        reset_tax = st.button("Reset GBIF taxonomy", type="primary")

        if reset_tax:
            reset_gbif_taxonomy(st.session_state["read_data_to_modify"])
            st.rerun()

    # if a proper column is selected create the dropdown and preview
    if not gbif_taxonomy_available and seq_data_col_names:
        st.write("Please select the column that holds the species names.")
        species_column = st.selectbox("Species column name", options=seq_data_col_names)

        if species_column:
            # generate a file preview
            preview = generate_file_preview(
                st.session_state["read_data_to_modify"], species_column
            )

            st.header("File preview")
            st.write(preview)

            start_query = st.button("Query GBIF API", type="primary")

            if start_query:
                # launch the API query
                with st.spinner("Querying GBIF API. Hold on.", show_time=True):
                    query_gbif(st.session_state["read_data_to_modify"], species_column)

                st.rerun()

    if not gbif_taxonomy_available and not sequence_metadata_available:
        st.write("**Please add proper sequence metadata first.**")


main()

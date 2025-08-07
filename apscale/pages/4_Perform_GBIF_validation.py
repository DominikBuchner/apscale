import duckdb, pyproj, time, requests
import pandas as pd
import pathlib
from joblib import Parallel, delayed
import streamlit as st
from pathlib import Path
from shapely.geometry import MultiPoint
from shapely.ops import transform
from shapely.geometry.polygon import orient
from shapely import wkt
from pygbif import occurrences as occ


def check_metadata(read_data_to_modify):
    read_data_to_modify = duckdb.connect(read_data_to_modify)
    # check for sequence metadata
    try:
        read_data_to_modify.execute("SELECT * FROM sequence_metadata LIMIT 1")
        sequence_metadata = True
    except duckdb.CatalogException:
        sequence_metadata = False

    # check for sample metadata
    try:
        read_data_to_modify.execute("SELECT * FROM sample_metadata LIMIT 1")
        sample_metadata = True
    except duckdb.CatalogException:
        sample_metadata = False

    # check for gbif taxonomy
    try:
        read_data_to_modify.execute(
            "SELECT gbif_taxonomy FROM sequence_metadata LIMIT 1"
        )
        gbif_taxonomy = True
    except duckdb.BinderException:
        gbif_taxonomy = False

    read_data_to_modify.close()

    return sequence_metadata, sample_metadata, gbif_taxonomy


def collect_sample_metadata_cols(read_data_to_modify: str) -> list:
    # create the connection
    read_data_to_modify = duckdb.connect(read_data_to_modify)
    info = read_data_to_modify.execute("PRAGMA table_info(sample_metadata)").df()
    read_data_to_modify.close()
    return info["name"].to_list()


def generate_file_preview(
    read_data_to_modify: str, lat_col: str, lon_col: str
) -> pd.DataFrame:
    # create the connection
    read_data_to_modify = duckdb.connect(read_data_to_modify)
    preview_frame = read_data_to_modify.execute(
        f"""
        SELECT 
            sample_idx,
            sample,
            {lat_col},
            {lon_col}
        FROM sample_metadata    
        LIMIT 20
        """
    ).df()

    return preview_frame


def compute_wkt_string(read_data_to_modify, sequence_idx, radius, lat_col, lon_col):
    # establish a new connection
    con = duckdb.connect(read_data_to_modify, read_only=True)
    data = con.execute(
        f"""SELECT
            "{lat_col}", 
            "{lon_col}"
            FROM temp_db.distribution_data
            WHERE sequence_idx = {sequence_idx}
        """
    ).df()

    # extract the coordinates
    coords = list(zip(data[lon_col], data[lat_col]))
    # compute the multiploint
    mp = MultiPoint(coords)
    hull = mp.convex_hull

    # change projection to expand by kilometers
    project = pyproj.Transformer.from_crs(
        "EPSG:4326", "EPSG:3857", always_xy=True
    ).transform
    hull_metric = transform(project, hull)
    buffered = hull_metric.buffer(radius * 1000)

    project_back = pyproj.Transformer.from_crs(
        "EPSG:3857", "EPSG:4326", always_xy=True
    ).transform
    buffered_geo = transform(project_back, buffered)
    buffered_geo = orient(buffered_geo, sign=1).wkt

    con.close()

    return buffered_geo


def compute_species_distributions(read_data_to_modify, lat_col, lon_col, radius):
    # create a temp path to save intermediate results
    temp_folder = read_data_to_modify.parent.parent.joinpath("temp")
    temp_folder.mkdir(exist_ok=True)
    connection_path = read_data_to_modify

    # connect to database
    read_data_to_modify = duckdb.connect(read_data_to_modify)
    temp_db = Path(temp_folder).joinpath("temp.duckdb")
    read_data_to_modify.execute(f"ATTACH '{temp_db}' as temp_db")

    # collect the read count data and join in the required sequence metadata
    read_data_to_modify.execute(
        f"""
        CREATE OR REPLACE TABLE temp_db.distribution_data AS 
        SELECT 
            srcd.sequence_idx,
            srcd.sample_idx, 
            smd.sequence_order,
            smd.gbif_taxonomy,
            samd."{lat_col}" AS lat,
            samd."{lon_col}" AS lon, 
        FROM main.sequence_read_count_data srcd
        LEFT JOIN main.sequence_metadata smd
            ON srcd.sequence_idx = smd.sequence_idx
        LEFT JOIN main.sample_metadata samd
            ON srcd.sample_idx = samd.sample_idx
        WHERE smd.gbif_taxonomy IS NOT NULL 
        """
    )

    # close connection that can be used for writing
    read_data_to_modify.close()

    # open a new read_only connection
    read_only = duckdb.connect(temp_db, read_only=True)

    # fetch all distinct sequence_idx values from there
    distinct_sequences = read_only.execute(
        "SELECT DISTINCT sequence_idx FROM temp_db.distribution_data"
    )
    chunk_count = 1

    # loop over in chunks
    while True:
        data_chunk = distinct_sequences.fetch_df_chunk(1)
        if data_chunk.empty:
            break
        else:
            # process each chunk seperately. data inside the chunk can be processed in parallel
            wkts = Parallel(n_jobs=-2)(
                delayed(compute_wkt_string)(temp_db, idx, radius, lat_col, lon_col)
                for idx in data_chunk["sequence_idx"]
            )
            # create a dataframe with idx wkt
            wkt_df = pd.DataFrame(
                data=zip(data_chunk["sequence_idx"], wkts),
                columns=["sequence_idx", "wkt_string"],
            )
            wkt_df["radius"] = radius
            # save to parquet for intermediate results
            output_file_name = temp_folder.joinpath(
                f"wkt_chunk_{chunk_count}.parquet.snappy"
            )
            chunk_count += 1
            wkt_df.to_parquet(output_file_name)
    # close the read only connection
    read_only.close()

    # ingest parquet files into duckdb for temporary display of data
    parquet_path = temp_folder.joinpath("wkt_*.parquet.snappy")

    # open a connection to the temp.db
    temp_db = duckdb.connect(temp_db)

    # add the wkt data to the temp db
    temp_db.execute(
        f"""
    CREATE OR REPLACE TABLE wkt_data AS
    SELECT * FROM read_parquet("{parquet_path}")
    """
    )

    # remove the temp parquet files
    for file in temp_folder.glob("wkt_*.parquet.snappy"):
        if file.is_file():
            file.unlink()

    # close the connection
    temp_db.close()


def wkt_calculated(temp_db):
    # create the path to the temp db
    if temp_db.is_file():
        temp_connection = duckdb.connect(temp_db)
    else:
        return False
    # try to read from the table
    try:
        wkt_data = temp_connection.execute("SELECT * FROM wkt_data LIMIT 1")
        return True
    except duckdb.CatalogException:
        return False


def generate_map_preview_table(temp_db) -> pd.DataFrame:
    # open the connection
    temp_db_con = duckdb.connect(temp_db)
    # create a view with thejoined tables
    temp_db_con.execute(
        f"""
        CREATE OR REPLACE TEMPORARY VIEW preview_table AS 
        SELECT 
            DISTINCT(wkt.sequence_idx),
            dd.gbif_taxonomy,
            wkt.radius
        FROM wkt_data AS wkt
        LEFT JOIN distribution_data AS dd
            ON wkt.sequence_idx = dd.sequence_idx
        ORDER BY dd.gbif_taxonomy, wkt.sequence_idx
        """
    )

    preview = temp_db_con.execute("SELECT * FROM preview_table").df()

    return preview, preview["sequence_idx"].to_list()


def select_map_data(temp_db, sequence_idx) -> pd.DataFrame:
    # create the duckdb connection
    temp_db_con = duckdb.connect(temp_db)

    occurence_data = temp_db_con.execute(
        f"""
        SELECT 
            lon, 
            lat
        FROM distribution_data AS dd
        WHERE dd.sequence_idx = {sequence_idx}
        """
    ).df()
    # add a label for the map
    occurence_data["label"] = "occurence"

    # extrakt the wkt from the wkt data
    wkt_string = (
        temp_db_con.execute(
            f"""
        SELECT * FROM wkt_data
        WHERE sequence_idx = {sequence_idx}
        """
        )
        .df()["wkt_string"]
        .item()
    )

    wkt_data = list(wkt.loads(wkt_string).exterior.coords)

    wkt_data = pd.DataFrame(data=wkt_data, columns=["lon", "lat"])
    wkt_data["label"] = "distribution"

    # concat both frames to get the map data
    map_data = pd.concat([occurence_data, wkt_data], ignore_index=True)
    # add a color column
    color_dict = {"occurence": "#800000", "distribution": "#008000"}
    map_data["color"] = map_data["label"].map(color_dict)

    return map_data


def validation_api_request(species_name, wkt_string) -> str:
    # perform the api call
    while True:
        try:
            validation_result = occ.search(
                scientificName=species_name, geometry=wkt_string, timeout=60
            )
            break
        except (
            requests.exceptions.HTTPError,
            requests.exceptions.ConnectionError,
        ) as e:
            if e.response.status_code == 429:
                time.sleep(1)
                continue

    if validation_result["count"] > 0:
        return "plausible"
    else:
        return "implausible"


def api_validation(temp_db, temp_folder):
    # establish the connection
    temp_db_con = duckdb.connect(temp_db)

    # collect the required data
    api_data = temp_db_con.execute(
        f"""
        SELECT 
            DISTINCT(wd.sequence_idx),
            dd.gbif_taxonomy,
            wd.wkt_string
        FROM wkt_data AS wd
        LEFT JOIN distribution_data AS dd
            ON wd.sequence_idx = dd.sequence_idx
        """
    )

    chunk_count = 1

    while True:
        api_data_chunk = api_data.fetch_df_chunk()
        if api_data_chunk.empty:
            break
        else:
            # create a filename for the chunk
            output_name = temp_folder.joinpath(
                f"gbif_validation_{chunk_count}.parquet.snappy"
            )
            chunk_args = zip(
                api_data_chunk["gbif_taxonomy"], api_data_chunk["wkt_string"]
            )
            validation_status = Parallel(n_jobs=-2)(
                delayed(validation_api_request)(species_name, wkt_string)
                for species_name, wkt_string in chunk_args
            )
            # some user output
            st.toast(f"Chunk {chunk_count} processed!")

            # save to parquet
            chunk_output = pd.DataFrame(
                zip(api_data_chunk["sequence_idx"], validation_status),
                columns=["sequence_idx", "gbif_validation"],
            )

            # write to parquet to later ingest into the main database
            chunk_output.to_parquet(output_name)

            # increase the chunk count
            chunk_count += 1

    temp_db_con.close()


def add_validation_to_metadata(read_data_to_modify, temp_folder):
    # connect to the main database
    read_data_to_modify_con = duckdb.connect(read_data_to_modify)
    # collect the parquet data
    parquet_files = temp_folder.joinpath("gbif_validation_*.parquet.snappy")

    # create a view over the parquet files
    read_data_to_modify_con.execute(
        f"""
        CREATE OR REPLACE TEMPORARY VIEW validation_parquet AS
        SELECT * 
        FROM read_parquet('{parquet_files}')
        """
    )

    # check if sequence gbif validation is in the database
    info = read_data_to_modify_con.execute("PRAGMA table_info(sequence_metadata)").df()
    info = info["name"].to_list()

    if "gbif_validation" not in info:
        read_data_to_modify_con.execute(
            "ALTER TABLE sequence_metadata ADD COLUMN gbif_validation TEXT"
        )

    # add the validation to the sequence metadata
    read_data_to_modify_con.execute(
        f"""
        UPDATE sequence_metadata AS smd
        SET gbif_validation = vp.gbif_validation
        FROM validation_parquet AS vp
        WHERE smd.sequence_idx = vp.sequence_idx
        """
    )

    # remove the parquet files
    for file in temp_folder.glob("gbif_validation_*.parquet.snappy"):
        if file.is_file():
            file.unlink()

    # infer gbif validation for groups from sequence metadata
    # create a view with joined in gbif validation
    read_data_to_modify_con.execute(
        f"""
        CREATE OR REPLACE TEMPORARY VIEW group_validation_temp AS
        (
            SELECT 
                gm.*,
                smd.gbif_validation
            FROM group_mapping AS gm
            LEFT JOIN sequence_metadata AS smd
                ON gm.sequence_idx = CAST(smd.sequence_idx AS HUGEINT)
            WHERE smd.gbif_validation IS NOT NULL
        )
        """
    )

    # group by group idx, cast plausible on output if any is plausible
    read_data_to_modify_con.execute(
        f"""
        CREATE OR REPLACE TEMPORARY VIEW group_validation AS
        SELECT
            group_idx,
            CASE
                WHEN BOOL_OR(gbif_validation = 'plausible') THEN 'plausible'
                ELSE 'implausible'
            END AS gbif_validation
        FROM group_validation_temp
        GROUP BY group_idx
        """
    )

    # add a new column to the group_metadata column and fill it with the gbif_validation
    # check if sequence gbif validation is in the database
    info = read_data_to_modify_con.execute("PRAGMA table_info(group_metadata)").df()
    info = info["name"].to_list()

    if "gbif_validation" not in info:
        read_data_to_modify_con.execute(
            "ALTER TABLE group_metadata ADD COLUMN gbif_validation TEXT"
        )

    # add the validation to the group metadata
    read_data_to_modify_con.execute(
        f"""
        UPDATE group_metadata AS gmd
        SET gbif_validation = gv.gbif_validation
        FROM group_validation AS gv
        WHERE gmd.sequence_idx = gv.group_idx
        """
    )

    # close the connection
    read_data_to_modify_con.close()


def gbif_validation_performed(read_data_to_modify: str) -> bool:
    # connect to database
    read_data_to_modify_con = duckdb.connect(read_data_to_modify)

    # get the table info
    info = read_data_to_modify_con.execute("PRAGMA table_info(sequence_metadata)").df()
    info = info["name"].to_list()

    read_data_to_modify_con.close()

    if "gbif_validation" in info:
        return True
    else:
        return False


def validation_preview(read_data_to_modify: str, n_rows: int) -> pd.DataFrame:
    # establish the connection
    read_data_to_modify_con = duckdb.connect(read_data_to_modify, read_only=True)

    # get the first n rows
    preview = read_data_to_modify_con.execute(
        f"""
        SELECT 
            sequence_idx,
            gbif_taxonomy,
            gbif_validation
        FROM sequence_metadata
        WHERE gbif_taxonomy IS NOT NULL
        LIMIT {n_rows}
        """
    ).df()

    read_data_to_modify_con.close()

    return preview


def reset_validation(read_data_to_modify: str) -> None:
    # establish the connection
    read_data_to_modify_con = duckdb.connect(read_data_to_modify)

    read_data_to_modify_con.execute(
        "ALTER TABLE sequence_metadata DROP COLUMN gbif_validation"
    )

    read_data_to_modify_con.close()


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

    # define the temp folder to look for intermediate saves
    temp_folder = st.session_state["read_data_to_modify"].parent.parent.joinpath("temp")
    temp_db = Path(temp_folder).joinpath("temp.duckdb")

    # header
    st.title("Validate species names via GBIF record search")
    st.write("This module needs sample metadata (lat, lon) for each sample.")
    st.write("This module needs harmonized GBIF taxonomy.")

    # check for the required data
    sequence_metadata, sample_metadata, gbif_taxonomy = check_metadata(
        st.session_state["read_data_to_modify"]
    )

    if not sequence_metadata:
        st.write("**Please add sequence metadata first.**")
    if not sample_metadata:
        st.write("**Please add sample metadata first.**")
    if not gbif_taxonomy:
        st.write("**Please harmonize taxonomy via GBIF.**")

    # check if wkt is calculated
    wkt_computed = wkt_calculated(temp_db)

    # check if the GBIF validation has been performed
    gbif_validated = gbif_validation_performed(st.session_state["read_data_to_modify"])

    # if all are true display a selector for lat lon
    if (
        sequence_metadata
        and sample_metadata
        and gbif_taxonomy
        and not wkt_computed
        and not gbif_validated
    ):
        sample_meta_cols = collect_sample_metadata_cols(
            st.session_state["read_data_to_modify"]
        )

        # display two columns
        lat, lon = st.columns(2)
        with lat:
            latitude = st.selectbox(
                label="Latitude",
                options=sample_meta_cols,
                help="Please select the latitude column.",
            )
        with lon:
            longitude = st.selectbox(
                label="Longitude",
                options=sample_meta_cols,
                help="Please select the longitude column.",
            )

        st.header("File preview")
        st.write(
            generate_file_preview(
                st.session_state["read_data_to_modify"], latitude, longitude
            )
        )

        # slider for the radius around sample points (growth of the polygon)
        radius = st.slider(
            label="Radius to check around sample",
            min_value=50,
            max_value=500,
            value=200,
            step=50,
        )

        compute_distributions = st.button(
            label="Compute species distributions", type="primary"
        )

        if compute_distributions:
            with st.spinner("Computing species distributions", show_time=True):
                compute_species_distributions(
                    st.session_state["read_data_to_modify"], latitude, longitude, radius
                )
                # rerun to update the interface
                st.rerun()

    # if the wkt strings are already computed, display the map
    if wkt_computed and not gbif_validated:
        st.write("Species distribution data detected!")

        # generate a preview
        preview, idx_selection = generate_map_preview_table(temp_db)
        st.write(preview)
        # select an idx to plot on map
        idx_to_plt = st.selectbox(
            label="Select any hash idx to plot on the map",
            options=idx_selection,
        )

        if idx_to_plt:
            map_data = select_map_data(temp_db, idx_to_plt)

            if not map_data.empty:
                st.map(map_data, latitude="lat", longitude="lon", color="color")

        # option to reset the species distribution data
        compute = st.button(label="Validate species via GBIF API", type="primary")
        reset = st.button(label="Reset species distributions", type="secondary")

        # reset wkt data
        if reset:
            temp_db.unlink()
            st.rerun()

        # perform the GBIF validation algorithm
        if compute:
            with st.spinner("Querying GBIF API. Hold on.", show_time=True):
                api_validation(temp_db, temp_folder)

            # ingest the data into the sequence metadata
            add_validation_to_metadata(
                st.session_state["read_data_to_modify"], temp_folder
            )
            # remove temp db and tempfolder
            temp_db.unlink()
            temp_folder.rmdir()
            st.rerun()

    # if validated data exists
    if gbif_validated:
        st.write("GBIF validation has already been performed.")
        # display preview rows
        preview_rows = st.number_input("Rows to display in preview", min_value=5)
        # display the preview
        st.write(
            validation_preview(
                st.session_state["read_data_to_modify"], n_rows=preview_rows
            )
        )

        reset_vali = st.button("Reset GBIF validation", type="primary")

        if reset_vali:
            reset_validation(st.session_state["read_data_to_modify"])
            st.rerun()


main()

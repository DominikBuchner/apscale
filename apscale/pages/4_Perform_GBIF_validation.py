import duckdb, pyproj, time, requests, random, zipfile, glob, shutil
import pandas as pd
import pathlib
from joblib import Parallel, delayed, load, dump
import streamlit as st
from pathlib import Path
from shapely.geometry import MultiPoint
from shapely.ops import transform
from shapely.geometry.polygon import orient
from shapely import wkt
from pygbif import occurrences as occ
from pygbif.gbifutils import NoResultException
from shapely.ops import unary_union
from streamlit_autorefresh import st_autorefresh


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
            smd.gbif_usage_key,
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
        data_chunk = distinct_sequences.fetch_df_chunk()
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
            dd.gbif_usage_key,
            wkt.radius,
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


def initialize_download(temp_db, temp_folder, username, password, email):
    # establish the connection to the duckdb database
    temp_db_con = duckdb.connect(temp_db, read_only=True)
    temp_db_con.execute("INSTALL spatial; LOAD spatial;")

    minx, miny, maxx, maxy = temp_db_con.execute(
        """
        SELECT
            MIN(ST_XMin(ST_GeomFromText(wkt_string))) AS minx,
            MIN(ST_YMin(ST_GeomFromText(wkt_string))) AS miny,
            MAX(ST_XMax(ST_GeomFromText(wkt_string))) AS maxx,
            MAX(ST_YMax(ST_GeomFromText(wkt_string))) AS maxy
        FROM wkt_data
    """
    ).fetchone()

    # bbox = f"POLYGON(({minx} {miny}, {minx} {maxy}, {maxx} {maxy}, {maxx} {miny}, {minx} {miny}))"

    # load all distinct usage keys
    usage_keys = (
        temp_db_con.execute(
            f"""
        SELECT
            DISTINCT(dd.gbif_usage_key)
        FROM distribution_data AS dd
        """
        )
        .df()["gbif_usage_key"]
        .to_list()
    )

    temp_db_con.close()

    download_predicate = {
        "type": "and",
        "predicates": [
            {
                "type": "in",
                "key": "TAXON_KEY",
                "values": [str(k) for k in usage_keys],
                "matchCase": False,
            },
            # Latitude min
            {
                "type": "greaterThanOrEquals",
                "key": "DECIMAL_LATITUDE",
                "value": miny,
            },
            # Latitude max
            {
                "type": "lessThanOrEquals",
                "key": "DECIMAL_LATITUDE",
                "value": maxy,
            },
            # Longitude min
            {
                "type": "greaterThanOrEquals",
                "key": "DECIMAL_LONGITUDE",
                "value": minx,
            },
            # Longitude max
            {
                "type": "lessThanOrEquals",
                "key": "DECIMAL_LONGITUDE",
                "value": maxx,
            },
        ],
    }

    # submit the download
    download_info = occ.download(
        queries=download_predicate,
        format="SIMPLE_PARQUET",
        user=username,
        pwd=password,
        email=email,
    )

    return download_info


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


def download_initialized(pickle_path):
    if pickle_path.is_file():
        return True
    else:
        return False


def download_info(pickle_path):
    if not download_initialized(pickle_path):
        st.info("There is no active download at the moment.")
        return False, False
    else:
        download_data = load(pickle_path)
        download_key = download_data[0]
        meta = occ.download_meta(download_key)
        created = meta["created"]
        status = meta["status"]
        st.info(
            f"""
                The download key is: {download_key}\n
                The download has been created at: {created}\n
                Current status: {status}\n
                This page will refresh every 20 seconds until the data is ready.
                """
        )

        return status, download_key


def download_gbif_data(download_key, temp_folder):
    # show some status updates
    status = st.status("Downloading and processing data.", expanded=True)
    # download the data to the temp folder, unpack it and load it into a duckdb database
    occ.download_get(download_key, path=temp_folder)
    status.write("Download completed.")
    # unzip and push data into duckdb
    status.write("Unzipping downloaded results.")
    download_file = temp_folder.joinpath(f"{download_key}.zip")

    with zipfile.ZipFile(download_file, "r") as in_stream:
        in_stream.extractall(temp_folder)
    # remove the zip file
    download_file.unlink()

    # Load data into duckDB
    status.write("Loading data into DuckDB.")

    # collect all files to load into duckdb (remove the empty file first)
    for file in temp_folder.joinpath("occurrence.parquet").iterdir():
        if file.stat().st_size == 0:
            file.unlink()

    # load into duckdb
    gbif_database = temp_folder.joinpath("gbif_db.duckdb")
    gbif_con = duckdb.connect(gbif_database)
    gbif_con.execute(
        f"""
        CREATE OR REPLACE TABLE 
            occ_data 
        AS SELECT 
            species,
            decimallatitude AS lat,
            decimallongitude AS lon,
            taxonkey
        FROM read_parquet('{temp_folder.joinpath("occurrence.parquet", ("*"))}')
        """
    )

    # close the connection
    gbif_con.close()

    status.write("Data successfully loaded.")

    # remove the downloaded parquet
    shutil.rmtree(temp_folder.joinpath("occurrence.parquet"))

    # return status widget to use further
    return status


def compute_validation(status):
    status.update(label="Finished", state="complete", expanded=False)


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
    download_pickle = temp_folder.joinpath("download.pkl")

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
            label="Select any sequence idx to plot on the map",
            options=idx_selection,
        )

        if idx_to_plt:
            map_data = select_map_data(temp_db, idx_to_plt)

            if not map_data.empty:
                st.map(map_data, latitude="lat", longitude="lon", color="color")

        st.info(
            """To submit a download to GBIF credentials are required. They will only be used to 
            initialize the download and never be stored. After sending the download the fields will
            be emptied automatically. You will receive an email once the download is finished.""",
            icon="ℹ️",
        )

        # give the input fields for username, password and email
        col1, col2, col3 = st.columns(3)

        with col1:
            user = st.text_input(label="GBIF username")
        with col2:
            pwd = st.text_input(label="GBIF password", type="password")
        with col3:
            email = st.text_input(label="Your mail adress")

        # option to initialize the download
        if user and pwd and email and not download_initialized(download_pickle):
            validate = st.button(
                label="Initialize download of validation data",
                type="primary",
                disabled=False,
            )
        else:
            validate = st.button(
                label="Initialize download of validation data",
                type="primary",
                disabled=True,
            )

        if validate:
            with st.spinner("Initializing download. Hold on.", show_time=True):
                pickle_data = initialize_download(
                    temp_db, temp_folder, user, pwd, email
                )
                # pickle the download data after the download has been requested
                dump(pickle_data, download_pickle)

                # rerun to clear all fields and disable the button
                st.rerun()

        # display download info
        status, download_key = download_info(download_pickle)

        if status == "SUCCEEDED":
            download_process = st.button(
                label="Download and process data", disabled=False, type="primary"
            )
        else:
            download_process = st.button(
                label="Download and process data", disabled=True
            )
            st_autorefresh(interval=20_000)

        # download and process the data, remove download files afterwards
        if download_process:
            # download the data and pass to duckdb
            status = download_gbif_data(download_key, temp_folder)

            # compute the actual validaton
            compute_validation(status)
        st.divider()

        reset = st.button(label="Reset species distributions", type="secondary")

        # reset wkt data
        if reset:
            temp_db.unlink()
            st.rerun()
            # # ingest the data into the sequence metadata
            # add_validation_to_metadata(
            #     st.session_state["read_data_to_modify"], temp_folder
            # )
            # # remove temp db and tempfolder
            # temp_db.unlink()
            # temp_folder.rmdir()
            # st.rerun()

    # if validated data exists
    # if gbif_validated:
    #     st.write("GBIF validation has already been performed.")
    #     # display preview rows
    #     preview_rows = st.number_input("Rows to display in preview", min_value=5)
    #     # display the preview
    #     st.write(
    #         validation_preview(
    #             st.session_state["read_data_to_modify"], n_rows=preview_rows
    #         )
    #     )

    #     reset_vali = st.button("Reset GBIF validation", type="primary")

    #     if reset_vali:
    #         reset_validation(st.session_state["read_data_to_modify"])
    #         st.rerun()


main()

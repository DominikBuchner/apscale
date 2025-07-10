import pyproj, json, time, requests
import dask.dataframe as dd
import streamlit as st
import pandas as pd
from pygbif import occurrences as occ
from pathlib import Path
from shapely.geometry import MultiPoint
from shapely.ops import transform
from shapely.geometry.polygon import orient
from shapely import wkt
import numpy as np


def check_metadata(read_data_to_modify: str) -> dict:
    """Function to check if the required metadata is available.

    Args:
        read_data_to_modify (str): Path to the modified read storage.

    Returns:
        dict: Returns a dict with {sample_metadata: bool, gbif_taxonomy: bool}
    """
    return_dict = {"sample_metadata": False, "gbif_taxonomy": False}

    # check if the required keys are in the read data to modify
    with pd.HDFStore(read_data_to_modify, "r") as hdf_store:
        if "/sample_metadata" in hdf_store.keys():
            return_dict["sample_metadata"] = True
        if "/gbif_taxonomy" in hdf_store.keys():
            return_dict["gbif_taxonomy"] = True

    return return_dict


def sample_metadata_columns(read_data_to_modify: str) -> bool:
    """Function to quickly extract the metadata columns from the read data.

    Returns:
        bool: True if available, else False
    """
    with pd.HDFStore(read_data_to_modify) as store:
        if "/sample_metadata" in store.keys():
            col_names = store.get_storer("sample_metadata").table.colnames
            col_names = [name for name in col_names if name != "index"]
            return col_names
        else:
            return []


def file_preview(read_data_to_modify: str, lat_col: str, lon_col: str) -> object:
    """Function to generate a file preview on the fly for the selected columns

    Args:
        read_data_to_modify (str): Read store to look up the data
        lat_col (str): selected column for latitude
        lon_col (str): selected column for longitude

    Returns:
        object: First few columns of the generated dataframe
    """
    with pd.HDFStore(read_data_to_modify, mode="r") as store:
        row_limit = min(store.get_storer("sample_metadata").nrows, 5)

    preview = pd.read_hdf(
        read_data_to_modify,
        key="sample_metadata",
        start=0,
        stop=row_limit,
        columns=["sample_name", lat_col, lon_col],
    )

    return preview


def compute_polygon_wkt(coordinates: list, radius: int) -> object:
    """Function to compute a polygon from coordinates and expand it by radius

    Args:
        coordinates (list): Coordinates to base the polygon on.
        radius (int): Radius to expand the polygon

    Returns:
        object:
    """
    # generate a multipoint object and compute the convex hull
    mp = MultiPoint(coordinates)
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

    return buffered_geo


def add_polygon_column(dataframe: object, radius: int):
    """Helper function to compute the polygons per row

    Args:
        dataframe (object): Dataframe to perform the operation on.
        radius (_type_): Radius to expand the polygon
    """
    dataframe["polygon_wkt"] = dataframe["lon_lat_list"].map(
        lambda coords: compute_polygon_wkt(coords, radius)
    )

    return dataframe


def compute_species_distributions(read_data_to_modify, lat_col, lon_col, radius):
    # read the gbif species
    gbif_species = dd.read_hdf(read_data_to_modify, key="gbif_taxonomy")
    gbif_species = gbif_species.dropna().reset_index()

    # load the readcount data
    read_count_data = dd.read_hdf(read_data_to_modify, key="read_count_data")

    # compute the readsum data that is needed for sorting in the end
    readsum_data = (
        read_count_data.groupby(by=["hash_idx"])
        .sum()
        .reset_index()[["hash_idx", "read_count"]]
    )

    # merge frames to only keep occurences of the gbif species
    read_count_data = read_count_data.merge(
        gbif_species, left_on="hash_idx", right_on="hash_idx", how="inner"
    )

    # add the sequence metadata to show hash values instead of idx, as idx is not sorted
    sequence_metadata = dd.read_hdf(
        read_data_to_modify, key="sequence_metadata", columns=["hash", "seq"]
    )
    sequence_metadata = sequence_metadata.reset_index()
    read_count_data = read_count_data.merge(
        sequence_metadata, left_on="hash_idx", right_on="hash_idx", how="left"
    )

    # add lon lat data
    sample_metadata = dd.read_hdf(
        read_data_to_modify, key="sample_metadata", columns=[lon_col, lat_col]
    )
    read_count_data = read_count_data.merge(
        sample_metadata, left_on="sample_idx", right_on="sample_idx", how="left"
    )

    read_count_data["lon_lat"] = read_count_data[[lon_col, lat_col]].apply(
        lambda row: (row[lon_col], row[lat_col]), axis=1, meta=("lon_lat", "object")
    )

    # group by hash_idx to collect the data for the polygon
    read_count_data = (
        read_count_data.groupby(
            by=["hash_idx", "gbif_taxonomy", "hash", "seq"], sort=False
        )["lon_lat"]
        .apply(lambda x: list(x), meta=("lon_lat_list", "object"))
        .reset_index()
    )

    meta = {
        "hash_idx": "int64",
        "gbif_taxonomy": "object",
        "hash": "object",
        "seq": "object",
        "read_count": "int64",
        "lon_lat_list": "object",
    }

    meta = read_count_data.head(1).iloc[0:0]
    read_count_data = read_count_data.map_partitions(lambda x: x, meta=meta)

    # add the readsum data to sort by it
    read_count_data = read_count_data.merge(
        readsum_data, left_on="hash_idx", right_on="hash_idx", how="left"
    ).sort_values(by="read_count", axis=0, ascending=False)

    # define the meta for the output
    sample = read_count_data.head(1)
    sample["polygon_wkt"] = ""
    meta = sample.iloc[0:0]

    # add the polygon data
    read_count_data = read_count_data.map_partitions(
        add_polygon_column, radius, meta=meta
    )

    # only strings can be streamed to hdf
    read_count_data["gbif_taxonomy"] = read_count_data["gbif_taxonomy"].astype(str)
    read_count_data["hash"] = read_count_data["hash"].astype(str)
    read_count_data["seq"] = read_count_data["seq"].astype(str)
    read_count_data["lon_lat_list"] = read_count_data["lon_lat_list"].apply(
        json.dumps, meta=("lon_lat_list", "str")
    )
    read_count_data["radius"] = radius

    # save to hdf to use later for the validation of records
    read_count_data.to_hdf(
        read_data_to_modify,
        key="gbif_taxonomy_distribution",
        mode="a",
        format="table",
    )


def extract_map_data(read_data_to_modify: object, hash_idx: int) -> object:
    """Function to extract map data for a given hash_idx to perform some basic plotting

    Args:
        read_data_to_modify (object): Read store to work in.
        hash_idx (int): The hash_idx to return a dataframe for.

    Returns:
        object: A dataframe that can be used for plotting
    """
    hash_idx_data = dd.read_hdf(
        read_data_to_modify,
        key="gbif_taxonomy_distribution",
        columns=["hash_idx", "lon_lat_list", "polygon_wkt"],
    )
    hash_idx_data = hash_idx_data[hash_idx_data["hash_idx"] == hash_idx].compute()

    # compute the dataframe and then transform it correctly
    hash_idx_data["lon_lat_list"] = hash_idx_data["lon_lat_list"].apply(
        lambda x: json.loads(x)
    )

    # restore the points from the wkt polygon
    hash_idx_data["polygon_wkt"] = hash_idx_data["polygon_wkt"].apply(
        lambda x: list(wkt.loads(x).exterior.coords)
    )

    # handle lon lad and polygon data seperatly
    occurence_data = hash_idx_data[["hash_idx", "lon_lat_list"]].explode("lon_lat_list")
    occurence_data["label"] = "occurence"
    occurence_data[["lon", "lat"]] = pd.DataFrame(
        occurence_data["lon_lat_list"].tolist(), index=occurence_data.index
    )
    occurence_data = occurence_data.drop(columns="lon_lat_list")

    distribution_data = hash_idx_data[["hash_idx", "polygon_wkt"]].explode(
        "polygon_wkt"
    )
    distribution_data["label"] = "distribution"
    distribution_data[["lon", "lat"]] = pd.DataFrame(
        distribution_data["polygon_wkt"].tolist(), index=distribution_data.index
    )
    distribution_data = distribution_data.drop(columns="polygon_wkt")

    # concat the map data, assign colors
    map_data = pd.concat([occurence_data, distribution_data], ignore_index=True)
    color_dict = {"occurence": "#800000", "distribution": "#008000"}
    map_data["color"] = map_data["label"].map(color_dict)

    return map_data


def perform_gbif_validation(read_data_to_modify: object):
    """Function to perform the GBIF validation for all records.

    Args:
        read_data_to_modify (object): Path to the HDF read store.
    """
    gbif_distribution_data = dd.read_hdf(
        read_data_to_modify,
        key="gbif_taxonomy_distribution",
        chunksize=1,
        columns=["hash_idx", "gbif_taxonomy", "polygon_wkt"],
    )

    # display the progress bar on the main page
    pbar = st.progress(0, text="Validating via GBIF.")
    nr_of_queries = gbif_distribution_data.npartitions

    results = []

    for i, line in enumerate(gbif_distribution_data.to_delayed()):
        # compute the values for that line
        line = line.compute()
        hash_idx, species, polygon = (
            line["hash_idx"].item(),
            line["gbif_taxonomy"].item(),
            line["polygon_wkt"].item(),
        )
        # perform the api call
        while True:
            try:
                validation_result = occ.search(
                    scientificName=species, geometry=polygon, timeout=60
                )
                break
            except requests.exceptions.HTTPError as e:
                if e.response.status_code == 429:
                    time.sleep(1)
                    continue

        results.append(
            (hash_idx, "True" if validation_result["count"] > 0 else "False")
        )
        progress = int((i + 1) / nr_of_queries * 100)
        pbar.progress(progress, text=f"Processing {hash_idx}: {species}")

    # transform to dataframe, save to hdf
    results = pd.DataFrame(results, columns=["hash_idx", "gbif_validation"])
    results.to_hdf(
        read_data_to_modify,
        key="gbif_validation_species",
        mode="a",
        format="table",
    )

    st.success("GBIF validation results successfully saved.")

    # map the results also for the species groups if there are any
    with pd.HDFStore(read_data_to_modify, mode="r") as store:
        groups_present = False
        if "/sequence_group_data" in store.keys():
            groups_present = True

    # give user output
    if groups_present:
        st.warning(
            "Species groups detected in the dataset. Deriving GBIF validation for species groups.",
            icon="ℹ️",
        )

        # derive for sequence groups, load the group mappings first
        group_mapping_path = read_data_to_modify.parents[2].joinpath(
            "11_read_table", "data", "group_mapping.csv"
        )

        # load gbif validation species and groups, add group annotation
        gbif_validation_species = dd.read_hdf(
            read_data_to_modify, key="gbif_validation_species"
        )
        gbif_validation_species["gbif_validation"] = gbif_validation_species[
            "gbif_validation"
        ].astype(bool)

        group_mapping = dd.read_csv(
            group_mapping_path,
            sep="\t",
            names=["hash_idx", "hash", "seq", "group_idx"],
            usecols=["hash_idx", "group_idx"],
        )

        # merge the two frames to have the groups in the gbif validation
        gbif_validation_sequence_groups = gbif_validation_species.merge(
            group_mapping, left_on="hash_idx", right_on="hash_idx", how="left"
        )

        gbif_validation_sequence_groups = (
            gbif_validation_sequence_groups[["group_idx", "gbif_validation"]]
            .groupby(by="group_idx")["gbif_validation"]
            .sum()
            .reset_index()
        )

        # transform back to boolean, rename column
        gbif_validation_sequence_groups = gbif_validation_sequence_groups.assign(
            gbif_validation=gbif_validation_sequence_groups["gbif_validation"] > 0
        )
        gbif_validation_sequence_groups["gbif_validation"] = (
            gbif_validation_sequence_groups["gbif_validation"].map(
                lambda x: "True" if x else "False", meta=("gbif_validation", "object")
            )
        )

        gbif_validation_sequence_groups = gbif_validation_sequence_groups.rename(
            columns={"group_idx": "hash_idx"}
        )

        gbif_validation_sequence_groups.to_hdf(
            read_data_to_modify,
            key="gbif_validation_species_groups",
            mode="a",
            format="table",
        )

        st.success("GBIF validation results for species groups successfully saved.")


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

    # header
    st.title("Validate species names via GBIF record search")
    st.write("This module needs sample metadata (lat, lon) for each sample.")
    st.write("This module needs harmonized GBIF taxonomy.")

    # look for metadata, look for harmonized species names
    return_dict = check_metadata(st.session_state["read_data_to_modify"])

    # give user feedback
    if not return_dict["sample_metadata"]:
        st.write("**Please add sample metadata first.**")
    if not return_dict["gbif_taxonomy"]:
        st.write("**Please add gbif taxonomy first.**")

    if return_dict["sample_metadata"] and return_dict["gbif_taxonomy"]:
        # extract column names from sample metadata
        col_names = sample_metadata_columns(st.session_state["read_data_to_modify"])

        # selectboxes for lattitude and longitude
        latitude, longitude = st.columns(2)
        with latitude:
            lat = st.selectbox(
                label="Latitude",
                options=col_names,
                help="Please select latitude column",
            )
        with longitude:
            lon = st.selectbox(
                label="Longitude",
                options=col_names,
                help="Please select longitude column",
            )

        # display file preview
        st.header("File preview")
        try:
            preview = file_preview(
                st.session_state["read_data_to_modify"], lat_col=lat, lon_col=lon
            )
            st.write(preview)
        except ValueError:
            st.write("Please select to distinct columns")

        # slider for the radius around sample points (growth of the polygon)
        radius = st.slider(
            label="Radius to check around sample",
            min_value=50,
            max_value=500,
            value=200,
            step=50,
        )

        compute_distributions = st.button(
            label="Compute / Update species distributions", type="primary"
        )

        if compute_distributions:
            with st.spinner("Computing species distributions. Hang on."):
                compute_species_distributions(
                    st.session_state["read_data_to_modify"], lat, lon, radius
                )
            st.success("Distributions have been computed and saved!")

        st.divider()
        # if gbif taxonomy distribution in keys, display this
        with pd.HDFStore(st.session_state["read_data_to_modify"]) as hdf_store:
            keys = hdf_store.keys()

        if "/gbif_taxonomy_distribution" in keys:
            gbif_taxonomy_distribution = dd.read_hdf(
                st.session_state["read_data_to_modify"],
                key="gbif_taxonomy_distribution",
                columns=["hash_idx", "radius", "gbif_taxonomy", "hash", "seq"],
            )

            species_list = gbif_taxonomy_distribution.compute()
            st.write("Detected data:")
            st.write(species_list)

            idx_to_plt = st.selectbox(
                label="Select any hash idx to plot on the map",
                options=species_list["hash_idx"],
            )

            # load the map data
            map_data = extract_map_data(
                st.session_state["read_data_to_modify"], hash_idx=idx_to_plt
            )

            if not map_data.empty:
                st.map(map_data, latitude="lat", longitude="lon", color="color")

            # perform the GBIF validation
            perform_validation = st.button(
                "Validate all species via GBIF API",
                type="primary",
                help="This function checks if there are any observations for the given species list within the selected radius. May take some time.",
            )

            if perform_validation:
                # function to perform the validation. Adds validation data to the read store
                perform_gbif_validation(st.session_state["read_data_to_modify"])


main()

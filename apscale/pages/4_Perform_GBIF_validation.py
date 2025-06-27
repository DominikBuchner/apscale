import pyproj, json
import dask.dataframe as dd
import streamlit as st
import pandas as pd
from pygbif import occurrences as occ
from shapely.geometry import MultiPoint
from shapely.ops import transform
from shapely.geometry.polygon import orient
from shapely import wkt
from shapely.geometry import Polygon


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

    # save to hdf to use later for the validation of records
    read_count_data.to_hdf(
        read_data_to_modify,
        key="gbif_taxonomy_distribution",
        mode="a",
        format="table",
    )

    # also create a duckdb database for lookups


def main():
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
        st.write("**Please add gbif taxonomy**")

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
            label="Compute species distributions", type="primary"
        )

        if compute_distributions:
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
                columns=["hash_idx", "gbif_taxonomy", "hash", "seq"],
            )

            species_list = gbif_taxonomy_distribution.compute()
            st.write(species_list)

            st.selectbox(
                label="Select any hash idx to plot on the map",
                options=species_list["hash_idx"],
            )


main()

#         string = "POLYGON((8.2 46.3, 8.8 46.3, 9.1 46.4, 7.1 48.6, 8.2 46.3))"
#         test = occ.search(scientificName=species, geometry=buffered_geo.wkt, timeout=30)
#         print(species, len(test["results"]))

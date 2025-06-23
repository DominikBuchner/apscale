from pygbif import occurrences as occ
import dask.dataframe as dd
from shapely.geometry import MultiPoint
import pyproj
from shapely.ops import transform
from shapely.geometry.polygon import orient
from shapely import wkt
from shapely.geometry import Polygon

gbif_species = dd.read_hdf(
    "C:\\Users\\Dominik\\Nextcloud\\data_storage_ag_leese\\dominik_buchner\\Apscale modular Update\\bmon24_2_filter_fwh_apscale\\12_analyze\\data\\read_data_storage_bmon24_2_filter_fwh.h5.lz",
    key="gbif_taxonomy",
)

# only keep valid gbif taxonomy
gbif_species = gbif_species.dropna().reset_index()

# load the readcount data
read_count_data = dd.read_hdf(
    "C:\\Users\\Dominik\\Nextcloud\\data_storage_ag_leese\\dominik_buchner\\Apscale modular Update\\bmon24_2_filter_fwh_apscale\\12_analyze\\data\\read_data_storage_bmon24_2_filter_fwh.h5.lz",
    key="read_count_data",
)

# merge the two frames to only look up samples that have occurence data for that species
read_count_data = read_count_data.merge(
    gbif_species, left_on="hash_idx", right_on="hash_idx", how="inner"
)

# add lon / lat data
sample_metadata = dd.read_hdf(
    "C:\\Users\\Dominik\\Nextcloud\\data_storage_ag_leese\\dominik_buchner\\Apscale modular Update\\bmon24_2_filter_fwh_apscale\\12_analyze\\data\\read_data_storage_bmon24_2_filter_fwh.h5.lz",
    key="sample_metadata",
)

read_count_data = read_count_data.merge(
    sample_metadata[["longitude", "latitude"]],
    left_on="sample_idx",
    right_on="sample_idx",
    how="left",
)


read_count_data["lon_lat"] = read_count_data[["longitude", "latitude"]].apply(
    lambda row: (row["longitude"], row["latitude"]),
    axis=1,
    meta=("lon_lat", "object"),
)

result = read_count_data.groupby(by=["hash_idx", "gbif_taxonomy"])["lon_lat"].apply(
    lambda x: list(x), meta=("lon_lat_list", "object")
)

result = result.reset_index()
print(result.compute())
for part in result.to_delayed():
    part = part.compute()

    for _, row in part.iterrows():
        test_hash, test_coords, species = (
            row["hash_idx"],
            row["lon_lat_list"],
            row["gbif_taxonomy"],
        )

        mp = MultiPoint(test_coords)
        hull = mp.convex_hull

        project = pyproj.Transformer.from_crs(
            "EPSG:4326", "EPSG:3857", always_xy=True
        ).transform
        hull_metric = transform(project, hull)
        buffered = hull_metric.buffer(500 * 1000)

        project_back = pyproj.Transformer.from_crs(
            "EPSG:3857", "EPSG:4326", always_xy=True
        ).transform
        buffered_geo = transform(project_back, buffered)
        buffered_geo = orient(buffered_geo, sign=1)

        string = "POLYGON((8.2 46.3, 8.8 46.3, 9.1 46.4, 7.1 48.6, 8.2 46.3))"
        test = occ.search(scientificName=species, geometry=buffered_geo.wkt, timeout=30)
        print(species, len(test["results"]))

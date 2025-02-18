import pandas as pd
from pathlib import Path
import numpy as np


def merge_replicates(
    esv_table: object, replicate_del: str, minimum_presence: int
) -> object:
    """Function to merge replicates based on the replicate delimiter and a minimum number of occurences

    Args:
        esv_table (object): ESV table as a pandas dataframe.
        replicate_del (str): Delimiter that identifies the replicate in the end of a sample name (e.g. _ for Sample_A)
        mininum_presence (int): Minimum number of replicates with reads present to accept as true signal (e.g. 2)

    Returns:
        object: Returns a pandas dataframe with merge replicates.
    """
    # preserve sequences and unique ids
    seqs, unique_ids = esv_table.pop("Seq"), esv_table.pop("unique_ID")

    # set the temporary id as index to transpose the dataframe
    # save a copy of the pre merged table for calculating stats
    pre_merged_table = esv_table.copy()

    # generate the merged table
    merged_table = esv_table.copy().set_index("temporary_ID").T

    # replace 0 values with np.nan, so they are not counted by groupby
    merged_table = merged_table.replace(0, np.nan)

    # remove the replicate names delimiter, so replicates have the same name
    merged_table.index = (
        merged_table.index.str.split(replicate_del).str[:-1].str.join("_")
    )

    # groupby index
    merged_table = merged_table.groupby(merged_table.index, sort=False).sum(
        min_count=minimum_presence
    )

    # rearange the dataframe, reset index, add 0 instead of np.nan again
    merged_table = merged_table.T.reset_index().replace(np.nan, 0)

    # extract the column names
    sample_names = merged_table.columns[1:]

    # add unique ids and seqs again
    merged_table = pd.concat(
        [merged_table.iloc[:, :1], unique_ids, merged_table.iloc[:, 1:], seqs], axis=1
    )

    # calculate removed ESVs and removed reads
    for sample_name in sample_names:
        # sum of reads pre and post filtered
        pre_merged_reads = pre_merged_table.filter(regex=sample_name).to_numpy().sum()
        post_merged_reads = merged_table.filter(regex=sample_name).to_numpy().sum()
        print(sample_name)
        print(post_merged_reads / pre_merged_reads)

def main(project=Path.cwd()) -> None:
    """Main function to merge replicates and remove negative controls from the dataset

    Args:
        project (str, optional): Path to the apscale project. Defaults to Path.cwd().
    """
    # collect the settings from the settings file
    settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="7_replicate_negative_controls",
    )

    # define the settings
    replicate_merging, replicate_del, minimum_presence, substract_ncs, nc_prefix = (
        bool(settings["merge replicates"].item()),
        settings["replicate delimiter"].item(),
        int(settings["minimum replicate presence"].item()),
        bool(settings["substract negative controls"].item()),
        settings["negative control prefix"].item(),
    )

    # read the esv table from the step before from parquet file
    esv_table = pd.read_parquet(
        Path(project).joinpath(
            "8_esv_table",
            "{}_ESV_table.parquet.snappy".format(
                Path(project).stem.replace("_apscale", "")
            ),
        )
    )

    # merge replicates first (optional)
    if replicate_merging:
        merge_replicates(esv_table, replicate_del, minimum_presence)

    # remove negative controls afterward (optional)

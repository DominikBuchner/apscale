import datetime
import pandas as pd
from pathlib import Path
import numpy as np


def merge_replicates(
    esv_table: object, replicate_del: str, minimum_presence: int
) -> tuple:
    """Function to merge replicates based on the replicate delimiter and a minimum number of occurences

    Args:
        esv_table (object): ESV table as a pandas dataframe.
        replicate_del (str): Delimiter that identifies the replicate in the end of a sample name (e.g. _ for Sample_A)
        mininum_presence (int): Minimum number of replicates with reads present to accept as true signal (e.g. 2)

    Returns:
        tuple: Returns a pandas dataframe with merged replicates and a dataframe with merging stats.
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

    # define merging stats
    merging_stats = []

    # calculate removed ESVs and removed reads
    for sample_name in sample_names:
        # sum of reads pre and post filtered
        pre_merged_reads = pre_merged_table.filter(regex=sample_name).to_numpy().sum()
        post_merged_reads = merged_table.filter(regex=sample_name).to_numpy().sum()

        # calculate pre merged esvs and post merged esvs
        pre_merged_esvs = np.count_nonzero(
            pre_merged_table.filter(regex=sample_name).sum(axis=1)
        )
        post_merged_esvs = np.count_nonzero(
            merged_table.filter(regex=sample_name).sum(axis=1)
        )
        merging_stats.append(
            (
                sample_name,
                pre_merged_reads,
                post_merged_reads,
                pre_merged_esvs,
                post_merged_esvs,
            )
        )

        # calculate removed reads
        if pre_merged_reads != 0:
            retained_reads = (
                (pre_merged_reads - post_merged_reads) * 100 / pre_merged_reads
            )
        else:
            retained_reads = 100

        if pre_merged_esvs != 0:
            retained_esvs = (pre_merged_esvs - post_merged_esvs) * 100 / pre_merged_esvs
        else:
            retained_esvs = 100
        # give user output
        print(
            "{}: {}: {:.2f} % of reads removed. {:.2f} % of ESVs removed.".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                sample_name,
                retained_reads,
                retained_esvs,
            )
        )

    # transform merging stats to dataframe
    merging_stats = pd.DataFrame(
        data=merging_stats,
        columns=[
            "sample name",
            "pre merging reads",
            "post merging reads",
            "pre merging esvs",
            "post merging esvs",
        ],
    )

    return merged_table, merging_stats


def substract_neg_controls(esv_table: object, nc_prefix: str) -> tuple:
    """Function to substract the maximum read in all negative controls.

    Args:
        esv_table (object): ESV table as pandas dataframe.
        nc_prefix (str): Prefix to detect the negative controls (e.g. NC_)

    Returns:
        tuple: Returns the ESV table with removed negative controls and the stats for removal.
    """
    # extract the negative controls from the esv table
    ncs = pd.concat(
        [esv_table.pop(col) for col in esv_table.columns if col.startswith(nc_prefix)],
        axis=1,
    )
    # calculate the maximum number of reads in the negative controls
    ncs = ncs.max(axis=1)
    nc_rem_stats = []

    # give user output
    for esv, nc_count in zip(esv_table["temporary_ID"], ncs):
        if nc_count > 0:
            print(
                "{}: {} reads substracted from {}.".format(
                    datetime.datetime.now().strftime("%H:%M:%S"), nc_count, esv
                )
            )
            nc_rem_stats.append((esv, nc_count))

    # generate the nc removal stats table
    nc_rem_stats = pd.DataFrame(data=nc_rem_stats, columns=["esv", "substracted reads"])

    # update the esv table
    esv_table = pd.concat(
        [
            esv_table.iloc[:, :2],
            esv_table.iloc[:, 2:-1].sub(ncs, axis=0).clip(0),
            esv_table.iloc[:, -1],
        ],
        axis=1,
    )

    # return updated esv table and statistics for project report
    return esv_table, nc_rem_stats


def remove_empty_rows(esv_table: object) -> object:
    """Function to remove empty rows from esv table after merging or negative control substraction

    Args:
        esv_table (object): Dataframe holding the esv table.

    Returns:
        object: Dataframe with empty rows removed
    """
    pre_esv_count = len(esv_table.index)
    esv_table["readsum"] = esv_table.iloc[:, 2:-1].sum(axis=1)
    esv_table = esv_table.loc[esv_table["readsum"] > 0]
    esv_table_clean = esv_table.drop("readsum", axis=1)

    print(
        "{}: {} empty lines removed from ESV table.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            pre_esv_count - len(esv_table_clean.index),
        )
    )

    return esv_table_clean


def save_outputs(esv_table: object, project:str, merging_stats=None, nc_rem_stats=None) -> None:
    """Function to write all output of the merging and negative removal step.

    Args:
        esv_table (object): Dataframe holding the updated ESV table.
        project (str): Path to the apscale project
        merging_stats (object, optional): Dataframe with the merging statistics. Defaults to None.
        nc_rem_stats (object, optional): Dataframe with the negative control removal statistics. Defaults to None.
    """
    # write the updated esv table as parquet and excel

    # write the updated fasta file

    # add the merging and nc stats to the project report
    pass

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
        # give user output
        print(
            "{}: Starting replicate merging.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

        esv_table, merging_stats = merge_replicates(
            esv_table, replicate_del, minimum_presence
        )

    # remove negative controls afterward (optional)
    if substract_ncs:
        # give user output
        print(
            "{}: Substracting reads found in negative controls.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

        esv_table, nc_rem_stats = substract_neg_controls(esv_table, nc_prefix)

    # remove empty rows and columns from the esv table
    if replicate_merging or substract_ncs:
        esv_table = remove_empty_rows(esv_table)

        # write outputs and project report

    else:
        print(
            "{}: Replicate merging and negative control substraction disabled. Skipping step.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )

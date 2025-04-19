import zarr, glob, gzip, pickle, shutil, dask, datetime, sys
import dask.dataframe as dd
from zict import File, Buffer, LRU, Func
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd
from apscale.a_create_project import choose_input

def main(project=Path.cwd()):
    """Main function to initially create the read table data.

    Args:
        project (str, optional): Path to the apscale project to work in. Defaults to Path.cwd().
    """
    # collect the settings from the excel sheet
    settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="11_read_table",
    )
    perform = settings["generate read table"].item()
    to_excel = settings["to excel"].item()
    to_parquet = settings["to parquet"].item()
    
    # find out which dataset to load
    prior_step = choose_input(project, "11_generate_read_table")
    
    # check that at least one hashing step has been performed
    if prior_step == "06_dereplication":
            ## give user output
        print(
        '{}: Please perform at least on data aggregation step before creating the read table.'.format(
            datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
        return None
    
    # collect the input for read table generation
    input = glob.glob(str(Path(project).joinpath(prior_step, "data", "*.fasta.gz")))
    
    # index all the data and save as hdf store
    
    # sort the hdf store to generate the fasta file
    
    # create parquet and excel outputs from the hdf
    
    

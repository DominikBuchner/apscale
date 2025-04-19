import dask.dataframe as dd
import zarr, glob, gzip, pickle, shutil, dask
from zict import File, Buffer, LRU, Func
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import dask.array
import pandas as pd

def main(project=Path.cwd()):
    """Main function to initially create the read table data.

    Args:
        project (str, optional): Path to the apscale project to work in. Defaults to Path.cwd().
    """
    settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="11_read_table",
    )
    perform = settings["generate read table"].item()
    to_excel = settings["to excel"].item()
    to_parquet = settings["to parquet"].item()
    
    print(perform, to_excel, to_parquet)

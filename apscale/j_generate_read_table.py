import zarr, glob, gzip, pickle, shutil, dask, datetime, sys, os
import dask.dataframe as dd
from zict import File, Buffer, LRU, Func
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import pandas as pd
from apscale.a_create_project import choose_input


def parse_fasta_data(input_path: str):
    """Generator function that takes an input path to a size annotated fasta file
    and yield the samples name, sequence hash, sequence and the sequence count.
    Data is yielded sequences wise

    Args:
        input_path (str): Input path to the fasta file.

    Yields:
        tuple: Yields a tuple with the sample name, sequence has, sequence and the sequence count.
    """
    # extract the sample name
    sample_name = Path(input_path).with_suffix("").with_suffix("").name
    # extract the data directly from the fasta itself
    with gzip.open(input_path, "rt") as data_stream:
        for header, seq in SimpleFastaParser(data_stream):
            header_data = header.split(";size=")
            seq_hash, seq_count = header_data[0], int(header_data[1])
            seq = seq.upper().strip()
            yield sample_name, seq_hash, seq, seq_count


def generate_buffered_dict(save_path: str, buffer_size: int) -> dict:
    """Function to generate a disk backed dict that can (infinetly) grow.

    Args:
        save_path (str): Path to save the dict to.
        buffer_size (int): How many values to keep in memory before pushing to disk.

    Returns:
        dict: Dict that is backed up if needed
    """
    # define cold and hot storage
    cold = Func(pickle.dumps, pickle.loads, File(Path(save_path)))
    hot = LRU(n=buffer_size, d={})
    buffered_dict = Buffer(hot, cold, n=int(buffer_size * 0.99))

    return buffered_dict


def save_data_lines_to_hdf(save_path: str, data_lines: dict) -> dict:
    """Function to append data lines to the hdf store

    Args:
        save_path (str): Savepath of the output hdf file.
        data_lines (dict): Dict that's holding the data lines

    Returns: Returns an empty dict to free memory space
    """
    # convert the data liens dict to dataframe
    data_lines = pd.DataFrame.from_dict(
        data_lines, orient="index", columns=["hash_idx", "sample_idx", "read_count"]
    )

    # write the data to hdf format
    with pd.HDFStore(
        save_path, mode="a", complib="blosc:blosclz", complevel=9
    ) as hdf_output:
        hdf_output.append(
            key="read_count_data", value=data_lines, format="table", data_columns=True
        )

    return {}


def index_data_to_hdf(project: str, input_files: str, buffer_size: int) -> str:
    """Function to index the dataset with buffered dicts and saving all data as dataframes to hdf.

    Args:
        project (str): Apscale project to work in.
        input_files (str): List of all input files to process
        buffer_size (int): How many values to keep in memory before flowing data to disk.

    Returns:
        str: Savepath of the hdf storage for easier access later on.
    """
    # generate savepaths for the buffered dicts first
    savepath_seq_hash_to_idx = Path(project).joinpath(
        "11_read_table", "temp", "seq_hash_to_idx"
    )
    savepath_seq_hash_to_seq = Path(project).joinpath(
        "11_read_table", "temp", "seq_hash_to_seq"
    )
    savepath_sample_to_idx = Path(project).joinpath(
        "11_read_table", "temp", "sample_to_idx"
    )

    # generate the buffered dicts
    seq_hash_to_idx = generate_buffered_dict(savepath_seq_hash_to_idx, buffer_size)
    seq_hash_to_seq = generate_buffered_dict(savepath_seq_hash_to_seq, buffer_size)
    sample_to_idx = generate_buffered_dict(savepath_sample_to_idx, buffer_size)

    # define all indeces prior to looping over the files
    hash_idx = 0
    sample_idx = 0
    data_idx = 0

    # spill over to hdf if too full, no lookups needed
    full_data = {}

    # define the hdf savepath, remove data from previous run to not append to it
    hdf_savename = Path(project).joinpath(
        "11_read_table",
        "data",
        "read_data_storage_{}.h5.lz".format(Path(project).stem.replace("_apscale", "")),
    )
    os.remove(hdf_savename)

    # loop over all input files
    for file in input_files:
        # loop over the data in the files line by line
        for sample_name, hash, seq, read_count in parse_fasta_data(file):
            # new hashes and samples only have to be added once a new one is found
            if hash not in seq_hash_to_idx:
                seq_hash_to_idx[hash] = hash_idx
                seq_hash_to_seq[hash] = seq
                # increase the hash idx
                hash_idx += 1
            if sample_name not in sample_to_idx:
                sample_to_idx[sample_name] = sample_idx
                sample_idx += 1
            # data has to be added at each iteration
            full_data[data_idx] = (
                seq_hash_to_idx[hash],
                sample_to_idx[sample_name],
                read_count,
            )
            data_idx += 1
            # as soon as buffer size data rows are collected, stream data to hdf
            if data_idx % buffer_size == 0:
                # add data to the hdf store
                save_data_lines_to_hdf(hdf_savename, full_data)
                full_data = {}
    else:
        save_data_lines_to_hdf(hdf_savename, full_data)

    # save the buffered dicts to the hdf store in chunks

    # return to hdf savename for reuse / ordering / reading
    return hdf_savename


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
            "{}: Please perform at least on data aggregation step before creating the read table.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
        return None

    # collect the input for read table generation
    # generate a data folder to save the hdf store
    input = glob.glob(str(Path(project).joinpath(prior_step, "data", "*.fasta.gz")))

    # generate a temp folder for saving the buffered dicts
    ## create temporal output folder
    for folder in ["temp", "data"]:
        try:
            os.mkdir(Path(project).joinpath("11_read_table", folder))
        except FileExistsError:
            pass

    # index all the data and save as hdf store
    hdf_savename = index_data_to_hdf(project, input, 100_000)

    # sort the hdf store to generate the fasta file

    # create parquet and excel outputs from the hdf

    # remove temporary savefiles
    shutil.rmtree(Path(project).joinpath("11_read_table", "temp"))

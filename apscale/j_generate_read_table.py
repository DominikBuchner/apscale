import glob, os, datetime, gzip
import pandas as pd
from joblib import Parallel, delayed
from pathlib import Path
from apscale.a_create_project import choose_input
from Bio.SeqIO.FastaIO import SimpleFastaParser


def fasta_to_parquet(fasta_path: str, output_folder: str) -> None:
    """Function to stream fasta to parquet as tables in the form of
    [header, sequence, count]

    Args:
        fasta_path (str): Path to the fasta file to transform.
        output_folder (str): Output folder to write to
    """
    # extract the sample name from the fasta_path
    fasta_name = Path(fasta_path).with_suffix("").stem

    # read the fasta with the simple fasta parser to generate rows
    all_rows = []

    with gzip.open(fasta_path, "rt") as data_stream:
        for header, seq in SimpleFastaParser(data_stream):
            header_data = header.split(";size=")
            seq_hash, seq_count = header_data[0], int(header_data[1])
            seq = seq.upper().strip()
            all_rows.append([fasta_name, seq_hash, seq, seq_count])

    # put into pandas dataframe
    fasta_data = pd.DataFrame(
        data=all_rows, columns=["sample", "hash", "sequence", "read_count"]
    )

    # add correct datatype annotations
    fasta_data = fasta_data.astype(
        {
            "sample": "string",
            "hash": "string",
            "sequence": "string",
            "read_count": "int64",
        }
    )

    # create an output path
    output_path = output_folder.joinpath(f"{fasta_name}.parquet.snappy")
    fasta_data.to_parquet(output_path)


def main(project=Path.cwd()):
    """Main function to initially create the read table data.

    Args:
        project (str, optional): Path to the apscale project to work in. Defaults to Path.cwd().
    """
    # collect the general settings
    gen_settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="0_general_settings",
    )
    cores = gen_settings["cores to use"].item()

    # collect the settings from the excel sheet
    settings = pd.read_excel(
        Path(project).joinpath(
            f"Settings_{Path(project).name.replace('_apscale', '')}.xlsx"
        ),
        sheet_name="11_read_table",
    )
    perform = settings["generate read table"].item()
    to_excel = settings["to excel"].item()
    to_parquet = settings["to parquet"].item()
    group_threshold = settings["sequence group threshold"].item()

    # find out which dataset to load
    prior_step = choose_input(project, "11_generate_read_table")

    # collect the input for read table generation
    # generate a data folder to save the hdf store
    input_files = glob.glob(
        str(Path(project).joinpath(prior_step, "data", "*.fasta.gz"))
    )
    # generate a temp folder for saving the buffered dicts
    ## create temporal output folder
    for folder in ["temp", "data"]:
        try:
            os.mkdir(Path(project).joinpath("11_read_table", folder))
        except FileExistsError:
            pass

    data_path = Path(project).joinpath("11_read_table", "data")
    temp_path = Path(project).joinpath("11_read_table", "temp")

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Building read data storage."
    )

    # transform all fasta files in input to parquet to ingest into duckdb
    Parallel(n_jobs=cores)(
        delayed(fasta_to_parquet)(fasta_file, temp_path) for fasta_file in input_files
    )

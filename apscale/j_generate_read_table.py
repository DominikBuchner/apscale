import glob, os, datetime, gzip, duckdb, hashlib
import pandas as pd
from joblib import Parallel, delayed
from pathlib import Path
from apscale.a_create_project import choose_input
from Bio.SeqIO.FastaIO import SimpleFastaParser


def fasta_to_parquet(fasta_path: str, output_folder: str, perform_hash: bool) -> None:
    """Function to stream fasta to parquet as tables in the form of
    [header, sequence, count]

    Args:
        fasta_path (str): Path to the fasta file to transform.
        output_folder (str): Output folder to write to
        perform_hash (bool): Wether to hash the sequences prior to parquet generation. Needed if no data aggregation step is used!
    """
    # extract the sample name from the fasta_path
    fasta_name = Path(fasta_path).with_suffix("").stem

    # read the fasta with the simple fasta parser to generate rows
    all_rows = []

    # add a dummy seq to each file to handle empty files greacefully
    all_rows.append([fasta_name, "dummy_hash", "dummy_seq", 0])

    with gzip.open(fasta_path, "rt") as data_stream:
        for header, seq in SimpleFastaParser(data_stream):
            header_data = header.split(";size=")
            seq_hash, seq_count = header_data[0], int(header_data[1])

            # only compute hash if needed
            if perform_hash:
                seq_hash = seq.encode("ascii")
                seq_hash = hashlib.sha3_256(seq_hash).hexdigest()

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

    # create an output path, save data to parquet to ingest into duckdb later
    output_path = output_folder.joinpath(f"{fasta_name}.parquet.snappy")
    fasta_data.to_parquet(output_path)


def build_read_store(input_folder, database_path):
    # fetch the parquet files
    parquet_path = input_folder.joinpath("*.parquet.snappy")

    # connect to the database
    read_data_store = duckdb.connect(database_path)

    # attach a temp database to the read_data_store
    temp_database = input_folder.joinpath("sequence_read_count_data.duckdb")
    read_data_store.execute(f"ATTACH '{temp_database}' AS temp_db")

    # add the readcounts from the parquet file to the table sequence_read_count_data
    read_data_store.execute(
        f"""
        CREATE TABLE temp_db.sequence_read_count_data AS
        SELECT * FROM read_parquet('{parquet_path}')
        """
    )

    # add the sample data table
    read_data_store.execute(
        """
        CREATE TABLE main.sample_data AS
        SELECT 
            row_number() OVER () AS sample_idx,
            sample,  
        FROM (
            SELECT DISTINCT sample
            FROM temp_db.sequence_read_count_data
        )
        """
    )

    # add the sequence data table
    read_data_store.execute(
        """
        CREATE TABLE main.sequence_data AS
        SELECT 
            row_number() OVER () AS hash_idx,
            hash,
            sequence,
        FROM (
            SELECT DISTINCT hash, sequence
            FROM temp_db.sequence_read_count_data
        )
        """
    )

    # create the indexed read_count_data, drop the raw read_count_data
    read_data_store.execute(
        """
        CREATE TABLE main.sequence_read_count_data AS
        SELECT
            sd.sample_idx,
            seqd.hash_idx AS sequence_idx,
            srd.read_count 
        FROM temp_db.sequence_read_count_data srd
        JOIN main.sample_data sd ON srd.sample = sd.sample
        JOIN main.sequence_data seqd ON srd.hash = seqd.hash AND srd.sequence = seqd.sequence
    """
    )

    # add the counts for the fasta data

    # give some user output

    read_data_store.close()


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
    project_name = Path(project).stem.replace("_apscale", "")
    database_path = data_path.joinpath(f"{project_name}_read_data_storage.duckdb")

    # check if there is an old db, remove it if it is
    if database_path.is_file():
        database_path.unlink()

    # check if there are old parquet files, remove them
    for file in temp_path.glob("*.parquet.snappy"):
        if file.is_file():
            file.unlink()

    # check if there are old parquet files and remove them first

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Building read data storage."
    )

    # only perform hashing if no aggregation step has been used before
    if prior_step == "06_dereplication":
        # transform all fasta files in input to parquet to ingest into duckdb
        Parallel(n_jobs=cores)(
            delayed(fasta_to_parquet)(fasta_file, temp_path, True)
            for fasta_file in input_files
        )
    else:
        # transform all fasta files in input to parquet to ingest into duckdb
        Parallel(n_jobs=cores)(
            delayed(fasta_to_parquet)(fasta_file, temp_path, False)
            for fasta_file in input_files
        )

    # build the duckdb database
    build_read_store(temp_path, database_path)

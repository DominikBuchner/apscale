import glob, os, datetime, gzip, duckdb, hashlib, more_itertools
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
            row_number() OVER () AS sequence_idx,
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
            seqd.sequence_idx,
            srd.read_count 
        FROM temp_db.sequence_read_count_data srd
        JOIN main.sample_data sd ON srd.sample = sd.sample
        JOIN main.sequence_data seqd ON srd.hash = seqd.hash AND srd.sequence = seqd.sequence
    """
    )

    # detach and remove the temp db
    read_data_store.execute("DETACH temp_db")
    temp_database.unlink()

    # add the counts and order for the fasta data
    read_data_store.execute(
        "ALTER TABLE sequence_data ADD COLUMN sequence_order BIGINT"
    )
    read_data_store.execute("ALTER TABLE sequence_data ADD COLUMN read_sum BIGINT")

    read_data_store.execute(
        """
        UPDATE sequence_data AS sd
        SET
            sequence_order = r.sequence_order,
            read_sum = r.read_sum
        FROM (
            SELECT
                row_number() OVER (ORDER BY sum(read_count) DESC) AS sequence_order, 
                sequence_idx, 
                sum(read_count) AS read_sum
            FROM sequence_read_count_data
            GROUP BY sequence_idx  
        ) r
        WHERE sd.sequence_idx = r.sequence_idx
    """
    )

    # give some user output
    nr_samples = read_data_store.execute(
        "SELECT COUNT(sample_idx) FROM sample_data"
    ).fetchone()[0]
    nr_sequences = (
        read_data_store.execute(
            "SELECT COUNT(sequence_idx) from sequence_data"
        ).fetchone()[0]
        - 1  # remove the empty seq
    )
    nr_reads = read_data_store.execute(
        "SELECT SUM(read_sum) FROM sequence_data"
    ).fetchone()[0]

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Dataset contains {nr_samples:,} samples, {nr_sequences:,} distinct sequences and {nr_reads:,} reads."
    )

    # close the read data store
    read_data_store.close()

    # remove the parquet files
    for file in Path(input_folder).glob("*.parquet.snappy"):
        if file.is_file():
            file.unlink()


def generate_fasta(
    project_name: str,
    output_folder: str,
    database_path: str,
    sequence_type: str,
) -> None:
    """Function to read the database and create a fasta file either for sequences or sequence groups

    Args:
        project_name (str): Name of the current apscale project.
        output_folder (str): Output folder to write to.
        database_path (str): Path to the database where the data is stored.
        sequence_type (str): Wether to write sequence or group
    """
    # establish the connection to the database
    read_data_store = duckdb.connect(database_path)
    read_data_store_cursor = read_data_store.execute(
        f"SELECT hash, sequence FROM {sequence_type}_data WHERE sequence != 'dummy_seq' ORDER BY {sequence_type}_order"
    )
    output_fasta = Path(output_folder).joinpath(
        f"0_{project_name}_{sequence_type}s.fasta"
    )

    with open(output_fasta, "w") as out_stream:
        while True:
            fasta_batch = read_data_store_cursor.fetchmany(100_000)
            if not fasta_batch:
                break
            lines = [f">{hash}\n{sequence}\n" for hash, sequence in fasta_batch]
            out_stream.writelines(lines)

    read_data_store.close()


def generate_read_table(
    project_name: str,
    output_folder: str,
    database_path: str,
    sequence_type: str,
    temp_folder: str,
    to_excel: bool,
) -> None:
    """Function to generate a read table from the database. Will be split into multiple tables if data is too large.

    Args:
        project_name (str): Name of the apscale project.
        output_folder (str): Output folder to write to.
        database_path (str): Path to the database to select the data from.
        sequence_type (str): Wether to write sequence or group data.
        temp_folder (str): Path to a temp folder for intermediate writes. Will be deleted at end of function.
        to_excel: Wether to write output to excel or not.
    """
    # connect the database
    read_data_store = duckdb.connect(database_path)
    # define the temporal database for writing intermediate steps
    temp_db = Path(temp_folder).joinpath("temp.duckdb")
    read_data_store.execute(f"ATTACH '{temp_db}' as temp_db")

    # define the select query. Build a table with hash_idx, sequence_idx, sequence_order, sample, read_count
    query_string = f"""
    SELECT
        seqd.sequence_idx AS sequence_idx,
        seqd.sequence_idx AS hash_idx,
        seqd.sequence_order,
        sd.sample,
        COALESCE(rcd.read_count, 0) AS read_count
    FROM main.{sequence_type}_data seqd
    CROSS JOIN main.sample_data sd
    LEFT JOIN main.{sequence_type}_read_count_data rcd
        ON seqd.sequence_idx = rcd.sequence_idx
        AND sd.sample_idx = rcd.sample_idx
    ORDER BY seqd.sequence_order ASC
    """

    # write the result to temp to perform multiple pivots if subsampling is needed
    read_data_store.execute(
        f"""CREATE TABLE temp_db.read_table_data AS
    {query_string}
    """
    )

    # extract the number of sequences
    sequence_count = read_data_store.execute(
        "SELECT COUNT(DISTINCT hash_idx) FROM temp_db.read_table_data"
    ).fetchone()[0]

    # extract the sample names as they are needed regardles of output format
    sample_names = read_data_store.execute(
        "SELECT DISTINCT sample FROM temp_db.read_table_data ORDER BY sample ASC"
    ).fetchall()
    sample_names = [row[0] for row in sample_names]

    # # if excel is used as output format compute the maximum samples per excel
    if to_excel:
        if sequence_count > 1_000_000:
            print(
                f"{datetime.datetime.now().strftime('%H:%M:%S')}: Dataset with {sequence_count} sequences is too large for excel. Only using parquet."
            )
        else:
            # compute the maximum number of samples per file
            chunksize = 10_000_000 // sequence_count

            # fetch the sample names from read store
            sample_names_iter = enumerate(
                more_itertools.chunked(sample_names, n=chunksize)
            )

            # for each chunk create a pivot table and stream it directly to excel
            for idx, chunk in sample_names_iter:
                sample_selector = ", ".join(f"'{sample}'" for sample in chunk)
                sample_columns = ", ".join(f'"{sample}"' for sample in chunk)
                # pivot the read_table data
                read_data_store.execute(
                    f"""
                    CREATE OR REPLACE TABLE temp_db.pivot AS
                        SELECT hash_idx, sequence_idx, {sample_columns} FROM (
                            SELECT * FROM temp_db.read_table_data
                            PIVOT (SUM(read_count) FOR sample IN ({sample_selector}))
                            ORDER BY sequence_order
                        )
                    """
                )

                # extract the read table with correct column names
                read_table = read_data_store.execute(
                    f"""
                    SELECT
                        sd.hash,
                        sd.sequence,
                        {sample_columns}
                    FROM main.sequence_data AS sd
                    LEFT JOIN temp_db.pivot AS pt
                        ON sd.sequence_idx = pt.sequence_idx
                    WHERE sd.sequence != 'dummy_seq'
                    ORDER BY sd.sequence_order
                    """
                ).df()

                # define the output path
                output_file_name_xlsx = Path(
                    output_folder.joinpath(
                        f"0_{project_name}_{sequence_type}_read_table_part_{idx}.xlsx"
                    )
                )

                # write the output, give user output
                read_table.to_excel(
                    output_file_name_xlsx, index=False, engine="xlsxwriter"
                )

    # always write output to parquet
    sample_selector = ", ".join(f"'{sample}'" for sample in sample_names)
    sample_columns = ", ".join(f'"{sample}"' for sample in sample_names)

    # define the output file name
    output_file_name_parquet = Path(
        output_folder.joinpath(
            f"0_{project_name}_{sequence_type}_read_table.parquet.snappy"
        )
    )

    read_data_store.execute(
        f"""
        CREATE OR REPLACE TABLE temp_db.pivot AS
            SELECT hash_idx, sequence_idx, {sample_columns} FROM (
                SELECT * FROM temp_db.read_table_data
                PIVOT (SUM(read_count) FOR sample IN ({sample_selector}))
                ORDER BY sequence_order
            )
        """
    )

    # extract the read table with correct column names
    read_table = read_data_store.execute(
        f"""
        COPY (
            SELECT
                sd.hash,
                sd.sequence,
                {sample_columns}
            FROM main.sequence_data AS sd
            LEFT JOIN temp_db.pivot AS pt
                ON sd.sequence_idx = pt.sequence_idx
            WHERE sd.sequence != 'dummy_seq'
            ORDER BY sd.sequence_order
        ) TO '{output_file_name_parquet}' (FORMAT 'parquet')
        """
    )

    # close the connection
    read_data_store.close()
    # remove the temp database
    temp_db.unlink()


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
    group_threshold = settings["sequence group threshold"].item()

    # find out which dataset to load
    prior_step = choose_input(project, "11_generate_read_table")

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Last analysis step was {prior_step}. Using data from this step."
    )

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

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Read data storage build successfully."
    )

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Creating fasta file for sequences."
    )

    # create the fasta file for the sequences
    generate_fasta(project_name, data_path, database_path, "sequence")

    print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Creating read table(s).")

    # create the read table for the sequences
    generate_read_table(
        project_name, data_path, database_path, "sequence", temp_path, to_excel
    )

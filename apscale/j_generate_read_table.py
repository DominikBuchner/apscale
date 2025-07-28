import glob, os, datetime, gzip, duckdb, hashlib, more_itertools, sys, subprocess
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

    # gather some basic stats about the read store
    # extract the number of sequences, includes the dummy_seq as it is removed after pivot
    sequence_count = read_data_store.execute(
        f"SELECT COUNT(hash) FROM main.{sequence_type}_data"
    ).fetchone()[0]

    # extract the sample names as they are needed regardles of output format
    sample_names = read_data_store.execute(
        "SELECT DISTINCT sample FROM main.sample_data ORDER BY sample ASC"
    ).fetchall()
    sample_names = [row[0] for row in sample_names]

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Building cross product of samples and sequences."
    )

    # build the cross product first
    read_data_store.execute(
        f"""
        CREATE OR REPLACE TABLE temp_db.cross_product AS
            SELECT
                seqd.sequence_idx,
                sd.sample_idx,
                sd.sample,
            FROM
                main.{sequence_type}_data AS seqd
            CROSS JOIN main.sample_data AS sd                    
        """
    )

    print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Joining in read counts.")

    read_data_store.execute(
        f"""
    CREATE OR REPLACE TABLE temp_db.read_table_data AS
        SELECT
            cp.sequence_idx,
            cp.sample_idx,
            cp.sample,
            COALESCE(rcd.read_count, 0) as read_count
        FROM temp_db.cross_product AS cp
        LEFT JOIN main.{sequence_type}_read_count_data AS rcd
            ON cp.sequence_idx = rcd.sequence_idx
            AND cp.sample_idx = rcd.sample_idx
    """
    )

    # if excel is used as output format compute the maximum samples per excel
    if sequence_count > 1_000_000:
        print(
            f"{datetime.datetime.now().strftime('%H:%M:%S')}: Dataset with {sequence_count - 1:,} sequences is too large for excel. Only using parquet."
        )
    else:
        print(
            f"{datetime.datetime.now().strftime('%H:%M:%S')}: Generating table(s) for Excel output."
        )

        # compute the maximum number of samples per file
        chunksize = 10_000_000 // sequence_count

        # fetch the sample names from read store
        sample_names_iter = enumerate(more_itertools.chunked(sample_names, n=chunksize))

        # for each chunk create a pivot table and stream it directly to excel
        for idx, chunk in sample_names_iter:
            sample_selector = ", ".join(f"'{sample}'" for sample in chunk)
            sample_columns = ", ".join(f'"{sample}"' for sample in sample_names)

            # select the subset of samples we are currently looking at
            read_data_store.execute(
                f"""
                CREATE OR REPLACE TABLE temp_db.filtered AS  
                    SELECT * FROM temp_db.read_table_data rtd
                    WHERE rtd.sample IN ({sample_selector})
                """
            )

            print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Pivoting table.")

            read_data_store.execute(
                f"""
                CREATE OR REPLACE TABLE temp_db.pivot AS
                    SELECT sequence_idx, {sample_columns}
                    FROM (
                        SELECT sequence_idx, sample, read_count
                        FROM temp_db.filtered
                    )
                    PIVOT (SUM(read_count) FOR sample IN ({sample_columns}))
                """
            )

            print(
                f"{datetime.datetime.now().strftime('%H:%M:%S')}: Adding hashes, sequences and ordering table."
            )

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

            print(
                f"{datetime.datetime.now().strftime('%H:%M:%S')}: Writing Excel output file."
            )

            # define the output path
            output_file_name_xlsx = Path(
                output_folder.joinpath(
                    f"0_{project_name}_{sequence_type}_read_table_part_{idx}.xlsx"
                )
            )

            # write the output, give user output
            read_table.to_excel(output_file_name_xlsx, index=False, engine="xlsxwriter")

    # always write output to parquet
    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Generating table for parquet output."
    )

    # define the output file name
    output_file_name_parquet = Path(
        output_folder.joinpath(
            f"0_{project_name}_{sequence_type}_read_table.parquet.snappy"
        )
    )

    # create a sql friendly sample columns variable
    sample_columns = ", ".join(f'"{sample}"' for sample in sample_names)
    sample_selector = ", ".join(f"'{sample}'" for sample in sample_names)

    print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Pivoting table.")

    # chunk the data
    max_pivot_cells = 25_000_000  # reasonable for a 32 Gb machine
    max_rows_per_chunk = max(1000, max_pivot_cells // len(sample_names))

    # collect the sequence ids
    sequence_ids = read_data_store.execute(
        "SELECT sequence_idx FROM main.sequence_data"
    ).fetchall()
    sequence_ids = [row[0] for row in sequence_ids]
    sequence_chunks = more_itertools.chunked(sequence_ids, max_rows_per_chunk)

    # write output to parquet in chunks
    for i, chunk in enumerate(sequence_chunks, start=1):
        min_idx, max_idx = min(chunk), max(chunk)

        # create a temp filename for output
        temp_filename = Path(temp_folder.joinpath(f"pivot_chunk_{i}.parquet.snappy"))

        # perform the pivot for this chunk
        read_data_store.execute(
            f"""
            COPY (
                SELECT sequence_idx, {sample_columns}
                FROM (
                    SELECT sequence_idx, sample, read_count
                    FROM temp_db.read_table_data
                    WHERE sequence_idx BETWEEN {min_idx} AND {max_idx}
                )
                PIVOT (SUM(read_count) FOR sample IN ({sample_selector}))
            )
            TO '{temp_filename}' (FORMAT PARQUET)
            """
        )

        # give user output
        print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Chunk {i} processed.")

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Adding hashes, sequences and ordering table and writing parquet output."
    )

    # create the parquet path to read the chunks from
    parquet_path = Path(temp_folder.joinpath("pivot_chunk_*.parquet.snappy"))

    # extract the read table with correct column names
    read_data_store.execute(
        f"""
        CREATE OR REPLACE TABLE temp_db.pivot AS 
            SELECT * FROM read_parquet('{parquet_path}')
        """
    )

    # remove the parquet files
    for file in temp_folder.glob("pivot_chunk_*.parquet.snappy"):
        if file.is_file():
            file.unlink()

    # build the final read_table
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


def generate_sequence_groups(
    project_name, data_path, database_path, temp_path, group_threshold
):
    # create the savename for the sequence fasta
    sequence_fasta = Path(data_path).joinpath(f"0_{project_name}_sequences.fasta")
    matchfile_path = Path(temp_path).joinpath(
        f"{project_name}_{group_threshold}_matchfile.tsv"
    )

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Computing matchfile with vsearch."
    )

    # compute the matchfile with vsearch
    subprocess.run(
        [
            "vsearch",
            "--usearch_global",
            sequence_fasta,
            "--db",
            sequence_fasta,
            "--id",
            str(group_threshold),
            "--userout",
            matchfile_path,
            "--quiet",
            "--userfields",
            "query+target+id",
            "--maxaccepts",
            str(0),
            "--self",
        ]
    )

    # ingest into duck db temp file to then join in data
    read_data_store = duckdb.connect(database_path)
    temp_db = Path(temp_path).joinpath("temp.duckdb")
    read_data_store.execute(f"ATTACH '{temp_db}' as temp_db")
    read_data_store.execute(
        f"""
        CREATE OR REPLACE TABLE temp_db.matchfile_temp AS
        SELECT * FROM read_csv('{matchfile_path}',
            delim = '\\t',
            header = False,
            columns = {{
                'query': 'STRING',
                'target': 'STRING',
                'pct_id': 'DOUBLE'
            }})
        """
    )

    read_data_store.close()
    if temp_db.is_file():
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
    perform_read_table_calc = settings["generate read table"].item()
    group_threshold = settings["sequence group threshold"].item()

    # check for valid group thereshold
    if not 0 < group_threshold < 1:
        print(
            f"{datetime.datetime.now().strftime('%H:%M:%S')}: Group threshold needs to be between 0 and 1."
        )
        sys.exit()

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

    # create the read table for the sequences
    if perform_read_table_calc:
        print(
            f"{datetime.datetime.now().strftime('%H:%M:%S')}: Creating read table(s)."
        )
        generate_read_table(
            project_name, data_path, database_path, "sequence", temp_path
        )

    # sequence grouping
    print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Grouping sequences.")

    generate_sequence_groups(
        project_name, data_path, database_path, temp_path, group_threshold
    )
    # fasta for sequence groups

    # read table for sequence groups

import glob, os, datetime, gzip, duckdb, hashlib, more_itertools, sys, subprocess
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
from pathlib import Path
from apscale.a_create_project import choose_input
from Bio.SeqIO.FastaIO import SimpleFastaParser

# chunking of the parquet pivot, ~25M cells is reasonable for a 32 Gb machine
MAX_PIVOT_CELLS = 25_000_000
MIN_ROWS_PER_CHUNK = 1000

# rows per logfile, well below the Excel limit of 1,048,575 data rows
LOG_ROWS_PER_FILE = 200_000

# Excel limits: columns per sheet (hash and sequence take two of them) and cells per file (memory)
EXCEL_MAX_COLUMNS = 16_384
EXCEL_MAX_CELLS = 10_000_000


# names that would collide with the columns of the read tables
RESERVED_COLUMN_NAMES = {"hash", "sequence", "sequence_idx"}


def quote_literal(value) -> str:
    """Quote a string (e.g. a path or sample name) as SQL string literal."""
    return "'" + str(value).replace("'", "''") + "'"


def quote_ident(name: str) -> str:
    """Quote a name as SQL identifier (e.g. a sample name used as column name)."""
    return '"' + str(name).replace('"', '""') + '"'


def check_sample_names(sample_names: list) -> None:
    """Sample names become column names in the read tables. DuckDB column names are case insensitive and
    must not collide with the fixed columns, so fail early with a readable message.

    Args:
        sample_names (list): Sample names to check.

    Raises:
        ValueError: If names differ only in case or use a reserved column name.
    """
    reserved = sorted(n for n in sample_names if str(n).lower() in RESERVED_COLUMN_NAMES)
    if reserved:
        raise ValueError(
            f"Sample names {reserved} are reserved column names ({sorted(RESERVED_COLUMN_NAMES)}, any case). Please rename the files."
        )

    seen = {}
    for name in sample_names:
        seen.setdefault(str(name).lower(), []).append(name)
    clashes = sorted(names for names in seen.values() if len(names) > 1)
    if clashes:
        raise ValueError(
            f"Sample names that differ only in case cannot be told apart in the read tables (case insensitive): {clashes}. Please rename the files."
        )


def remove_old_read_tables(
    output_folder: str, project_name: str, sequence_type: str = None
) -> None:
    """Remove read tables (Excel parts and parquet) of earlier runs, so they are not mixed up with new ones.

    Args:
        output_folder (str): Folder the read tables are written to.
        project_name (str): Name of the apscale project.
        sequence_type (str, optional): Only remove "sequence" or "group" tables. Defaults to both.
    """
    sequence_or_group = {"sequence": 0, "group": 1}
    for seq_type in [sequence_type] if sequence_type else sequence_or_group:
        base = f"{sequence_or_group[seq_type]}_{glob.escape(project_name)}_{seq_type}_read_table"
        for pattern in (f"{base}_part_*.xlsx", f"{base}.parquet.snappy"):
            for old_file in Path(output_folder).glob(pattern):
                old_file.unlink()


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
            seq_hash, seq_count = header_data[0], int(header_data[1].rstrip(";"))

            # normalize first, so that the hash does not depend on case or whitespace
            seq = seq.upper().strip()

            # only compute hash if needed
            if perform_hash:
                seq_hash = hashlib.sha3_256(seq.encode("ascii")).hexdigest()

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
    read_data_store.execute(f"ATTACH {quote_literal(temp_database)} AS temp_db")

    # add the readcounts from the parquet file to the table sequence_read_count_data
    read_data_store.execute(f"""
        CREATE TABLE temp_db.sequence_read_count_data AS
        SELECT
            CAST(sample AS VARCHAR) AS sample,
            CAST(hash AS VARCHAR) AS hash,
            CAST(sequence AS VARCHAR) AS sequence,
            read_count
        FROM read_parquet({quote_literal(parquet_path)})
        """)

    # add the sample data table
    read_data_store.execute("""
        CREATE TABLE main.sample_data AS
        SELECT 
            row_number() OVER (ORDER BY sample) AS sample_idx,
            sample,  
        FROM (
            SELECT DISTINCT sample
            FROM temp_db.sequence_read_count_data
        )
        """)

    # add the sequence data table
    read_data_store.execute("""
        CREATE TABLE main.sequence_data AS
        SELECT 
            row_number() OVER (ORDER BY hash, sequence) AS sequence_idx,
            hash,
            sequence,
        FROM (
            SELECT DISTINCT hash, sequence
            FROM temp_db.sequence_read_count_data
        )
        """)

    # create the indexed read_count_data, drop the raw read_count_data
    read_data_store.execute("""
        CREATE TABLE main.sequence_read_count_data AS
        SELECT
            sd.sample_idx,
            seqd.sequence_idx,
            srd.read_count 
        FROM temp_db.sequence_read_count_data srd
        JOIN main.sample_data sd ON srd.sample = sd.sample
        JOIN main.sequence_data seqd ON srd.hash = seqd.hash AND srd.sequence = seqd.sequence
    """)

    # detach and remove the temp db
    read_data_store.execute("DETACH temp_db")
    temp_database.unlink()

    # add the counts and order for the fasta data
    read_data_store.execute(
        "ALTER TABLE sequence_data ADD COLUMN sequence_order BIGINT"
    )
    read_data_store.execute("ALTER TABLE sequence_data ADD COLUMN read_sum BIGINT")

    read_data_store.execute("""
        UPDATE sequence_data AS sd
        SET
            sequence_order = r.sequence_order,
            read_sum = r.read_sum
        FROM (
            SELECT
                row_number() OVER (
                    ORDER BY sum(srcd.read_count) DESC, seqd.hash ASC, seqd.sequence ASC
                ) AS sequence_order, 
                srcd.sequence_idx, 
                sum(srcd.read_count)::BIGINT AS read_sum
            FROM sequence_read_count_data AS srcd
            JOIN sequence_data AS seqd ON srcd.sequence_idx = seqd.sequence_idx
            GROUP BY srcd.sequence_idx, seqd.hash, seqd.sequence
        ) r
        WHERE sd.sequence_idx = r.sequence_idx
    """)

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
        f"SELECT hash, sequence FROM {sequence_type}_data WHERE sequence != 'dummy_seq' ORDER BY sequence_order"
    )

    sequence_or_group = {"sequence": 0, "group": 1}

    output_fasta = Path(output_folder).joinpath(
        f"{sequence_or_group[sequence_type]}_{project_name}_{sequence_type}s.fasta"
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
    read_data_store.execute(f"ATTACH {quote_literal(temp_db)} as temp_db")
    sequence_or_group = {"sequence": 0, "group": 1}

    # remove read tables of earlier runs, a run with fewer parts would leave old ones behind
    remove_old_read_tables(output_folder, project_name, sequence_type)

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
    check_sample_names(sample_names)

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Building cross product of samples and sequences."
    )

    # build the cross product first
    read_data_store.execute(f"""
        CREATE OR REPLACE TABLE temp_db.cross_product AS
            SELECT
                seqd.sequence_idx,
                sd.sample_idx,
                sd.sample,
            FROM
                main.{sequence_type}_data AS seqd
            CROSS JOIN main.sample_data AS sd                    
        """)

    print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Joining in read counts.")

    read_data_store.execute(f"""
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
    """)

    # if excel is used as output format compute the maximum samples per excel
    if sequence_count > 1_000_000:
        print(
            f"{datetime.datetime.now().strftime('%H:%M:%S')}: Dataset with {sequence_count - 1:,} sequences is too large for excel. Only using parquet."
        )
    else:
        print(
            f"{datetime.datetime.now().strftime('%H:%M:%S')}: Generating table(s) for Excel output."
        )

        # compute the maximum number of samples per file, respect the column limit of Excel
        chunksize = max(
            1, min(EXCEL_MAX_CELLS // sequence_count, EXCEL_MAX_COLUMNS - 2)
        )

        # fetch the sample names from read store
        sample_names_iter = enumerate(more_itertools.chunked(sample_names, n=chunksize))

        # for each chunk create a pivot table and stream it directly to excel
        for idx, chunk in sample_names_iter:
            sample_selector = ", ".join(quote_literal(sample) for sample in chunk)
            sample_columns = ", ".join(quote_ident(sample) for sample in chunk)

            # select the subset of samples we are currently looking at
            read_data_store.execute(f"""
                CREATE OR REPLACE TABLE temp_db.filtered AS  
                    SELECT * FROM temp_db.read_table_data rtd
                    WHERE rtd.sample IN ({sample_selector})
                """)

            print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Pivoting table.")

            read_data_store.execute(f"""
                CREATE OR REPLACE TABLE temp_db.pivot AS
                    SELECT sequence_idx, {sample_columns}
                    FROM (
                        SELECT sequence_idx, sample, read_count
                        FROM temp_db.filtered
                    )
                    PIVOT (SUM(read_count)::BIGINT FOR sample IN ({sample_columns}))
                """)

            print(
                f"{datetime.datetime.now().strftime('%H:%M:%S')}: Adding hashes, sequences and ordering table."
            )

            read_table = read_data_store.execute(f"""
                SELECT
                    sd.hash,
                    sd.sequence,
                    {sample_columns}
                FROM main.{sequence_type}_data AS sd
                LEFT JOIN temp_db.pivot AS pt
                    ON sd.sequence_idx = pt.sequence_idx
                WHERE sd.sequence != 'dummy_seq'
                ORDER BY sd.sequence_order
                """).df()

            print(
                f"{datetime.datetime.now().strftime('%H:%M:%S')}: Writing Excel output file."
            )

            # define the output path
            output_file_name_xlsx = Path(
                output_folder.joinpath(
                    f"{sequence_or_group[sequence_type]}_{project_name}_{sequence_type}_read_table_part_{idx}.xlsx"
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
            f"{sequence_or_group[sequence_type]}_{project_name}_{sequence_type}_read_table.parquet.snappy"
        )
    )

    # create a sql friendly sample columns variable
    sample_columns = ", ".join(quote_ident(sample) for sample in sample_names)
    sample_selector = ", ".join(quote_literal(sample) for sample in sample_names)

    print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Pivoting table.")

    # chunk the data
    max_rows_per_chunk = max(MIN_ROWS_PER_CHUNK, MAX_PIVOT_CELLS // len(sample_names))

    # split the sequence_idx range into non-overlapping chunks of max_rows_per_chunk rows,
    # independent of the physical row order of the table
    chunk_ranges = read_data_store.execute(f"""
        SELECT min(sequence_idx), max(sequence_idx)
        FROM (
            SELECT
                sequence_idx,
                (row_number() OVER (ORDER BY sequence_idx) - 1) // {max_rows_per_chunk} AS chunk
            FROM main.{sequence_type}_data
        )
        GROUP BY chunk
        ORDER BY chunk
        """).fetchall()

    # write output to parquet in chunks
    for i, (min_idx, max_idx) in enumerate(chunk_ranges, start=1):

        # create a temp filename for output
        temp_filename = Path(temp_folder.joinpath(f"pivot_chunk_{i}.parquet.snappy"))

        # perform the pivot for this chunk
        read_data_store.execute(f"""
            COPY (
                SELECT sequence_idx, {sample_columns}
                FROM (
                    SELECT sequence_idx, sample, read_count
                    FROM temp_db.read_table_data
                    WHERE sequence_idx BETWEEN {min_idx} AND {max_idx}
                )
                PIVOT (SUM(read_count)::BIGINT FOR sample IN ({sample_selector}))
            )
            TO {quote_literal(temp_filename)} (FORMAT PARQUET)
            """)

        # give user output
        print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Chunk {i} processed.")

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Adding hashes and sequences, ordering table and writing parquet output."
    )

    # create the parquet path to read the chunks from
    parquet_path = Path(temp_folder.joinpath("pivot_chunk_*.parquet.snappy"))

    # extract the read table with correct column names
    read_data_store.execute(f"""
        CREATE OR REPLACE TABLE temp_db.pivot AS 
            SELECT * FROM read_parquet({quote_literal(parquet_path)})
        """)

    # remove the parquet files
    for file in temp_folder.glob("pivot_chunk_*.parquet.snappy"):
        if file.is_file():
            file.unlink()

    # build the final read_table
    read_table = read_data_store.execute(f"""
        COPY (
            SELECT
                sd.hash,
                sd.sequence,
                {sample_columns}
            FROM main.{sequence_type}_data AS sd
            LEFT JOIN temp_db.pivot AS pt
                ON sd.sequence_idx = pt.sequence_idx
            WHERE sd.sequence != 'dummy_seq'
            ORDER BY sd.sequence_order
        ) TO {quote_literal(output_file_name_parquet)} (FORMAT 'parquet')
        """)

    # close the connection
    read_data_store.close()
    # remove the temp database
    temp_db.unlink()


def generate_sequence_groups(
    project_name: str,
    data_path: str,
    database_path: str,
    temp_path: str,
    group_threshold: float,
):
    """Function to cluster sequences into groups with a fixed threshold.

    Args:
        project_name (str): Name of the Apscale project.
        data_path (str): Path to save the final data.
        database_path (str): Path to the duckdb database.
        temp_path (str): Path to save temporary products.
        group_threshold (float): Threshold to perform the grouping with.
    """
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
    read_data_store.execute(f"ATTACH {quote_literal(temp_db)} as temp_db")
    read_data_store.execute(f"""
        CREATE OR REPLACE TABLE temp_db.matchfile_temp AS
        SELECT * FROM read_csv({quote_literal(matchfile_path)},
            delim = '\\t',
            header = False,
            columns = {{
                'query': 'STRING',
                'target': 'STRING',
                'pct_id': 'DOUBLE'
            }})
        """)

    # add sequence_idx and order to the table and save in main table
    read_data_store.execute(f"""
    CREATE OR REPLACE TABLE temp_db.group_lookup_table AS
        SELECT 
            mt.query,
            mt.target,
            mt.pct_id,
            sd_query.sequence_idx AS query_idx,
            sd_query.sequence_order AS query_order,
            sd_target.sequence_idx AS target_idx,
            sd_target.sequence_order AS target_order
        FROM temp_db.matchfile_temp AS mt
        LEFT JOIN main.sequence_data AS sd_query
            ON mt.query = sd_query.hash
        LEFT JOIN main.sequence_data AS sd_target
            ON mt.target = sd_target.hash
    """)

    # fetch the number of rows in the sequence data
    nr_of_sequences = (
        read_data_store.execute("SELECT COUNT(*) FROM main.sequence_data").fetchone()[0]
        - 1
    )

    # fetch the ordered sequence data to loop over
    sequence_data = read_data_store.execute(f"""
        SELECT * FROM main.sequence_data AS sd
        WHERE sd.hash != 'dummy_hash'
        ORDER BY sd.sequence_order
        """)

    # flag to use the first sequence always as the first group
    first = True

    # open a second connection for the grouping operations
    grouping_connection = duckdb.connect(database_path)

    # create a temporary csv for pushing out the group mappings
    mapping_csv_path = temp_path.joinpath(f"{project_name}_group_mapping.csv")

    with open(mapping_csv_path, "w") as mapping_csv:
        # give user output
        sequence_counter = 0
        with tqdm(total=nr_of_sequences, desc="Grouping sequences") as pbar:
            while True:
                chunk = sequence_data.fetch_df_chunk()
                # stop when finished
                if chunk.empty:
                    break
                # perform the sequence grouping
                for (
                    sequence_idx,
                    hash,
                    sequence,
                    sequence_order,
                    read_sum,
                ) in chunk.itertuples(index=False):
                    # extract the matches for this sequence_idx
                    all_hash_matches = grouping_connection.execute(
                        f"SELECT * FROM temp_db.group_lookup_table WHERE query_idx = {sequence_idx}"
                    ).df()

                    # if it is the first hit, write it to the groups found table
                    if first:
                        first = False
                        # add this to the groups found table since it is by definition the first group found
                        grouping_connection.execute(f"""
                            CREATE OR REPLACE TABLE temp_db.groups_found AS
                            SELECT
                                CAST(query AS VARCHAR) AS query,
                                CAST(target AS VARCHAR) AS target,
                                pct_id,
                                query_idx,
                                query_order,
                                target_idx,
                                target_order
                            FROM all_hash_matches
                            """)

                        mapping_csv.write(f"{sequence_idx}\t{sequence_idx}\n")
                        sequence_counter += 1
                        pbar.update(1)
                        continue

                    # perform the actual grouping, select all potential matches for the current hash
                    target_matches = grouping_connection.execute(f"""
                        SELECT * FROM temp_db.groups_found sgf 
                        WHERE sgf.target_idx = {sequence_idx}
                        ORDER BY sgf.query_order LIMIT 1
                        """).df()
                    # if there is not matching group, add the sequence and all target to the groups found --> it forms a new group
                    if target_matches.empty:
                        grouping_connection.execute(f"""
                            INSERT INTO temp_db.groups_found SELECT * FROM all_hash_matches
                            """)
                        mapping_csv.write(f"{sequence_idx}\t{sequence_idx}\n")
                        sequence_counter += 1
                        pbar.update(1)
                    # if there is a matching group it can be added to this group
                    else:
                        group_idx = target_matches["query_idx"].item()
                        mapping_csv.write(f"{sequence_idx}\t{group_idx}\n")
                        sequence_counter += 1
                        pbar.update(1)

    # add the group mapping to the main db
    read_data_store.execute(f"""
        CREATE OR REPLACE TABLE main.group_mapping AS
        SELECT * FROM read_csv({quote_literal(mapping_csv_path)},
        delim = '\\t',
        header = False,
        columns = {{
            'sequence_idx': 'BIGINT',
            'group_idx': 'BIGINT'
        }})                  
        """)

    # add readcounts, merge the groups
    read_data_store.execute(f"""
        CREATE OR REPLACE TABLE temp_db.group_data_temp AS
        SELECT 
            gmt.group_idx AS sequence_idx,
            SUM(sd.read_sum)::BIGINT AS read_sum
        FROM main.group_mapping AS gmt
        LEFT JOIN main.sequence_data sd
            ON gmt.sequence_idx = sd.sequence_idx
        GROUP BY gmt.group_idx
        """)

    # add hash and sequence
    read_data_store.execute(f"""
        CREATE OR REPLACE TABLE main.group_data AS
        SELECT
            gdt.sequence_idx,
            sd.hash,
            sd.sequence,
            row_number() OVER (ORDER BY gdt.read_sum DESC, sd.hash ASC) AS sequence_order,
            gdt.read_sum
        FROM temp_db.group_data_temp AS gdt
        LEFT JOIN main.sequence_data AS sd
            ON gdt.sequence_idx = sd.sequence_idx
        """)

    nr_of_groups = read_data_store.execute(
        "SELECT COUNT(*) FROM main.group_data"
    ).fetchone()[0]

    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Found {nr_of_groups} sequence groups. Updating read data storage."
    )

    # update the sequence_idx column in sequence_read_count_data with group_idx
    # coalesce to include the dummy hash / dummy seq
    read_data_store.execute(f"""
        CREATE OR REPLACE TABLE temp_db.group_data_temp AS
        SELECT 
            srcd.sample_idx,
            COALESCE(gm.group_idx, srcd.sequence_idx) AS sequence_idx,
            srcd.read_count
        FROM main.sequence_read_count_data AS srcd
        LEFT JOIN main.group_mapping AS gm
            ON srcd.sequence_idx = gm.sequence_idx
        """)

    # groupby sample_idx (any), sequence_idx (any), sum(read_count)
    read_data_store.execute(f"""
    CREATE OR REPLACE TABLE main.group_read_count_data AS
        SELECT
            gdt.sample_idx,
            gdt.sequence_idx,
            SUM(gdt.read_count)::BIGINT AS read_count
        FROM temp_db.group_data_temp AS gdt
        GROUP BY gdt.sample_idx, gdt.sequence_idx
    """)

    read_data_store.close()
    grouping_connection.close()

    # remove files no longer needed
    if temp_db.is_file():
        temp_db.unlink()

    if matchfile_path.is_file():
        matchfile_path.unlink()

    if mapping_csv_path.is_file():
        mapping_csv_path.unlink()


def write_log(logfile_base_path: str, database_path: str):
    """Function to log the group mapping. Writes one Excel file per LOG_ROWS_PER_FILE rows.

    Args:
        logfile_base_path (str): Path to the folder to write the logfiles to.
        database_path (str): Path to the duckdb database
    """
    # remove logs of earlier runs, they would be mixed up with the new ones
    for old_log in Path(logfile_base_path).glob("Logfile_11_read_table_part_*.xlsx"):
        old_log.unlink()

    # connect to the read data store
    read_data_store = duckdb.connect(database_path)

    logfile_cursor = read_data_store.execute("""
        SELECT
            gd.hash AS group_hash,
            sd.hash AS sequence_hash,
            sd.read_sum
        FROM group_mapping AS gm
        LEFT JOIN sequence_data AS sd
            ON gm.sequence_idx = sd.sequence_idx
        LEFT JOIN group_data AS gd
            ON gm.group_idx = gd.sequence_idx
        ORDER BY gd.sequence_order, sd.read_sum DESC, sd.hash ASC
        """)
    columns = [column[0] for column in logfile_cursor.description]

    print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Writing logs.")

    file_nr = 1
    while True:
        rows = logfile_cursor.fetchmany(LOG_ROWS_PER_FILE)
        if not rows:
            break

        # generate a logfile name, one file per chunk
        logfile_name = Path(logfile_base_path).joinpath(
            f"Logfile_11_read_table_part_{file_nr}.xlsx"
        )
        pd.DataFrame(rows, columns=columns).to_excel(
            logfile_name, engine="xlsxwriter", index=False
        )
        file_nr += 1

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
    perform_read_table_calc = settings["generate read table"].item()
    group_threshold = settings["sequence group threshold"].item()

    # check for valid group thereshold
    if not 0 < group_threshold <= 1:
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
    # fail early, before any heavy computation, if sample names cannot be used as column names
    check_sample_names([Path(file).with_suffix("").stem for file in input_files])
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
    logfile_base_path = Path(project).joinpath("11_read_table")

    # check if there is an old db, remove it if it is
    if database_path.is_file():
        database_path.unlink()

    # remove read tables of earlier runs, they would be stale if read table generation is switched off
    remove_old_read_tables(data_path, project_name)

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
    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Starting sequence grouping with a threshold of {group_threshold}."
    )

    generate_sequence_groups(
        project_name, data_path, database_path, temp_path, group_threshold
    )

    # sequence grouping
    print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: Grouping finished.")

    # fasta for sequence groups
    print(
        f"{datetime.datetime.now().strftime('%H:%M:%S')}: Creating fasta file for sequence groups."
    )
    generate_fasta(project_name, data_path, database_path, "group")

    # create the read table for the groups
    if perform_read_table_calc:
        print(
            f"{datetime.datetime.now().strftime('%H:%M:%S')}: Creating read table(s)."
        )
        generate_read_table(project_name, data_path, database_path, "group", temp_path)

    # finally write a human readable log for the group mapping
    write_log(logfile_base_path, database_path)

    # remove the temporary folder
    if temp_path.is_dir():
        temp_path.rmdir()

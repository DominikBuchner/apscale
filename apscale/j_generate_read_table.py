import glob, datetime, duckdb, more_itertools, psutil, shutil, subprocess
import pandas as pd
from joblib import Parallel, delayed
from joblib.externals.loky import get_reusable_executor
from pathlib import Path
from apscale.a_create_project import choose_input
from apscale.j_fasta_to_parquet import fasta_to_parquet

# chunking of the parquet pivot, ~25M cells is reasonable for a 32 Gb machine
MAX_PIVOT_CELLS = 25_000_000
MIN_ROWS_PER_CHUNK = 1000

# rows per logfile, well below the Excel limit of 1,048,575 data rows
LOG_ROWS_PER_FILE = 200_000

# Excel limits: columns per sheet (hash and sequence take two of them) and cells per file (memory)
EXCEL_MAX_COLUMNS = 16_384
EXCEL_MAX_CELLS = 10_000_000


# additional vsearch options for the sequence grouping, e.g. for experiments
VSEARCH_EXTRA_ARGS = []

# resources, set from the general settings ("cores to use") in main. None: use all cores / the default share of the RAM
CORES = None
DUCKDB_MEMORY_LIMIT = None
DUCKDB_RAM_SHARE = 0.5  # DuckDB would use 80 % of the RAM by default, leave room for vsearch and other programs
DUCKDB_MB_PER_THREAD = 256  # working memory per thread; measured: 1 GB works with 4 threads, fails with 15 (Excel export of 639 sample columns)

# names that would collide with the columns of the read tables
RESERVED_COLUMN_NAMES = {"hash", "sequence", "sequence_idx"}


def duckdb_memory_limit_mb(connection) -> float:
    """Memory limit of a DuckDB connection in MB, as reported by DuckDB (e.g. '953.6 MiB')."""
    value, unit = connection.execute("SELECT current_setting('memory_limit')").fetchone()[0].split()
    factors = {"KiB": 2**10, "MiB": 2**20, "GiB": 2**30, "TiB": 2**40, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12, "bytes": 1}
    return float(value) * factors[unit] / 2**20


def connect_database(database_path: Path, temp_directory: Path = None):
    """Connect to a duckdb database with limited resources. Data that does not fit into the memory limit is
    spilled to the temp directory instead of using more memory.

    Args:
        database_path (Path): Path to the duckdb database.
        temp_directory (Path, optional): Folder for spilled data. Defaults to a folder next to the database.

    Returns:
        The duckdb connection.
    """
    connection = duckdb.connect(str(database_path))
    memory_limit = DUCKDB_MEMORY_LIMIT or (
        f"{int(psutil.virtual_memory().total * DUCKDB_RAM_SHARE / 1e6)}MB"
    )
    connection.execute(f"SET memory_limit = {quote_literal(memory_limit)}")

    # more threads than the memory limit can feed run out of memory instead of spilling to disk
    limit_mb = duckdb_memory_limit_mb(connection)
    threads = min(CORES or psutil.cpu_count(), max(1, int(limit_mb // DUCKDB_MB_PER_THREAD)))
    connection.execute(f"SET threads = {threads}")
    if temp_directory is not None:
        connection.execute(f"SET temp_directory = {quote_literal(temp_directory)}")
    return connection


def log(message: str) -> None:
    """Print a message with the current time."""
    print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: {message}")


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


def write_excel(connection, query: str, output_file: str) -> None:
    """Stream the result of a query into an Excel file with the DuckDB excel extension (single sheet with header).

    Args:
        connection: DuckDB connection to run the query on.
        query (str): SELECT statement producing the table.
        output_file (str): Path of the xlsx file to write.
    """
    try:
        connection.execute("LOAD excel")
    except duckdb.Error:
        try:
            connection.execute("INSTALL excel")
            connection.execute("LOAD excel")
        except duckdb.Error as error:
            raise RuntimeError(
                "The DuckDB excel extension is needed to write Excel files and could not be installed "
                f"(is there an internet connection for the first use?): {error}"
            ) from None
    connection.execute(
        f"COPY ({query}) TO {quote_literal(output_file)} (FORMAT xlsx, HEADER true)"
    )


def pivot_query(
    sequence_type: str, samples: list, min_idx: int = None, max_idx: int = None
) -> str:
    """SQL to pivot the sparse read counts to one column per sample. Sequences without reads in any of the
    samples have no row, missing counts of other sequences are NULL.

    Args:
        sequence_type (str): Wether to pivot sequence or group data.
        samples (list): Sample names, one column each.
        min_idx (int, optional): First sequence_idx to include. Defaults to all.
        max_idx (int, optional): Last sequence_idx to include. Defaults to all.

    Returns:
        str: SELECT statement returning sequence_idx and the sample columns.
    """
    selector = ", ".join(quote_literal(sample) for sample in samples)
    columns = ", ".join(quote_ident(sample) for sample in samples)
    idx_filter = (
        f"AND rcd.sequence_idx BETWEEN {min_idx} AND {max_idx}"
        if min_idx is not None
        else ""
    )
    return f"""
        SELECT sequence_idx, {columns}
        FROM (
            SELECT rcd.sequence_idx, sd.sample, rcd.read_count
            FROM main.{sequence_type}_read_count_data AS rcd
            JOIN main.sample_data AS sd
                ON rcd.sample_idx = sd.sample_idx
            WHERE sd.sample IN ({selector}) {idx_filter}
        )
        PIVOT (SUM(read_count)::BIGINT FOR sample IN ({selector}))
        """


def zero_filled_columns(samples: list, table_alias: str = "pt") -> str:
    """SELECT list with one column per sample where NULL (no reads) becomes 0."""
    return ", ".join(
        f"COALESCE({table_alias}.{quote_ident(sample)}, 0)::BIGINT AS {quote_ident(sample)}"
        for sample in samples
    )


def remove_old_read_tables(
    output_folder: Path, project_name: str, sequence_type: str = None
) -> None:
    """Remove read tables (Excel parts and parquet) of earlier runs, so they are not mixed up with new ones.

    Args:
        output_folder (Path): Folder the read tables are written to.
        project_name (str): Name of the apscale project.
        sequence_type (str, optional): Only remove "sequence" or "group" tables. Defaults to both.
    """
    sequence_or_group = {"sequence": 0, "group": 1}
    for seq_type in [sequence_type] if sequence_type else sequence_or_group:
        base = f"{sequence_or_group[seq_type]}_{glob.escape(project_name)}_{seq_type}_read_table"
        for pattern in (f"{base}_part_*.xlsx", f"{base}.parquet.snappy"):
            for old_file in Path(output_folder).glob(pattern):
                old_file.unlink()


def convert_fastas_to_parquet(
    input_files: list, output_folder: Path, perform_hash: bool, cores: int
) -> None:
    """Transform all fasta files to parquet files to ingest into duckdb, using several processes.

    The worker processes are shut down afterwards. joblib would keep them alive (with their memory) for the
    rest of the run, although they are only needed here.

    Args:
        input_files (list): Paths of the fasta.gz files.
        output_folder (Path): Folder to write the parquet files to.
        perform_hash (bool): Wether to hash the sequences, see fasta_to_parquet.
        cores (int): Number of processes to use.
    """
    try:
        Parallel(n_jobs=cores)(
            delayed(fasta_to_parquet)(fasta_file, output_folder, perform_hash)
            for fasta_file in input_files
        )
    finally:
        get_reusable_executor().shutdown(wait=True)


def build_read_store(input_folder: Path, database_path: Path) -> None:
    """Build the read data storage (sample, sequence and read count tables) from the parquet files of a folder.

    Args:
        input_folder (Path): Folder with the parquet files of the samples, they are removed afterwards.
        database_path (Path): Path of the duckdb database to create.
    """
    input_folder = Path(input_folder)
    # fetch the parquet files
    parquet_path = input_folder.joinpath("*.parquet.snappy")

    # connect to the database
    read_data_store = connect_database(database_path, input_folder)

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

    # one hash has to stand for exactly one sequence, the hash identifies the sequence from here on
    ambiguous = read_data_store.execute("""
        SELECT hash, COUNT(DISTINCT sequence) AS nr_sequences
        FROM temp_db.sequence_read_count_data
        GROUP BY hash
        HAVING COUNT(DISTINCT sequence) > 1
        LIMIT 5
        """).fetchall()
    if ambiguous:
        raise ValueError(
            f"Sequence hashes have to be unique, but these hashes belong to more than one sequence: {ambiguous}"
        )

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

    # add the sequence data table with its read sums and order (most reads first, ties by hash)
    read_data_store.execute("""
        CREATE TABLE main.sequence_data AS
        SELECT
            row_number() OVER (ORDER BY hash) AS sequence_idx,
            hash,
            any_value(sequence) AS sequence,
            row_number() OVER (ORDER BY sum(read_count) DESC, hash ASC) AS sequence_order,
            sum(read_count)::BIGINT AS read_sum
        FROM temp_db.sequence_read_count_data
        GROUP BY hash
        """)

    # create the indexed read_count_data
    read_data_store.execute("""
        CREATE TABLE main.sequence_read_count_data AS
        SELECT
            sd.sample_idx,
            seqd.sequence_idx,
            srd.read_count 
        FROM temp_db.sequence_read_count_data srd
        JOIN main.sample_data sd ON srd.sample = sd.sample
        JOIN main.sequence_data seqd ON srd.hash = seqd.hash
    """)

    # detach and remove the temp db
    read_data_store.execute("DETACH temp_db")
    temp_database.unlink()

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

    log(f"Dataset contains {nr_samples:,} samples, {nr_sequences:,} distinct sequences and {nr_reads:,} reads.")

    # close the read data store
    read_data_store.close()

    # remove the parquet files
    for file in Path(input_folder).glob("*.parquet.snappy"):
        if file.is_file():
            file.unlink()


def generate_fasta(
    project_name: str,
    output_folder: Path,
    database_path: Path,
    sequence_type: str,
) -> None:
    """Function to read the database and create a fasta file either for sequences or sequence groups

    Args:
        project_name (str): Name of the current apscale project.
        output_folder (Path): Output folder to write to.
        database_path (Path): Path to the database where the data is stored.
        sequence_type (str): Wether to write sequence or group
    """
    # establish the connection to the database
    read_data_store = connect_database(database_path)
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
    output_folder: Path,
    database_path: Path,
    sequence_type: str,
    temp_folder: Path,
) -> None:
    """Function to generate a read table from the database. Will be split into multiple tables if data is too large.

    Args:
        project_name (str): Name of the apscale project.
        output_folder (Path): Output folder to write to.
        database_path (Path): Path to the database to select the data from.
        sequence_type (str): Wether to write sequence or group data.
        temp_folder (Path): Path to a temp folder for intermediate writes.
    """
    output_folder, temp_folder = Path(output_folder), Path(temp_folder)
    # connect the database
    read_data_store = connect_database(database_path, temp_folder)
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

    # if excel is used as output format compute the maximum samples per excel
    if sequence_count > 1_000_000:
        log(f"Dataset with {sequence_count - 1:,} sequences is too large for excel. Only using parquet.")
    else:
        log(f"Generating table(s) for Excel output.")

        # compute the maximum number of samples per file, respect the column limit of Excel
        chunksize = max(
            1, min(EXCEL_MAX_CELLS // sequence_count, EXCEL_MAX_COLUMNS - 2)
        )

        # fetch the sample names from read store
        sample_names_iter = enumerate(
            more_itertools.chunked(sample_names, n=chunksize), start=1
        )

        # for each chunk create a pivot table and stream it directly to excel
        for idx, chunk in sample_names_iter:
            log(f"Pivoting table, adding hashes and sequences.")

            # pivot the sparse read counts of this chunk of samples, missing counts are zeros
            read_table_query = f"""
                SELECT
                    sd.hash,
                    sd.sequence,
                    {zero_filled_columns(chunk)}
                FROM main.{sequence_type}_data AS sd
                LEFT JOIN ({pivot_query(sequence_type, chunk)}) AS pt
                    ON sd.sequence_idx = pt.sequence_idx
                WHERE sd.sequence != 'dummy_seq'
                ORDER BY sd.sequence_order
                """

            log(f"Writing Excel output file.")

            # define the output path
            output_file_name_xlsx = Path(
                output_folder.joinpath(
                    f"{sequence_or_group[sequence_type]}_{project_name}_{sequence_type}_read_table_part_{idx}.xlsx"
                )
            )

            # let DuckDB stream the result into the Excel file
            write_excel(read_data_store, read_table_query, output_file_name_xlsx)

    # always write output to parquet
    log(f"Generating table for parquet output.")

    # define the output file name
    output_file_name_parquet = Path(
        output_folder.joinpath(
            f"{sequence_or_group[sequence_type]}_{project_name}_{sequence_type}_read_table.parquet.snappy"
        )
    )

    log(f"Pivoting table.")

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
            COPY ({pivot_query(sequence_type, sample_names, min_idx, max_idx)})
            TO {quote_literal(temp_filename)} (FORMAT PARQUET)
            """)

        # give user output
        log(f"Chunk {i} processed.")

    log(f"Adding hashes and sequences, ordering table and writing parquet output.")

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
    read_data_store.execute(f"""
        COPY (
            SELECT
                sd.hash,
                sd.sequence,
                {zero_filled_columns(sample_names)}
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
    data_path: Path,
    database_path: Path,
    temp_path: Path,
    group_threshold: float,
):
    """Function to cluster sequences into groups with a fixed threshold.

    Args:
        project_name (str): Name of the Apscale project.
        data_path (Path): Path to save the final data.
        database_path (Path): Path to the duckdb database.
        temp_path (Path): Path to save temporary products.
        group_threshold (float): Threshold to perform the grouping with.
    """
    cluster_fasta = Path(temp_path).joinpath(f"{project_name}_cluster_input.fasta")
    uc_path = Path(temp_path).joinpath(f"{project_name}_{group_threshold}_clusters.uc")

    read_data_store = connect_database(database_path, temp_path)
    temp_db = Path(temp_path).joinpath("temp.duckdb")
    read_data_store.execute(f"ATTACH {quote_literal(temp_db)} as temp_db")

    # the hash is the label of the sequence in vsearch, it has to be unique
    nr_of_sequences, nr_of_hashes = read_data_store.execute(
        "SELECT COUNT(*), COUNT(DISTINCT hash) FROM main.sequence_data WHERE sequence != 'dummy_seq'"
    ).fetchone()
    if nr_of_sequences != nr_of_hashes:
        raise ValueError(
            f"Sequence hashes are not unique ({nr_of_sequences - nr_of_hashes} duplicates), sequences cannot be grouped."
        )

    if nr_of_sequences == 0:
        # nothing to group, keep the table layout
        read_data_store.execute(
            "CREATE OR REPLACE TABLE main.group_mapping (sequence_idx BIGINT, group_idx BIGINT)"
        )
    else:
        # vsearch sorts by size (descending) and label, the same order as sequence_order. The sizes have to be the summed reads
        cursor = read_data_store.execute(
            "SELECT hash, sequence, read_sum FROM main.sequence_data WHERE sequence != 'dummy_seq' ORDER BY sequence_order"
        )
        with open(cluster_fasta, "w") as out_stream:
            while True:
                batch = cursor.fetchmany(100_000)
                if not batch:
                    break
                out_stream.writelines(
                    [f">{hash};size={read_sum}\n{sequence}\n" for hash, sequence, read_sum in batch]
                )

        log(f"Clustering sequences with vsearch.")

        # greedy, abundance based clustering (AGC): a sequence joins the most abundant centroid it is similar to
        try:
            result = subprocess.run(
                [
                    "vsearch",
                    "--cluster_size",
                    str(cluster_fasta),
                    "--id",
                    str(group_threshold),
                    "--sizeorder",
                    "--maxaccepts",
                    str(0),
                    "--minseqlength",
                    str(1),
                    "--uc",
                    str(uc_path),
                    "--quiet",
                ]
                + (["--threads", str(CORES)] if CORES else [])
                + VSEARCH_EXTRA_ARGS,
                capture_output=True,
                text=True,
            )
        except FileNotFoundError:
            raise RuntimeError(
                "vsearch was not found. Please install vsearch and make sure it is available in the PATH."
            ) from None
        if result.returncode != 0:
            raise RuntimeError(
                f"vsearch failed with exit code {result.returncode}: {result.stderr.strip()}"
            )

        # the uc file has one S (centroid) or H (member of a centroid) line per sequence, labels are hash;size=n
        label = "regexp_replace({}, ';size=[0-9]+;?$', '')"
        read_data_store.execute(f"""
            CREATE OR REPLACE TABLE main.group_mapping AS
            SELECT
                member.sequence_idx AS sequence_idx,
                centroid.sequence_idx AS group_idx
            FROM (
                SELECT
                    {label.format("c8")} AS member_hash,
                    {label.format("CASE WHEN c0 = 'S' THEN c8 ELSE c9 END")} AS centroid_hash
                FROM read_csv({quote_literal(uc_path)},
                    delim = '\t',
                    header = False,
                    quote = '',
                    columns = {{
                        'c0': 'VARCHAR', 'c1': 'VARCHAR', 'c2': 'VARCHAR', 'c3': 'VARCHAR', 'c4': 'VARCHAR',
                        'c5': 'VARCHAR', 'c6': 'VARCHAR', 'c7': 'VARCHAR', 'c8': 'VARCHAR', 'c9': 'VARCHAR'
                    }})
                WHERE c0 IN ('S', 'H')
            ) AS uc
            JOIN main.sequence_data AS member ON uc.member_hash = member.hash
            JOIN main.sequence_data AS centroid ON uc.centroid_hash = centroid.hash
            ORDER BY member.sequence_idx
            """)

        nr_mapped = read_data_store.execute(
            "SELECT COUNT(*) FROM main.group_mapping"
        ).fetchone()[0]
        if nr_mapped != nr_of_sequences:
            raise RuntimeError(
                f"vsearch assigned {nr_mapped} of {nr_of_sequences} sequences to a group."
            )

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

    log(f"Found {nr_of_groups} sequence groups. Updating read data storage.")

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

    # remove files no longer needed
    for file in (temp_db, cluster_fasta, uc_path):
        if file.is_file():
            file.unlink()


def write_log(logfile_base_path: Path, database_path: Path):
    """Function to log the group mapping. Writes one Excel file per LOG_ROWS_PER_FILE rows.

    Args:
        logfile_base_path (Path): Path to the folder to write the logfiles to.
        database_path (Path): Path to the duckdb database
    """
    # remove logs of earlier runs, they would be mixed up with the new ones
    for old_log in Path(logfile_base_path).glob("Logfile_11_read_table_part_*.xlsx"):
        old_log.unlink()

    # connect to the read data store
    read_data_store = connect_database(database_path)

    # number the log rows once, so that the files can be cut into fixed ranges
    read_data_store.execute("""
        CREATE OR REPLACE TEMP TABLE log_rows AS
        SELECT
            row_number() OVER (
                ORDER BY gd.sequence_order, sd.read_sum DESC, sd.hash ASC
            ) AS row_nr,
            gd.hash AS group_hash,
            sd.hash AS sequence_hash,
            sd.read_sum
        FROM group_mapping AS gm
        LEFT JOIN sequence_data AS sd
            ON gm.sequence_idx = sd.sequence_idx
        LEFT JOIN group_data AS gd
            ON gm.group_idx = gd.sequence_idx
        """)
    nr_of_rows = read_data_store.execute("SELECT COUNT(*) FROM log_rows").fetchone()[0]

    log(f"Writing logs.")

    for file_nr, first_row in enumerate(range(1, nr_of_rows + 1, LOG_ROWS_PER_FILE), start=1):
        # one file per range of rows
        logfile_name = Path(logfile_base_path).joinpath(
            f"Logfile_11_read_table_part_{file_nr}.xlsx"
        )
        write_excel(
            read_data_store,
            f"""
            SELECT group_hash, sequence_hash, read_sum FROM log_rows
            WHERE row_nr BETWEEN {first_row} AND {first_row + LOG_ROWS_PER_FILE - 1}
            ORDER BY row_nr
            """,
            logfile_name,
        )

    read_data_store.close()


def main(project: Path = None):
    """Main function to initially create the read table data.

    Args:
        project (Path, optional): Path to the apscale project to work in. Defaults to the current working directory.
    """
    project = Path.cwd() if project is None else Path(project)
    # collect the general settings
    gen_settings = pd.read_excel(
        Path(project).joinpath(
            "Settings_{}.xlsx".format(Path(project).name.replace("_apscale", ""))
        ),
        sheet_name="0_general_settings",
    )
    cores = int(gen_settings["cores to use"].item())

    # DuckDB and vsearch use the same number of cores as the file conversion
    global CORES
    CORES = cores

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
        raise ValueError(
            f"Group threshold needs to be between 0 and 1 (0 excluded), got {group_threshold}."
        )

    # find out which dataset to load
    prior_step = choose_input(project, "11_generate_read_table")
    if prior_step is None:
        raise ValueError(
            "No previous processing step is enabled in the settings, there is no data to build read tables from."
        )

    log(f"Last analysis step was {prior_step}. Using data from this step.")

    # collect the input for read table generation
    input_files = glob.glob(
        str(Path(project).joinpath(prior_step, "data", "*.fasta.gz"))
    )
    if not input_files:
        raise FileNotFoundError(
            f"No *.fasta.gz files found in {Path(project).joinpath(prior_step, 'data')}."
        )

    # fail early, before any heavy computation, if sample names cannot be used as column names
    check_sample_names([Path(file).with_suffix("").stem for file in input_files])
    data_path = Path(project).joinpath("11_read_table", "data")
    temp_path = Path(project).joinpath("11_read_table", "temp")

    # the temp folder belongs to this step only, start from scratch to get rid of leftovers of crashed runs
    shutil.rmtree(temp_path, ignore_errors=True)
    for folder in (temp_path, data_path):
        folder.mkdir(parents=True, exist_ok=True)

    project_name = Path(project).stem.replace("_apscale", "")
    database_path = data_path.joinpath(f"{project_name}_read_data_storage.duckdb")
    logfile_base_path = Path(project).joinpath("11_read_table")

    # check if there is an old db, remove it if it is
    if database_path.is_file():
        database_path.unlink()

    # remove read tables of earlier runs, they would be stale if read table generation is switched off
    remove_old_read_tables(data_path, project_name)

    log(f"Building read data storage.")

    # only perform hashing if no aggregation step has been used before
    convert_fastas_to_parquet(
        input_files, temp_path, prior_step == "06_dereplication", cores
    )

    # build the duckdb database
    build_read_store(temp_path, database_path)

    log(f"Read data storage build successfully.")

    log(f"Creating fasta file for sequences.")

    # create the fasta file for the sequences
    generate_fasta(project_name, data_path, database_path, "sequence")

    # create the read table for the sequences
    if perform_read_table_calc:
        log(f"Creating read table(s).")
        generate_read_table(
            project_name, data_path, database_path, "sequence", temp_path
        )

    # sequence grouping
    log(f"Starting sequence grouping with a threshold of {group_threshold}.")

    generate_sequence_groups(
        project_name, data_path, database_path, temp_path, group_threshold
    )

    # sequence grouping
    log(f"Grouping finished.")

    # fasta for sequence groups
    log(f"Creating fasta file for sequence groups.")
    generate_fasta(project_name, data_path, database_path, "group")

    # create the read table for the groups
    if perform_read_table_calc:
        log(f"Creating read table(s).")
        generate_read_table(project_name, data_path, database_path, "group", temp_path)

    # finally write a human readable log for the group mapping
    write_log(logfile_base_path, database_path)

    # remove the temporary folder
    shutil.rmtree(temp_path, ignore_errors=True)

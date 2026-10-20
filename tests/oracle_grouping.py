"""Test oracle: the grouping of apscale before it was replaced by `vsearch --cluster_size`.

All-vs-all search with vsearch (`--usearch_global`), followed by a sequence-by-sequence assignment in Python.
It is run with `--maxrejects 0` (exhaustive search), which is the reference result the new implementation must
reproduce. Slow (one DuckDB query per sequence), only meant for small test datasets.
"""
import datetime, subprocess
from pathlib import Path

import duckdb
from tqdm import tqdm

from apscale.j_generate_read_table import quote_literal


def oracle_generate_sequence_groups(
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
    try:
        result = subprocess.run(
            [
                "vsearch",
                "--usearch_global",
                str(sequence_fasta),
                "--db",
                str(sequence_fasta),
                "--id",
                str(group_threshold),
                "--userout",
                str(matchfile_path),
                "--quiet",
                "--userfields",
                "query+target+id",
                "--maxaccepts",
                str(0),
                "--self",
            ]
            + ["--maxrejects", "0"],  # exhaustive search: the reference result
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

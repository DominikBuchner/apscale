import streamlit as st
import duckdb
import pandas as pd
import datetime
from pathlib import Path
import more_itertools


def param_to_literal(val):
    if isinstance(val, str):
        return "'" + val.replace("'", "''") + "'"
    elif isinstance(val, (int, float)):
        return str(val)
    elif isinstance(val, (datetime.date, datetime.datetime)):
        return "'" + val.strftime("%Y-%m-%d %H:%M:%S") + "'"
    elif isinstance(val, bool):
        return "TRUE" if val else "FALSE"
    else:
        raise ValueError(f"Unsupported type for literal conversion: {type(val)}")


def expand_sql_clause(clause, params):
    # This will replace each '?' in the clause by the literal param in order.
    parts = clause.split("?")
    if len(parts) - 1 != len(params):
        raise ValueError("Number of parameters does not match number of placeholders")

    expanded = []
    for part, param in zip(parts, params + [None]):
        expanded.append(part)
        if param is not None:
            expanded.append(param_to_literal(param))
    return "".join(expanded)


def metadata_available(read_data_to_modify: str) -> tuple:
    # connect to read data store
    read_data_to_modify_con = duckdb.connect(read_data_to_modify, read_only=True)

    # extract the tables
    tables = read_data_to_modify_con.execute("SHOW TABLES").fetchall()
    tables = [table[0] for table in tables]

    read_data_to_modify_con.close()

    sequence_metadata_available, sample_metadata_available = False, False

    if "sequence_metadata" in tables:
        sequence_metadata_available = True
    if "sample_metadata" in tables:
        sample_metadata_available = True

    return sequence_metadata_available, sample_metadata_available


def collect_unique_values(read_data_to_modify, column, table):
    read_data_to_modify_con = duckdb.connect(read_data_to_modify, read_only=True)

    # get the unique column values from the selected column
    unique_values = (
        read_data_to_modify_con.execute(
            f"""
        SELECT DISTINCT("{column}")
        FROM {table}
        WHERE {table}."{column}" IS NOT NULL
        ORDER BY {table}."{column}" ASC
        """
        )
        .df()[column]
        .to_list()
    )

    return unique_values


def generate_dynamic_selection(read_data_to_modify, table, type):
    # extract the column names from duckdb database
    read_data_to_modify_con = duckdb.connect(read_data_to_modify, read_only=True)

    # extract table info
    table_info = read_data_to_modify_con.execute(f"PRAGMA table_info({table})").df()

    read_data_to_modify_con.close()

    # store the columns and datatypes here
    columns_to_type = dict(zip(table_info["name"], table_info["type"]))

    # let the user select columns to filter
    columns_to_filter = st.multiselect(
        "Apply filter to:", columns_to_type.keys(), key=type
    )

    # show a dropdown for each selected filter
    filters = {}

    for column in columns_to_filter:
        values = collect_unique_values(read_data_to_modify, column, table)
        dtype = columns_to_type[column]

        # all categorial values
        if dtype == "VARCHAR" or dtype == "BOOLEAN":
            selected_values = st.multiselect(
                f"Keep {type} where {column} matches:",
                options=values,
                default=values[0],
            )
            if selected_values:
                filters[column] = selected_values

        # all numeric values
        if dtype == "DOUBLE" or dtype == "BIGINT":
            min_value = st.number_input(
                label=f"Keep {type} of {column} which are greater or equal than:",
                min_value=min(values),
                max_value=max(values),
                value=min(values),
            )

            max_value = st.number_input(
                label=f"Keep {type} of {column} which are smaller or equal than:",
                min_value=min(values),
                max_value=max(values),
                value=max(values),
            )
            if min_value and max_value:
                filters[column] = (min_value, max_value)

        if dtype == "TIMESTAMP_NS":
            date_range = st.date_input(
                f"Keep {type} in selected date range:",
                value=(values[0], values[1]),
                min_value=values[0],
                max_value=values[-1],
            )

            if len(date_range) == 2:
                start_date, end_date = pd.to_datetime(date_range)
                filters[column] = (start_date, end_date)
    # create an SQL filter string
    clauses = []
    params = []

    for column, value in filters.items():
        dtype = columns_to_type[column]

        # categoric selections
        if dtype in ("VARCHAR", "BOOLEAN"):
            placeholders = ", ".join(["?"] * len(value))
            clauses.append(f'"{column}" IN ({placeholders})')
            params.extend(value)

        elif dtype in ("DOUBLE", "BIGINT", "TIMESTAMP_NS"):
            clauses.append(f'"{column}" BETWEEN ? AND ?')
            params.extend(value)

    where_clause = " AND ".join(clauses) if clauses else "1=1"

    return where_clause, params, filters


def generate_sample_data_splits(
    read_data_to_modify: str, sample_where_clause, sample_params
) -> list:
    # connect to the database
    read_data_to_modify_con = duckdb.connect(read_data_to_modify, read_only=True)

    # create a selector to perform the filtering
    selector = expand_sql_clause(sample_where_clause, sample_params)

    # get the table info
    info = read_data_to_modify_con.execute("PRAGMA table_info(sample_metadata)").df()

    # only get categorial values
    info = info.loc[
        (info["type"].str.upper() == "VARCHAR")
        | (info["type"].str.upper() == "BOOLEAN")
    ]

    selectbox_options = info["name"].to_list()

    # create a multi select
    selection = st.multiselect(
        label="Split the output table on these columns:", options=selectbox_options
    )

    if selection:
        # create a column selector
        column_selector = ", ".join(f'"{col}"' for col in selection)

        # filter the table to get unique value combinations
        sample_data_splits = read_data_to_modify_con.execute(
            f"""
            SELECT DISTINCT {column_selector} 
            FROM sample_metadata
            WHERE {selector}
            """
        )

        sample_data_splits = sample_data_splits.df().to_dict(orient="records")

        # close the connection
        read_data_to_modify_con.close()

        return sample_data_splits
    else:
        return []


def get_export_columns(read_data_to_modify, seq_metadata_available):
    # connect to database
    read_data_to_modify_con = duckdb.connect(read_data_to_modify, read_only=True)

    if seq_metadata_available:
        info = read_data_to_modify_con.execute(
            "PRAGMA table_info(sequence_metadata)"
        ).df()
    else:
        info = read_data_to_modify_con.execute("PRAGMA table_info(sequence_data)").df()

    read_data_to_modify_con.close()

    return info["name"].to_list()


def prepare_export_data(
    read_data_to_modify: str,
    sequence_meta_available: bool,
    sample_meta_available: bool,
    seq_where_clause: str,
    seq_params: list,
    sample_where_clause: str,
    sample_params: list,
) -> None:

    # create a temp folder for exports
    temp_folder = read_data_to_modify.parent.parent.joinpath("temp")
    temp_folder.mkdir(exist_ok=True)

    # create a connection to the database and attach the temp connection
    read_data_to_modify_con = duckdb.connect(read_data_to_modify)
    temp_db = Path(temp_folder).joinpath("temp.duckdb")
    read_data_to_modify_con.execute(f"ATTACH '{temp_db}' as temp_db")

    if sequence_meta_available:
        # filter sequence metadata
        where_filter = expand_sql_clause(seq_where_clause, seq_params)

        read_data_to_modify_con.execute(
            f"""
        CREATE OR REPLACE TABLE temp_db.seq_data AS
        (
            SELECT * FROM main.sequence_metadata
            WHERE {where_filter}
        )
        """
        )
        # filter group metadata
        read_data_to_modify_con.execute(
            f"""
        CREATE OR REPLACE TABLE temp_db.group_data AS
        (
            SELECT * FROM main.group_metadata
            WHERE {where_filter}
        )
        """
        )
    else:
        read_data_to_modify_con.execute(
            f"""
        CREATE OR REPLACE TABLE temp_db.seq_data AS
            SELECT * FROM main.sequence_data
        """
        )
        # filter group metadata
        read_data_to_modify_con.execute(
            f"""
        CREATE OR REPLACE TABLE temp_db.group_data AS
            SELECT * FROM main.group_data
        """
        )
    if sample_meta_available:
        # filter sequence metadata
        where_filter = expand_sql_clause(sample_where_clause, sample_params)

        read_data_to_modify_con.execute(
            f"""
        CREATE OR REPLACE TABLE temp_db.sample_data AS
        (
            SELECT * FROM main.sample_metadata
            WHERE {where_filter}
        )
        """
        )
    else:
        read_data_to_modify_con.execute(
            f"""
        CREATE OR REPLACE TABLE temp_db.sample_data AS
            SELECT * FROM main.sample_data
        """
        )

    # close the connection in the end
    read_data_to_modify_con.close()

    return temp_folder, temp_db


def build_where_clause(filter: dict):
    if filter == {"1=1": "1=1"}:
        return "1=1", []

    conditions = []
    params = []

    for col, val in filter.items():
        conditions.append(f"{col} = ?")
        params.append(val)

    where_clause = " AND ".join(conditions)
    return where_clause, params


def build_output_folder_name(filter: dict):
    if filter == {"1=1": "1=1"}:
        return "read_table"
    output_string = [f"{key}-{value}" for key, value in filter.items()]
    output_string = "_".join(output_string)
    return output_string


def export_read_tables(
    read_data_to_modify: str,
    temp_folder: str,
    temp_db: str,
    sample_data_splits: list,
    display_columns: list,
    export_path: str,
):
    # connect to the database and attach the temp_db
    read_data_to_modify_con = duckdb.connect(read_data_to_modify)
    read_data_to_modify_con.execute(f"ATTACH '{temp_db}' as temp_db")

    # if the data splits are empty create a dummy selection that is always true
    if not sample_data_splits:
        sample_data_splits = [{"1=1": "1=1"}]

    # loop over the data splits and dynamically create an sql selection for each subtable
    for selection in sample_data_splits:
        where_clause, params = build_where_clause(selection)
        where_string = expand_sql_clause(where_clause, params)

        # create a temporary view for that dataset
        read_data_to_modify_con.execute(
            f"""
            CREATE OR REPLACE VIEW temp_db.sample_data_selection AS
            (
            SELECT * FROM temp_db.sample_data WHERE {where_string}
            )
            """
        )

        # generate an output folder name, create the folder
        output_folder = build_output_folder_name(selection)
        output_folder = export_path.joinpath(output_folder)
        output_folder.mkdir(exist_ok=True)

        # extract the sample names
        sample_names = read_data_to_modify_con.execute(
            "SELECT DISTINCT sample FROM temp_db.sample_data_selection ORDER BY sample ASC"
        ).fetchall()
        sample_names = [row[0] for row in sample_names]

        # build the cross product of samples for sequences and groups
        # SEQUENCES
        read_data_to_modify_con.execute(
            f"""
            CREATE OR REPLACE TABLE temp_db.seq_cross_product AS
                SELECT
                    seqd.sequence_idx,
                    sds.sample_idx,
                    sds.sample
                FROM
                    temp_db.seq_data seqd
                CROSS JOIN temp_db.sample_data_selection AS sds
            """
        )

        # GROUPS
        read_data_to_modify_con.execute(
            f"""
            CREATE OR REPLACE TABLE temp_db.group_cross_product AS
                SELECT
                    gd.sequence_idx,
                    sds.sample_idx,
                    sds.sample
                FROM
                    temp_db.group_data gd
                CROSS JOIN temp_db.sample_data_selection AS sds
            """
        )

        # join in the read counts
        # SEQUENCES
        read_data_to_modify_con.execute(
            f"""
            CREATE OR REPLACE TABLE temp_db.seq_read_table_data AS
            SELECT 
                scp.sequence_idx,
                scp.sample_idx,
                scp.sample,
                COALESCE(rcd.read_count, 0) AS read_count
            FROM temp_db.seq_cross_product AS scp
            LEFT JOIN main.sequence_read_count_data AS rcd
                ON scp.sequence_idx = rcd.sequence_idx
                AND scp.sample_idx = rcd.sample_idx
            """
        )

        # GROUPS
        read_data_to_modify_con.execute(
            f"""
            CREATE OR REPLACE TABLE temp_db.group_read_table_data AS
            SELECT 
                gcp.sequence_idx,
                gcp.sample_idx,
                gcp.sample,
                COALESCE(rcd.read_count, 0) AS read_count
            FROM temp_db.group_cross_product AS gcp
            LEFT JOIN main.sequence_read_count_data AS rcd
                ON gcp.sequence_idx = rcd.sequence_idx
                AND gcp.sample_idx = rcd.sample_idx
            """
        )

        # write the outputs
        output_names = {}

        # define output files
        for type in ["sequence", "group"]:
            for file_ext in ["parquet.snappy", "csv"]:
                output_names[(type, file_ext)] = f"{type}s_read_table.{file_ext}"

        # create a sql friendly sample columns variable
        sample_columns = ", ".join(f'"{sample}"' for sample in sample_names)
        sample_selector = ", ".join(f"'{sample}'" for sample in sample_names)

        # chunk the data
        max_pivot_cells = 25_000_000  # reasonable for a 32 Gb machine
        max_rows_per_chunk = max(1000, max_pivot_cells // len(sample_names))

        # collect the sequence_ids
        sequence_ids = read_data_to_modify_con.execute(
            "SELECT sequence_idx FROM temp_db.seq_data"
        ).fetchall()
        sequence_ids = [row[0] for row in sequence_ids]
        sequence_chunks = more_itertools.chunked(sequence_ids, max_rows_per_chunk)

        # write the output to parquet in chunks
        # SEQUENCES
        for i, chunk in enumerate(sequence_chunks, start=1):
            min_idx, max_idx = min(chunk), max(chunk)

            # create a temp filename for output
            temp_filename = Path(
                output_folder.joinpath(f"pivot_chunk_{i}.parquet.snappy")
            )

            # perform the pivot for this chunk
            read_data_to_modify_con.execute(
                f"""
            COPY (
                SELECT sequence_idx, {sample_columns}
                FROM (
                    SELECT sequence_idx, sample, read_count
                    FROM temp_db.seq_read_table_data
                    WHERE sequence_idx BETWEEN {min_idx} AND {max_idx}
                )
                PIVOT (SUM(read_count) FOR sample IN ({sample_selector}))
            )
            TO '{temp_filename}' (FORMAT PARQUET)
            """
            )

        # create the path to read the pivoted parquet from
        input_parquet_path = Path(output_folder).joinpath(
            "pivot_chunk_*.parquet.snappy"
        )

        # ingest the chunks into the database
        read_data_to_modify_con.execute(
            f"""
            CREATE OR REPLACE TABLE temp_db.pivot AS
                SELECT * FROM read_parquet('{input_parquet_path}')
            """
        )

        # remove the parquet chunks
        for file in output_folder.glob("pivot_chunk_*.parquet.snappy"):
            if file.is_file():
                file.unlink()

        sample_columns = ", ".join(
            [f'sd."{col}"' for col in display_columns]
            + [f'pt."{col}"' for col in sample_names]
        )

        output_name = output_folder.joinpath(
            output_names[("sequence", "parquet.snappy")]
        )

        # build the final read_table to parquet
        read_data_to_modify_con.execute(
            f"""
            COPY (
                SELECT
                    {sample_columns}
                FROM temp_db.seq_data AS sd
                LEFT JOIN temp_db.pivot AS pt
                    ON sd.sequence_idx = pt.sequence_idx
                WHERE sd.sequence != 'dummy_seq'
                ORDER BY sd.sequence_order
            ) TO '{output_name}' (FORMAT 'parquet')
            """
        )

        # build the final read_table to csv
        output_name = output_folder.joinpath(output_names[("sequence", "csv")])

        read_data_to_modify_con.execute(
            f"""
            COPY (
                SELECT
                    {sample_columns}
                FROM temp_db.seq_data AS sd
                LEFT JOIN temp_db.pivot AS pt
                    ON sd.sequence_idx = pt.sequence_idx
                WHERE sd.sequence != 'dummy_seq'
                ORDER BY sd.sequence_order
            ) TO '{output_name}' (FORMAT 'csv', DELIMITER '\t')
            """
        )

        # collect the sequence_ids
        sequence_ids = read_data_to_modify_con.execute(
            "SELECT sequence_idx FROM temp_db.group_data"
        ).fetchall()
        sequence_ids = [row[0] for row in sequence_ids]
        sequence_chunks = more_itertools.chunked(sequence_ids, max_rows_per_chunk)

        # create a sql friendly sample columns variable
        sample_columns = ", ".join(f'"{sample}"' for sample in sample_names)
        sample_selector = ", ".join(f"'{sample}'" for sample in sample_names)

        # write the output to parquet in chunks
        # GROUPS
        for i, chunk in enumerate(sequence_chunks, start=1):
            min_idx, max_idx = min(chunk), max(chunk)

            # create a temp filename for output
            temp_filename = Path(
                output_folder.joinpath(f"pivot_chunk_{i}.parquet.snappy")
            )

            # perform the pivot for this chunk
            read_data_to_modify_con.execute(
                f"""
            COPY (
                SELECT sequence_idx, {sample_columns}
                FROM (
                    SELECT sequence_idx, sample, read_count
                    FROM temp_db.group_read_table_data
                    WHERE sequence_idx BETWEEN {min_idx} AND {max_idx}
                )
                PIVOT (SUM(read_count) FOR sample IN ({sample_selector}))
            )
            TO '{temp_filename}' (FORMAT PARQUET)
            """
            )

        # create the path to read the pivoted parquet from
        input_parquet_path = Path(output_folder).joinpath(
            "pivot_chunk_*.parquet.snappy"
        )

        # ingest the chunks into the database
        read_data_to_modify_con.execute(
            f"""
            CREATE OR REPLACE TABLE temp_db.pivot AS
                SELECT * FROM read_parquet('{input_parquet_path}')
            """
        )

        # remove the parquet chunks
        for file in output_folder.glob("pivot_chunk_*.parquet.snappy"):
            if file.is_file():
                file.unlink()

        sample_columns = ", ".join(
            [f'gd."{col}"' for col in display_columns]
            + [f'pt."{col}"' for col in sample_names]
        )

        output_name = output_folder.joinpath(output_names[("group", "parquet.snappy")])

        # build the final read_table to parquet
        read_data_to_modify_con.execute(
            f"""
            COPY (
                SELECT
                    {sample_columns}
                FROM temp_db.group_data AS gd
                LEFT JOIN temp_db.pivot AS pt
                    ON gd.sequence_idx = pt.sequence_idx
                WHERE gd.sequence != 'dummy_seq'
                ORDER BY gd.sequence_order
            ) TO '{output_name}' (FORMAT 'parquet')
            """
        )

        # build the final read_table to csv
        output_name = output_folder.joinpath(output_names[("group", "csv")])

        read_data_to_modify_con.execute(
            f"""
            COPY (
                SELECT
                    {sample_columns}
                FROM temp_db.group_data AS gd
                LEFT JOIN temp_db.pivot AS pt
                    ON gd.sequence_idx = pt.sequence_idx
                WHERE gd.sequence != 'dummy_seq'
                ORDER BY gd.sequence_order            
                ) TO '{output_name}' (FORMAT 'csv', DELIMITER '\t')
            """
        )

    # close the connection
    read_data_to_modify_con.close()


def filters_to_log_df(seq_filters: dict, samp_filters: dict) -> pd.DataFrame:
    records = []

    for col, values in seq_filters.items():
        for val in values:
            records.append({"filter_type": "sequence", "column": col, "value": val})

    for col, values in samp_filters.items():
        for val in values:
            records.append({"filter_type": "sample", "column": col, "value": val})

    return pd.DataFrame(records)


def main():
    # prevent page from scroling up on click
    st.markdown(
        """
    <style>
        * {
        overflow-anchor: none !important;
        }
    </style>""",
        unsafe_allow_html=True,
    )

    # add helptext
    help = """This module can be used to export read tables in either parquet or excel format.
    The read tables can be seperated by metadata fields and will contain the GBIF taxonomy and GBIF
    validation column if those steps where performed. The metadata to be displayed and exported can
    also be selected."""

    # define the title
    st.title("Export read tables", help=help)

    # check if sample and sequence metadata are available
    seq_meta_available, sample_meta_available = metadata_available(
        st.session_state["read_data_to_modify"]
    )

    # initialize sample and sequence filters
    sample_filters, sequence_filters, sample_data_splits = {}, {}, []

    # sequence filtering
    st.header("Select filters for sequence metadata")
    sequence_where_clause, sequence_params = "", []
    if seq_meta_available:
        sequence_where_clause, sequence_params, sequence_filters = (
            generate_dynamic_selection(
                st.session_state["read_data_to_modify"],
                "sequence_metadata",
                "sequences",
            )
        )
    else:
        st.write("Sequence metadata is needed for this module.")
    st.divider()

    st.header("Select filters for sample metadata")
    sample_where_clause, sample_params = "", []
    if sample_meta_available:
        sample_where_clause, sample_params, sample_filters = generate_dynamic_selection(
            st.session_state["read_data_to_modify"], "sample_metadata", "samples"
        )
    else:
        st.write("Sample metadata is needed for this module.")
    st.divider()

    st.header("Split output tables by sample metadata")
    if sample_meta_available:
        sample_data_splits = generate_sample_data_splits(
            st.session_state["read_data_to_modify"], sample_where_clause, sample_params
        )
    else:
        st.write("Sample metadata is needed for this module.")

    st.divider()
    st.header("Export read tables")

    # get the columns to display in the multiselect
    display_columns = get_export_columns(
        st.session_state["read_data_to_modify"], seq_meta_available
    )

    cols_to_include = st.multiselect(
        label="Select columns to include in exported table",
        options=display_columns,
        default=display_columns,
    )

    if len(cols_to_include) == 0:
        export_tables = st.button(
            label="Export read tables", type="primary", disabled=True
        )
    else:
        export_tables = st.button(
            label="Export read tables", type="primary", disabled=False
        )

    if export_tables:
        # create an export folder with a timestamp
        data_folder = st.session_state["read_data_to_modify"].parent
        export_folder = (
            f"export_{datetime.datetime.now().strftime('%d-%m-%Y_%H-%M-%S')}"
        )
        export_path = Path(data_folder.joinpath(export_folder))
        export_path.mkdir()

        # write the filtering log to that export path
        log = filters_to_log_df(sequence_filters, sample_filters)
        log.to_excel(export_path.joinpath("filter_log.xlsx"), index=False)

        with st.spinner("Building read tables. Hang on!", show_time=True):
            # prepare temp data in first step
            temp_folder, temp_db = prepare_export_data(
                st.session_state["read_data_to_modify"],
                seq_meta_available,
                sample_meta_available,
                sequence_where_clause,
                sequence_params,
                sample_where_clause,
                sample_params,
            )

            # export read tables then
            export_read_tables(
                st.session_state["read_data_to_modify"],
                temp_folder,
                temp_db,
                sample_data_splits,
                cols_to_include,
                export_path,
            )

            st.toast("Read tables created!")

        # remove the temp db and temp folder
        temp_db.unlink()
        temp_folder.rmdir()


main()

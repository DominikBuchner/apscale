import streamlit as st
import duckdb
import pandas as pd
import datetime


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

    return where_clause, params


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

        return sample_data_splits

    else:
        return []


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

    # sequence filtering
    st.header("Select filters for sequence metadata")
    if seq_meta_available:
        sequence_where_clause, sequence_params = generate_dynamic_selection(
            st.session_state["read_data_to_modify"], "sequence_metadata", "sequences"
        )
    else:
        st.write("Sequence metadata is needed for this module.")
    st.divider()

    st.header("Select filters for sample metadata")
    if sample_meta_available:
        sample_where_clause, sample_params = generate_dynamic_selection(
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


main()

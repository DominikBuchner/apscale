"""Sample names (taken from file names) end up in SQL: odd characters must work, ambiguous names must fail early."""
import gzip

import duckdb
import openpyxl
import pytest

from apscale import j_generate_read_table as j

ODD_NAMES = ["O'Brien", 'say "hi"', "with space", "dot.name", "semi;colon", "mix'\"ed", "ünï", "[brackets]*"]


def _excel_columns(data):
    columns = []
    for part in sorted(data.glob("0_t_sequence_read_table_part_*.xlsx")):
        wb = openpyxl.load_workbook(part, read_only=True)
        header = next(wb.active.iter_rows(min_row=1, max_row=1, values_only=True))
        columns += [c for c in header if c not in ("hash", "sequence")]
        wb.close()
    return columns


def test_odd_sample_names_in_read_tables(tmp_path, make_storage):
    # the project folder contains an apostrophe as well, it ends up in SQL file paths
    db, data, temp = make_storage(tmp_path / "it's a project", ODD_NAMES)
    j.generate_read_table("t", data, db, "sequence", temp)

    assert sorted(_excel_columns(data)) == sorted(ODD_NAMES)
    parquet = next(data.glob("0_t_sequence_read_table.parquet.snappy"))
    parquet_columns = duckdb.sql(f"SELECT * FROM read_parquet({j.quote_literal(parquet)}) LIMIT 0").columns
    assert sorted(parquet_columns) == sorted(["hash", "sequence"] + ODD_NAMES)
    # counts end up in the right column
    row = duckdb.sql(f"""SELECT "O'Brien", "ünï" FROM read_parquet({j.quote_literal(parquet)}) WHERE hash = 'hash0'""").fetchone()
    assert row[0] == 1 + (0 + 0) % 5 and row[1] == 1 + (ODD_NAMES.index("ünï") + 0) % 5


def test_odd_file_names_through_fasta_to_parquet(tmp_path):
    folder = tmp_path / "it's here"
    folder.mkdir()
    fasta = folder / "O'Brien.fasta.gz"
    with gzip.open(fasta, "wt") as f:
        f.write(">h1;size=4\nACGT\n")
    j.fasta_to_parquet(str(fasta), folder, False)
    db = folder / "db.duckdb"
    j.build_read_store(folder, db)
    assert duckdb.connect(str(db)).execute("SELECT sample FROM sample_data").fetchall() == [("O'Brien",)]


@pytest.mark.parametrize(
    "names, message",
    [
        (["Sample1", "sample1"], "case"),
        (["hash", "other"], "reserved"),
        (["Sequence", "other"], "reserved"),
        (["sequence_idx", "other"], "reserved"),
    ],
)
def test_ambiguous_sample_names_fail_early(names, message):
    with pytest.raises(ValueError, match=message):
        j.check_sample_names(names)


def test_ambiguous_names_stop_read_table_generation(tmp_path, make_storage):
    db, data, temp = make_storage(tmp_path, ["Sample1", "sample1"])
    with pytest.raises(ValueError, match="case"):
        j.generate_read_table("t", data, db, "sequence", temp)


def test_hash_with_two_sequences_is_rejected(tmp_path):
    import pandas as pd

    temp = tmp_path / "temp"
    temp.mkdir()
    pd.DataFrame(
        [("s1", "h1", "ACGT", 3), ("s2", "h1", "ACGA", 4), ("s1", "dummy_hash", "dummy_seq", 0), ("s2", "dummy_hash", "dummy_seq", 0)],
        columns=["sample", "hash", "sequence", "read_count"],
    ).to_parquet(temp / "all.parquet.snappy")
    with pytest.raises(ValueError, match="unique"):
        j.build_read_store(temp, tmp_path / "db.duckdb")

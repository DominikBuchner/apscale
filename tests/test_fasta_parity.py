"""The lightweight pyarrow converter must write the same tables as the former pandas based one."""
import gzip
from pathlib import Path

import duckdb
import pytest

from apscale.j_fasta_to_parquet import fasta_to_parquet

from oracle_fasta_to_parquet import oracle_fasta_to_parquet

SMALL = Path("D:/apscale_test_data/small_dataset")

CASES = {
    "plain": ">a;size=3\nACGT\n>b;size=5;\nGGCC\n",
    "multiline": ">a;size=3\nACGT\nTTGA\nCC\n>b;size=1\nAAAA\n",
    "lowercase_and_spaces": ">a;size=2\nacgt \n>b;size=2\n  ACGT\n",
    "blank_lines": ">a;size=2\nACGT\n\n>b;size=4\nTTTT\n\n",
    "windows_line_endings": ">a;size=2\r\nACGT\r\n>b;size=4\r\nTTTT\r\n",
    "empty": "",
}


def _table(parquet):
    con = duckdb.connect()
    cols = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{parquet}')").fetchall()
    rows = con.execute(f"SELECT * FROM read_parquet('{parquet}') ORDER BY hash, sequence, read_count").fetchall()
    return [(c[0], c[1]) for c in cols], rows


@pytest.mark.parametrize("perform_hash", [False, True])
@pytest.mark.parametrize("name", list(CASES))
def test_same_table_as_pandas_version(tmp_path, name, perform_hash):
    fasta = tmp_path / "in" / "sample.one.fasta.gz"
    fasta.parent.mkdir()
    with gzip.open(fasta, "wt", newline="") as f:
        f.write(CASES[name])
    new_dir, old_dir = tmp_path / "new", tmp_path / "old"
    new_dir.mkdir(), old_dir.mkdir()
    fasta_to_parquet(str(fasta), new_dir, perform_hash)
    oracle_fasta_to_parquet(str(fasta), old_dir, perform_hash)

    new_file, old_file = next(new_dir.glob("*.parquet.snappy")), next(old_dir.glob("*.parquet.snappy"))
    assert new_file.name == old_file.name == "sample.one.parquet.snappy"
    assert _table(new_file) == _table(old_file)


@pytest.mark.skipif(not SMALL.is_dir(), reason="real test data not available")
@pytest.mark.parametrize("perform_hash", [False, True])
def test_same_table_on_real_data(tmp_path, perform_hash):
    new_dir, old_dir = tmp_path / "new", tmp_path / "old"
    new_dir.mkdir(), old_dir.mkdir()
    files = sorted(SMALL.glob("*.fasta.gz"))
    assert len(files) == 400
    for f in files:
        fasta_to_parquet(str(f), new_dir, perform_hash)
        oracle_fasta_to_parquet(str(f), old_dir, perform_hash)
    for f in files:
        name = f"{f.name[:-len('.fasta.gz')]}.parquet.snappy"
        assert _table(new_dir / name) == _table(old_dir / name), name

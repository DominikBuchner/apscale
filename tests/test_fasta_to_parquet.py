import gzip, hashlib, sys
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from apscale.j_generate_read_table import fasta_to_parquet


def _run(tmp_path: Path, records: str, perform_hash: bool):
    fasta = tmp_path / "sample_1.fasta.gz"
    with gzip.open(fasta, "wt") as f:
        f.write(records)
    fasta_to_parquet(str(fasta), tmp_path, perform_hash)
    return duckdb.sql(
        f"SELECT hash, sequence, read_count FROM read_parquet('{tmp_path}/sample_1.parquet.snappy') "
        "WHERE sequence != 'dummy_seq' ORDER BY read_count"
    ).fetchall()


def test_hash_independent_of_case_and_whitespace(tmp_path):
    rows = _run(tmp_path, ">a;size=3;\nacgtACGT\n>b;size=5;\nACGTACGT\n", perform_hash=True)
    expected = hashlib.sha3_256(b"ACGTACGT").hexdigest()
    assert rows == [(expected, "ACGTACGT", 3), (expected, "ACGTACGT", 5)]


def test_header_hash_is_kept_without_hashing(tmp_path):
    rows = _run(tmp_path, ">abc123;size=7\nacgt\n", perform_hash=False)
    assert rows == [("abc123", "ACGT", 7)]

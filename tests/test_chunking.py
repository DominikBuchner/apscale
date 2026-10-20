"""Chunked parquet pivot must neither duplicate nor lose rows, whatever the physical table order is.

Builds a read data storage from the real small dataset (skipped if it is not available).
"""
import shutil, sys
from pathlib import Path

import duckdb
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from apscale import j_generate_read_table as j

SMALL = Path("D:/apscale_test_data/small_dataset")
pytestmark = pytest.mark.skipif(not SMALL.is_dir(), reason="real test data not available")


@pytest.fixture(scope="module")
def storage(tmp_path_factory):
    root = tmp_path_factory.mktemp("storage")
    temp, data = root / "temp", root / "data"
    temp.mkdir(), data.mkdir()
    for f in SMALL.glob("*.fasta.gz"):
        j.fasta_to_parquet(str(f), temp, False)
    db = data / "read_data_storage.duckdb"
    j.build_read_store(temp, db)
    j.generate_fasta("t", data, db, "sequence")
    j.generate_sequence_groups("t", data, db, temp, 0.97)
    return db


@pytest.mark.parametrize("sequence_type", ["sequence", "group"])
@pytest.mark.parametrize("shuffle", [False, True])
def test_chunked_parquet_complete_and_ordered(storage, tmp_path, monkeypatch, sequence_type, shuffle):
    db = tmp_path / "db.duckdb"
    shutil.copy(storage, db)
    if shuffle:  # same content, physical row order unrelated to sequence_idx
        con = duckdb.connect(str(db))
        for t in ("sequence_data", "group_data"):
            con.execute(f"CREATE OR REPLACE TABLE {t} AS SELECT * FROM {t} ORDER BY random()")
        con.close()
    monkeypatch.setattr(j, "MIN_ROWS_PER_CHUNK", 7)  # force many small chunks
    monkeypatch.setattr(j, "MAX_PIVOT_CELLS", 1)
    out, temp = tmp_path / "out", tmp_path / "temp"
    out.mkdir(), temp.mkdir()
    j.generate_read_table("t", out, db, sequence_type, temp)

    parquet = next(out.glob("*_read_table.parquet.snappy"))
    n_expected = duckdb.connect(str(db)).execute(
        f"SELECT count(*) FROM {sequence_type}_data WHERE sequence != 'dummy_seq'"
    ).fetchone()[0]
    rows = duckdb.sql(f"SELECT count(*), count(DISTINCT hash) FROM read_parquet('{parquet}')").fetchone()
    assert rows == (n_expected, n_expected)
    hashes = [r[0] for r in duckdb.sql(f"SELECT hash FROM read_parquet('{parquet}')").fetchall()]
    expected = [r[0] for r in duckdb.connect(str(db)).execute(
        f"SELECT hash FROM {sequence_type}_data WHERE sequence != 'dummy_seq' ORDER BY sequence_order").fetchall()]
    assert hashes == expected

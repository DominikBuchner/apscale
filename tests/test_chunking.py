"""Chunked parquet pivot must neither duplicate nor lose rows, whatever the physical table order is.

Uses the `storage` fixture (conftest.py).
"""
import shutil

import duckdb
import pytest

from apscale import j_generate_read_table as j


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

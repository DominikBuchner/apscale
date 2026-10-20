"""The group mapping log must be split into numbered files without losing rows."""
import shutil

import duckdb
import pandas as pd

from apscale import j_generate_read_table as j


def test_log_split_into_numbered_files(storage, tmp_path, monkeypatch):
    db = tmp_path / "db.duckdb"
    shutil.copy(storage, db)
    monkeypatch.setattr(j, "LOG_ROWS_PER_FILE", 100)
    stale = tmp_path / "Logfile_11_read_table_part_99.xlsx"  # left over from an earlier run
    stale.write_text("old")

    j.write_log(tmp_path, db)

    files = sorted(tmp_path.glob("Logfile_11_read_table_part_*.xlsx"), key=lambda f: int(f.stem.split("_")[-1]))
    n_rows = duckdb.connect(str(db)).execute("SELECT count(*) FROM group_mapping").fetchone()[0]
    assert [f.name[-6:-5] for f in files] == [str(i) for i in range(1, -(-n_rows // 100) + 1)]
    assert not stale.exists()

    log = pd.concat([pd.read_excel(f) for f in files], ignore_index=True)
    assert list(log.columns) == ["group_hash", "sequence_hash", "read_sum"]
    assert len(log) == n_rows
    # every sequence appears exactly once and belongs to the group stored in the database
    expected = duckdb.connect(str(db)).execute("""
        SELECT gd.hash, sd.hash FROM group_mapping gm
        JOIN sequence_data sd ON gm.sequence_idx = sd.sequence_idx
        JOIN group_data gd ON gm.group_idx = gd.sequence_idx""").fetchall()
    assert sorted(zip(log["group_hash"], log["sequence_hash"])) == sorted(expected)
    # groups in order of group_data.sequence_order, so the first group is the most abundant one
    first_group = duckdb.connect(str(db)).execute(
        "SELECT hash FROM group_data ORDER BY sequence_order LIMIT 1").fetchone()[0]
    assert log["group_hash"].iloc[0] == first_group

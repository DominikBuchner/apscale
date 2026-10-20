"""Sequence grouping with `vsearch --cluster_size` must reproduce the exhaustive all-vs-all reference."""
import shutil

import duckdb
import pytest

from apscale import j_generate_read_table as j

from oracle_grouping import oracle_generate_sequence_groups

pytestmark = pytest.mark.skipif(shutil.which("vsearch") is None, reason="vsearch not installed")


def _mapping(db):
    con = duckdb.connect(str(db), read_only=True)
    rows = con.execute("""
        SELECT m.hash, g.hash FROM group_mapping gm
        JOIN sequence_data m ON gm.sequence_idx = m.sequence_idx
        JOIN sequence_data g ON gm.group_idx = g.sequence_idx ORDER BY 1""").fetchall()
    con.close()
    return rows


def _fresh_copy(storage, tmp_path, name):
    folder = tmp_path / name
    (folder / "temp").mkdir(parents=True)
    db = folder / "db.duckdb"
    shutil.copy(storage, db)
    return folder, db


@pytest.mark.parametrize("threshold", [0.97, 0.9])
def test_mapping_equals_exhaustive_reference(storage, tmp_path, threshold):
    ref_folder, ref_db = _fresh_copy(storage, tmp_path, "reference")
    j.generate_fasta("t", ref_folder, ref_db, "sequence")
    oracle_generate_sequence_groups("t", ref_folder, ref_db, ref_folder / "temp", threshold)

    new_folder, new_db = _fresh_copy(storage, tmp_path, "new")
    j.generate_sequence_groups("t", new_folder, new_db, new_folder / "temp", threshold)

    reference, new = _mapping(ref_db), _mapping(new_db)
    assert len(new) == duckdb.connect(str(new_db)).execute(
        "SELECT count(*) FROM sequence_data WHERE sequence != 'dummy_seq'").fetchone()[0]
    assert new == reference
    assert len({g for _, g in new}) < len(new)  # something was actually grouped


def test_mapping_independent_of_threads(storage, tmp_path, monkeypatch):
    results = []
    for threads in (1, 3):
        folder, db = _fresh_copy(storage, tmp_path, f"t{threads}")
        monkeypatch.setattr(j, "VSEARCH_EXTRA_ARGS", ["--threads", str(threads)])
        j.generate_sequence_groups("t", folder, db, folder / "temp", 0.97)
        results.append(_mapping(db))
    assert results[0] == results[1]


def test_group_tables_consistent(storage, tmp_path):
    folder, db = _fresh_copy(storage, tmp_path, "tables")
    j.generate_sequence_groups("t", folder, db, folder / "temp", 0.97)
    con = duckdb.connect(str(db), read_only=True)
    # group read sums == sum of member read sums, total reads unchanged
    assert con.execute("SELECT sum(read_sum) FROM group_data").fetchone() == con.execute(
        "SELECT sum(read_sum) FROM sequence_data WHERE sequence != 'dummy_seq'").fetchone()
    assert con.execute("SELECT sum(read_count) FROM group_read_count_data").fetchone() == con.execute(
        "SELECT sum(read_count) FROM sequence_read_count_data").fetchone()
    # every group representative is a member of its own group and the most abundant one (lowest order)
    bad = con.execute("""
        SELECT count(*) FROM group_mapping gm
        JOIN sequence_data m ON gm.sequence_idx = m.sequence_idx
        JOIN sequence_data g ON gm.group_idx = g.sequence_idx
        WHERE g.sequence_order > m.sequence_order""").fetchone()[0]
    assert bad == 0
    assert con.execute("""SELECT count(*) FROM group_data gd JOIN group_mapping gm
        ON gd.sequence_idx = gm.sequence_idx AND gm.group_idx != gm.sequence_idx""").fetchone()[0] == 0


def test_single_sequence_and_empty_input(tmp_path, make_storage):
    db, data, temp = make_storage(tmp_path / "one", ["s1", "s2"], n_sequences=1)
    j.generate_sequence_groups("t", data, db, temp, 0.97)
    assert _mapping(db) == [("hash0", "hash0")]

    db, data, temp = make_storage(tmp_path / "none", ["s1", "s2"], n_sequences=0)
    j.generate_sequence_groups("t", data, db, temp, 0.97)
    assert duckdb.connect(str(db)).execute("SELECT count(*) FROM group_mapping").fetchone()[0] == 0

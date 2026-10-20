"""Outputs of earlier runs must not be mixed up with the ones of the current run."""
from apscale import j_generate_read_table as j

SAMPLES = [f"s{i}" for i in range(10)]


def _parts(folder):
    return sorted(p.name for p in folder.glob("0_t_sequence_read_table_part_*.xlsx"))


def test_old_excel_parts_are_removed(tmp_path, make_storage, monkeypatch):
    db, data, temp = make_storage(tmp_path, SAMPLES)

    monkeypatch.setattr(j, "EXCEL_MAX_CELLS", 12)  # 4 sequences incl. dummy -> 3 samples per part
    j.generate_read_table("t", data, db, "sequence", temp)
    assert len(_parts(data)) == 4

    monkeypatch.setattr(j, "EXCEL_MAX_CELLS", 10_000_000)  # everything fits into one part
    j.generate_read_table("t", data, db, "sequence", temp)
    assert _parts(data) == ["0_t_sequence_read_table_part_0.xlsx"]


def test_other_files_are_kept(tmp_path, make_storage):
    db, data, temp = make_storage(tmp_path, SAMPLES)
    keep = [data / "1_t_group_read_table_part_0.xlsx", data / "0_other_sequence_read_table_part_0.xlsx", data / "notes.xlsx"]
    for f in keep:
        f.write_text("x")
    j.generate_read_table("t", data, db, "sequence", temp)
    assert all(f.exists() for f in keep)


def test_remove_old_read_tables_covers_parquet_and_both_types(tmp_path):
    names = [
        "0_t_sequence_read_table_part_0.xlsx",
        "0_t_sequence_read_table.parquet.snappy",
        "1_t_group_read_table_part_3.xlsx",
        "1_t_group_read_table.parquet.snappy",
    ]
    for n in names:
        (tmp_path / n).write_text("x")
    (tmp_path / "0_t_sequences.fasta").write_text("x")
    j.remove_old_read_tables(tmp_path, "t")
    assert [p.name for p in tmp_path.iterdir()] == ["0_t_sequences.fasta"]

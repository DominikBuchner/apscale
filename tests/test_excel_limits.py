"""Excel files must stay within the sheet limit of 16,384 columns, however many samples there are."""
import duckdb
import openpyxl
import pandas as pd

from apscale import j_generate_read_table as j

N_SAMPLES, N_SEQUENCES = 17_000, 3  # more samples than fit into one sheet


def _storage(root):
    temp, data = root / "temp", root / "data"
    temp.mkdir(), data.mkdir()
    samples = [f"sample_{i:05d}" for i in range(N_SAMPLES)]
    rows = [(s, f"hash{k}", "ACGT" * (k + 1), 1 + (i + k) % 5) for i, s in enumerate(samples) for k in range(N_SEQUENCES)]
    rows += [(s, "dummy_hash", "dummy_seq", 0) for s in samples]
    pd.DataFrame(rows, columns=["sample", "hash", "sequence", "read_count"]).to_parquet(temp / "all.parquet.snappy")
    db = data / "db.duckdb"
    j.build_read_store(temp, db)
    return db, data, temp, samples


def test_excel_parts_respect_column_limit(tmp_path):
    db, data, temp, samples = _storage(tmp_path)
    j.generate_read_table("t", data, db, "sequence", temp)

    parts = sorted(data.glob("0_t_sequence_read_table_part_*.xlsx"))
    columns = []
    for part in parts:
        wb = openpyxl.load_workbook(part, read_only=True)
        ws = wb.active
        header = next(ws.iter_rows(min_row=1, max_row=1, values_only=True))
        assert len(header) <= 16_384
        columns += [c for c in header if c not in ("hash", "sequence")]
        wb.close()
    assert len(parts) >= 2
    assert sorted(columns) == samples

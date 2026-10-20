"""Shared fixtures: a read data storage built from the real small dataset (skipped if unavailable)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from apscale import j_generate_read_table as j

SMALL = Path("D:/apscale_test_data/small_dataset")


@pytest.fixture(scope="session")
def storage(tmp_path_factory):
    if not SMALL.is_dir():
        pytest.skip("real test data not available")
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


@pytest.fixture
def make_storage():
    """Factory for small synthetic read data storages: make_storage(root, sample_names, n_sequences=3)."""
    import pandas as pd

    def _make(root: Path, samples: list, n_sequences: int = 3):
        temp, data = root / "temp", root / "data"
        temp.mkdir(parents=True), data.mkdir(parents=True)
        rows = [
            (s, f"hash{k}", "ACGT" * (k + 1), 1 + (i + k) % 5)
            for i, s in enumerate(samples)
            for k in range(n_sequences)
        ]
        rows += [(s, "dummy_hash", "dummy_seq", 0) for s in samples]
        pd.DataFrame(rows, columns=["sample", "hash", "sequence", "read_count"]).to_parquet(
            temp / "all.parquet.snappy"
        )
        db = data / "db.duckdb"
        j.build_read_store(temp, db)
        return db, data, temp

    return _make

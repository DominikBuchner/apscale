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

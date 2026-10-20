"""The DuckDB connection must respect the configured resources."""
import psutil

from apscale import j_generate_read_table as j


def _setting(connection, name):
    return connection.execute(f"SELECT current_setting('{name}')").fetchone()[0]


def test_memory_limit_and_threads_follow_settings(tmp_path, monkeypatch):
    monkeypatch.setattr(j, "CORES", 6)
    monkeypatch.setattr(j, "DUCKDB_MEMORY_LIMIT", "4GB")
    con = j.connect_database(tmp_path / "db.duckdb", tmp_path / "spill")
    assert 3.5 < j.duckdb_memory_limit_mb(con) / 1024 < 4.0  # 4 GB
    assert int(_setting(con, "threads")) == 6
    assert str(tmp_path / "spill") in _setting(con, "temp_directory") or "spill" in _setting(con, "temp_directory")
    con.close()


def test_threads_are_capped_by_memory_limit(tmp_path, monkeypatch):
    monkeypatch.setattr(j, "CORES", 22)
    monkeypatch.setattr(j, "DUCKDB_MEMORY_LIMIT", "512MB")
    con = j.connect_database(tmp_path / "db.duckdb")
    assert int(_setting(con, "threads")) == int(j.duckdb_memory_limit_mb(con) // j.DUCKDB_MB_PER_THREAD) < 22
    con.close()


def test_default_limit_is_a_share_of_the_ram(tmp_path, monkeypatch):
    monkeypatch.setattr(j, "CORES", None)
    monkeypatch.setattr(j, "DUCKDB_MEMORY_LIMIT", None)
    con = j.connect_database(tmp_path / "db.duckdb")
    expected_mb = psutil.virtual_memory().total * j.DUCKDB_RAM_SHARE / 2**20
    assert abs(j.duckdb_memory_limit_mb(con) - expected_mb) / expected_mb < 0.05
    assert 1 <= int(_setting(con, "threads")) <= psutil.cpu_count()
    con.close()

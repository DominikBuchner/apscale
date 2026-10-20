"""Input validation and cleanup behaviour of j_generate_read_table.main()."""
import gzip, shutil

import openpyxl
import pytest

from apscale import j_generate_read_table as j
from apscale.a_create_project import create_project

needs_vsearch = pytest.mark.skipif(shutil.which("vsearch") is None, reason="vsearch not installed")


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    create_project("demo")
    proj = tmp_path / "demo_apscale"
    settings = proj / "Settings_demo.xlsx"
    wb = openpyxl.load_workbook(settings)
    wb["0_general_settings"]["A2"] = 2
    wb.save(settings)
    return proj


def _add_samples(proj, n=3):
    for i in range(n):
        with gzip.open(proj / "10_nc_removal" / "data" / f"sample_{i}.fasta.gz", "wt") as f:
            f.write(f">h{i};size={10 + i}\n{'ACGT' * 10}\n>shared;size=5\n{'TTGCA' * 8}\n")


def _set_threshold(proj, value):
    settings = proj / "Settings_demo.xlsx"
    wb = openpyxl.load_workbook(settings)
    wb["11_read_table"]["B2"] = value
    wb.save(settings)


@pytest.mark.parametrize("threshold", [0, 1.5, -0.2])
def test_invalid_group_threshold_raises(project, threshold):
    _add_samples(project)
    _set_threshold(project, threshold)
    with pytest.raises(ValueError, match="between 0 and 1"):
        j.main(project)


def test_no_input_files_raises_clear_error(project):
    with pytest.raises(FileNotFoundError, match=r"No \*\.fasta\.gz files"):
        j.main(project)


@needs_vsearch
def test_leftovers_in_temp_folder_do_not_break_run(project):
    _add_samples(project)
    temp = project / "11_read_table" / "temp"
    temp.mkdir(parents=True)
    for leftover in ("temp.duckdb", "sequence_read_count_data.duckdb", "pivot_chunk_9.parquet.snappy", "x.tsv"):
        (temp / leftover).write_text("left over from a crashed run")
    j.main(project)
    assert not temp.exists()


@needs_vsearch
def test_missing_output_folder_is_created(project):
    _add_samples(project)
    shutil.rmtree(project / "11_read_table")
    j.main(project)
    assert list((project / "11_read_table" / "data").glob("*_read_table.parquet.snappy"))


def test_missing_vsearch_gives_clear_error(storage, tmp_path, monkeypatch):
    db = tmp_path / "db.duckdb"
    shutil.copy(storage, db)
    monkeypatch.setenv("PATH", str(tmp_path))  # no vsearch here
    with pytest.raises(RuntimeError, match="vsearch was not found"):
        j.generate_sequence_groups("t", tmp_path, db, tmp_path, 0.97)


@needs_vsearch
def test_failing_vsearch_reports_its_message(storage, tmp_path, monkeypatch):
    db = tmp_path / "db.duckdb"
    shutil.copy(storage, db)
    monkeypatch.setattr(j, "VSEARCH_EXTRA_ARGS", ["--no-such-option"])
    with pytest.raises(RuntimeError, match="vsearch failed"):
        j.generate_sequence_groups("t", tmp_path, db, tmp_path, 0.97)

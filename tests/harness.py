"""Test and benchmark harness for apscale/j_generate_read_table.py.

Usage:
    python tests/harness.py run <dataset_dir> <out_dir> [--label NAME] [--work-dir DIR]

Builds a fresh apscale project from a folder of *.fasta.gz files (the folder is only read,
never modified), runs j_generate_read_table.main() on it, measures wall time and peak RSS
(main process + all children), and stores a canonical snapshot of the results in <out_dir>.
Snapshots are keyed by sample name / hash (never by internal indices), so they can be
compared between code versions.
"""

import argparse, json, os, shutil, subprocess, sys, threading, time
from pathlib import Path

import duckdb
import psutil

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))


class PeakRss:
    """Samples the RSS of this process plus all children in a background thread."""

    def __init__(self, interval: float = 0.1):
        self.interval, self.peak, self._stop = interval, 0, threading.Event()
        self._proc = psutil.Process()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        while not self._stop.is_set():
            total = 0
            try:
                for p in [self._proc] + self._proc.children(recursive=True):
                    try:
                        total += p.memory_info().rss
                    except psutil.Error:
                        pass
            except psutil.Error:
                pass
            self.peak = max(self.peak, total)
            time.sleep(self.interval)

    def __enter__(self):
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._thread.join()


def make_project(dataset_dir: Path, work_dir: Path, name: str) -> Path:
    """Create an apscale project whose last step is 10_nc_removal (hashes come from the headers)."""
    from apscale.a_create_project import create_project

    work_dir.mkdir(parents=True, exist_ok=True)
    project = work_dir / f"{name}_apscale"
    if project.exists():
        shutil.rmtree(project)
    cwd = os.getcwd()
    os.chdir(work_dir)
    try:
        create_project(name)
    finally:
        os.chdir(cwd)
    target = project / "10_nc_removal" / "data"
    for f in sorted(Path(dataset_dir).glob("*.fasta.gz")):
        shutil.copy2(f, target / f.name)
    return project


def snapshot(project: Path, out_dir: Path) -> dict:
    """Write canonical, index-free parquet tables from the read data storage."""
    out_dir.mkdir(parents=True, exist_ok=True)
    db_files = list((project / "11_read_table" / "data").glob("*_read_data_storage.duckdb"))
    assert len(db_files) == 1, db_files
    con = duckdb.connect(str(db_files[0]), read_only=True)
    queries = {
        # sequences with their order and read sum (dummy removed)
        "sequences": """
            SELECT hash, sequence, sequence_order, read_sum FROM sequence_data
            WHERE sequence != 'dummy_seq' ORDER BY hash""",
        # sparse read counts per sample and hash
        "read_counts": """
            SELECT s.sample, q.hash, r.read_count
            FROM sequence_read_count_data r
            JOIN sample_data s USING (sample_idx)
            JOIN sequence_data q USING (sequence_idx)
            WHERE q.sequence != 'dummy_seq' AND r.read_count > 0
            ORDER BY s.sample, q.hash""",
        # member hash -> group hash
        "group_mapping": """
            SELECT m.hash AS member_hash, g.hash AS group_hash
            FROM group_mapping gm
            JOIN sequence_data m ON gm.sequence_idx = m.sequence_idx
            JOIN sequence_data g ON gm.group_idx = g.sequence_idx
            ORDER BY member_hash""",
        "groups": """
            SELECT hash, sequence, sequence_order, read_sum FROM group_data ORDER BY hash""",
        "group_read_counts": """
            SELECT s.sample, q.hash, r.read_count
            FROM group_read_count_data r
            JOIN sample_data s USING (sample_idx)
            JOIN sequence_data q USING (sequence_idx)
            WHERE q.sequence != 'dummy_seq' AND r.read_count > 0
            ORDER BY s.sample, q.hash""",
    }
    info = {}
    for name, query in queries.items():
        path = out_dir / f"{name}.parquet"
        con.execute(f"COPY ({query}) TO '{path}' (FORMAT parquet)")
        info[name] = con.execute(f"SELECT count(*) FROM read_parquet('{path}')").fetchone()[0]
    con.close()
    return info


def run(dataset_dir: Path, out_dir: Path, label: str, work_dir: Path) -> None:
    from apscale import j_generate_read_table

    name = Path(dataset_dir).name
    project = make_project(dataset_dir, work_dir, name)
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = out_dir / "stdout.log"
    start = time.perf_counter()
    with PeakRss() as rss, open(log_path, "w", encoding="utf-8") as log:
        old_stdout, sys.stdout = sys.stdout, log
        try:
            j_generate_read_table.main(project)
        finally:
            sys.stdout = old_stdout
    wall = time.perf_counter() - start

    # copy the produced files, snapshot the canonical tables
    produced = out_dir / "produced"
    if produced.exists():
        shutil.rmtree(produced)
    shutil.copytree(project / "11_read_table", produced)
    for db in produced.rglob("*.duckdb"):
        db.unlink()  # large and not comparable byte-wise; canonical snapshot below
    counts = snapshot(project, out_dir / "snapshot")

    commit = subprocess.run(
        ["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"], capture_output=True, text=True
    ).stdout.strip()
    metrics = {
        "label": label,
        "dataset": name,
        "git_commit": commit,
        "wall_seconds": round(wall, 2),
        "peak_rss_mb": round(rss.peak / 1e6, 1),
        "snapshot_rows": counts,
        "produced_files": sorted(str(p.relative_to(produced)) for p in produced.rglob("*") if p.is_file()),
    }
    (out_dir / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("dataset_dir", type=Path)
    r.add_argument("out_dir", type=Path)
    r.add_argument("--label", default="unlabelled")
    r.add_argument("--work-dir", type=Path, default=Path("D:/apscale_test_data/work"))
    args = ap.parse_args()
    run(args.dataset_dir, args.out_dir, args.label, args.work_dir)

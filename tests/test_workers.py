"""Worker processes of the FASTA conversion must not stay alive for the rest of the run."""
import gzip
import time

import psutil

from apscale import j_generate_read_table as j


def test_no_worker_processes_left_after_conversion(tmp_path):
    files = []
    for i in range(4):
        fasta = tmp_path / f"s{i}.fasta.gz"
        with gzip.open(fasta, "wt") as f:
            f.write(f">h{i};size=3\nACGT\n")
        files.append(str(fasta))

    j.convert_fastas_to_parquet(files, tmp_path, perform_hash=False, cores=2)

    assert len(list(tmp_path.glob("*.parquet.snappy"))) == 4
    # joblib's small resource tracker process may stay, the workers must not
    def workers():
        return [
            c for c in psutil.Process().children(recursive=True)
            if "resource_tracker" not in " ".join(c.cmdline())
        ]

    deadline = time.time() + 10  # workers need a moment to exit
    while workers() and time.time() < deadline:
        time.sleep(0.2)
    assert workers() == []

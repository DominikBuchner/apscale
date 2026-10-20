"""Peak memory (RSS) and time of vsearch for the old and the new grouping call.

Usage: python tests/experiments/vsearch_memory.py <fasta with ;size= headers> [threshold]
"""
import subprocess, sys, tempfile, time
from pathlib import Path

import psutil

fasta = Path(sys.argv[1])
thr = sys.argv[2] if len(sys.argv) > 2 else "0.97"
out = Path(tempfile.mkdtemp())
n_seq = sum(1 for line in open(fasta) if line[0] == ">")
print(f"{n_seq:,} sequences, {psutil.cpu_count()} logical cores")

CALLS = {
    "new: cluster_size, maxaccepts 0": ["--cluster_size", str(fasta), "--sizeorder", "--maxaccepts", "0", "--uc", str(out / "a.uc")],
    "new: cluster_size, maxaccepts 100": ["--cluster_size", str(fasta), "--sizeorder", "--maxaccepts", "100", "--uc", str(out / "b.uc")],
    "new: cluster_size, maxaccepts 0, 1 thread": ["--cluster_size", str(fasta), "--sizeorder", "--maxaccepts", "0", "--threads", "1", "--uc", str(out / "c.uc")],
    "old: usearch_global all-vs-all, maxaccepts 0": ["--usearch_global", str(fasta), "--db", str(fasta), "--maxaccepts", "0", "--self",
                                                     "--userout", str(out / "d.tsv"), "--userfields", "query+target+id"],
}
only = sys.argv[3:] if len(sys.argv) > 3 else list(CALLS)
print(f"{'call':48}{'seconds':>9}{'peak RSS MB':>13}")
for name, args in CALLS.items():
    if not any(o in name for o in only) and only != list(CALLS):
        continue
    t = time.perf_counter()
    proc = subprocess.Popen(["vsearch", "--id", thr, "--quiet", *args])
    p, peak = psutil.Process(proc.pid), 0
    while proc.poll() is None:
        try:
            peak = max(peak, p.memory_info().rss)
        except psutil.Error:
            pass
        time.sleep(0.05)
    print(f"{name:48}{time.perf_counter() - t:9.1f}{peak / 1e6:13.0f}", flush=True)

"""Spread a scaled fasta (labels hash;size=n) over many samples as fasta.gz files, one per sample.

Usage: python tests/experiments/make_scaled_samples.py <scaled.fasta> <n_samples> <out_dir> [seed]
Each sequence appears in 1-6 random samples, its reads are split between them.
"""
import gzip, sys
from pathlib import Path

import numpy as np

src, n_samples, out = sys.argv[1], int(sys.argv[2]), Path(sys.argv[3])
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
out.mkdir(parents=True, exist_ok=True)
files = [gzip.open(out / f"sample_{i:04d}.fasta.gz", "wt", compresslevel=1) for i in range(n_samples)]
records = 0
with open(src) as f:
    while header := f.readline():
        seq = f.readline()
        label, size = header[1:].strip().split(";size=")
        k = min(int(rng.integers(1, 7)), n_samples)
        for sample, share in zip(rng.choice(n_samples, size=k, replace=False), rng.multinomial(int(size) + k, [1 / k] * k)):
            files[sample].write(f">{label};size={max(int(share), 1)}\n{seq}")
            records += 1
for f in files:
    f.close()
print(records, "records in", n_samples, "samples")

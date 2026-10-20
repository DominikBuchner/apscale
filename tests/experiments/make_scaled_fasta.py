"""Blow up the real sequences of a run to a larger synthetic dataset for memory/time scaling tests.

Usage: python tests/experiments/make_scaled_fasta.py <sequences.parquet (hash, sequence)> <n_sequences> <out.fasta> [seed]
Each new sequence is a copy of a random real sequence with 0-40 random substitutions (0-13 % divergence), so a
mix of near duplicates and new centroids appears. Read sums follow a heavy-tailed distribution. Labels are hashes.
"""
import hashlib, sys

import duckdb
import numpy as np

src, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
base = [r[0] for r in duckdb.sql(f"SELECT sequence FROM read_parquet('{src}')").fetchall()]
seen, written = set(), 0
with open(out, "w") as f:
    while written < n:
        seq = list(base[rng.integers(len(base))])
        for pos in rng.integers(len(seq), size=int(rng.integers(0, 41))):
            seq[pos] = "ACGT"[rng.integers(4)]
        seq = "".join(seq)
        if seq in seen:
            continue
        seen.add(seq)
        size = int(rng.pareto(1.2) * 3) + 1
        f.write(f">{hashlib.sha3_256(seq.encode()).hexdigest()};size={size}\n{seq}\n")
        written += 1
print(written, "sequences written to", out)

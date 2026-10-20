"""Experiment: do vsearch defaults (maxrejects 32, dust masking) change hits and sequence groups?

Usage: python tests/experiments/vsearch_params.py <data_dir_of_a_finished_run> <scratch_dir> [threshold]
<data_dir> holds <project>_read_data_storage.duckdb and 0_<project>_sequences.fasta (11_read_table/data of a run).
Nothing in <data_dir> is modified, every variant works on a copy.
"""
import shutil, subprocess, sys, time
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from apscale import j_generate_read_table as j

data_dir, scratch = Path(sys.argv[1]), Path(sys.argv[2])
threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.97
db_src = next(data_dir.glob("*_read_data_storage.duckdb"))
fasta_src = next(data_dir.glob("0_*_sequences.fasta"))
project = fasta_src.name[2:-len("_sequences.fasta")]

VARIANTS = {
    "default": [],
    "maxrejects0": ["--maxrejects", "0"],
    "nomask": ["--qmask", "none", "--dbmask", "none"],
    "maxrejects0+nomask": ["--maxrejects", "0", "--qmask", "none", "--dbmask", "none"],
}
shutil.rmtree(scratch, ignore_errors=True)
results = {}
for name, extra in VARIANTS.items():
    folder = scratch / name.replace("+", "_")
    (folder / "temp").mkdir(parents=True)
    shutil.copy(db_src, folder / db_src.name)
    shutil.copy(fasta_src, folder / fasta_src.name)

    # raw hits, to see what the search finds at all
    t = time.perf_counter()
    hits = folder / "hits.tsv"
    subprocess.run(
        ["vsearch", "--usearch_global", str(fasta_src), "--db", str(fasta_src), "--id", str(threshold), "--userout", str(hits),
         "--quiet", "--userfields", "query+target+id", "--maxaccepts", "0", "--self"] + extra,
        check=True,
    )
    t_vsearch = time.perf_counter() - t

    # grouping with the real function
    j.VSEARCH_EXTRA_ARGS = extra
    db = folder / db_src.name
    j.generate_sequence_groups(project, folder, db, folder / "temp", threshold)
    con = duckdb.connect(str(db), read_only=True)
    mapping = folder / "mapping.parquet"
    con.execute(f"""COPY (SELECT m.hash AS member, g.hash AS grp FROM group_mapping gm
        JOIN sequence_data m ON gm.sequence_idx = m.sequence_idx
        JOIN sequence_data g ON gm.group_idx = g.sequence_idx ORDER BY member) TO '{mapping}' (FORMAT parquet)""")
    con.close()
    results[name] = (hits, mapping, t_vsearch)

Q = lambda f: f"read_csv('{f}', delim='\t', header=false, columns={{'q':'VARCHAR','t':'VARCHAR','id':'DOUBLE'}})"
base_hits, base_map, _ = results["default"]
print(f"\n{'variant':22} {'vsearch s':>9} {'hit pairs':>10} {'only here':>10} {'only default':>12} {'groups':>7} {'members regrouped':>18}")
for name, (hits, mapping, t_v) in results.items():
    n = duckdb.sql(f"select count(*) from {Q(hits)}").fetchone()[0]
    only_here = duckdb.sql(f"select count(*) from (select q,t from {Q(hits)} except select q,t from {Q(base_hits)})").fetchone()[0]
    only_base = duckdb.sql(f"select count(*) from (select q,t from {Q(base_hits)} except select q,t from {Q(hits)})").fetchone()[0]
    groups = duckdb.sql(f"select count(distinct grp) from read_parquet('{mapping}')").fetchone()[0]
    regrouped = duckdb.sql(f"select count(*) from read_parquet('{mapping}') a join read_parquet('{base_map}') b using (member) where a.grp != b.grp").fetchone()[0]
    print(f"{name:22} {t_v:9.1f} {n:10,} {only_here:10,} {only_base:12,} {groups:7,} {regrouped:18,}")

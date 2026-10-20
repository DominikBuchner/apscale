"""Compare the canonical snapshots of two harness runs: python tests/compare_runs.py <snapshot_a> <snapshot_b>"""
import sys, duckdb

a, b = sys.argv[1], sys.argv[2]
P = lambda d, t: f"read_parquet('{d}/{t}.parquet')"
checks = {
    "sequence_order verschieden": f"SELECT count(*) FROM {P(a,'sequences')} x JOIN {P(b,'sequences')} y USING (hash) WHERE x.sequence_order != y.sequence_order",
    "sequences nur in einem Lauf": f"SELECT (SELECT count(*) FROM {P(a,'sequences')} WHERE hash NOT IN (SELECT hash FROM {P(b,'sequences')})) + (SELECT count(*) FROM {P(b,'sequences')} WHERE hash NOT IN (SELECT hash FROM {P(a,'sequences')}))",
    "Gruppenzuordnung verschieden": f"SELECT count(*) FROM {P(a,'group_mapping')} x JOIN {P(b,'group_mapping')} y USING (member_hash) WHERE x.group_hash != y.group_hash",
    "Gruppenrepraesentanten nur in einem Lauf": f"SELECT (SELECT count(*) FROM {P(a,'groups')} WHERE hash NOT IN (SELECT hash FROM {P(b,'groups')})) + (SELECT count(*) FROM {P(b,'groups')} WHERE hash NOT IN (SELECT hash FROM {P(a,'groups')}))",
    "read_counts verschieden": f"SELECT count(*) FROM (SELECT * FROM {P(a,'read_counts')} EXCEPT SELECT * FROM {P(b,'read_counts')})",
    "group_read_counts verschieden": f"SELECT count(*) FROM (SELECT * FROM {P(a,'group_read_counts')} EXCEPT SELECT * FROM {P(b,'group_read_counts')})",
}
for name, q in checks.items():
    print(f"{name}: {duckdb.sql(q).fetchone()[0]}")

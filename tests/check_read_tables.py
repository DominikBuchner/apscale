"""Check the parquet read tables of a harness run against its canonical snapshot.

Usage: python tests/check_read_tables.py <run_dir>     (run_dir contains produced/ and snapshot/)
Exit code 1 if any check fails.
"""
import sys
from pathlib import Path

import duckdb

run = Path(sys.argv[1])
snap, data = run / "snapshot", next((run / "produced" / "data").glob("0_*_sequence_read_table.parquet.snappy")).parent
failed = False


def check(name: str, ok: bool, detail: str = "") -> None:
    global failed
    failed |= not ok
    print(f"{'OK  ' if ok else 'FAIL'} {name} {detail}")


for kind, table, counts, prefix in (
    ("sequence", "sequences", "read_counts", "0_"),
    ("group", "groups", "group_read_counts", "1_"),
):
    parquet = next(data.glob(f"{prefix}*_read_table.parquet.snappy"))
    P = f"read_parquet('{parquet}')"
    S = lambda t: f"read_parquet('{snap}/{t}.parquet')"
    n_rows = duckdb.sql(f"SELECT count(*) FROM {P}").fetchone()[0]
    n_expected = duckdb.sql(f"SELECT count(*) FROM {S(table)}").fetchone()[0]
    check(f"{kind}: rows", n_rows == n_expected, f"{n_rows} vs {n_expected}")

    wrong_order = duckdb.sql(f"""
        SELECT count(*) FROM (SELECT hash, row_number() OVER () AS pos FROM {P}) p
        JOIN {S(table)} s USING (hash) WHERE p.pos != s.sequence_order""").fetchone()[0]
    check(f"{kind}: row order == sequence_order", wrong_order == 0, f"{wrong_order} rows differ")

    cols = [c for c in duckdb.sql(f"SELECT * FROM {P} LIMIT 0").columns if c not in ("hash", "sequence")]
    nulls = duckdb.sql(
        f"SELECT count(*) FROM {P} WHERE " + " OR ".join(f'"{c}" IS NULL' for c in cols)
    ).fetchone()[0]
    check(f"{kind}: no NULL counts", nulls == 0, f"{nulls} rows with NULL")

    # per-sample sums, parquet vs. sparse snapshot
    unpivoted = duckdb.sql(f"SELECT hash, UNPIVOT_COLS FROM {P}".replace("UNPIVOT_COLS", ", ".join(f'"{c}"' for c in cols)))
    diff = duckdb.sql(f"""
        WITH a AS (SELECT sample, sum(read_count)::BIGINT AS s FROM (
                       UNPIVOT unpivoted ON COLUMNS(* EXCLUDE hash) INTO NAME sample VALUE read_count) GROUP BY sample),
             b AS (SELECT sample, sum(read_count)::BIGINT AS s FROM {S(counts)} GROUP BY sample)
        SELECT count(*) FROM a FULL JOIN b USING (sample) WHERE COALESCE(a.s, 0) != COALESCE(b.s, 0)""").fetchone()[0]
    check(f"{kind}: per-sample read sums", diff == 0, f"{diff} samples differ")

sys.exit(1 if failed else 0)

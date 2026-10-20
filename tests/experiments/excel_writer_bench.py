"""Prototype benchmark for step 3.1: time and peak memory of three ways to write a read table to Excel.

Usage: python tests/experiments/excel_writer_bench.py <out_dir> [case ...]
Cases mirror the largest file apscale writes (about 10M cells):
  long: 800,000 rows x (hash, sequence, 10 count columns)
  wide: 10,000 rows x (hash, sequence, 1,000 count columns)
Every measurement runs in its own process, memory is the peak RSS of that process.
"""
import subprocess, sys, time
from pathlib import Path

import psutil

CASES = {"long": (800_000, 10), "wide": (10_000, 1_000)}
METHODS = ["pandas", "duckdb_xlsx", "xlsxwriter_stream"]


def query(rows: int, cols: int) -> str:
    counts = ", ".join(f'((i * {c + 1}) % 977)::BIGINT AS "sample_{c + 1}"' for c in range(cols))
    return (f"SELECT sha256(i::VARCHAR) AS hash, substr(repeat(md5(i::VARCHAR), 10), 1, 313) AS sequence, {counts} "
            f"FROM range({rows}) t(i)")


def child(case: str, method: str, out: str) -> None:
    import duckdb

    rows, cols = CASES[case]
    con = duckdb.connect()
    q = query(rows, cols)
    if method == "pandas":
        con.execute(q).df().to_excel(out, index=False, engine="xlsxwriter")
    elif method == "duckdb_xlsx":
        con.execute("LOAD excel")
        con.execute(f"COPY ({q}) TO '{out}' (FORMAT xlsx, HEADER true)")
    elif method == "xlsxwriter_stream":
        import xlsxwriter

        cursor = con.execute(q)
        wb = xlsxwriter.Workbook(out, {"constant_memory": True})
        ws = wb.add_worksheet()
        ws.write_row(0, 0, [c[0] for c in cursor.description])
        row = 1
        while batch := cursor.fetchmany(20_000):
            for r in batch:
                ws.write_row(row, 0, r)
                row += 1
        wb.close()


if __name__ == "__main__":
    if sys.argv[1] == "--child":
        child(*sys.argv[2:5])
        sys.exit()
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"{'case':6}{'method':20}{'seconds':>9}{'peak RSS MB':>13}{'file MB':>9}")
    for case in sys.argv[2:] or CASES:
        for method in METHODS:
            out = out_dir / f"{case}_{method}.xlsx"
            t = time.perf_counter()
            proc = subprocess.Popen([sys.executable, __file__, "--child", case, method, str(out)])
            p, peak = psutil.Process(proc.pid), 0
            while proc.poll() is None:
                try:
                    peak = max(peak, p.memory_info().rss)
                except psutil.Error:
                    pass
                time.sleep(0.05)
            size = out.stat().st_size / 1e6 if out.exists() else float("nan")
            print(f"{case:6}{method:20}{time.perf_counter() - t:9.1f}{peak / 1e6:13.0f}{size:9.1f}   rc={proc.returncode}", flush=True)

"""Peak memory per phase of a harness run: python tests/phase_memory.py <run_dir>
Aligns memory_timeline.csv with the HH:MM:SS timestamps of stdout.log."""
import csv, re, sys, time
from pathlib import Path

run = Path(sys.argv[1])
samples = [(float(r["epoch"]), int(r["total_rss"]), int(r["main_rss"]), int(r["children"])) for r in csv.DictReader(open(run / "memory_timeline.csv"))]
day0 = time.mktime(time.strptime(time.strftime("%Y-%m-%d", time.localtime(samples[0][0])), "%Y-%m-%d"))
phases = []
for line in open(run / "stdout.log", encoding="utf-8", errors="replace"):
    m = re.match(r"(\d\d):(\d\d):(\d\d): (.*)", line.strip())
    if m:
        h, mi, s, text = m.groups()
        phases.append((day0 + int(h) * 3600 + int(mi) * 60 + int(s), text[:70]))
print(f"{'start':>8} {'sec':>5} {'peak total MB':>14} {'peak main MB':>13} {'max children':>13}  phase")
for i, (start, text) in enumerate(phases):
    end = phases[i + 1][0] if i + 1 < len(phases) else samples[-1][0]
    window = [s for s in samples if start - 0.5 <= s[0] <= end + 0.5]
    if window:
        print(f"{time.strftime('%H:%M:%S', time.localtime(start)):>8} {end - start:5.0f} {max(w[1] for w in window) / 1e6:14.0f} {max(w[2] for w in window) / 1e6:13.0f} {max(w[3] for w in window):13d}  {text}")

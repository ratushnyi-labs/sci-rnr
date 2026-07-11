#!/usr/bin/env python3
"""bench_mem.py -- the memory-decoupling proof (the core claim).

Two families of isolated peak-RSS probes (each a fresh mem_probe.py process):

  A. DECOUPLING FROM FILE SIZE.  The hashed model's peak RSS is flat across
     growing input, because its tables are pre-allocated to a fixed
     2^HBITS * (W+1) size and never grow.  The exact-dictionary baselines
     (rnr1 NGramPredictor, cm CMPredictor) allocate one entry per distinct
     context, so their RSS climbs ~linearly with bytes and would OOM at the
     64 MB / 256+ MB scale.  This is the whole point of the design.

  B. BOUNDED GROWTH IN W.  At fixed HBITS the hashed peak RSS grows only
     linearly in (W+1) tables (and is otherwise file-size-independent);
     table_mb = 3 * (W+1) * 2^HBITS is reported so the bound is explicit.
     At a small HBITS the interpreter base dominates and W=16 stays within
     ~1.5x of W=3 (the run_checks memory gate).

Writes results/mem_size.csv and results/mem_w.csv.  Serialized, one process
at a time.  Run with -u.

    python -u bench_mem.py
"""

import json
import subprocess
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
PY = sys.executable
PROBE = str(_HERE / "mem_probe.py")


def probe(model, W, HBITS, size_mb, corpus="text"):
    p = subprocess.run(
        [PY, "-u", PROBE, model, "--W", str(W), "--HBITS", str(HBITS),
         "--size-mb", str(size_mb), "--corpus", corpus],
        capture_output=True, text=True, timeout=7200)
    if p.returncode != 0:
        raise RuntimeError("probe failed: " + p.stderr.strip()[:300])
    return json.loads(p.stdout.strip().splitlines()[-1])


def write_csv(path, rows, fields):
    import csv
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def part_a_decoupling_from_size():
    """Hashed (flat) vs exact-dict baselines (growing) across file size."""
    rows = []
    W, HBITS = 6, 22
    sizes = [0.25, 0.5, 1.0]
    print("== A. peak RSS vs FILE SIZE (W=%d, HBITS=%d) ==" % (W, HBITS),
          flush=True)
    for model in ("hashed", "dict"):
        for s in sizes:
            r = probe(model, W, HBITS, s)
            rows.append(r)
            print("  %-6s size=%4.2fMB  peakRSS=%7.1fMB  table=%s  (%.0fs)"
                  % (model, s, r["peak_rss_mb"],
                     ("%dMB" % r["table_mb"]) if r["table_mb"] else "-",
                     r["wall_s"]), flush=True)
        write_csv(_HERE / "results" / "mem_size.csv", rows,
                  ["model", "corpus", "W", "HBITS", "size_mb",
                   "peak_rss_mb", "table_mb", "wall_s"])
    return rows


def part_b_bounded_in_w():
    """Hashed peak RSS across W at fixed HBITS (bounded, ~linear in W+1)."""
    rows = []
    HBITS = 22
    size_mb = 0.25
    print("== B. peak RSS vs W (HBITS=%d, size=%.2fMB) ==" % (HBITS, size_mb),
          flush=True)
    for W in (3, 6, 9, 12, 16):
        r = probe("hashed", W, HBITS, size_mb)
        rows.append(r)
        print("  W=%2d  peakRSS=%7.1fMB  table=%dMB  (%.0fs)"
              % (W, r["peak_rss_mb"], r["table_mb"], r["wall_s"]), flush=True)
        write_csv(_HERE / "results" / "mem_w.csv", rows,
                  ["model", "W", "HBITS", "size_mb", "peak_rss_mb",
                   "table_mb", "wall_s"])
    if len(rows) >= 2:
        ratio = rows[-1]["peak_rss_mb"] / rows[0]["peak_rss_mb"]
        print("  W=16/W=3 peak-RSS ratio at HBITS=%d: %.2fx "
              "(bounded; grows with the per-order table count)" % (HBITS,
              ratio), flush=True)
    return rows


if __name__ == "__main__":
    part_a_decoupling_from_size()
    part_b_bounded_in_w()
    print("wrote results/mem_size.csv, results/mem_w.csv", flush=True)

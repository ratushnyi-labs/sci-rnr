#!/usr/bin/env python3
"""Mass-repetition harness for timing statistics.

Rates are BIT-DETERMINISTIC: every repetition of a (tier, type, coder)
cell must produce a byte-identical archive, so compression numbers need
one run.  What repetition estimates is the TIMING / RSS distribution.
This helper runs N repetitions of selected cells efficiently (corpus
materialized once, cells appended with increasing run_idx, journal-
checkpointed so an interrupted 100000-run campaign resumes), asserts the
determinism witness (identical compressed_sha256 across repetitions),
and prints distribution summaries.

Feasibility guide (wall time PER CELL, this host):
    tier 500  x brotli-q11 : ~20 min/rep  -> 100000 reps infeasible
    tier 500  x zstd-19    : ~75 s/rep    -> 1000 reps ~ 21 h
    tier smoke (16 MiB) x zstd-19 : ~2.5 s/rep -> 100000 reps ~ 3 days
    tier smoke x gzip/zstd-fast lanes: seconds -> 100000 reps feasible
Pick the tier to match the repetition count; the record schema is
identical across tiers, so distributions compose in analysis.

Usage:
    python repeat_stats.py --tier smoke --types text --coders zstd-19 \
        --reps 1000 [--batch 50]
Appends to out/<tier>/results.jsonl (same store as the campaign) and
re-exports out/results.csv including run_idx, so downstream statistics
scripts see every repetition as a row.
"""
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable
OUT = HERE / "out"


def existing_runs(tier: str) -> dict:
    """cell prefix -> (max run_idx seen, {compressed_sha256 set})"""
    seen = defaultdict(lambda: [-1, set()])
    p = OUT / tier / "results.jsonl"
    if p.exists():
        with open(p) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if r.get("status") != "done":
                    continue
                key = (r["ctype"], r["coder"])
                seen[key][0] = max(seen[key][0], r.get("run_idx", 0))
                sha = (r.get("measured") or {}).get("compressed_sha256")
                if sha:
                    seen[key][1].add(sha)
    return seen


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tier", required=True)
    ap.add_argument("--types", required=True)
    ap.add_argument("--coders", required=True)
    ap.add_argument("--reps", type=int, required=True,
                    help="TOTAL repetitions per cell (existing runs count)")
    ap.add_argument("--batch", type=int, default=100,
                    help="campaign invocations are chunked so a kill only "
                         "loses the in-flight cell")
    args = ap.parse_args(argv)

    types = args.types.split(",")
    coders = args.coders.split(",")

    # drive run_campaign with increasing --runs; its journal skips done cells
    target = args.reps
    step = max(1, args.batch)
    reached = 0
    while reached < target:
        reached = min(target, reached + step)
        cmd = [PY, "-u", str(HERE / "run_campaign.py"),
               "--tier", args.tier, "--types", args.types,
               "--coders", args.coders, "--runs", str(reached)]
        rc = subprocess.call(cmd)
        if rc != 0:
            print(f"[repeat] campaign exited rc={rc} at runs={reached}; "
                  f"journal preserved -- rerun to resume")
            return rc

    # determinism assertion + timing summary
    ok = True
    seen = existing_runs(args.tier)
    p = OUT / args.tier / "results.jsonl"
    series = defaultdict(lambda: {"enc": [], "dec": []})
    with open(p) as f:
        for line in f:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get("status") != "done":
                continue
            key = (r["ctype"], r["coder"])
            if r["ctype"] in types and r["coder"] in coders:
                m = r["measured"]
                series[key]["enc"].append(m["enc_s"])
                series[key]["dec"].append(m["dec_s"])
    print(f"{'cell':<40}{'reps':>6}{'shas':>6}"
          f"{'enc med':>10}{'enc p95':>10}{'dec med':>10}")
    for key, s in sorted(series.items()):
        shas = seen[key][1]
        n = len(s["enc"])
        det = len(shas) <= 1
        ok &= det
        enc, dec = sorted(s["enc"]), sorted(s["dec"])
        med = statistics.median(enc)
        p95 = enc[min(n - 1, int(0.95 * n))]
        dmed = statistics.median(dec)
        flag = "" if det else "  << NONDETERMINISTIC OUTPUT"
        print(f"{'/'.join(key):<40}{n:>6}{len(shas):>6}"
              f"{med:>10.3f}{p95:>10.3f}{dmed:>10.3f}{flag}")
    subprocess.call([PY, str(HERE / "to_csv.py"),
                     "--tiers", f"500,1500,5000,{args.tier}"])
    if not ok:
        print("[repeat] FAIL: some cell produced differing archives across "
              "repetitions -- investigate before using timing statistics")
        return 1
    print("[repeat] determinism witness OK: one archive sha per cell")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

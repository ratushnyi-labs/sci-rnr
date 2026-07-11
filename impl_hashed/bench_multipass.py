#!/usr/bin/env python3
"""bench_multipass.py -- does epoch-training the mixer weights help the rate?

Measures, at a small (reaper-safe) slice, the coded-pass cross-entropy for
passes = 1..P at fixed (W, HBITS), plus the effective bpb once the shipped
weight blob is charged.  The coded CE is the realizable rate; the weight
overhead is (W+1)*256*4 bytes, negligible at scale and dominant on tiny
files.  Resumable; writes results/multipass.csv.

    python -u bench_multipass.py [--size-mb 0.5] [--W 12] [--HBITS 22]
"""

import argparse
import csv
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "data"))
import bench_util as bu  # noqa: E402
from hashed_predictor import HashedCMPredictor  # noqa: E402
from hashed_coder import train_weights  # noqa: E402

FIELDS = ["corpus", "size_mb", "W", "HBITS", "passes",
          "coded_ce_bpb", "weight_overhead_bpb", "effective_bpb",
          "train_time_s", "per_pass_time_s"]


def load_done(out):
    rows, done = [], set()
    if out.exists():
        for r in csv.DictReader(open(out)):
            rows.append(r)
            done.add((r["corpus"], float(r["size_mb"]), int(r["W"]),
                      int(r["HBITS"]), int(r["passes"])))
    return rows, done


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--size-mb", type=float, default=0.5)
    ap.add_argument("--W", type=int, default=12)
    ap.add_argument("--HBITS", type=int, default=22)
    ap.add_argument("--passes", default="1,2,3,4")
    ap.add_argument("--corpora", default="text,code")
    args = ap.parse_args(argv)
    out = _HERE / "results" / "multipass.csv"
    rows, done = load_done(out)
    for corpus in args.corpora.split(","):
        data = bu.load_slice(corpus, int(args.size_mb * bu.MB))
        n = len(data)
        for P in [int(x) for x in args.passes.split(",")]:
            key = (corpus, args.size_mb, args.W, args.HBITS, P)
            if key in done:
                print("(cached) %s P=%d" % (corpus, P), flush=True)
                continue
            per_pass, train_s = [], 0.0
            if P <= 1:
                p = HashedCMPredictor(args.W, HBITS=args.HBITS)
                w_over = 0.0
            else:
                w, per_pass = train_weights(data, args.W, args.HBITS, P)
                train_s = sum(per_pass)
                p = HashedCMPredictor(args.W, HBITS=args.HBITS)
                p.load_weights(w)
                p.freeze = True
                p.reset_counts()
                w_over = 8.0 * (w.size * 4) / n
            t0 = time.perf_counter()
            coded = bu.model_ce_bpb(p, data)
            rows.append({
                "corpus": corpus, "size_mb": args.size_mb, "W": args.W,
                "HBITS": args.HBITS, "passes": P, "coded_ce_bpb": coded,
                "weight_overhead_bpb": w_over, "effective_bpb": coded + w_over,
                "train_time_s": train_s,
                "per_pass_time_s": ";".join("%.1f" % t for t in per_pass)})
            bu.csv_write(out, rows, FIELDS)
            print("[%s %.2gMB] passes=%d  coded=%.4f  eff=%.4f bpb  "
                  "(train %.0fs, code %.0fs)" % (corpus, args.size_mb, P,
                  coded, coded + w_over, train_s, time.perf_counter() - t0),
                  flush=True)
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()

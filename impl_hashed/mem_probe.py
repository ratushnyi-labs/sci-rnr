#!/usr/bin/env python3
"""mem_probe.py -- isolated peak-RSS probe for one predictor over one slice.

Runs a single predictor (hashed / exact-dict n-gram / exact-dict CM) over
the first size-MB of a corpus in a FRESH process and prints one JSON line
with peak RSS (ru_maxrss, bytes on macOS).  Isolation gives a clean peak per
config and keeps the host to one process at a time.

    python -u mem_probe.py <hashed|dict|cmdict> --W --HBITS --size-mb --corpus
"""

import argparse
import json
import resource
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "impl"))
sys.path.insert(0, str(_HERE.parent / "impl_cm"))
sys.path.insert(0, str(_HERE.parent / "data"))

import bench_util as bu  # noqa: E402


def build(model, W, HBITS):
    if model == "hashed":
        from hashed_predictor import HashedCMPredictor
        return HashedCMPredictor(W, HBITS=HBITS)
    if model == "dict":
        import rnr1
        return rnr1.NGramPredictor(W)          # exact-dict counting baseline
    if model == "cmdict":
        from cm_predictor import CMPredictor
        return CMPredictor(W)                   # exact-dict CM baseline
    raise ValueError(model)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("model", choices=["hashed", "dict", "cmdict"])
    ap.add_argument("--W", type=int, default=6)
    ap.add_argument("--HBITS", type=int, default=22)
    ap.add_argument("--size-mb", type=float, default=0.5)
    ap.add_argument("--corpus", default="text")
    args = ap.parse_args(argv)

    data = bu.load_slice(args.corpus, int(args.size_mb * bu.MB))
    p = build(args.model, args.W, args.HBITS)
    t0 = time.perf_counter()
    for x in data:
        p.dist()
        p.update(x)
    wall = time.perf_counter() - t0
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    table_mb = (p.table_bytes() / 1e6) if hasattr(p, "table_bytes") else None
    print(json.dumps({
        "model": args.model, "W": args.W, "HBITS": args.HBITS,
        "size_mb": args.size_mb, "corpus": args.corpus,
        "peak_rss_mb": peak, "wall_s": wall, "table_mb": table_mb,
    }), flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""bench_rate.py -- rate (bpb) vs W, HBITS, and passes for the hashed CM
predictor, with external-compressor references on the SAME slice.

Rate is measured as model cross-entropy (ideal code length of the
predictor's integer distribution); the reference coder's redundancy over
this is < 1e-4 bpb (impl_cm gate G5), so archive bpb == CE bpb to 4 places.
Multi-pass configs report BOTH the coded-pass CE and the effective bpb that
includes the shipped weight blob (negligible at scale, dominant on 1 MB).

Writes results/rate.csv incrementally (one row per config) so partial runs
are still usable.  Slow by design (maximum-rate, not speed) -- run with -u.

    python -u bench_rate.py [--size-mb 1] [--corpora text,code]
"""

import argparse
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "impl"))
sys.path.insert(0, str(_HERE.parent / "data"))

import bench_util as bu  # noqa: E402
from hashed_predictor import HashedCMPredictor  # noqa: E402
from hashed_coder import train_weights  # noqa: E402

FIELDS = ["corpus", "size_mb", "W", "HBITS", "passes", "K",
          "coded_ce_bpb", "weight_overhead_bpb", "effective_bpb",
          "table_mb", "train_time_s", "code_ce_time_s", "per_pass_time_s"]


def config_ce(data, W, HBITS, passes):
    n = len(data)
    per_pass, train_s = [], 0.0
    if passes <= 1:
        p = HashedCMPredictor(W, HBITS=HBITS)
        w_overhead = 0.0
    else:
        w, per_pass = train_weights(data, W, HBITS, passes)
        train_s = sum(per_pass)
        p = HashedCMPredictor(W, HBITS=HBITS)
        p.load_weights(w)
        p.freeze = True
        p.reset_counts()
        w_overhead = 8.0 * (w.size * 4) / n     # int32 weight blob shipped
    t0 = time.perf_counter()
    coded = bu.model_ce_bpb(p, data)
    code_s = time.perf_counter() - t0
    return {
        "coded_ce_bpb": coded, "weight_overhead_bpb": w_overhead,
        "effective_bpb": coded + w_overhead,
        "table_mb": p.table_bytes() / 1e6,
        "train_time_s": train_s, "code_ce_time_s": code_s,
        "per_pass_time_s": ";".join("%.1f" % t for t in per_pass),
    }


def _load_done(out):
    """Resume support: return existing rows + a set of done config keys."""
    rows, done = [], set()
    if out.exists():
        import csv
        with open(out) as f:
            for r in csv.DictReader(f):
                rows.append(r)
                done.add((r["corpus"], float(r["size_mb"]), int(r["W"]),
                          int(r["HBITS"]), int(r["passes"])))
    return rows, done


def sweep(size_mb, corpora):
    size = int(size_mb * bu.MB)
    out = _HERE / "results" / "rate.csv"
    rows, done = _load_done(out)
    # phase A: W plateau (HBITS 22, single pass); phase B: HBITS plateau
    # (W 16); phase C: multi-pass benefit (W 16, HBITS 22).
    # Kept per-config short (<~6 min) so a periodic python reaper cannot
    # starve progress: this file is resumable (skips cached configs) and is
    # driven by a shell retry loop.  HBITS=24 and multi-pass at 1 MB are too
    # long to survive; they are characterized separately at a smaller size.
    plan = {
        "text": ([("A", W, 22, 1) for W in (4, 8, 12, 16, 20)] +
                 [("B", 16, HB, 1) for HB in (18, 20)]),
        "code": ([("A", W, 22, 1) for W in (4, 8, 12, 16)]),
    }
    for corpus in corpora:
        data = bu.load_slice(corpus, size)
        for _phase, W, HB, P in plan.get(corpus, []):
            if (corpus, float(size_mb), W, HB, P) in done:
                print("[%s %2.0fMB] W=%2d HBITS=%2d P=%d  (cached)"
                      % (corpus, size_mb, W, HB, P), flush=True)
                continue
            t0 = time.perf_counter()
            res = config_ce(data, W, HB, P)
            res.update({"corpus": corpus, "size_mb": size_mb, "W": W,
                        "HBITS": HB, "passes": P, "K": len(data)})
            rows.append(res)
            bu.csv_write(out, rows, FIELDS)
            print("[%s %2.0fMB] W=%2d HBITS=%2d P=%d  coded=%.4f eff=%.4f "
                  "bpb  table=%dMB  (%.0fs)"
                  % (corpus, size_mb, W, HB, P, res["coded_ce_bpb"],
                     res["effective_bpb"], res["table_mb"],
                     time.perf_counter() - t0), flush=True)
    # external references on the same slices
    ext_rows = []
    import tempfile
    for corpus in corpora:
        data = bu.load_slice(corpus, size)
        with tempfile.TemporaryDirectory() as td:
            ip = Path(td) / "in.bin"
            ip.write_bytes(data)
            for name in ("gzip-9", "bz2-9", "xz-6", "xz-9e", "zstd-19",
                         "zstd-22u", "brotli-11"):
                r = bu.run_external(name, ip, td)
                if r is None:
                    continue
                ext_rows.append({"corpus": corpus, "tool": name,
                                 "bpb": r["bpb"], "enc_mb_s": r["enc_mb_s"],
                                 "dec_mb_s": r["dec_mb_s"]})
                print("[%s %2.0fMB] %-10s bpb=%.4f" % (corpus, size_mb, name,
                                                       r["bpb"]), flush=True)
    bu.csv_write(_HERE / "results" / "rate_external.csv", ext_rows,
                 ["corpus", "tool", "bpb", "enc_mb_s", "dec_mb_s"])
    print("wrote", out, "and rate_external.csv", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--size-mb", type=float, default=1.0)
    ap.add_argument("--corpora", default="text,code")
    args = ap.parse_args(argv)
    sweep(args.size_mb, args.corpora.split(","))


if __name__ == "__main__":
    main()

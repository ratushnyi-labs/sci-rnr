#!/usr/bin/env python3
"""hashed_worker.py -- single-config subprocess worker for the hashed coder.

Each invocation performs exactly ONE operation (encode or decode) in a fresh
process and prints one JSON line of metrics -- including peak RSS from
resource.getrusage(RUSAGE_SELF).ru_maxrss (bytes on macOS) -- so the driver
gets a clean, isolated peak for encode vs decode and keeps the host to one
process at a time.

    python -u hashed_worker.py enc <in> <archive> --W --HBITS --passes --K
    python -u hashed_worker.py dec <archive> <out> [--orig <in>]
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

from hashed_coder import pack_hashed_mp, unpack_hashed_mp  # noqa: E402


def peak_rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="op", required=True)
    e = sub.add_parser("enc")
    e.add_argument("infile"); e.add_argument("archive")
    e.add_argument("--W", type=int, default=12)
    e.add_argument("--HBITS", type=int, default=22)
    e.add_argument("--passes", type=int, default=1)
    e.add_argument("--K", type=int, default=0, help="0 = warm whole-block")
    d = sub.add_parser("dec")
    d.add_argument("archive"); d.add_argument("outfile")
    d.add_argument("--orig", default=None)
    args = ap.parse_args(argv)

    if args.op == "enc":
        data = Path(args.infile).read_bytes()
        K = None if args.K <= 0 else args.K
        t0 = time.perf_counter()
        raw, meta = pack_hashed_mp(data, W=args.W, HBITS=args.HBITS,
                                   passes=args.passes, K=K,
                                   log=lambda s: print(s, file=sys.stderr,
                                                       flush=True))
        wall = time.perf_counter() - t0
        Path(args.archive).write_bytes(raw)
        n = len(data)
        out = {
            "op": "enc", "W": args.W, "HBITS": args.HBITS,
            "passes": args.passes, "K": (K if K else n),
            "bytes": n, "comp_bytes": len(raw),
            "bpb": 8.0 * len(raw) / max(1, n),
            "enc_wall_s": wall, "enc_mb_s": n / 1e6 / max(wall, 1e-9),
            "train_time_s": meta["train_time_s"],
            "encode_time_s": meta["encode_time_s"],
            "per_pass_time_s": meta["per_pass_time_s"],
            "n_passes": meta["n_passes"],
            "enc_peak_rss_mb": peak_rss_mb(),
        }
    else:
        raw = Path(args.archive).read_bytes()
        t0 = time.perf_counter()
        data = unpack_hashed_mp(raw, verify=True)
        wall = time.perf_counter() - t0
        Path(args.outfile).write_bytes(data)
        ok = True
        if args.orig:
            ok = data == Path(args.orig).read_bytes()
        n = len(data)
        out = {
            "op": "dec", "bytes": n,
            "dec_wall_s": wall, "dec_mb_s": n / 1e6 / max(wall, 1e-9),
            "dec_peak_rss_mb": peak_rss_mb(), "ok": bool(ok),
        }
    print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()

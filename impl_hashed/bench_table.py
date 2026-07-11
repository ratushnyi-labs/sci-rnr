#!/usr/bin/env python3
"""bench_table.py -- feasible-size compression table (first 1..64 MB slice).

Columns per (corpus, size):
  externals   gzip-9, bz2-9, xz-6, xz-9e, zstd-19, zstd-22u, brotli-11
  ours        rnr1-count-W3   (fast C engine impl_fast, exact-dict n-gram)
              rnr1-hashed-CM  (this work, best (W,HBITS,passes))

Every row records the full schema: bpb, enc/dec MB/s, enc/dec peak RSS,
plus for the hashed coder W/HBITS/passes/K, train_time_s, per_pass_time_s,
and a `measured` flag.  Sizes 1/4/16/64 MB are the feasible single-process
slice; 256/500/1500 MB rows are PROJECTED (time/bpb extrapolated from the
measured points; peak RSS is the SAME fixed table -- that is the point) and
need the C port to be run for real (next stage).

    python -u bench_table.py [--sizes 1,4] [--corpora text,code] \
        [--W 16 --HBITS 22 --passes 1] [--proj 16,64,256,500,1500]
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
PY = sys.executable
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "data"))
import bench_util as bu  # noqa: E402

FIELDS = ["corpus", "size_mb", "coder", "measured", "bpb",
          "enc_mb_s", "dec_mb_s", "enc_peak_rss_mb", "dec_peak_rss_mb",
          "W", "HBITS", "passes", "K", "train_time_s", "per_pass_time_s",
          "note"]

FAST_ENGINE = _HERE.parent / "impl_fast" / "rnr1fast.py"


def _time_l(argv, out_path):
    """Run argv under /usr/bin/time -l, stdout->out_path; parse peak RSS."""
    t0 = time.perf_counter()
    with open(out_path, "wb") as fo:
        p = subprocess.run(["/usr/bin/time", "-l"] + argv,
                           stdout=fo, stderr=subprocess.PIPE)
    wall = time.perf_counter() - t0
    m = re.search(rb"(\d+)\s+maximum resident set size", p.stderr)
    peak = int(m.group(1)) / 1e6 if m else float("nan")
    return wall, peak, p.returncode


def run_fast_ngram(in_path, tmpdir, W=3, K=65536):
    """rnr1-count-W3 via the fast C engine (impl_fast).  Measured metrics."""
    if not FAST_ENGINE.exists():
        return None
    in_path = Path(in_path)
    n = in_path.stat().st_size
    arc = Path(tmpdir) / "ng.rnr"
    dec = Path(tmpdir) / "ng.dec"
    enc = _time_l([PY, "-u", str(FAST_ENGINE), "pack", str(in_path),
                   str(arc), "--W", str(W), "--K", str(K)], Path(tmpdir)/"e.log")
    if enc[2] != 0:
        return None
    dc = _time_l([PY, "-u", str(FAST_ENGINE), "unpack", str(arc), str(dec)],
                 Path(tmpdir) / "d.log")
    if dc[2] != 0:
        return None
    ok = dec.read_bytes() == in_path.read_bytes()
    cb = arc.stat().st_size
    return {"bpb": 8.0 * cb / max(1, n),
            "enc_mb_s": n / 1e6 / max(enc[0], 1e-9), "enc_peak_rss_mb": enc[1],
            "dec_mb_s": n / 1e6 / max(dc[0], 1e-9), "dec_peak_rss_mb": dc[1],
            "ok": ok, "W": W, "K": K}


def run_hashed(in_path, tmpdir, W, HBITS, passes):
    """rnr1-hashed-CM via the isolated worker (enc then dec).  Measured."""
    worker = str(_HERE / "hashed_worker.py")
    arc = Path(tmpdir) / "hc.rnr"
    dec = Path(tmpdir) / "hc.dec"
    e = subprocess.run([PY, "-u", worker, "enc", str(in_path), str(arc),
                        "--W", str(W), "--HBITS", str(HBITS),
                        "--passes", str(passes)],
                       capture_output=True, text=True, timeout=6 * 3600)
    if e.returncode != 0:
        raise RuntimeError("hashed enc failed: " + e.stderr.strip()[-400:])
    em = json.loads(e.stdout.strip().splitlines()[-1])
    d = subprocess.run([PY, "-u", worker, "dec", str(arc), str(dec),
                        "--orig", str(in_path)],
                       capture_output=True, text=True, timeout=6 * 3600)
    if d.returncode != 0:
        raise RuntimeError("hashed dec failed: " + d.stderr.strip()[-400:])
    dm = json.loads(d.stdout.strip().splitlines()[-1])
    return em, dm


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="1,4,16,64",
                    help="sizes (MB) for externals + fast C n-gram")
    ap.add_argument("--hashed-sizes", default="1,4",
                    help="sizes (MB) at which the slow hashed CM is MEASURED")
    ap.add_argument("--corpora", default="text,code")
    ap.add_argument("--W", type=int, default=16)
    ap.add_argument("--HBITS", type=int, default=22)
    ap.add_argument("--passes", type=int, default=1)
    ap.add_argument("--proj", default="16,64,256,500,1500",
                    help="hashed-CM projected-only sizes (MB)")
    ap.add_argument("--no-hashed", action="store_true")
    args = ap.parse_args(argv)

    sizes = [float(s) for s in args.sizes.split(",")]
    hashed_sizes = set(float(s) for s in args.hashed_sizes.split(",") if s)
    proj_sizes = [float(s) for s in args.proj.split(",") if s]
    corpora = args.corpora.split(",")
    out = _HERE / "results" / "table.csv"
    rows, done = [], set()
    if out.exists():
        import csv
        with open(out) as f:
            for r in csv.DictReader(f):
                rows.append(r)
                done.add((r["corpus"], float(r["size_mb"]), r["coder"],
                          str(r["measured"])))

    def is_done(corpus, size_mb, coder, measured=True):
        return (corpus, float(size_mb), coder, str(measured)) in done

    def flush():
        bu.csv_write(out, rows, FIELDS)

    for corpus in corpora:
        # keep per-corpus anchor for projection
        anchor = None
        for size_mb in sizes:
            data = bu.load_slice(corpus, int(size_mb * bu.MB))
            with tempfile.TemporaryDirectory() as td:
                ip = Path(td) / "in.bin"
                ip.write_bytes(data)
                # externals
                for name in ("gzip-9", "bz2-9", "xz-6", "xz-9e", "zstd-19",
                             "zstd-22u", "brotli-11"):
                    if is_done(corpus, size_mb, name):
                        continue
                    r = bu.run_external(name, ip, td)
                    if r is None:
                        continue
                    rows.append({"corpus": corpus, "size_mb": size_mb,
                                 "coder": name, "measured": True,
                                 "bpb": r["bpb"], "enc_mb_s": r["enc_mb_s"],
                                 "dec_mb_s": r["dec_mb_s"],
                                 "enc_peak_rss_mb": r["enc_peak_mb"],
                                 "dec_peak_rss_mb": r["dec_peak_mb"],
                                 "note": "" if r["ok"] else "ROUNDTRIP_FAIL"})
                    print("[%s %gMB] %-12s bpb=%.4f" % (corpus, size_mb, name,
                          r["bpb"]), flush=True)
                    flush()
                # rnr1 counting W3 (fast C engine)
                ng = None if is_done(corpus, size_mb, "rnr1-count-W3") \
                    else run_fast_ngram(ip, td)
                if ng:
                    rows.append({"corpus": corpus, "size_mb": size_mb,
                                 "coder": "rnr1-count-W3", "measured": True,
                                 "bpb": ng["bpb"], "enc_mb_s": ng["enc_mb_s"],
                                 "dec_mb_s": ng["dec_mb_s"],
                                 "enc_peak_rss_mb": ng["enc_peak_rss_mb"],
                                 "dec_peak_rss_mb": ng["dec_peak_rss_mb"],
                                 "W": 3, "K": 65536,
                                 "note": "" if ng["ok"] else "ROUNDTRIP_FAIL"})
                    print("[%s %gMB] %-12s bpb=%.4f (C engine)"
                          % (corpus, size_mb, "rnr1-count-W3", ng["bpb"]),
                          flush=True)
                    flush()
                # rnr1 hashed CM (this work) -- measured only at small sizes
                if (not args.no_hashed and size_mb in hashed_sizes
                        and not is_done(corpus, size_mb, "rnr1-hashed-CM")):
                    em, dm = run_hashed(ip, td, args.W, args.HBITS,
                                        args.passes)
                    row = {"corpus": corpus, "size_mb": size_mb,
                           "coder": "rnr1-hashed-CM", "measured": True,
                           "bpb": em["bpb"], "enc_mb_s": em["enc_mb_s"],
                           "dec_mb_s": dm["dec_mb_s"],
                           "enc_peak_rss_mb": em["enc_peak_rss_mb"],
                           "dec_peak_rss_mb": dm["dec_peak_rss_mb"],
                           "W": args.W, "HBITS": args.HBITS,
                           "passes": args.passes, "K": em["K"],
                           "train_time_s": em["train_time_s"],
                           "per_pass_time_s": ";".join(
                               "%.1f" % t for t in em["per_pass_time_s"]),
                           "note": "" if dm["ok"] else "ROUNDTRIP_FAIL"}
                    rows.append(row)
                    anchor = em  # remember throughput/bpb for projection
                    print("[%s %gMB] %-12s bpb=%.4f enc=%.4fMB/s dec=%.4fMB/s "
                          "encRSS=%.0fMB" % (corpus, size_mb,
                          "rnr1-hashed-CM", em["bpb"], em["enc_mb_s"],
                          dm["dec_mb_s"], em["enc_peak_rss_mb"]), flush=True)
                    flush()
        # recover anchor from cached rows if the hashed row was skipped
        if anchor is None and not args.no_hashed:
            best = None
            for r in rows:
                if (r["corpus"] == corpus and r["coder"] == "rnr1-hashed-CM"
                        and str(r["measured"]) == "True"):
                    if best is None or float(r["size_mb"]) >= float(best["size_mb"]):
                        best = r
            if best:
                anchor = {"enc_mb_s": float(best["enc_mb_s"]),
                          "bpb": float(best["bpb"]),
                          "enc_peak_rss_mb": float(best["enc_peak_rss_mb"])}
        # projected hashed rows for the larger sizes (C-port stage)
        if anchor and not args.no_hashed:
            enc_mb_s = anchor["enc_mb_s"]
            dec_ref = None
            # use measured dec throughput if present in rows
            for r in rows:
                if (r["corpus"] == corpus and r["coder"] == "rnr1-hashed-CM"
                        and r.get("dec_mb_s")):
                    dec_ref = r["dec_mb_s"]
            for size_mb in proj_sizes:
                if is_done(corpus, size_mb, "rnr1-hashed-CM", measured=False):
                    continue
                rows.append({
                    "corpus": corpus, "size_mb": size_mb,
                    "coder": "rnr1-hashed-CM", "measured": False,
                    "bpb": anchor["bpb"],
                    "enc_mb_s": enc_mb_s, "dec_mb_s": dec_ref,
                    "enc_peak_rss_mb": anchor["enc_peak_rss_mb"],
                    "dec_peak_rss_mb": anchor["enc_peak_rss_mb"],
                    "W": args.W, "HBITS": args.HBITS, "passes": args.passes,
                    "note": "PROJECTED: bpb~=largest measured (improves with "
                            "warm-up); time=size/enc_mb_s; peakRSS FIXED "
                            "(file-size-independent); needs C port"})
            flush()
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()

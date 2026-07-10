#!/usr/bin/env python3
"""Per-run generated test sets: the mass-repetition statistics engine.

DESIGN (differs from the fixed-prefix tier campaign):
  On EACH run a fresh set of test files is GENERATED according to the
  content types, and all tests are executed against that generated set.
  Generation is seeded by (MASTER_SEED, run_idx), so any run is exactly
  reproducible from its recorded run_idx, yet different runs sample
  different data.  Across N runs the recorded rates/times are therefore
  distributions over BOTH corpus draws and machine timing -- the design
  intended for very large N (e.g. 100000 script invocations).

  Per type, run r generates:
    white_noise    fresh SHAKE-256 stream, domain-separated by run seed
    text/video/images/archive_mix/precompressed
                   a seeded random contiguous WINDOW of the pinned
                   data/big master (window offset recorded; the master
                   is immutable and sha-pinned, so window = derivable)
  Every record carries: run_idx, generation seed, window offset, the
  generated file's sha256, and the archive's sha256 (repeat-determinism
  witness: re-running the same run_idx must reproduce both).

USAGE
  python run_testset.py --size-mb 64 --runs 10            # runs next 10
  python run_testset.py --size-mb 64 --run-idx 7 --verify-repeat
  python run_testset.py --size-mb 16 --types text,video --coders zstd-19
Appends to out/runset<size>/results.jsonl (journal-checkpointed per
(run,type,coder) cell) and re-exports out/runset<size>/results.csv.
NOTE: do not run concurrently with a tier campaign on the same host if
timing statistics matter (CPU contention).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import time
from pathlib import Path

import coders
import corpora

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
MASTER_SEED = 20260710
_CHUNK = 8 * 1024 * 1024

CSV_COLUMNS = [
    "run_idx", "content_type", "coder", "size_bytes", "gen_seed",
    "window_offset", "corpus_sha256", "compressed_bytes",
    "compressed_sha256", "bpb", "ratio", "enc_s", "dec_s",
    "enc_peak_rss_mb", "dec_peak_rss_mb", "roundtrip_ok", "timestamp",
]


def _seed_for(run_idx: int, ctype: str) -> int:
    h = hashlib.sha256(
        f"rnr-runset:{MASTER_SEED}:{run_idx}:{ctype}".encode()).digest()
    return int.from_bytes(h[:8], "big")


def _gen_noise(path: Path, nbytes: int, seed: int) -> None:
    left, i = nbytes, 0
    with open(path, "wb") as f:
        while left > 0:
            x = hashlib.shake_256(
                f"rnr-runset-noise:{seed}:chunk={i}".encode()
            ).digest(min(_CHUNK, left))
            f.write(x)
            left -= len(x)
            i += 1


def _gen_window(path: Path, nbytes: int, seed: int, ctype: str) -> int:
    ent = corpora._big_entries()[ctype]
    master = corpora.BIG_DIR / ent["file"]
    msize = master.stat().st_size
    if msize < nbytes:
        raise RuntimeError(f"{ctype} master smaller than requested size")
    span = msize - nbytes
    offset = (seed % (span + 1)) if span > 0 else 0
    left = nbytes
    with open(master, "rb") as src, open(path, "wb") as out:
        src.seek(offset)
        while left > 0:
            chunk = src.read(min(_CHUNK, left))
            if not chunk:
                raise RuntimeError("master shorter than expected")
            out.write(chunk)
            left -= len(chunk)
    return offset


def generate(run_idx: int, ctype: str, nbytes: int, workdir: Path):
    """Generate this run's file for ctype; return (path, seed, offset, sha)."""
    seed = _seed_for(run_idx, ctype)
    path = workdir / f"r{run_idx}_{ctype}.bin"
    if ctype == "white_noise":
        _gen_noise(path, nbytes, seed)
        offset = None
    else:
        offset = _gen_window(path, nbytes, seed, ctype)
    return path, seed, offset, corpora._sha256_file(path)


def load_journal(p: Path) -> set:
    done = set()
    if p.exists():
        with open(p) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2 and parts[1] == "done":
                    done.add(parts[0])
    return done


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--size-mb", type=int, default=64)
    ap.add_argument("--types", default=",".join(corpora.CONTENT_TYPES))
    ap.add_argument("--coders", default=None,
                    help="default: every registered lane except rnr1-ref")
    ap.add_argument("--runs", type=int, default=1,
                    help="how many NEW runs to execute (run_idx continues "
                         "from the journal)")
    ap.add_argument("--run-idx", type=int, default=None,
                    help="execute exactly this run_idx (for reproduction)")
    ap.add_argument("--verify-repeat", action="store_true",
                    help="with --run-idx: assert regenerated files and "
                         "archives match previously recorded sha256s")
    ap.add_argument("--cap-min", type=float, default=30.0)
    args = ap.parse_args(argv)

    nbytes = args.size_mb * 1_000_000
    types = [t for t in args.types.split(",") if t]
    lanes = (args.coders.split(",") if args.coders else
             [c for c in coders.REGISTRY if c != "rnr1-ref"])
    out_dir = OUT / f"runset{args.size_mb}"
    out_dir.mkdir(parents=True, exist_ok=True)
    journal_p = out_dir / "cells.journal"
    results_p = out_dir / "results.jsonl"
    work = out_dir / "work"
    work.mkdir(exist_ok=True)
    done = load_journal(journal_p)

    prior = {}
    if args.verify_repeat and results_p.exists():
        with open(results_p) as f:
            for line in f:
                r = json.loads(line)
                prior[(r["run_idx"], r["content_type"], r["coder"])] = r

    if args.run_idx is not None:
        run_indices = [args.run_idx]
    else:
        used = set()
        if results_p.exists():
            with open(results_p) as f:
                used = {json.loads(l)["run_idx"] for l in f if l.strip()}
        start = (max(used) + 1) if used else 0
        run_indices = list(range(start, start + args.runs))

    for ri in run_indices:
        for ctype in types:
            if ctype != "white_noise" and not corpora.available(
                    "500", ctype):
                continue
            path, seed, offset, csha = generate(ri, ctype, nbytes, work)
            try:
                for lane in lanes:
                    cid = f"r{ri}/{ctype}/{lane}"
                    if cid in done and not args.verify_repeat:
                        continue
                    m = coders.run_cell(lane, path, csha, work,
                                        cap_s=args.cap_min * 60.0)
                    rec = {
                        "run_idx": ri, "content_type": ctype,
                        "coder": lane, "size_bytes": nbytes,
                        "gen_seed": seed, "window_offset": offset,
                        "corpus_sha256": csha, "measured": m,
                        "timestamp": time.strftime(
                            "%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    }
                    if args.verify_repeat:
                        p = prior.get((ri, ctype, lane))
                        if p:
                            same_c = p["corpus_sha256"] == csha
                            same_a = (p["measured"].get("compressed_sha256")
                                      == m.get("compressed_sha256"))
                            status = "REPRO-OK" if (same_c and same_a) \
                                else "REPRO-FAIL"
                            print(f"[runset] {cid}: corpus "
                                  f"{'==' if same_c else '!='} archive "
                                  f"{'==' if same_a else '!='} -> {status}")
                            if not (same_c and same_a):
                                return 1
                            continue
                    with open(results_p, "a") as f:
                        f.write(json.dumps(rec, sort_keys=True) + "\n")
                    with open(journal_p, "a") as f:
                        f.write(f"{cid}\tdone\n")
                    done.add(cid)
                    print(f"[runset] {cid}: bpb={m['bpb']:.4f} "
                          f"enc={m['enc_s']:.2f}s")
            finally:
                path.unlink(missing_ok=True)

    # CSV export
    rows = []
    with open(results_p) as f:
        for line in f:
            r = json.loads(line)
            m = r["measured"]
            rows.append({
                "run_idx": r["run_idx"], "content_type": r["content_type"],
                "coder": r["coder"], "size_bytes": r["size_bytes"],
                "gen_seed": r["gen_seed"],
                "window_offset": r["window_offset"],
                "corpus_sha256": r["corpus_sha256"],
                "compressed_bytes": m["compressed_size"],
                "compressed_sha256": m.get("compressed_sha256"),
                "bpb": m["bpb"], "ratio": m["ratio"],
                "enc_s": m["enc_s"], "dec_s": m["dec_s"],
                "enc_peak_rss_mb": (round(m["enc_peak_rss"] / 1e6, 1)
                                    if m.get("enc_peak_rss") else None),
                "dec_peak_rss_mb": (round(m["dec_peak_rss"] / 1e6, 1)
                                    if m.get("dec_peak_rss") else None),
                "roundtrip_ok": m["roundtrip_ok"],
                "timestamp": r["timestamp"],
            })
    rows.sort(key=lambda r: (r["run_idx"], r["content_type"], r["coder"]))
    with open(out_dir / "results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"[runset] {len(rows)} rows -> {out_dir/'results.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Deferred 16-64 MB smoke campaign for the two-pass coding mode.

NOT run by run_checks.py -- this is the measurement the gate defers.
Grid (defaults; all overridable):

    content type   payload (byte-prefix)                tier caps
    text           data/payloads/enwik8                 16 / 64 MB
    code           data/payloads/rnr_scripts_src.tar    capped at 3.4 MB
    structured     data/payloads/sqlite_synth.db        capped at 12.4 MB

  x  K in {16 KiB, 64 KiB, 256 KiB}
  x  coder in {rnr1-single, rnr1-2pass}

Tiers are decimal MB (10^6 bytes) prefixes of the payload, matching the
convention of experiments/type1_scale.  Where the payload is shorter
than the tier the ACTUAL byte count is recorded (size_bytes column) and
a "capped" marker is set -- the repo currently has no 16-64 MB code or
structured corpus, which is an honest limitation of this campaign, not
hidden by repetition tricks.

Output: one JSONL record per completed cell (append-only, resumable)
plus a flat CSV in the style of experiments/type1_scale/to_csv.py.

Costs (measured on this host at 2 MiB, run_checks T6): PASS-1 training
~1 MB/s, PASS-2 frozen encode ~0.5 MB/s, frozen decode ~0.4 MB/s, all
pure Python single-core; the single-pass side uses the byte-identity-
gated C engine (impl_fast) when built (~7-9 MB/s), else the pure
reference (~0.06 MB/s -- the default grid then takes DAYS; build the C
engine first).  With the C engine the full default grid is roughly
30-60 minutes single-core.  --no-decode replaces the full two-pass
decode with 32 verified random-access spot reads per cell if wall
clock matters more than the full round-trip column.

Memory caveat: PASS-1 keeps sparse Python dict tables for every
distinct context of orders 0..W.  At the 64 MB text tier with W=3 that
is millions of contexts, i.e. an expected 2-4 GB peak RSS.  Use
--tiers 16 (or a smaller W) on memory-constrained hosts.

Usage:
    measure_smoke.py --plan          # print the grid and estimates, exit
    measure_smoke.py --run [...]     # actually measure
Options:
    --types text,code,structured  --tiers 16,64  --Ks 16384,65536,262144
    --W 3  --out-dir impl_2pass/out  --no-decode  --redo
"""

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "impl"))

import two_pass as tp  # noqa: E402
import rnr1  # noqa: E402

MB = 10 ** 6

PAYLOADS = {
    "text": os.path.join(REPO, "data", "payloads", "enwik8"),
    "code": os.path.join(REPO, "data", "payloads", "rnr_scripts_src.tar"),
    "structured": os.path.join(REPO, "data", "payloads", "sqlite_synth.db"),
}

COLUMNS = [
    "tier", "size_bytes", "size_mb", "content_type", "coder", "family",
    "coder_version", "run_idx", "compressed_bytes", "bpb", "ratio",
    "enc_s", "dec_s", "enc_mb_s", "dec_mb_s",
    "roundtrip_ok", "censored", "status", "threads",
    "K", "W", "mode_chosen", "model_bytes", "model_bpb_share",
    "train_s", "capped", "corpus_sha256", "seed", "timestamp",
]

CODER_VERSIONS = {
    "rnr1-single": "1.0.0",
    "rnr1-2pass": "1.1.0+prune-v1",
}


def load_fast():
    """Byte-identity-gated fast single-pass engine, if built."""
    try:
        sys.path.insert(0, os.path.join(REPO, "impl_fast"))
        import rnr1fast
        # spot identity check before trusting it (cheap, 64 KiB)
        probe = open(PAYLOADS["text"], "rb").read(65536) \
            if os.path.exists(PAYLOADS["text"]) else b"probe" * 13107
        if rnr1fast.pack(probe, W=3, K=16384) != rnr1.pack(probe, W=3, K=16384):
            print("[smoke] WARNING: fast engine NOT byte-identical on probe; "
                  "using pure reference", flush=True)
            return None
        return rnr1fast
    except Exception as e:
        print("[smoke] fast engine unavailable (%s); pure reference used "
              "for the single-pass side (much slower)" % e, flush=True)
        return None


def sha256_hex(data):
    return hashlib.sha256(data).hexdigest()


def utcnow():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def cell_key(r):
    return (r["content_type"], r["tier"], r["K"], r["coder"])


def spot_reads(arc, data, rng_seed=20260711, n_reads=32):
    """Verified random-access spot reads (used with --no-decode)."""
    import random
    rng = random.Random(rng_seed)
    for _ in range(n_reads):
        p = rng.randrange(len(data))
        k = rng.randrange(1, 16384)
        got, _ = arc.read(p, k)
        if got != data[p: p + k]:
            return False
    return True


def run_cell(ctype, tier, K, coder, data, capped, sha, W, fast, frozen,
             train_s, no_decode):
    n = len(data)
    rec = {
        "tier": tier, "size_bytes": n, "size_mb": round(n / 1e6, 3),
        "content_type": ctype, "coder": coder,
        "family": ("c-fast" if (coder == "rnr1-single" and fast) else "py-ref"),
        "coder_version": CODER_VERSIONS[coder], "run_idx": 0,
        "threads": 1, "censored": False, "status": "done",
        "K": K, "W": W, "capped": capped, "corpus_sha256": sha,
        "seed": "", "timestamp": utcnow(),
        "mode_chosen": None, "model_bytes": None, "model_bpb_share": None,
        "train_s": None,
    }
    if coder == "rnr1-single":
        t0 = time.time()
        raw = fast.pack(data, W=W, K=K) if fast else rnr1.pack(data, W=W, K=K)
        rec["enc_s"] = round(time.time() - t0, 3)
        t0 = time.time()
        out = fast.unpack(raw) if fast else rnr1.unpack(raw)
        rec["dec_s"] = round(time.time() - t0, 3)
        rec["roundtrip_ok"] = out == data
    else:
        sraw = fast.pack(data, W=W, K=K) if fast else rnr1.pack(data, W=W, K=K)
        t0 = time.time()
        raw, info = tp.pack2(data, W=W, K=K, single_raw=sraw, frozen=frozen,
                             return_info=True)
        rec["enc_s"] = round(time.time() - t0, 3)  # PASS 2 (+ compare) only
        rec["train_s"] = round(train_s, 3)         # PASS 1, shared across K
        rec["mode_chosen"] = info["mode"]
        rec["model_bytes"] = info["model_bytes"]
        rec["model_bpb_share"] = (round(8.0 * info["model_bytes"] / n, 6)
                                  if info["model_bytes"] else None)
        t0 = time.time()
        if no_decode:
            arc = tp.open_archive(raw)
            rec["roundtrip_ok"] = spot_reads(arc, data)
            rec["dec_s"] = None
            rec["status"] = "done-spot-reads"
        else:
            out = tp.unpack2(raw)
            rec["dec_s"] = round(time.time() - t0, 3)
            rec["roundtrip_ok"] = out == data
    rec["compressed_bytes"] = len(raw)
    rec["bpb"] = round(8.0 * len(raw) / n, 6)
    rec["ratio"] = round(n / len(raw), 4)
    rec["enc_mb_s"] = (round(n / 1e6 / rec["enc_s"], 3)
                       if rec.get("enc_s") else None)
    rec["dec_mb_s"] = (round(n / 1e6 / rec["dec_s"], 3)
                       if rec.get("dec_s") else None)
    return rec


def write_csv(records, path):
    import csv
    records = sorted(records, key=lambda r: (str(r["content_type"]),
                                             int(r["tier"]), int(r["K"]),
                                             str(r["coder"])))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", action="store_true",
                    help="actually measure (otherwise only the plan prints)")
    ap.add_argument("--plan", action="store_true", help="print plan and exit")
    ap.add_argument("--types", default="text,code,structured")
    ap.add_argument("--tiers", default="16,64", help="decimal MB prefixes")
    ap.add_argument("--Ks", default="16384,65536,262144")
    ap.add_argument("--W", type=int, default=3)
    ap.add_argument("--out-dir", default=os.path.join(HERE, "out"))
    ap.add_argument("--no-decode", action="store_true",
                    help="two-pass cells: verified spot reads instead of a "
                         "full sequential decode (faster, weaker evidence)")
    ap.add_argument("--redo", action="store_true",
                    help="re-measure cells already present in results.jsonl")
    args = ap.parse_args(argv)

    types = [t.strip() for t in args.types.split(",") if t.strip()]
    tiers = [int(t) for t in args.tiers.split(",")]
    Ks = [int(k) for k in args.Ks.split(",")]

    plan = []
    for ctype in types:
        path = PAYLOADS.get(ctype)
        if path is None or not os.path.exists(path):
            print("[smoke] SKIP type %r: payload missing (%s)" % (ctype, path))
            continue
        fsize = os.path.getsize(path)
        for tier in tiers:
            want = tier * MB
            actual = min(want, fsize)
            capped = actual < want
            if capped and (tier != tiers[0]) and actual <= (tiers[0] * MB):
                # a fully-capped larger tier would duplicate the smaller one
                print("[smoke] SKIP %s tier %d MB: payload only %.1f MB "
                      "(already covered by the smaller tier)"
                      % (ctype, tier, fsize / 1e6))
                continue
            for K in Ks:
                for coder in ("rnr1-single", "rnr1-2pass"):
                    plan.append((ctype, tier, actual, capped, K, coder))

    print("[smoke] plan: %d cells" % len(plan))
    for ctype, tier, actual, capped, K, coder in plan:
        print("    %-11s tier=%-3d MB actual=%8.1f MB%s K=%-7d %s"
              % (ctype, tier, actual / 1e6, " (capped)" if capped else "",
                 K, coder))
    two_mb = sum(a for _, _, a, _, _, c in plan if c == "rnr1-2pass") / 1e6
    train_mb = sum({(t, a): a for t, _, a, _, _, c in plan
                    if c == "rnr1-2pass"}.values()) / 1e6
    est_min = (train_mb / 1.0 + two_mb / 0.5
               + (0 if args.no_decode else two_mb / 0.4)) / 60
    print("[smoke] rough pure-Python cost at measured rates: ~%.0f min "
          "(train %.0f MB + encode%s %.0f MB two-pass); single-pass side "
          "adds minutes with the C engine built, DAYS without it"
          % (est_min, train_mb,
             "" if args.no_decode else "+decode", two_mb))
    if not args.run:
        print("[smoke] dry run only -- pass --run to measure "
              "(deliberately not run by run_checks.py)")
        return 0

    os.makedirs(args.out_dir, exist_ok=True)
    jsonl_path = os.path.join(args.out_dir, "smoke_results.jsonl")
    csv_path = os.path.join(args.out_dir, "smoke_results.csv")
    done = {}
    if os.path.exists(jsonl_path) and not args.redo:
        with open(jsonl_path) as f:
            for line in f:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    done[cell_key(r)] = r

    fast = load_fast()
    records = list(done.values())
    cache = {}   # (ctype, actual) -> (data, sha, frozen, train_s)
    for ctype, tier, actual, capped, K, coder in plan:
        key = (ctype, tier, K, coder)
        if key in done:
            print("[smoke] skip (done): %s" % (key,))
            continue
        ck = (ctype, actual)
        if ck not in cache:
            with open(PAYLOADS[ctype], "rb") as f:
                data = f.read(actual)
            sha = sha256_hex(data)
            t0 = time.time()
            model, blob, _ = tp.build_frozen_model(data, args.W)
            train_s = time.time() - t0
            cache[ck] = (data, sha, (model, blob), train_s)
            # keep at most one corpus in memory
            for old in [c for c in cache if c != ck]:
                del cache[old]
        data, sha, frozen, train_s = cache[ck]
        print("[smoke] cell %s ..." % (key,), flush=True)
        t0 = time.time()
        rec = run_cell(ctype, tier, K, coder, data, capped, sha, args.W,
                       fast, frozen, train_s, args.no_decode)
        print("[smoke]   -> %.4f bpb (%s) in %.1f s"
              % (rec["bpb"], rec.get("mode_chosen") or coder,
                 time.time() - t0), flush=True)
        with open(jsonl_path, "a") as f:
            f.write(json.dumps(rec, sort_keys=True) + "\n")
        records.append(rec)
        write_csv(records, csv_path)   # rewrite after every cell (resumable)
    print("[smoke] %d records -> %s" % (len(records), csv_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())

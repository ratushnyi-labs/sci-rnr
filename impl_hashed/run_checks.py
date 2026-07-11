#!/usr/bin/env python3
"""run_checks.py -- PASS/FAIL gate for the impl_hashed deliverable
(memory-bounded hashed context-mixing predictor for the RNR Type-I coder).

Gates (tiny inputs; a few minutes total)
----------------------------------------
  G0  predictor self-test        splitmix determinism + learns a pattern
  G1  round-trip bit-exact       hashed pack -> unpack == input on text /
                                 mixed / incompressible (store-mode) slices,
                                 W beyond the reference cap, plain + multipass,
                                 + random-access (reset-per-sync) variant
  G2  two-process determinism    the worker packs the same input to the
                                 byte-identical archive in two fresh
                                 processes and matches the in-process archive
  G3  memory-bound property      peak RSS at W=16 within ~1.5x of W=3 at the
                                 same (small) HBITS -- the tables do not blow
                                 up as the order grows
  G4  never-worse-than-counting  hashed-CM model cross-entropy <= the exact
                                 n-gram baseline's at the same W on text

Exit 0 iff every gate PASSes.   python -u run_checks.py
"""

import hashlib
import json
import math
import subprocess
import sys
import tempfile
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
PY = sys.executable
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "impl"))
sys.path.insert(0, str(_HERE.parent / "impl_cm"))
sys.path.insert(0, str(_HERE.parent / "data"))

import loader  # noqa: E402
import rnr1  # noqa: E402
import hashed_predictor as hp  # noqa: E402
from hashed_coder import (pack_hashed, unpack_hashed, read_hashed,  # noqa: E402
                          pack_hashed_mp, unpack_hashed_mp)

RESULTS = []


def gate(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", name,
                         (" -- " + detail) if detail else ""), flush=True)
    return bool(ok)


def g0_selftest():
    try:
        bpb = hp._self_test()
        return gate("G0 predictor self-test", bpb < 3.0,
                    "repeating-text %.3f bpb" % bpb)
    except Exception as exc:  # noqa: BLE001
        return gate("G0 predictor self-test", False, str(exc)[:200])


def g1_roundtrip():
    text = loader.read_range("enwik8", 1 << 20, 48 * 1024)
    mixed = loader.read_range("calgary", 100 * 1024, 32 * 1024)
    rand = loader.read_range("ctrl_urandom", 0, 16 * 1024)
    ok = True
    for label, data in (("text", text), ("mixed", mixed),
                        ("incompressible", rand)):
        for W, HB in ((6, 18), (12, 20), (16, 18)):
            raw = pack_hashed(data, W=W, HBITS=HB)
            good = unpack_hashed(raw, HBITS=HB) == data
            ok &= good
            if not good:
                print("    MISMATCH %s W=%d HB=%d" % (label, W, HB), flush=True)
    # store-mode fail-safe engages on incompressible input (K forces blocks)
    raw = pack_hashed(rand, W=6, HBITS=16, K=4096)
    ok &= len(raw) <= len(rand) + 4096
    # multipass, decodable
    for P in (2, 3):
        raw = pack_hashed_mp(text, W=8, HBITS=18, passes=P)[0]
        ok &= unpack_hashed_mp(raw) == text
    # random-access reset-per-sync variant
    raw = pack_hashed(text, W=8, HBITS=18, K=8192)
    piece, st = read_hashed(raw, 20000, 500, HBITS=18)
    ra_ok = piece == text[20000:20500]
    ok &= ra_ok
    print("    RA warmup=%d decoded=%d ok=%s"
          % (st["warmup_bytes"], st["decoded_bytes"], ra_ok), flush=True)
    return gate("G1 round-trip bit-exact (+multipass, +RA, +store-mode)", ok)


def g2_determinism():
    data = loader.read_range("enwik8", 2 << 20, 24 * 1024)
    worker = str(_HERE / "hashed_worker.py")
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.bin"
        src.write_bytes(data)
        digests = []
        for i in (1, 2):
            arc = Path(td) / ("o%d.rnr" % i)
            p = subprocess.run([PY, "-u", worker, "enc", str(src), str(arc),
                                "--W", "10", "--HBITS", "18", "--passes", "1"],
                               capture_output=True, text=True, timeout=1200)
            if p.returncode != 0:
                return gate("G2 two-process determinism", False,
                            p.stderr.strip()[-200:])
            digests.append(hashlib.sha256(arc.read_bytes()).hexdigest())
        inproc = hashlib.sha256(pack_hashed(data, W=10, HBITS=18)).hexdigest()
        ok = digests[0] == digests[1] == inproc
        return gate("G2 two-process determinism", ok, "sha256 " + digests[0][:16])


def g3_memory_bound():
    probe = str(_HERE / "mem_probe.py")
    HB = 16
    peaks = {}
    for W in (3, 16):
        p = subprocess.run([PY, "-u", probe, "hashed", "--W", str(W),
                            "--HBITS", str(HB), "--size-mb", "0.1",
                            "--corpus", "text"],
                           capture_output=True, text=True, timeout=1200)
        if p.returncode != 0:
            return gate("G3 memory-bound (W16 <= 1.5x W3)", False,
                        p.stderr.strip()[-200:])
        peaks[W] = json.loads(p.stdout.strip().splitlines()[-1])["peak_rss_mb"]
    ratio = peaks[16] / peaks[3]
    return gate("G3 memory-bound (W16 <= 1.5x W3 @ HBITS=%d)" % HB,
                ratio <= 1.5,
                "peakRSS W3=%.1fMB W16=%.1fMB ratio=%.2fx"
                % (peaks[3], peaks[16], ratio))


def g4_beats_counting():
    data = loader.read_range("enwik8", 3 << 20, 96 * 1024)
    W = 4
    ph = hp.HashedCMPredictor(W, HBITS=20)
    bits_h = 0.0
    for x in data:
        f, t = ph.dist()
        bits_h += -math.log2(f[x] / t)
        ph.update(x)
    ce_h = bits_h / len(data)
    pn = rnr1.NGramPredictor(W)
    bits_n = 0.0
    for x in data:
        f, t = pn.dist()
        bits_n += -math.log2(f[x] / t)
        pn.update(x)
    ce_n = bits_n / len(data)
    return gate("G4 hashed-CM never worse than counting (text, W=%d)" % W,
                ce_h <= ce_n + 1e-9,
                "hashed CE=%.4f  n-gram CE=%.4f bpb" % (ce_h, ce_n))


def main():
    t0 = time.perf_counter()
    g0_selftest()
    g1_roundtrip()
    g2_determinism()
    g3_memory_bound()
    g4_beats_counting()
    ok = all(r for _, r in RESULTS)
    print("=" * 60, flush=True)
    for name, r in RESULTS:
        print("  %-54s %s" % (name, "PASS" if r else "FAIL"))
    print("OVERALL: %s (%.1f s)" % ("PASS" if ok else "FAIL",
                                    time.perf_counter() - t0), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Executable checks for the reference RNR Type-I coder (impl/rnr1.py).

Prints one PASS/FAIL line per check and an OVERALL line, in the style of
scripts/verify/.  Floating-point arithmetic appears only here (offline
analysis, R-3.6); the coder under test is integer-only.

Checks
------
C1  E12-A structure: injectivity, parity equations, minimum nonzero mask
    weight 3 (classes 1/2 structurally empty under byte-aligned prediction).
C2  Round-trip bit-exactness on generated smoke data (random, repetitive,
    zeros, empty, 1 byte, Markov) + tamper detection via H_v / sub-block
    hashes.
C3  Round-trip bit-exactness on a real file (tex/rnr_coding.tex).
C4  Determinism: two in-process encodes identical; two separate-process
    CLI encodes identical; decode-of-encode across two separate processes
    identical.
C5  Random access: 200 random (p, k) reads equal sequential-decode slices;
    warm-up discard <= K-1; decoded window <= K + k (Theorem 10.3);
    verified-read path; short reads at end of source.
C6  Rate sanity: on an order-2 Markov source with known entropy rate h,
    measured repair-stream bpb within 5% of h + predicted adaptive/sync
    overhead at K = 64 KiB; arithmetic-coder redundancy negligible.
C7  Sync overhead: measured bpb(K=64KiB) - bpb(no sync) matches the
    Monte-Carlo prediction and stays below 5% of h.
C8  Adversarial input: uniform-random input codes at <= 8.2 bpb total
    (fail-safe store mode, Theorem 7.22 semantics) and random access
    works on store-mode sub-blocks.

Run:  /Users/para/.venvs/rnr/bin/python impl/run_checks.py
"""

import math
import os
import random
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rnr1  # noqa: E402

import numpy as np  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PYTHON = sys.executable
RNR1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rnr1.py")

RESULTS = []


def record(name, ok, detail=""):
    line = "PASS" if ok else "FAIL"
    print("%s: %s%s" % (line, name, (" -- " + detail) if detail else ""))
    RESULTS.append((name, ok))
    return ok


# ---------------------------------------------------------------------------
# Order-2 Markov source with known entropy rate (4 symbols over bytes a-d)
# ---------------------------------------------------------------------------

ALPHABET = [97, 98, 99, 100]  # 'a'..'d'
BASE_ROW = [0.60, 0.22, 0.13, 0.05]


def markov_tensor(seed=7):
    """Order-2 transition tensor P[s1][s2] -> row over 4 next symbols."""
    rng = random.Random(seed)
    tensor = {}
    for s1 in range(4):
        for s2 in range(4):
            row = BASE_ROW[:]
            rng.shuffle(row)
            tensor[(s1, s2)] = row
    return tensor


def markov_entropy_rate(tensor):
    """Exact entropy rate in bits/symbol via the stationary pair chain."""
    # Pair chain: (s1, s2) -> (s2, s3) with prob tensor[(s1,s2)][s3].
    T = np.zeros((16, 16))
    for s1 in range(4):
        for s2 in range(4):
            for s3 in range(4):
                T[4 * s1 + s2, 4 * s2 + s3] = tensor[(s1, s2)][s3]
    pi = np.full(16, 1.0 / 16)
    for _ in range(20000):
        nxt = pi @ T
        if np.abs(nxt - pi).max() < 1e-15:
            pi = nxt
            break
        pi = nxt
    h = 0.0
    for s1 in range(4):
        for s2 in range(4):
            row = tensor[(s1, s2)]
            h += pi[4 * s1 + s2] * (-sum(p * math.log2(p) for p in row if p > 0))
    return h, pi


def markov_sample(tensor, n, seed):
    rng = random.Random(seed)
    s1, s2 = rng.randrange(4), rng.randrange(4)
    out = bytearray()
    symbols = list(range(4))
    for _ in range(n):
        row = tensor[(s1, s2)]
        s3 = rng.choices(symbols, weights=row)[0]
        out.append(ALPHABET[s3])
        s1, s2 = s2, s3
    return bytes(out)


def assignment_bpb(data, W, K):
    """Ideal (pre-arithmetic-coding) code length of the predictor's
    sequential probability assignment, in bits per byte, with sync-point
    resets every K bytes.  This is the theory-side cost model: the coder
    should match it up to arithmetic-coder redundancy."""
    total_bits = 0.0
    n = len(data)
    for start in range(0, n, K):
        block = data[start: start + K]
        pred = rnr1.NGramPredictor(W)
        for x in block:
            freq, total = pred.dist()
            total_bits += math.log2(total) - math.log2(int(freq[x]))
            pred.update(x)
    return total_bits / max(1, n)


def repair_bpb(archive_bytes, n):
    arc = rnr1.Archive(archive_bytes)
    return 8.0 * len(arc.repair) / max(1, n)


# ---------------------------------------------------------------------------
# C1: E12-A structure
# ---------------------------------------------------------------------------

def check_c1():
    ok = True
    vals = set(int(v) for v in rnr1.E12)
    ok &= len(vals) == 256  # injective
    # Recompute parity equations independently for a few bytes.
    for b in [0, 1, 5, 0x55, 0xAA, 0xFF, 200]:
        w = int(rnr1.E12[b])
        c = [0] + [(w >> i) & 1 for i in range(12)]  # 1-based
        bits = [(b >> j) & 1 for j in range(8)]
        ok &= [c[3], c[5], c[6], c[7], c[9], c[10], c[11], c[12]] == bits
        ok &= c[1] == c[3] ^ c[5] ^ c[7] ^ c[9] ^ c[11]
        ok &= c[2] == c[3] ^ c[6] ^ c[7] ^ c[10] ^ c[11]
        ok &= c[4] == c[5] ^ c[6] ^ c[7] ^ c[12]
        ok &= c[8] == c[9] ^ c[10] ^ c[11] ^ c[12]
    nz = rnr1.WT[rnr1.WT > 0]
    min_wt = int(nz.min())
    ok &= min_wt == 3  # linear code, min distance 3 -> classes 1/2 empty
    return record(
        "C1 E12-A injective, parity equations match Section 4.2, "
        "min nonzero mask weight = 3",
        ok,
        "classes 1/2 structurally empty as documented",
    )


# ---------------------------------------------------------------------------
# C2: round-trip on generated smoke data + tamper detection
# ---------------------------------------------------------------------------

def check_c2(markov_data):
    rng = random.Random(123)
    cases = {
        "random-8KiB": bytes(rng.randrange(256) for _ in range(8192)),
        "repetitive-8KiB": (b"the quick brown fox " * 500)[:8192],
        "zeros-4KiB": bytes(4096),
        "empty": b"",
        "one-byte": b"Q",
        "markov-64KiB": markov_data[:65536],
    }
    ok = True
    details = []
    for name, data in cases.items():
        raw = rnr1.pack(data, W=2, K=16384)
        out = rnr1.unpack(raw)
        good = out == data
        ok &= good
        details.append("%s %.3f bpb%s" % (name, 8.0 * len(raw) / max(1, len(data)),
                                          "" if good else " MISMATCH"))
    # Tamper detection: flip one byte in the repair stream area.
    data = cases["repetitive-8KiB"]
    raw = bytearray(rnr1.pack(data, W=2, K=16384))
    raw[-10] ^= 0x40
    try:
        rnr1.unpack(bytes(raw))
        ok = False
        details.append("tampered archive accepted (BAD)")
    except ValueError:
        details.append("tamper rejected")
    return record("C2 round-trip bit-exact on smoke data + tamper rejection",
                  ok, "; ".join(details))


# ---------------------------------------------------------------------------
# C3: round-trip on a real file
# ---------------------------------------------------------------------------

def check_c3():
    path = os.path.join(REPO, "tex", "rnr_coding.tex")
    with open(path, "rb") as f:
        data = f.read()
    t0 = time.time()
    raw = rnr1.pack(data, W=3, K=65536)
    t1 = time.time()
    out = rnr1.unpack(raw)
    t2 = time.time()
    ok = out == data
    bpb = 8.0 * len(raw) / len(data)
    return record(
        "C3 round-trip bit-exact on real file tex/rnr_coding.tex",
        ok,
        "%d bytes, total %.4f bpb (W=3, K=64KiB), enc %.1fs dec %.1fs"
        % (len(data), bpb, t1 - t0, t2 - t1),
    )


# ---------------------------------------------------------------------------
# C4: determinism (in-process and across separate processes)
# ---------------------------------------------------------------------------

def check_c4(markov_data):
    data = markov_data[:65536]
    ok = True
    details = []

    raw_a = rnr1.pack(data, W=2, K=16384)
    raw_b = rnr1.pack(data, W=2, K=16384)
    ok &= raw_a == raw_b
    details.append("in-process encode x2 %s" % ("identical" if raw_a == raw_b else "DIFFER"))

    with tempfile.TemporaryDirectory(prefix="rnr1_checks_") as td:
        src = os.path.join(td, "src.bin")
        with open(src, "wb") as f:
            f.write(data)
        arcs, outs = [], []
        for i in range(2):
            arc = os.path.join(td, "a%d.rnr" % i)
            out = os.path.join(td, "o%d.bin" % i)
            subprocess.run([PYTHON, RNR1, "pack", src, arc, "--W", "2",
                            "--K", "16384"], check=True, capture_output=True)
            subprocess.run([PYTHON, RNR1, "unpack", arc, out],
                           check=True, capture_output=True)
            with open(arc, "rb") as f:
                arcs.append(f.read())
            with open(out, "rb") as f:
                outs.append(f.read())
        ok &= arcs[0] == arcs[1]
        details.append("2-process encode %s" % ("identical" if arcs[0] == arcs[1] else "DIFFER"))
        ok &= arcs[0] == raw_a
        details.append("CLI == API archive %s" % ("yes" if arcs[0] == raw_a else "NO"))
        ok &= outs[0] == outs[1] == data
        details.append("2-process decode %s" %
                       ("identical & correct" if outs[0] == outs[1] == data else "DIFFER"))
    return record("C4 determinism: encode/decode byte-identical across "
                  "repeated runs and separate processes", ok, "; ".join(details))


# ---------------------------------------------------------------------------
# C5: random access equals sequential-decode slices (Theorem 10.3 pattern)
# ---------------------------------------------------------------------------

def check_c5(markov_data, markov_archive):
    K = 65536
    arc = rnr1.Archive(markov_archive)
    n = arc.n
    ok = True
    rng = random.Random(2026)
    max_warm = 0
    bad = 0
    for q in range(200):
        p = rng.randrange(n)
        # k spans 1 byte .. ~3 sub-blocks, log-uniform-ish
        k = min(1 << rng.randrange(0, 18), 3 * K)
        got, stats = arc.read(p, k)
        want = markov_data[p: p + k]
        if got != want:
            bad += 1
            ok = False
        if stats["warmup_bytes"] > K - 1:
            ok = False
        if stats["decoded_bytes"] > K + k:
            ok = False
        max_warm = max(max_warm, stats["warmup_bytes"])
    # Verified-read path (full covering sub-blocks, R-11.1.1 hashes)
    for q in range(5):
        p = rng.randrange(n - 4096)
        got, stats = arc.read(p, 4096, verify=True)
        ok &= got == markov_data[p: p + 4096]
        ok &= stats["decoded_bytes"] <= 2 * K + 4096
    # Short read at end + empty read at n
    got, _ = arc.read(n - 100, 5000)
    ok &= got == markov_data[n - 100:]
    got, _ = arc.read(n, 10)
    ok &= got == b""
    return record(
        "C5 random access: 200 random (p,k) reads == sequential slices; "
        "warm-up <= K-1; decode window <= K+k",
        ok,
        "mismatches=%d, max warm-up %d bytes (K=%d)" % (bad, max_warm, K),
    )


# ---------------------------------------------------------------------------
# C6/C7: rate sanity on order-2 Markov source, sync overhead prediction
# ---------------------------------------------------------------------------

def check_c6_c7(tensor, h, markov_data, markov_archive):
    W, K, n = 2, 65536, len(markov_data)

    bpb_meas = repair_bpb(markov_archive, n)

    # Arithmetic-coder redundancy: measured bits vs the assignment cost of
    # the very same data (difference = AC quantization + flush + padding).
    assign_actual = assignment_bpb(markov_data, W, K)
    ac_red = bpb_meas - assign_actual

    # Theory-side prediction, independent of the coded sample: Monte-Carlo
    # of the predictor's sequential probability assignment on fresh chains,
    # plus a flush allowance (<= 3 bytes per sub-block).
    m = (n + K - 1) // K
    reps = [assignment_bpb(markov_sample(tensor, n, seed=1000 + r), W, K)
            for r in range(2)]
    flush_allow = 8.0 * 3 * m / n
    pred_bpb = sum(reps) / len(reps) + flush_allow
    red_pred = pred_bpb - h

    ok6 = True
    ok6 &= abs(bpb_meas - (h + red_pred)) <= 0.05 * h
    ok6 &= 0 <= ac_red <= 0.01
    ok6 &= red_pred < 0.05 * h  # the band is not inflated by the overhead term
    record(
        "C6 rate sanity on order-2 Markov source (W=2, K=64KiB)",
        ok6,
        "h=%.4f bpb, predicted overhead=%.4f, target=%.4f, measured=%.4f "
        "(gap %.4f, band +/-%.4f); AC redundancy=%.5f bpb"
        % (h, red_pred, h + red_pred, bpb_meas,
           bpb_meas - (h + red_pred), 0.05 * h, ac_red),
    )

    # C7: sync overhead = bpb(K=64KiB) - bpb(no sync), vs prediction.
    raw_nosync = rnr1.pack(markov_data, W=W, K=n)
    bpb_nosync = repair_bpb(raw_nosync, n)
    over_meas = bpb_meas - bpb_nosync
    reps_ns = [assignment_bpb(markov_sample(tensor, n, seed=1000 + r), W, n)
               for r in range(2)]
    pred_nosync = sum(reps_ns) / len(reps_ns) + 8.0 * 3 / n
    over_pred = pred_bpb - pred_nosync
    ok7 = True
    ok7 &= rnr1.unpack(raw_nosync) == markov_data
    ok7 &= over_meas >= 0
    ok7 &= over_meas <= 0.05 * h
    ok7 &= abs(over_meas - over_pred) <= 0.005 + 0.5 * over_pred
    record(
        "C7 sync-point overhead at K=64KiB matches prediction and is < 5% of h",
        ok7,
        "measured=%.5f bpb, predicted=%.5f bpb (no-sync bpb=%.4f)"
        % (over_meas, over_pred, bpb_nosync),
    )
    return ok6 and ok7


# ---------------------------------------------------------------------------
# C8: adversarial floor (uniform random input near store mode)
# ---------------------------------------------------------------------------

def check_c8():
    rng = random.Random(99)
    data = bytes(rng.randrange(256) for _ in range(65536))
    raw = rnr1.pack(data, W=2, K=16384)
    ok = rnr1.unpack(raw) == data
    bpb = 8.0 * len(raw) / len(data)
    ok &= bpb <= 8.2
    arc = rnr1.Archive(raw)
    n_raw_mode = sum(1 for e in arc.entries if e[3] == rnr1.MODE_RAW)
    ok &= n_raw_mode == arc.m  # every sub-block should fall back to raw
    # Random access must work on store-mode sub-blocks too.
    for _ in range(20):
        p = rng.randrange(len(data))
        k = rng.randrange(1, 40000)
        got, stats = arc.read(p, k)
        ok &= got == data[p: p + k]
        ok &= stats["decoded_bytes"] <= 16384 + k
    return record("C8 adversarial input: store-mode fail-safe keeps "
                  "uniform-random input near 8 bpb; random access on raw "
                  "sub-blocks correct",
                  ok, "total %.4f bpb (<= 8.2), %d/%d sub-blocks raw"
                  % (bpb, n_raw_mode, arc.m))


# ---------------------------------------------------------------------------

def main():
    t_start = time.time()
    print("=" * 72)
    print("run_checks: reference RNR Type-I coder (impl/rnr1.py)")
    print("=" * 72)

    tensor = markov_tensor()
    h, _ = markov_entropy_rate(tensor)
    n = 8 * 65536  # 512 KiB -> 8 sub-blocks at K = 64 KiB
    markov_data = markov_sample(tensor, n, seed=555)
    markov_archive = rnr1.pack(markov_data, W=2, K=65536)

    check_c1()
    check_c2(markov_data)
    check_c3()
    check_c4(markov_data)
    check_c5(markov_data, markov_archive)
    check_c6_c7(tensor, h, markov_data, markov_archive)
    check_c8()

    n_pass = sum(1 for _, ok in RESULTS if ok)
    n_all = len(RESULTS)
    overall = n_pass == n_all
    print("-" * 72)
    print("OVERALL: %s (%d/%d checks passed, %.1f s)"
          % ("PASS" if overall else "FAIL", n_pass, n_all, time.time() - t_start))
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())

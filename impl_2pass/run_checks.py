#!/usr/bin/env python3
"""
Executable checks for the two-pass (train-then-freeze) coding mode
(impl_2pass/two_pass.py).  Tiny inputs only (<= 4 MB total per check);
the 16-64 MB campaign is deferred to measure_smoke.sh, which this gate
deliberately does NOT run.

Prints one PASS/FAIL line per check and an OVERALL line, in the style
of impl/run_checks.py.  Floats appear only here (offline analysis); the
coder under test is integer-only on the coding path (the encoder-side
prune estimate uses float log2, documented in FORMAT2.md).

Checks
------
T0  Fast-engine byte-identity spot check (impl_fast vs impl/rnr1.py) on
    text and noise; decides whether later checks may use the fast
    engine for single-pass baselines (pure reference otherwise).
T1  Round-trip bit-exactness on smoke data (random, repetitive, zeros,
    empty, 1 byte, Markov, real tex), both auto mode and forced
    two-pass; tamper rejection in the repair stream AND in the model
    section of a two-pass archive.
T2  Random access on a two-pass archive: 80 random (p, k) reads equal
    sequential slices; warm-up <= K-1; decode window <= K + k
    (Theorem 10.3 pattern is unchanged: frozen model = no cross-block
    state); verified-read path; end-of-source edge reads.
T3  Determinism: two in-process pack2 calls identical; two separate-
    process CLI encodes identical to each other and to the in-process
    archive; decodes identical.
T4  Accounting identity (two-pass): repair-stream bits of coded
    sub-blocks == sum of -log2(model probability) with the frozen
    model, up to coder redundancy in [0, 1e-3) bpb at K = 64 KiB.
    Single-pass redundancy reported for reference.
T5  Fallback on white noise: two-pass total (model + stream) >= single-
    pass size, so pack2 emits the single-pass archive (bit-identical to
    the reference), <= 8.2 bpb, still decodable via unpack2.
T6  FIRST SIZE SIGNAL (2 MiB of enwik8, 2 MiB of scripts-src tar;
    K in {16 KiB, 64 KiB}; W=3): bpb single-pass vs two-pass, model-
    header bpb share; gate = two-pass never loses (chosen <= single)
    and spot random-access reads on the two-pass archives are correct.
    The predicted trend (two-pass advantage grows as K shrinks) is
    reported, not gated.

Run:  /Users/para/.venvs/rnr/bin/python impl_2pass/run_checks.py [--quick]
      (--quick skips T6, the ~2-3 minute size-signal part)
"""

import argparse
import math
import os
import random
import struct
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "impl"))

import two_pass as tp  # noqa: E402
import rnr1  # noqa: E402

PYTHON = sys.executable
TWO_PASS_CLI = os.path.join(HERE, "two_pass.py")

RESULTS = []


def record(name, ok, detail=""):
    line = "PASS" if ok else "FAIL"
    print("%s: %s%s" % (line, name, (" -- " + detail) if detail else ""))
    RESULTS.append((name, ok))
    return ok


def note(msg):
    print("      " + msg)


# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

def markov_sample(n, seed=555):
    """Order-2 Markov chain over 4 byte values (compressible, non-text)."""
    rng = random.Random(7)
    base = [0.60, 0.22, 0.13, 0.05]
    tensor = {}
    for s1 in range(4):
        for s2 in range(4):
            row = base[:]
            rng.shuffle(row)
            tensor[(s1, s2)] = row
    rng = random.Random(seed)
    s1, s2 = rng.randrange(4), rng.randrange(4)
    out = bytearray()
    symbols = list(range(4))
    alphabet = [97, 98, 99, 100]
    for _ in range(n):
        row = tensor[(s1, s2)]
        s3 = rng.choices(symbols, weights=row)[0]
        out.append(alphabet[s3])
        s1, s2 = s2, s3
    return bytes(out)


def tex_bytes(n):
    with open(os.path.join(REPO, "tex", "rnr_coding.tex"), "rb") as f:
        return f.read()[:n]


# ---------------------------------------------------------------------------
# T0: fast single-pass engine identity (lets T5/T6 stay CPU-cheap)
# ---------------------------------------------------------------------------

FAST = None  # module or None


def check_t0():
    global FAST
    detail = ""
    ok = True
    try:
        sys.path.insert(0, os.path.join(REPO, "impl_fast"))
        import rnr1fast  # noqa: F401
        text = tex_bytes(65536)
        noise = bytes(random.Random(4).randrange(256) for _ in range(32768))
        id1 = rnr1fast.pack(text, W=3, K=16384) == rnr1.pack(text, W=3, K=16384)
        id2 = rnr1fast.pack(noise, W=2, K=8192) == rnr1.pack(noise, W=2, K=8192)
        ok = id1 and id2
        if ok:
            FAST = rnr1fast
            detail = "text and noise archives byte-identical; fast engine " \
                     "used for single-pass baselines below"
        else:
            detail = "MISMATCH (text=%s noise=%s); falling back to the " \
                     "pure reference everywhere" % (id1, id2)
    except Exception as e:  # library missing/unbuilt
        detail = "fast engine unavailable (%s); pure reference used " \
                 "(slower, same result)" % e
    return record("T0 fast-engine byte-identity spot check", ok, detail)


def single_pack(data, W, K):
    """Reference single-pass archive bytes (fast engine if T0 passed)."""
    if FAST is not None:
        return FAST.pack(data, W=W, K=K)
    return rnr1.pack(data, W=W, K=K)


# ---------------------------------------------------------------------------
# T1: round trips + tamper rejection
# ---------------------------------------------------------------------------

def check_t1():
    rng = random.Random(123)
    cases = {
        "random-8KiB": bytes(rng.randrange(256) for _ in range(8192)),
        "repetitive-8KiB": (b"the quick brown fox " * 500)[:8192],
        "zeros-4KiB": bytes(4096),
        "empty": b"",
        "one-byte": b"Q",
        "markov-64KiB": markov_sample(65536),
        "tex-128KiB": tex_bytes(131072),
    }
    ok = True
    details = []
    for name, data in cases.items():
        raw, info = tp.pack2(data, W=2, K=16384, return_info=True)
        good = tp.unpack2(raw) == data
        # never-lose invariant of the emitted archive
        good &= info["single_bytes"] is None \
            or info["archive_bytes"] <= info["single_bytes"]
        if data:  # forced two-pass path must round-trip as well
            raw2 = tp.pack2(data, W=2, K=16384, force="two")
            good &= tp.unpack2(raw2) == data
        ok &= good
        details.append("%s %.3f bpb %s%s"
                       % (name, 8.0 * len(raw) / max(1, len(data)),
                          info["mode"], "" if good else " MISMATCH"))
    # Tamper rejection on a two-pass archive: repair stream, then model.
    data = cases["repetitive-8KiB"]
    raw = tp.pack2(data, W=2, K=16384, force="two")
    t1 = bytearray(raw)
    t1[-10] ^= 0x40  # repair stream area
    try:
        tp.unpack2(bytes(t1))
        ok = False
        details.append("tampered stream accepted (BAD)")
    except (ValueError, AssertionError):
        details.append("stream tamper rejected")
    t2 = bytearray(raw)
    t2[rnr1.HEADER_SIZE + 4 + 8] ^= 0x01  # inside the model blob
    try:
        tp.unpack2(bytes(t2))
        ok = False
        details.append("tampered model accepted (BAD)")
    except (ValueError, AssertionError):
        details.append("model tamper rejected")
    return record("T1 round-trip bit-exact (auto + forced two-pass) + "
                  "tamper rejection", ok, "; ".join(details))


# ---------------------------------------------------------------------------
# T2: random access on a two-pass archive
# ---------------------------------------------------------------------------

def check_t2():
    K = 16384
    data = tex_bytes(262144)
    raw = tp.pack2(data, W=3, K=K, force="two")
    arc = tp.open_archive(raw)
    ok = isinstance(arc, tp.TwoPassArchive)
    n = arc.n
    rng = random.Random(2026)
    max_warm = 0
    bad = 0
    for _ in range(80):
        p = rng.randrange(n)
        k = min(1 << rng.randrange(0, 16), 3 * K)
        got, stats = arc.read(p, k)
        if got != data[p: p + k]:
            bad += 1
            ok = False
        ok &= stats["warmup_bytes"] <= K - 1
        ok &= stats["decoded_bytes"] <= K + k
        max_warm = max(max_warm, stats["warmup_bytes"])
    for _ in range(5):  # verified-read path (full covering sub-blocks)
        p = rng.randrange(n - 4096)
        got, stats = arc.read(p, 4096, verify=True)
        ok &= got == data[p: p + 4096]
        ok &= stats["decoded_bytes"] <= 2 * K + 4096
    got, _ = arc.read(n - 100, 5000)   # short read at end
    ok &= got == data[n - 100:]
    got, _ = arc.read(n, 10)           # empty read at n
    ok &= got == b""
    return record(
        "T2 random access on two-pass archive: 80 reads == slices; "
        "warm-up <= K-1; window <= K+k",
        ok, "mismatches=%d, max warm-up %d bytes (K=%d)" % (bad, max_warm, K))


# ---------------------------------------------------------------------------
# T3: determinism (in-process, and across separate processes)
# ---------------------------------------------------------------------------

def check_t3():
    data = markov_sample(65536)
    ok = True
    details = []
    raw_a = tp.pack2(data, W=2, K=16384)
    raw_b = tp.pack2(data, W=2, K=16384)
    ok &= raw_a == raw_b
    details.append("in-process encode x2 %s"
                   % ("identical" if raw_a == raw_b else "DIFFER"))
    with tempfile.TemporaryDirectory(prefix="rnr_2pass_checks_") as td:
        src = os.path.join(td, "src.bin")
        with open(src, "wb") as f:
            f.write(data)
        arcs, outs = [], []
        for i in range(2):
            arc = os.path.join(td, "a%d.rnr" % i)
            out = os.path.join(td, "o%d.bin" % i)
            subprocess.run([PYTHON, TWO_PASS_CLI, "pack", src, arc,
                            "--W", "2", "--K", "16384"],
                           check=True, capture_output=True)
            subprocess.run([PYTHON, TWO_PASS_CLI, "unpack", arc, out],
                           check=True, capture_output=True)
            with open(arc, "rb") as f:
                arcs.append(f.read())
            with open(out, "rb") as f:
                outs.append(f.read())
        ok &= arcs[0] == arcs[1]
        details.append("2-process encode %s"
                       % ("identical" if arcs[0] == arcs[1] else "DIFFER"))
        ok &= arcs[0] == raw_a
        details.append("CLI == API archive %s" % ("yes" if arcs[0] == raw_a else "NO"))
        ok &= outs[0] == outs[1] == data
        details.append("2-process decode %s"
                       % ("identical & correct" if outs[0] == outs[1] == data
                          else "DIFFER"))
    return record("T3 determinism: archives byte-identical across repeated "
                  "runs and separate processes", ok, "; ".join(details))


# ---------------------------------------------------------------------------
# T4: accounting identity (stream bits vs frozen-model probabilities)
# ---------------------------------------------------------------------------

def frozen_ideal_bits(block, model):
    """Sum of -log2(model prob) over one sub-block under the frozen
    model (offline float analysis of the exact integer distributions)."""
    bits = 0.0
    hk = 0
    avail = 0
    W = model.W
    mask = model.mask
    for x in block:
        c = model.ctx(hk, avail)
        bits += math.log2(c.total) - math.log2(int(c.freq[x]))
        hk = ((hk << 8) | x) & mask
        if avail < W:
            avail += 1
    return bits


def adaptive_ideal_bits(block, W):
    """Same for the reference adaptive path (per-block cold start)."""
    pred = rnr1.NGramPredictor(W)
    bits = 0.0
    for x in block:
        freq, total = pred.dist()
        bits += math.log2(total) - math.log2(int(freq[x]))
        pred.update(x)
    return bits


def check_t4():
    W, K = 3, 65536
    data = tex_bytes(262144)
    raw = tp.pack2(data, W=W, K=K, force="two")
    arc = tp.open_archive(raw)
    stream_bits = ideal_bits = coded_bytes = 0
    for j in range(arc.m):
        if arc.entries[j][3] != rnr1.MODE_CODED:
            continue
        blk = data[j * K: (j + 1) * K]
        stream_bits += 8 * len(arc._block_blob(j))
        ideal_bits += frozen_ideal_bits(blk, arc.model)
        coded_bytes += len(blk)
    red2 = (stream_bits - ideal_bits) / coded_bytes
    ok = coded_bytes == len(data)  # expect all sub-blocks coded on tex
    ok &= 0.0 <= red2 < 1e-3

    # single-pass redundancy on the same data, for reference
    sraw = single_pack(data, W, K)
    sarc = rnr1.Archive(sraw)
    s_stream = s_ideal = s_bytes = 0
    for j in range(sarc.m):
        if sarc.entries[j][3] != rnr1.MODE_CODED:
            continue
        blk = data[j * K: (j + 1) * K]
        s_stream += 8 * len(sarc._block_blob(j))
        s_ideal += adaptive_ideal_bits(blk, W)
        s_bytes += len(blk)
    red1 = (s_stream - s_ideal) / max(1, s_bytes)
    return record(
        "T4 accounting identity: two-pass stream bits == sum -log2 p "
        "(frozen model) + redundancy in [0, 1e-3) bpb",
        ok,
        "two-pass redundancy %.6f bpb over %d coded bytes (K=64KiB); "
        "single-pass reference redundancy %.6f bpb"
        % (red2, coded_bytes, red1))


# ---------------------------------------------------------------------------
# T5: fallback triggers on white noise
# ---------------------------------------------------------------------------

def check_t5():
    rng = random.Random(99)
    data = bytes(rng.randrange(256) for _ in range(262144))
    sraw = single_pack(data, 3, 65536)
    raw, info = tp.pack2(data, W=3, K=65536, single_raw=sraw,
                         return_info=True)
    ok = info["mode"] == "single-pass"
    ok &= info["two_pass_bytes"] is not None \
        and info["two_pass_bytes"] >= info["single_bytes"]
    ok &= raw == sraw  # emitted bytes are exactly the single-pass archive
    ok &= tp.unpack2(raw) == data
    bpb = 8.0 * len(raw) / len(data)
    ok &= bpb <= 8.2
    return record(
        "T5 fallback on white noise: two-pass (model+stream) >= single, "
        "single-pass archive emitted, <= 8.2 bpb",
        ok,
        "single %d B, two-pass %d B (model %d B), emitted %.4f bpb"
        % (info["single_bytes"], info["two_pass_bytes"],
           info["model_bytes"], bpb))


# ---------------------------------------------------------------------------
# T6: first size signal (2 MiB enwik8 + scripts-src, K in {16, 64} KiB)
# ---------------------------------------------------------------------------

SIZE_SIGNAL_MIB = 2  # 2-4 MB window requested; 2 MiB keeps CPU minutes small


def check_t6():
    W = 3
    nbytes = SIZE_SIGNAL_MIB << 20
    corpora = {
        "enwik8": os.path.join(REPO, "data", "payloads", "enwik8"),
        "scripts-src": os.path.join(REPO, "data", "payloads",
                                    "rnr_scripts_src.tar"),
    }
    ok = True
    rows = []
    rng = random.Random(31337)
    for cname, path in corpora.items():
        if not os.path.exists(path):
            ok = False
            rows.append((cname, None))
            note("corpus missing: %s" % path)
            continue
        with open(path, "rb") as f:
            data = f.read(nbytes)
        t0 = time.time()
        frozen = tp.build_frozen_model(data, W)[:2]
        train_s = time.time() - t0
        deltas = {}
        for K in (16384, 65536):
            sraw = single_pack(data, W, K)
            t1 = time.time()
            raw, info = tp.pack2(data, W=W, K=K, single_raw=sraw,
                                 frozen=frozen, return_info=True)
            enc_s = time.time() - t1
            sbpb = 8.0 * info["single_bytes"] / len(data)
            tbpb = 8.0 * info["two_pass_bytes"] / len(data)
            mshare = 8.0 * info["model_bytes"] / len(data)
            deltas[K] = sbpb - tbpb
            ok &= info["archive_bytes"] <= info["single_bytes"]  # never lose
            # spot random-access reads when the two-pass side won
            ra_ok = True
            if info["mode"] == "two-pass":
                arc = tp.open_archive(raw)
                for _ in range(8):
                    p = rng.randrange(len(data))
                    k = rng.randrange(1, 8192)
                    got, st = arc.read(p, k)
                    ra_ok &= got == data[p: p + k]
                    ra_ok &= st["decoded_bytes"] <= K + k
            ok &= ra_ok
            rows.append((cname, dict(K=K, single_bpb=sbpb, two_bpb=tbpb,
                                     model_bpb=mshare, mode=info["mode"],
                                     ra_ok=ra_ok, enc_s=enc_s,
                                     train_s=train_s)))
        if 16384 in deltas and 65536 in deltas:
            trend = "confirmed" if deltas[16384] >= deltas[65536] else \
                    "NOT confirmed on this input"
            note("%s: two-pass advantage %.4f bpb at K=16KiB vs %.4f at "
                 "K=64KiB -- prediction (wins grow as K shrinks) %s"
                 % (cname, deltas[16384], deltas[65536], trend))
    for cname, r in rows:
        if r is None:
            continue
        note("%s K=%-6d single %.4f bpb | two-pass %.4f bpb "
             "(model %.4f bpb) | mode=%s | train %.1fs + enc %.1fs"
             % (cname, r["K"], r["single_bpb"], r["two_bpb"],
                r["model_bpb"], r["mode"], r["train_s"], r["enc_s"]))
    return record(
        "T6 size signal on %d MiB enwik8 + scripts-src, K in {16,64} KiB: "
        "two-pass never loses; spot random access correct" % SIZE_SIGNAL_MIB,
        ok, "%d cells measured" % sum(1 for _, r in rows if r))


# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="skip T6 (the ~2-3 min size-signal measurement)")
    args = ap.parse_args(argv)

    t_start = time.time()
    print("=" * 72)
    print("run_checks: two-pass (train-then-freeze) mode "
          "(impl_2pass/two_pass.py)")
    print("=" * 72)

    check_t0()
    check_t1()
    check_t2()
    check_t3()
    check_t4()
    check_t5()
    if not args.quick:
        check_t6()
    else:
        print("SKIP: T6 size signal (--quick)")

    n_pass = sum(1 for _, ok in RESULTS if ok)
    n_all = len(RESULTS)
    overall = n_pass == n_all
    print("-" * 72)
    print("OVERALL: %s (%d/%d checks passed, %.1f s)"
          % ("PASS" if overall else "FAIL", n_pass, n_all,
             time.time() - t_start))
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())

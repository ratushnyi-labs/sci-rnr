#!/usr/bin/env python3
"""
Executable checks for the fast C RNR Type-I coder (impl_fast/rnr1_fast.c).

Prints one PASS/FAIL line per check and an OVERALL line, in the style of
scripts/verify/.  The Python reference impl/rnr1.py stays NORMATIVE; the C
implementation is accepted only through bit-identity against it.

Checks
------
F0  Build: librnr1fast.dylib compiles from rnr1_fast.c with -O3.
F1  Parallel-reference trust gate: ref_driver.parallel_pack/parallel_unpack
    (process pool over independent sub-blocks, container assembled by the
    reference's own code paths) is byte-identical to sequential rnr1.pack /
    rnr1.unpack on multi-sub-block inputs.  Only after this gate passes are
    parallel reference results used for the corpus matrix.
F2  Smoke-class bit-identity matrix at W in {2,3}: the generated classes of
    impl/run_checks.py C2 (random/repetitive/zeros/empty/one-byte/Markov)
    plus a real text slice, at K = 64 KiB and K = 16 KiB, checking
    (a) fast pack == reference pack (byte-identical),
    (b) fast unpack(reference archive) == source,
    (c) reference unpack(fast archive) == source,
    (d) fast unpack(fast archive) == source.
F3  Corpus bit-identity matrix at W in {2,3}, K = 64 KiB: EVERY non-deferred
    manifest entry (enwik8, calgary, canterbury, scripts_src, sqlite_synth,
    telemetry_grid, ctrl_urandom, ctrl_zstd, ctrl_base64): fast pack ==
    reference pack, byte for byte; the multi-threaded pack (8 threads over
    independent sub-blocks) is additionally required to equal the single-
    threaded archive byte for byte (F3b within the same matrix).
F4  Corpus cross-decode matrix: fast decodes every reference archive of F3
    and the REFERENCE decodes every fast archive of F3 (both with full
    verification: sub-block hashes + H_v), plus fast-decodes-fast and
    multi-threaded fast decode (F4b).
F5  Error-path conformance: tampered repair stream rejected; truncated
    header rejected; wrong magic rejected; model-hash mismatch rejected.
F6  Determinism: two separate-process CLI encodes of the same input are
    byte-identical, and equal to the in-process API result.
F7  Throughput: C encode speed on enwik8 (W=3, K=64 KiB), target
    >= 10 MB/s; decode speed reported; peak RSS reported.

Reference archives for F3/F4 are produced by the parallel reference driver
(gated by F1) and cached in --refdir to make re-runs cheap.  A full cold run
regenerates them (~25 min on 8 cores; the reference codes ~70-160 KB/s/core).

Run:  /Users/para/.venvs/rnr/bin/python impl_fast/run_checks.py
      [--refdir DIR] [--skip-corpus] [--jobs J]
"""

import argparse
import hashlib
import os
import random
import resource
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "impl"))
sys.path.insert(0, os.path.join(REPO, "data"))

PYTHON = sys.executable
RESULTS = []


def record(name, ok, detail=""):
    line = "PASS" if ok else "FAIL"
    print("%s: %s%s" % (line, name, (" -- " + detail) if detail else ""), flush=True)
    RESULTS.append((name, ok))
    return ok


def markov_bytes(n, seed=555):
    """Same construction family as impl/run_checks.py (order-2 Markov)."""
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
    for _ in range(n):
        row = tensor[(s1, s2)]
        s3 = rng.choices(range(4), weights=row)[0]
        out.append(97 + s3)
        s1, s2 = s2, s3
    return bytes(out)


def check_f0():
    t0 = time.time()
    r = subprocess.run(["sh", os.path.join(HERE, "build.sh")],
                       capture_output=True, text=True)
    ok = r.returncode == 0 and os.path.exists(os.path.join(HERE, "librnr1fast.dylib"))
    return record("F0 build librnr1fast.dylib (clang -O3)", ok,
                  "%.1fs%s" % (time.time() - t0,
                               ("; " + r.stderr.strip()) if r.returncode else ""))


def check_f1(rnr1, ref_driver, jobs):
    data = open(os.path.join(REPO, "tex", "rnr_coding.tex"), "rb").read()[:300000]
    ok = True
    details = []
    for W, K in ((2, 65536), (3, 16384)):
        a = rnr1.pack(data, W=W, K=K)
        b = ref_driver.parallel_pack(data, W, K, jobs)
        ok &= a == b
        ok &= ref_driver.parallel_unpack(a, jobs) == data
        ok &= rnr1.unpack(b) == data
        details.append("W=%d K=%d %s" % (W, K, "identical" if a == b else "DIFFER"))
    return record("F1 parallel reference driver == sequential rnr1.pack/unpack",
                  ok, "; ".join(details))


def smoke_cases():
    rng = random.Random(123)
    return {
        "random-8KiB": bytes(rng.randrange(256) for _ in range(8192)),
        "repetitive-8KiB": (b"the quick brown fox " * 500)[:8192],
        "zeros-4KiB": bytes(4096),
        "empty": b"",
        "one-byte": b"Q",
        "markov-64KiB": markov_bytes(65536),
        "text-192KiB": open(os.path.join(REPO, "tex", "rnr_coding.tex"),
                            "rb").read()[:196608],
    }


def check_f2(rnr1, rnr1fast):
    ok = True
    n_pairs = 0
    bad = []
    for W in (2, 3):
        for K in (65536, 16384):
            for name, data in smoke_cases().items():
                a = rnr1.pack(data, W=W, K=K)
                b = rnr1fast.pack(data, W=W, K=K)
                good = (a == b
                        and rnr1fast.unpack(a) == data
                        and rnr1.unpack(b) == data
                        and rnr1fast.unpack(b) == data)
                ok &= good
                n_pairs += 1
                if not good:
                    bad.append("W=%d K=%d %s" % (W, K, name))
    return record("F2 smoke-class bit-identity + cross-decode (W in {2,3}, "
                  "K in {64,16} KiB)", ok,
                  "%d configurations%s" % (n_pairs,
                                           ("; BAD: " + ", ".join(bad)) if bad else ""))


def ensure_ref_archive(ref_driver, loader, refdir, name, W, K, jobs):
    dst = os.path.join(refdir, "%s.W%d.rnr1" % (name, W))
    if not os.path.exists(dst):
        with open(loader.payload_path(name), "rb") as f:
            data = f.read()
        raw = ref_driver.parallel_pack(data, W, K, jobs)
        with open(dst + ".tmp", "wb") as f:
            f.write(raw)
        os.rename(dst + ".tmp", dst)
    return dst


def check_f3_f4(rnr1, rnr1fast, ref_driver, loader, refdir, fastdir, jobs):
    K = 65536
    names = loader.corpus_names()
    names.sort(key=lambda n: loader.corpus_bytes(n))  # small first: fail fast
    ok3 = True
    ok4 = True
    lines = []
    os.makedirs(fastdir, exist_ok=True)
    for W in (2, 3):
        for name in names:
            t0 = time.time()
            with open(loader.payload_path(name), "rb") as f:
                data = f.read()
            ref_path = ensure_ref_archive(ref_driver, loader, refdir, name, W, K, jobs)
            with open(ref_path, "rb") as f:
                ref_arc = f.read()
            fast_arc = rnr1fast.pack(data, W=W, K=K)
            ident = fast_arc == ref_arc
            mt_ident = rnr1fast.pack(data, W=W, K=K, threads=8) == fast_arc
            ok3 &= ident and mt_ident
            fast_path = os.path.join(fastdir, "%s.W%d.rnr1" % (name, W))
            with open(fast_path, "wb") as f:
                f.write(fast_arc)
            # cross-decode: fast <- ref archive; reference <- fast archive
            d_fast_ref = rnr1fast.unpack(ref_arc) == data
            d_fast_fast = rnr1fast.unpack(fast_arc) == data
            d_fast_mt = rnr1fast.unpack(fast_arc, threads=8) == data
            with open(fast_path, "rb") as f:
                d_ref_fast = ref_driver.parallel_unpack(f.read(), jobs) == data
            dec_ok = d_fast_ref and d_ref_fast and d_fast_fast and d_fast_mt
            ok4 &= dec_ok
            lines.append("%-14s W=%d %s%s%s (%.0fs)"
                         % (name, W,
                            "identical" if ident else "DIFFER",
                            "" if mt_ident else " MT-DIFFER",
                            "" if dec_ok else " DECODE-BAD",
                            time.time() - t0))
            print("   f3/f4: " + lines[-1], flush=True)
            del data, ref_arc, fast_arc
    record("F3 corpus bit-identity: fast pack == reference pack, mt pack "
           "== st pack (all %d non-deferred entries x W in {2,3}, K=64KiB)"
           % len(names), ok3)
    record("F4 corpus cross-decode: fast<-ref, ref<-fast, fast<-fast, "
           "mt<-fast (full verification)", ok4)
    return ok3 and ok4


def check_f5(rnr1, rnr1fast):
    data = markov_bytes(65536)
    raw = bytearray(rnr1fast.pack(data, W=2, K=16384))
    ok = True
    details = []
    # tamper with repair stream
    t = bytearray(raw); t[-10] ^= 0x40
    for impl, tag in ((rnr1fast, "fast"), (rnr1, "ref")):
        try:
            impl.unpack(bytes(t))
            ok = False
            details.append("%s accepted tamper (BAD)" % tag)
        except ValueError:
            details.append("%s rejects tamper" % tag)
    # truncated header
    try:
        rnr1fast.unpack(bytes(raw[:20])); ok = False
    except ValueError:
        details.append("truncated rejected")
    # wrong magic
    t = bytearray(raw); t[0] ^= 0xFF
    try:
        rnr1fast.unpack(bytes(t)); ok = False
    except ValueError:
        details.append("bad magic rejected")
    # model hash mismatch
    t = bytearray(raw); t[24] ^= 0x01
    try:
        rnr1fast.unpack(bytes(t)); ok = False
    except ValueError:
        details.append("model-hash mismatch rejected")
    return record("F5 error paths: tamper/truncation/magic/model-hash all "
                  "rejected", ok, "; ".join(details))


def check_f6(rnr1fast):
    data = markov_bytes(262144, seed=777)
    ok = True
    details = []
    api = rnr1fast.pack(data, W=3, K=65536)
    with tempfile.TemporaryDirectory(prefix="rnr1fast_checks_") as td:
        src = os.path.join(td, "src.bin")
        with open(src, "wb") as f:
            f.write(data)
        arcs = []
        for i in range(2):
            arc = os.path.join(td, "a%d.rnr" % i)
            threads = "1" if i == 0 else "8"   # MT must equal ST across processes
            subprocess.run([PYTHON, os.path.join(HERE, "rnr1fast.py"), "pack",
                            src, arc, "--W", "3", "--K", "65536",
                            "--threads", threads],
                           check=True, capture_output=True)
            with open(arc, "rb") as f:
                arcs.append(f.read())
    ok &= arcs[0] == arcs[1]
    details.append("2-process encode (1 vs 8 threads) %s"
                   % ("identical" if arcs[0] == arcs[1] else "DIFFER"))
    ok &= arcs[0] == api
    details.append("CLI == API %s" % ("yes" if arcs[0] == api else "NO"))
    return record("F6 determinism: separate-process encodes byte-identical", ok,
                  "; ".join(details))


def check_f7(rnr1fast, runs=2):
    path = os.path.join(REPO, "data", "payloads", "enwik8")
    with open(path, "rb") as f:
        data = f.read()
    res = {}
    raw = None
    for threads in (1, 8):
        enc_t, dec_t = [], []
        for _ in range(runs):
            t0 = time.time()
            raw = rnr1fast.pack(data, W=3, K=65536, threads=threads)
            t1 = time.time()
            out = rnr1fast.unpack(raw, threads=threads)
            t2 = time.time()
            assert out == data
            enc_t.append(t1 - t0); dec_t.append(t2 - t1)
            del out
        res[threads] = (len(data) / 1e6 / min(enc_t),
                        len(data) / 1e6 / min(dec_t))
    rss_mib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1 << 20)
    ok = res[8][0] >= 10.0     # target met by the deliverable configuration
    return record("F7 throughput on enwik8 (W=3, K=64KiB): encode >= 10 MB/s",
                  ok,
                  "1-thread enc %.2f dec %.2f MB/s; 8-thread enc %.2f dec "
                  "%.2f MB/s (best of %d); %.4f bpb; peak RSS %.0f MiB "
                  "(process-wide, incl. corpus buffers)"
                  % (res[1][0], res[1][1], res[8][0], res[8][1],
                     runs, 8.0 * len(raw) / len(data), rss_mib))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refdir", default=os.environ.get(
        "RNR_FAST_REFDIR", os.path.join(HERE, "workdir", "refarc")),
        help="cache dir for parallel-reference archives (F3/F4)")
    ap.add_argument("--fastdir", default=os.path.join(HERE, "workdir", "fastarc"))
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--skip-corpus", action="store_true",
                    help="skip F3/F4 (full-corpus matrix; slow on cold cache)")
    args = ap.parse_args()

    t_start = time.time()
    print("=" * 72)
    print("run_checks: fast C RNR Type-I coder (impl_fast/rnr1_fast.c)")
    print("reference: impl/rnr1.py (normative)")
    print("=" * 72)

    check_f0()
    import rnr1
    import rnr1fast
    import ref_driver
    import loader

    check_f1(rnr1, ref_driver, args.jobs)
    check_f2(rnr1, rnr1fast)
    if args.skip_corpus:
        print("SKIP: F3/F4 corpus matrix (--skip-corpus)")
    else:
        os.makedirs(args.refdir, exist_ok=True)
        check_f3_f4(rnr1, rnr1fast, ref_driver, loader,
                    args.refdir, args.fastdir, args.jobs)
    check_f5(rnr1, rnr1fast)
    check_f6(rnr1fast)
    check_f7(rnr1fast)

    n_pass = sum(1 for _, ok in RESULTS if ok)
    n_all = len(RESULTS)
    overall = n_pass == n_all
    print("-" * 72)
    print("OVERALL: %s (%d/%d checks passed, %.1f s)"
          % ("PASS" if overall else "FAIL", n_pass, n_all, time.time() - t_start))
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())

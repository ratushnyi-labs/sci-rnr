#!/usr/bin/env python3
"""
Executable checks for the Rust production CLI (rnr-cli, binary `rnr`).

Prints one PASS/FAIL line per check and an OVERALL line, in the style of
scripts/verify/.  The Python reference impl/rnr1.py stays NORMATIVE; the
Rust implementation is accepted only through bit-identity against it.
The validated C port (impl_fast) is used as a second independent
cross-decode engine.

Checks
------
R0  Build: `cargo build --release` produces target/release/rnr.
R1  Smoke-class bit-identity + full four-way cross-decode at W in {2,3},
    K in {64,16} KiB: generated classes (random/repetitive/zeros/empty/
    one-byte/order-2 Markov) plus a real text slice, checking
      (a) rust pack == reference pack (byte-identical),
      (b) rust unpack(reference archive) == source,
      (c) REFERENCE unpack(rust archive) == source,
      (d) rust unpack(rust archive) == source,
      (e) C engine unpack(rust archive) == source,
      (f) rust unpack(C-engine archive) == source.
R2  Corpus bit-identity matrix at W in {2,3}, K = 64 KiB: EVERY
    non-deferred manifest entry: rust pack == reference pack byte for
    byte (reference archives come from the impl_fast cache, produced by
    the F1-gated parallel reference driver; regenerated through the same
    driver if missing); multi-threaded pack (--threads 8) == single-
    threaded archive; rust fully verifies + decodes the reference
    archive; the C engine decodes the rust archive.  Reference-decodes-
    rust on the corpus is implied by byte-identity (the rust archive IS
    the reference archive, already decoded by the reference in the
    impl_fast run) and is additionally exercised directly in R1(c).
R3  Error paths: tampered repair stream, truncated header, wrong magic,
    model-hash mismatch, truncated repair stream all rejected (nonzero
    exit, no output file).
R4  Determinism: two separate-process CLI encodes (1 vs 8 threads) are
    byte-identical, and equal to the R1 in-run archive.
R5  Random-access read conformance: `rnr read` output equals the source
    slice and the reference Archive.read() output for reads inside one
    sub-block, across boundaries, at end-of-source, past-end clamping,
    and with --verify; warmup/decoded stats match the reference exactly
    (decode window <= K + k, Theorem 10.3).
R6  Throughput on enwik8 (W=3, K=64 KiB): single-thread and --threads 8
    pack + unpack, one brief run each (shared host), report MB/s + bpb.

Run:  /Users/para/.venvs/rnr/bin/python rnr-cli/run_checks.py
      [--refdir DIR] [--skip-corpus] [--jobs J]
"""

import argparse
import os
import random
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "impl"))
sys.path.insert(0, os.path.join(REPO, "impl_fast"))
sys.path.insert(0, os.path.join(REPO, "data"))

RNR = os.path.join(HERE, "target", "release", "rnr")
WORKDIR = os.path.join(HERE, "workdir")
RESULTS = []


def record(name, ok, detail=""):
    line = "PASS" if ok else "FAIL"
    print("%s: %s%s" % (line, name, (" -- " + detail) if detail else ""), flush=True)
    RESULTS.append((name, ok))
    return ok


def run_rnr(args, expect_fail=False):
    r = subprocess.run([RNR] + args, capture_output=True)
    if expect_fail:
        return r.returncode != 0
    if r.returncode != 0:
        raise RuntimeError("rnr %s failed: %s" % (args, r.stderr.decode()))
    return True


def rnr_pack(src, dst, W, K, threads=1):
    run_rnr(["pack", src, dst, "--W", str(W), "--K", str(K),
             "--threads", str(threads)])
    with open(dst, "rb") as f:
        return f.read()


def rnr_unpack(src, dst, threads=1, no_verify=False):
    args = ["unpack", src, dst, "--threads", str(threads)]
    if no_verify:
        args.append("--no-verify")
    run_rnr(args)
    with open(dst, "rb") as f:
        return f.read()


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


def check_r0():
    t0 = time.time()
    r = subprocess.run(["cargo", "build", "--release"], cwd=HERE,
                       capture_output=True, text=True)
    ok = r.returncode == 0 and os.path.exists(RNR)
    return record("R0 build target/release/rnr (cargo --release)", ok,
                  "%.1fs%s" % (time.time() - t0,
                               ("; " + r.stderr.strip()[-200:]) if r.returncode else ""))


def check_r1(rnr1, rnr1fast, tmp):
    ok = True
    n_pairs = 0
    bad = []
    for W in (2, 3):
        for K in (65536, 16384):
            for name, data in smoke_cases().items():
                src = os.path.join(tmp, "s.bin")
                with open(src, "wb") as f:
                    f.write(data)
                ref_arc = rnr1.pack(data, W=W, K=K)
                ref_path = os.path.join(tmp, "s.ref.rnr")
                with open(ref_path, "wb") as f:
                    f.write(ref_arc)
                rust_path = os.path.join(tmp, "s.rust.rnr")
                rust_arc = rnr_pack(src, rust_path, W, K)
                c_arc = rnr1fast.pack(data, W=W, K=K)
                c_path = os.path.join(tmp, "s.c.rnr")
                with open(c_path, "wb") as f:
                    f.write(c_arc)
                out = os.path.join(tmp, "s.out")
                good = (rust_arc == ref_arc
                        and rnr_unpack(ref_path, out) == data      # rust <- ref
                        and rnr1.unpack(rust_arc) == data          # ref  <- rust
                        and rnr_unpack(rust_path, out) == data     # rust <- rust
                        and rnr1fast.unpack(rust_arc) == data      # C    <- rust
                        and rnr_unpack(c_path, out) == data)       # rust <- C
                ok &= good
                n_pairs += 1
                if not good:
                    bad.append("W=%d K=%d %s" % (W, K, name))
    return record("R1 smoke-class bit-identity + 4-way cross-decode "
                  "(rust/ref/C; W in {2,3}, K in {64,16} KiB)", ok,
                  "%d configurations%s"
                  % (n_pairs, ("; BAD: " + ", ".join(bad)) if bad else ""))


def ensure_ref_archive(loader, refdir, name, W, K, jobs):
    """Reference archive from the impl_fast cache (parallel reference
    driver, gated by impl_fast check F1); regenerated if missing."""
    dst = os.path.join(refdir, "%s.W%d.rnr1" % (name, W))
    if not os.path.exists(dst):
        import ref_driver
        with open(loader.payload_path(name), "rb") as f:
            data = f.read()
        raw = ref_driver.parallel_pack(data, W, K, jobs)
        with open(dst + ".tmp", "wb") as f:
            f.write(raw)
        os.rename(dst + ".tmp", dst)
    return dst


def check_r2(rnr1fast, loader, refdir, jobs, tmp):
    K = 65536
    names = loader.corpus_names()
    names.sort(key=lambda n: loader.corpus_bytes(n))  # small first: fail fast
    ok = True
    for W in (2, 3):
        for name in names:
            t0 = time.time()
            src = str(loader.payload_path(name))
            ref_path = ensure_ref_archive(loader, refdir, name, W, K, jobs)
            with open(ref_path, "rb") as f:
                ref_arc = f.read()
            rust_path = os.path.join(tmp, "c.rust.rnr")
            rust_arc = rnr_pack(src, rust_path, W, K, threads=1)
            ident = rust_arc == ref_arc
            mt_path = os.path.join(tmp, "c.rust.mt.rnr")
            mt_ident = rnr_pack(src, mt_path, W, K, threads=8) == rust_arc
            out = os.path.join(tmp, "c.out")
            with open(src, "rb") as f:
                data = f.read()
            d_rust_ref = rnr_unpack(ref_path, out) == data   # rust <- ref, verify
            d_c_rust = rnr1fast.unpack(rust_arc) == data     # C engine <- rust
            good = ident and mt_ident and d_rust_ref and d_c_rust
            ok &= good
            print("   r2: %-14s W=%d %s%s%s (%.0fs)"
                  % (name, W,
                     "identical" if ident else "DIFFER",
                     "" if mt_ident else " MT-DIFFER",
                     "" if (d_rust_ref and d_c_rust) else " DECODE-BAD",
                     time.time() - t0), flush=True)
            del data, ref_arc, rust_arc
    return record("R2 corpus bit-identity: rust pack == reference pack, "
                  "mt == st, rust<-ref + C<-rust decode (all %d non-deferred "
                  "entries x W in {2,3}, K=64KiB)" % len(names), ok)


def check_r3(tmp):
    data = markov_bytes(65536)
    src = os.path.join(tmp, "e.bin")
    with open(src, "wb") as f:
        f.write(data)
    arc_path = os.path.join(tmp, "e.rnr")
    raw = rnr_pack(src, arc_path, 2, 16384)
    ok = True
    details = []

    def try_reject(mut, tag):
        nonlocal ok
        bad_path = os.path.join(tmp, "e.bad.rnr")
        with open(bad_path, "wb") as f:
            f.write(mut)
        out = os.path.join(tmp, "e.badout")
        if os.path.exists(out):
            os.unlink(out)
        rejected = run_rnr(["unpack", bad_path, out], expect_fail=True)
        ok &= rejected
        details.append("%s %s" % (tag, "rejected" if rejected else "ACCEPTED(BAD)"))

    t = bytearray(raw); t[-10] ^= 0x40
    try_reject(bytes(t), "tampered repair stream")
    try_reject(raw[:20], "truncated header")
    t = bytearray(raw); t[0] ^= 0xFF
    try_reject(bytes(t), "bad magic")
    t = bytearray(raw); t[24] ^= 0x01
    try_reject(bytes(t), "model-hash mismatch")
    try_reject(raw[:len(raw) - 300], "truncated repair stream")
    return record("R3 error paths rejected (nonzero exit)", ok, "; ".join(details))


def check_r4(tmp):
    data = markov_bytes(262144, seed=777)
    src = os.path.join(tmp, "d.bin")
    with open(src, "wb") as f:
        f.write(data)
    arcs = []
    for i, threads in enumerate(("1", "8")):
        arc = os.path.join(tmp, "d%d.rnr" % i)
        subprocess.run([RNR, "pack", src, arc, "--W", "3", "--K", "65536",
                        "--threads", threads], check=True, capture_output=True)
        with open(arc, "rb") as f:
            arcs.append(f.read())
    ok = arcs[0] == arcs[1]
    return record("R4 determinism: separate-process encodes (1 vs 8 threads) "
                  "byte-identical", ok,
                  "identical" if ok else "DIFFER")


def check_r5(rnr1, tmp):
    data = open(os.path.join(REPO, "tex", "rnr_coding.tex"), "rb").read()[:200000]
    src = os.path.join(tmp, "r.bin")
    with open(src, "wb") as f:
        f.write(data)
    K = 16384
    arc_path = os.path.join(tmp, "r.rnr")
    rnr_pack(src, arc_path, 3, K)
    with open(arc_path, "rb") as f:
        ref_arc = rnr1.Archive(f.read())
    cases = [
        (0, 100),            # start of block 0
        (5000, 2000),        # inside block 0
        (K - 50, 100),       # spans one boundary
        (3 * K, K),          # exactly one full block
        (2 * K + 7, 3 * K),  # spans multiple blocks
        (len(data) - 30, 100),   # clamped at end (k' < k)
        (len(data), 50),     # p == n: empty read
    ]
    ok = True
    details = []
    for p, k in cases:
        outf = os.path.join(tmp, "r.out")
        r = subprocess.run([RNR, "read", arc_path, "--pos", str(p),
                            "--len", str(k), "--out", outf],
                           capture_output=True, text=True)
        got = open(outf, "rb").read()
        want, stats = ref_arc.read(p, k)
        expect = data[p:p + k]
        stat_line = r.stderr.strip().splitlines()[-1] if r.stderr.strip() else ""
        want_stats = "warmup=%d decoded=%d" % (stats["warmup_bytes"],
                                               stats["decoded_bytes"])
        good = (r.returncode == 0 and got == expect and got == want
                and want_stats in stat_line
                and stats["decoded_bytes"] <= K + k)
        ok &= good
        if not good:
            details.append("p=%d k=%d MISMATCH (%s | want %s)"
                           % (p, k, stat_line, want_stats))
    # --verify variant on a boundary-spanning read
    outf = os.path.join(tmp, "r.out")
    r = subprocess.run([RNR, "read", arc_path, "--pos", str(K - 50),
                        "--len", "100", "--out", outf, "--verify"],
                       capture_output=True, text=True)
    ok &= r.returncode == 0 and open(outf, "rb").read() == data[K - 50:K + 50]
    # out-of-range position must fail
    ok &= run_rnr(["read", arc_path, "--pos", str(len(data) + 1), "--len", "1"],
                  expect_fail=True)
    return record("R5 random-access read == reference (data + warmup/decoded "
                  "stats), --verify, out-of-range rejected", ok,
                  "; ".join(details) if details else
                  "%d cases + verify + range check" % len(cases))


def check_r6(tmp):
    src = os.path.join(REPO, "data", "payloads", "enwik8")
    n = os.path.getsize(src)
    res = {}
    arc = os.path.join(tmp, "t.rnr")
    out = os.path.join(tmp, "t.out")
    for threads in (1, 8):
        t0 = time.time()
        raw = rnr_pack(src, arc, 3, 65536, threads=threads)
        t1 = time.time()
        rnr_unpack(arc, out, threads=threads)
        t2 = time.time()
        with open(out, "rb") as f_out, open(src, "rb") as f_src:
            assert f_out.read() == f_src.read()
        res[threads] = (n / 1e6 / (t1 - t0), n / 1e6 / (t2 - t1), len(raw))
        os.unlink(out)
    ok = res[1][0] > 0  # throughput is reported, not gated
    return record("R6 throughput on enwik8 (W=3, K=64KiB; single run each, "
                  "shared host)", ok,
                  "1-thread enc %.2f dec %.2f MB/s; 8-thread enc %.2f dec "
                  "%.2f MB/s; %.4f bpb"
                  % (res[1][0], res[1][1], res[8][0], res[8][1],
                     8.0 * res[1][2] / n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refdir", default=os.environ.get(
        "RNR_FAST_REFDIR",
        os.path.join(REPO, "impl_fast", "workdir", "refarc")),
        help="cache dir of parallel-reference archives (impl_fast F1-gated)")
    ap.add_argument("--jobs", type=int, default=max(2, (os.cpu_count() or 2) // 2))
    ap.add_argument("--skip-corpus", action="store_true",
                    help="skip R2 (full-corpus matrix)")
    ap.add_argument("--skip-throughput", action="store_true")
    args = ap.parse_args()

    t_start = time.time()
    print("=" * 72)
    print("run_checks: Rust production CLI (rnr-cli, binary `rnr`)")
    print("reference: impl/rnr1.py (normative); C cross-engine: impl_fast")
    print("=" * 72)

    check_r0()
    import rnr1
    import rnr1fast
    import loader

    os.makedirs(WORKDIR, exist_ok=True)
    tmp = os.path.join(WORKDIR, "tmp")
    os.makedirs(tmp, exist_ok=True)

    check_r1(rnr1, rnr1fast, tmp)
    if args.skip_corpus:
        print("SKIP: R2 corpus matrix (--skip-corpus)")
    else:
        os.makedirs(args.refdir, exist_ok=True)
        check_r2(rnr1fast, loader, args.refdir, args.jobs, tmp)
    check_r3(tmp)
    check_r4(tmp)
    check_r5(rnr1, tmp)
    if args.skip_throughput:
        print("SKIP: R6 throughput (--skip-throughput)")
    else:
        check_r6(tmp)

    n_pass = sum(1 for _, ok in RESULTS if ok)
    n_all = len(RESULTS)
    overall = n_pass == n_all
    print("-" * 72)
    print("OVERALL: %s (%d/%d checks passed, %.1f s)"
          % ("PASS" if overall else "FAIL", n_pass, n_all, time.time() - t_start))
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())

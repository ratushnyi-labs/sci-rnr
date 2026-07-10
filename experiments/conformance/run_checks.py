#!/usr/bin/env python3
"""Determinism & conformance gate for the reference RNR Type-I coder.

Covers hypotheses H8/H9/H14 (exp-design sections 2.7, 2.10, 8; REPRO
family of bench/METHODS.md section 5: binary pass/fail conformance,
exact counts, no p-values) and the engineering-spec edge requirements
R-3.2/R-3.4/R-5.2 (overflow guard, reduction order), for Part I
Theorems 10.1/10.2 (bit-exact forward pass; compatible online
optimizer).

Legs
----
  A   native, separate process, PYTHONHASHSEED=0   (reference leg)
  A2  native, separate process, PYTHONHASHSEED=0   (across-process rerun)
  B   native, separate process, PYTHONHASHSEED=1   (hash-seed variation)
  C   native, separate process, PYTHONHASHSEED unset (randomized hashing)
  D   x86-64 Linux container via docker --platform linux/amd64
      (cross-ISA leg; arm64 host runs A/A2/B/C natively)

Every leg re-encodes the identical deterministic test-block corpus and
reports archive SHA-256s (H8), per-sub-block integer inference-trace
digests (H14), and the per-sub-block adaptive-state hash ladder (H9).
Leg D also cross-decodes leg A's archive files, and a native decode-only
leg E cross-decodes leg D's archive files (exp-design 8.3 steps 5-6).

Checks
------
K1  corpus deterministically built; spec hashes stable
K2  round-trip + sub-block/H_v verification on every block (leg A)
K3  same-process double-encode byte-identical (leg A, first blocks)
K4  across-process byte-identity: leg A == leg A2 (H8)
K5  PYTHONHASHSEED invariance: A == B == C (H8)
K6  cross-ISA archive byte-identity: leg A == leg D, file bytes compared
    directly, not only hashes (H8)                      [SKIP if no docker]
K7  cross-ISA integer inference-trace identity: every per-sub-block
    trace digest equal A vs D (H14)                     [SKIP if no docker]
K8  cross-ISA adaptive-state ladder identity: every rung equal A vs D,
    plus encoder-vs-decoder state agreement per leg (H9) [SKIP if no docker]
K9  cross-decode: D decodes A's archives, E (native) decodes D's; all
    outputs verify against the original block hashes    [SKIP if no docker]
K10 overflow-guard headroom: stress block drives per-context totals to
    the COUNT_CAP halving regime; max arithmetic-coder total observed
    equals 256 + 32*(COUNT_CAP-1) = 2097376, an 8.0x margin below the
    2^24 bound (R-3.2/R-5.2)
K11 overflow guard trips: coder rejects total == 2^24 and accepts
    2^24 - 1 at the exact boundary; fault-injected predictor scale
    (COUNT_SCALE -> 2^20) trips the dist() guard within a few bytes
K12 reduction-order independence probe: on sampled live predictor
    states, class partition sums / totals / prefix sums recomputed in
    randomized orders equal the fixed-order values exactly (R-3.4);
    dtype audit int64 (R-3.1)

R-3.4 coverage argument (documented here, verified by K6/K7/K12): the
coding path is integer-only, so every reduction is exact and its VALUE
is order-independent; the implementation additionally fixes the order
syntactically (np.cumsum / np.add.at / left-to-right loops).  Any
order-sensitive reduction could only arise from floating point, and
would produce different bits on at least one of the two ISAs (NEON
vs SSE/AVX kernels differ in blocking and fused ops); the per-position
cross-ISA trace-hash identity of K7 is therefore an end-to-end witness
for R-3.4 on real hardware, and K12 is the direct in-situ probe.

Run:
  /Users/para/.venvs/rnr/bin/python -u experiments/conformance/run_checks.py
        (reduced scale, < 10 min; runs the docker leg when available,
         SKIPs it -- never PASSes it -- when unavailable)
  ... run_checks.py --full          (100-block campaign, exp-design 8.1)
  ... run_checks.py --no-docker     (force-skip the cross-ISA leg)
  ... run_checks.py --no-record     (do not append to results JSONL)
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "impl"))
sys.path.insert(0, os.path.join(REPO, "bench"))

import numpy as np  # noqa: E402

import metrics  # noqa: E402
import results as results_store  # noqa: E402
import rnr1  # noqa: E402

PYTHON = sys.executable
HARNESS = os.path.join(HERE, "harness.py")
MAKE_CORPUS = os.path.join(HERE, "make_corpus.py")
IMAGE = "rnr-conformance-amd64"
BASE_IMAGE = "python:3.12-slim"
NUMPY_PIN = "2.4.6"
EXPERIMENT_ID = "conformance-h8-h9-h14"
EXPECTED_STRESS_MAX_TOTAL = 256 + 32 * (rnr1.COUNT_CAP - 1)  # 2097376

RESULTS = []


def record(name, ok, detail=""):
    print("%s: %s%s" % ("PASS" if ok else "FAIL", name,
                        (" -- " + detail) if detail else ""), flush=True)
    RESULTS.append((name, bool(ok)))
    return ok


def skip(name, why):
    print("SKIP: %s -- %s" % (name, why), flush=True)
    RESULTS.append((name, None))


def run_logged(cmd, log_path, env=None, timeout=1800):
    with open(log_path, "w", encoding="utf-8") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                              env=env, timeout=timeout)
    return proc.returncode


# ---------------------------------------------------------------------------
# Legs
# ---------------------------------------------------------------------------

def native_leg(workdir, leg, hashseed, cross=None, timeout=1800):
    """Run the harness in a fresh native process; returns report path."""
    outdir = os.path.join(workdir, leg)
    os.makedirs(outdir, exist_ok=True)
    report = os.path.join(outdir, "report.json")
    env = dict(os.environ)
    env.pop("PYTHONHASHSEED", None)
    if hashseed is not None:
        env["PYTHONHASHSEED"] = str(hashseed)
    cmd = [PYTHON, "-u", HARNESS,
           "--spec", os.path.join(workdir, "blocks", "spec.json"),
           "--out", report,
           "--archives-dir", os.path.join(outdir, "archives"),
           "--double-encode", "2"]
    if cross:
        cmd += ["--cross-decode-archives", cross[0],
                "--cross-decode-manifest", cross[1]]
    rc = run_logged(cmd, os.path.join(workdir, "logs", leg + ".log"),
                    env=env, timeout=timeout)
    if rc != 0:
        raise RuntimeError("native leg %s failed (rc=%d); see logs" % (leg, rc))
    return report


def native_decode_leg(workdir, leg, archives_dir, manifest, timeout=1800):
    """Decode-only native leg (cross-decodes another leg's archives)."""
    outdir = os.path.join(workdir, leg)
    os.makedirs(outdir, exist_ok=True)
    report = os.path.join(outdir, "report.json")
    cmd = [PYTHON, "-u", HARNESS, "--out", report,
           "--cross-decode-archives", archives_dir,
           "--cross-decode-manifest", manifest]
    rc = run_logged(cmd, os.path.join(workdir, "logs", leg + ".log"),
                    timeout=timeout)
    if rc != 0:
        raise RuntimeError("decode leg %s failed (rc=%d)" % (leg, rc))
    return report


def docker_available():
    try:
        rc = subprocess.run(["docker", "info"], capture_output=True,
                            timeout=60).returncode
        return rc == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def ensure_docker_image():
    """Return (image, runtime_pip: bool, digests: dict) or None."""
    def _inspect(img):
        p = subprocess.run(
            ["docker", "image", "inspect", img, "--format",
             "{{.Id}}|{{.Architecture}}|{{join .RepoDigests \",\"}}"],
            capture_output=True, text=True, timeout=60)
        if p.returncode != 0:
            return None
        ident, arch, repod = p.stdout.strip().split("|")
        return {"id": ident, "architecture": arch, "repo_digests": repod}

    info = _inspect(IMAGE)
    if info is None:
        rc = subprocess.run(
            ["docker", "build", "--platform", "linux/amd64",
             "-t", IMAGE, os.path.join(HERE, "docker")],
            capture_output=True, timeout=900).returncode
        info = _inspect(IMAGE) if rc == 0 else None
    if info is not None and info["architecture"] == "amd64":
        digests = {"image": IMAGE, **info}
        base = _inspect(BASE_IMAGE)
        if base:
            digests["base_image"] = {"image": BASE_IMAGE, **base}
        return IMAGE, False, digests
    # Fallback: pinned base image + runtime pip numpy (records digest).
    subprocess.run(["docker", "pull", "--platform", "linux/amd64",
                    BASE_IMAGE], capture_output=True, timeout=900)
    base = _inspect(BASE_IMAGE)
    if base is not None and base["architecture"] == "amd64":
        return BASE_IMAGE, True, {"image": BASE_IMAGE, **base,
                                  "runtime_pip_numpy": NUMPY_PIN}
    return None


def docker_leg(workdir, leg, image, runtime_pip, cross, timeout=3600):
    """Run the harness under docker --platform linux/amd64 (x86-64)."""
    outdir = os.path.join(workdir, leg)
    os.makedirs(outdir, exist_ok=True)
    report = os.path.join(outdir, "report.json")
    inner = ("python -u /repo/experiments/conformance/harness.py"
             " --spec /work/blocks/spec.json"
             " --out /work/%s/report.json"
             " --archives-dir /work/%s/archives"
             " --double-encode 2"
             " --cross-decode-archives /work/%s"
             " --cross-decode-manifest /work/%s"
             % (leg, leg,
                os.path.relpath(cross[0], workdir),
                os.path.relpath(cross[1], workdir)))
    if runtime_pip:
        inner = ("pip install -q --no-cache-dir numpy==%s && %s"
                 % (NUMPY_PIN, inner))
    cmd = ["docker", "run", "--rm", "--platform", "linux/amd64",
           "-e", "PYTHONHASHSEED=0",
           "-v", "%s:/repo:ro" % REPO,
           "-v", "%s:/work" % workdir,
           image, "sh", "-c", inner]
    rc = run_logged(cmd, os.path.join(workdir, "logs", leg + ".log"),
                    timeout=timeout)
    if rc != 0:
        raise RuntimeError("docker leg failed (rc=%d); see logs" % rc)
    return report


# ---------------------------------------------------------------------------
# Report comparison
# ---------------------------------------------------------------------------

def load_report(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def blocks_projection(rep):
    """Platform-invariant projection of a report's per-block results."""
    out = []
    for b in rep["blocks"]:
        b = dict(b)
        b.pop("double_encode_ok", None)  # leg-local probe, not identity data
        out.append(b)
    return out


def count_field_matches(ra, rb, field):
    a = {b["name"]: b[field] for b in ra["blocks"]}
    bb = {b["name"]: b[field] for b in rb["blocks"]}
    assert set(a) == set(bb)
    return sum(1 for k in a if a[k] == bb[k]), len(a)


def compare_archive_files(dir_a, dir_b, names):
    same = 0
    for name in names:
        with open(os.path.join(dir_a, name + ".rnr1"), "rb") as f:
            xa = f.read()
        with open(os.path.join(dir_b, name + ".rnr1"), "rb") as f:
            xb = f.read()
        same += int(xa == xb)
    return same, len(names)


# ---------------------------------------------------------------------------
# Edge checks (native, offline analysis -- floats permitted here, R-3.6)
# ---------------------------------------------------------------------------

def check_overflow_guard_trips():
    """K11: the 2^24 guard is live at the exact boundary and under fault."""
    if not __debug__:
        return None, "interpreter running with -O; assertions stripped"
    # (i) exact boundary on the arithmetic-coder guard.
    enc = rnr1.ArithmeticEncoder()
    enc.encode(0, 1, rnr1.PROB_TOTAL_BOUND - 1)  # must be accepted
    tripped_at_bound = False
    try:
        rnr1.ArithmeticEncoder().encode(0, 1, rnr1.PROB_TOTAL_BOUND)
    except AssertionError:
        tripped_at_bound = True
    # (ii) fault injection: inflate the predictor's count scale so that a
    # short repeated-byte input pushes dist() totals past 2^24; the
    # dist() invariant guard must trip (and coding must not proceed).
    old_scale = rnr1.COUNT_SCALE
    tripped_fault = False
    try:
        rnr1.COUNT_SCALE = 1 << 20
        try:
            rnr1.encode_subblock(b"\x00" * 64, 3)
        except AssertionError:
            tripped_fault = True
    finally:
        rnr1.COUNT_SCALE = old_scale
    # (iii) sanity after restore: normal encode still works.
    blob = rnr1.encode_subblock(b"hello world" * 30, 3)
    restored_ok = rnr1.decode_subblock(blob, 330, 3) == b"hello world" * 30
    ok = tripped_at_bound and tripped_fault and restored_ok
    detail = ("boundary accept 2^24-1 / reject 2^24: %s; fault trip: %s; "
              "restored: %s" % (tripped_at_bound, tripped_fault, restored_ok))
    return ok, detail


def check_reduction_order_probe(sample):
    """K12: randomized-order integer reductions match fixed-order exactly."""
    import random
    rng = random.Random(metrics.derive_seed(EXPERIMENT_ID, 0, tag="run"))
    pred = rnr1.NGramPredictor(3)
    probes = 0
    for i, x in enumerate(sample):
        freq, total = pred.dist()
        if i % 97 == 0:
            if freq.dtype != np.int64:
                return False, "freq dtype %s != int64 (R-3.1)" % freq.dtype
            bstar, cand, classid, class_freq = rnr1._partition(freq)
            cf_fixed = [int(v) for v in class_freq]
            tot_fixed = int(freq.sum())
            if tot_fixed != int(total):
                return False, "freq.sum() != incrementally tracked total"
            cum_fixed = np.cumsum(class_freq).tolist()
            for _ in range(3):
                idxs = list(range(256))
                rng.shuffle(idxs)
                acc = [0] * rnr1.NUM_CLASSES
                tot = 0
                for b in idxs:
                    acc[int(classid[b])] += int(freq[b])
                    tot += int(freq[b])
                if acc != cf_fixed or tot != tot_fixed:
                    return False, "shuffled reduction diverged at pos %d" % i
            run_sum, cum_manual = 0, []
            for v in cf_fixed:
                run_sum += v
                cum_manual.append(run_sum)
            if cum_manual != [int(v) for v in cum_fixed]:
                return False, "prefix-sum mismatch at pos %d" % i
            probes += 1
        pred.update(x)
    return True, "%d probe points x 3 shuffles, all exact (int64 path)" % probes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--full", action="store_true",
                    help="100-block campaign corpus (exp-design 8.1)")
    ap.add_argument("--no-docker", action="store_true",
                    help="skip the cross-ISA container leg")
    ap.add_argument("--no-record", action="store_true",
                    help="do not append to the results JSONL store")
    ap.add_argument("--workdir",
                    default=os.path.join(HERE, "workdir"))
    ap.add_argument("--resume", action="store_true",
                    help="crash recovery: reuse legs whose report.json "
                         "already exists in the workdir instead of "
                         "re-running them (never skips the comparisons)")
    args = ap.parse_args(argv)

    t_start = time.time()
    mode = "full" if args.full else "fast"
    workdir = os.path.abspath(args.workdir)
    if os.path.isdir(workdir) and not args.resume:
        shutil.rmtree(workdir)
    os.makedirs(os.path.join(workdir, "logs"), exist_ok=True)

    def leg_done(leg):
        return args.resume and os.path.exists(
            os.path.join(workdir, leg, "report.json"))

    def leg_path(leg):
        return os.path.join(workdir, leg, "report.json")

    # ---- K1: deterministic corpus -------------------------------------
    spec_path = os.path.join(workdir, "blocks", "spec.json")
    if args.resume and os.path.exists(spec_path):
        rc = 0
    else:
        cmd = [PYTHON, "-u", MAKE_CORPUS, "--workdir", workdir]
        if args.full:
            cmd.append("--full")
        rc = run_logged(cmd, os.path.join(workdir, "logs", "corpus.log"))
    ok = rc == 0 and os.path.exists(spec_path)
    with open(spec_path, "rb") as f:
        spec_bytes = f.read()
    spec_sha = hashlib.sha256(spec_bytes).hexdigest()
    spec = json.loads(spec_bytes)
    record("K1 corpus built (%s mode)" % mode, ok,
           "%d blocks, %d bytes, spec sha256=%s"
           % (len(spec), sum(s["n"] for s in spec), spec_sha[:16]))
    if not ok:
        return finish(args, mode, spec_sha, {}, t_start)

    # ---- native legs ---------------------------------------------------
    print("running native legs (A, A2, B, C)...", flush=True)
    rep_a = load_report(leg_path("native_a") if leg_done("native_a")
                        else native_leg(workdir, "native_a", 0))
    rep_a2 = load_report(leg_path("native_a2") if leg_done("native_a2")
                         else native_leg(workdir, "native_a2", 0))
    rep_b = load_report(leg_path("native_b") if leg_done("native_b")
                        else native_leg(workdir, "native_b", 1))
    rep_c = load_report(leg_path("native_c") if leg_done("native_c")
                        else native_leg(workdir, "native_c", None))

    names = [b["name"] for b in rep_a["blocks"]]
    n_blocks = len(names)
    n_positions = sum(b["n"] for b in rep_a["blocks"])
    n_rungs = sum(b["n_subblocks"] for b in rep_a["blocks"])

    # K2: round-trip + verification + encoder/decoder agreement on leg A.
    rt = all(b["roundtrip_ok"] for b in rep_a["blocks"])
    ed_t = all(b["enc_dec_trace_ok"] for b in rep_a["blocks"])
    ed_s = all(b["enc_dec_state_ok"] for b in rep_a["blocks"])
    record("K2 round-trip + H_v verify, all blocks (leg A)", rt and ed_t and ed_s,
           "roundtrip %d/%d; enc==dec trace %s, state %s"
           % (sum(b["roundtrip_ok"] for b in rep_a["blocks"]), n_blocks,
              ed_t, ed_s))

    # K3: same-process double-encode.
    dbl = [b.get("double_encode_ok") for b in rep_a["blocks"]
           if "double_encode_ok" in b]
    record("K3 same-process double-encode identical (H8)",
           len(dbl) > 0 and all(dbl), "%d/%d blocks probed" % (sum(map(bool, dbl)), len(dbl)))

    # K4: across-process identity (identical env).
    same_aa2 = blocks_projection(rep_a) == blocks_projection(rep_a2)
    record("K4 across-process byte-identity A==A2 (H8)", same_aa2,
           "archives chain %s" % rep_a["global"]["archives_chain"][:16])

    # K5: PYTHONHASHSEED invariance.
    same_b = blocks_projection(rep_a) == blocks_projection(rep_b)
    same_c = blocks_projection(rep_a) == blocks_projection(rep_c)
    record("K5 PYTHONHASHSEED invariance A==B==C (H8)", same_b and same_c,
           "seed 0 vs 1 vs unset")

    # ---- docker cross-ISA leg -------------------------------------------
    counts = {
        "n_blocks": n_blocks, "n_positions": n_positions, "n_rungs": n_rungs,
        "native_process_identity": same_aa2, "hashseed_invariance":
        same_b and same_c, "roundtrip": rt,
        "arch_native": rep_a["env"]["machine"],
    }
    docker_ok = None
    digests = {}
    if args.no_docker:
        why = "--no-docker given"
    elif not docker_available():
        why = "docker daemon unavailable"
    else:
        why = None
        img = ensure_docker_image()
        if img is None:
            why = "no amd64 python+numpy image obtainable"
    if why is not None:
        for k in ("K6 cross-ISA archive byte-identity (H8)",
                  "K7 cross-ISA inference-trace identity (H14)",
                  "K8 cross-ISA adaptive-state ladder identity (H9)",
                  "K9 cross-decode arm64<->x86-64 (H8)"):
            skip(k, why)
    else:
        image, runtime_pip, digests = img
        if leg_done("docker_d"):
            rep_d = load_report(leg_path("docker_d"))
        else:
            print("running docker leg D (image %s, runtime_pip=%s)..."
                  % (image, runtime_pip), flush=True)
            rep_d = load_report(docker_leg(
                workdir, "docker_d", image, runtime_pip,
                cross=(os.path.join(workdir, "native_a", "archives"),
                       os.path.join(workdir, "native_a", "report.json"))))
        rep_e = load_report(
            leg_path("native_e") if leg_done("native_e")
            else native_decode_leg(
                workdir, "native_e",
                os.path.join(workdir, "docker_d", "archives"),
                os.path.join(workdir, "docker_d", "report.json")))

        arch_pair = "%s->%s" % (rep_a["env"]["machine"], rep_d["env"]["machine"])
        counts["arch_docker"] = rep_d["env"]["machine"]
        counts["docker_env"] = rep_d["env"]

        # K6: archives byte-identical (compare actual file bytes).
        n_same_arc, _ = count_field_matches(rep_a, rep_d, "archive_sha256")
        n_same_bytes, n_tot = compare_archive_files(
            os.path.join(workdir, "native_a", "archives"),
            os.path.join(workdir, "docker_d", "archives"), names)
        counts["cross_isa_identical_archives"] = n_same_bytes
        record("K6 cross-ISA archive byte-identity (H8)",
               n_same_arc == n_tot and n_same_bytes == n_tot,
               "%s: %d/%d archives byte-identical (file compare)"
               % (arch_pair, n_same_bytes, n_tot))

        # K7: full integer inference trace identical (H14).
        n_same_tr, _ = count_field_matches(rep_a, rep_d, "enc_traces")
        chain_eq = (rep_a["global"]["traces_chain"]
                    == rep_d["global"]["traces_chain"])
        counts["cross_isa_identical_traces"] = n_same_tr
        record("K7 cross-ISA inference-trace identity (H14)",
               n_same_tr == n_tot and chain_eq,
               "%d/%d blocks, %d positions traced, chain %s"
               % (n_same_tr, n_tot, n_positions,
                  rep_d["global"]["traces_chain"][:16]))

        # K8: adaptive-state ladder identical (H9).
        n_same_lad, _ = count_field_matches(rep_a, rep_d, "enc_state_ladder")
        ed_t_d = all(b["enc_dec_trace_ok"] for b in rep_d["blocks"])
        ed_s_d = all(b["enc_dec_state_ok"] for b in rep_d["blocks"])
        full_proj = blocks_projection(rep_a) == blocks_projection(rep_d)
        counts["cross_isa_identical_ladders"] = n_same_lad
        counts["cross_isa_full_projection_equal"] = full_proj
        record("K8 cross-ISA adaptive-state ladder identity (H9)",
               n_same_lad == n_tot and ed_t_d and ed_s_d and full_proj,
               "%d/%d blocks, %d ladder rungs; enc==dec on x86-64: %s"
               % (n_same_lad, n_tot, n_rungs, ed_t_d and ed_s_d))

        # K9: cross-decode both directions.
        xd = rep_d.get("cross_decode", {"n": 0, "n_ok": -1})
        xe = rep_e.get("cross_decode", {"n": 0, "n_ok": -1})
        counts["cross_decode_d"] = (xd["n_ok"], xd["n"])
        counts["cross_decode_e"] = (xe["n_ok"], xe["n"])
        record("K9 cross-decode arm64<->x86-64 (H8)",
               xd["n_ok"] == xd["n"] == n_tot and xe["n_ok"] == xe["n"] == n_tot,
               "x86-64 decodes native: %d/%d; native decodes x86-64: %d/%d"
               % (xd["n_ok"], xd["n"], xe["n_ok"], xe["n"]))
        docker_ok = all(r[1] for r in RESULTS[-4:])

    # ---- K10/K11/K12 edge checks ---------------------------------------
    stress = next(b for b in rep_a["blocks"] if b["name"] == "stress_cap_repeat")
    headroom = rnr1.PROB_TOTAL_BOUND / stress["max_total"]
    counts["stress_max_total"] = stress["max_total"]
    record("K10 overflow-guard headroom at COUNT_CAP regime (R-3.2)",
           stress["max_total"] == EXPECTED_STRESS_MAX_TOTAL
           and stress["max_total"] < rnr1.PROB_TOTAL_BOUND,
           "max_total=%d == 256+32*(2^16-1)=%d; bound 2^24=%d; headroom %.2fx"
           % (stress["max_total"], EXPECTED_STRESS_MAX_TOTAL,
              rnr1.PROB_TOTAL_BOUND, headroom))

    trip = check_overflow_guard_trips()
    if trip[0] is None:
        skip("K11 overflow guard trips at 2^24 (R-3.2/R-5.2)", trip[1])
    else:
        record("K11 overflow guard trips at 2^24 (R-3.2/R-5.2)", trip[0], trip[1])
        counts["overflow_guard_trips"] = trip[0]

    sample_item = next(s for s in spec if s["source"].startswith("enwik8"))
    with open(os.path.join(workdir, "blocks", sample_item["file"]), "rb") as f:
        sample = f.read()[:2048]
    red = check_reduction_order_probe(sample)
    record("K12 reduction-order independence probe (R-3.4)", red[0], red[1])
    counts["reduction_probe"] = red[0]

    return finish(args, mode, spec_sha, counts, t_start,
                  docker_ok=docker_ok, digests=digests)


def finish(args, mode, spec_sha, counts, t_start, docker_ok=None, digests=None):
    n_pass = sum(1 for _, ok in RESULTS if ok is True)
    n_fail = sum(1 for _, ok in RESULTS if ok is False)
    n_skip = sum(1 for _, ok in RESULTS if ok is None)
    all_ok = n_fail == 0
    print("OVERALL: %s -- %d pass, %d fail, %d skip (%.1fs, %s mode)"
          % ("PASS" if all_ok else "FAIL", n_pass, n_fail, n_skip,
             time.time() - t_start, mode), flush=True)

    if not args.no_record and counts:
        write_records(mode, spec_sha, counts, docker_ok, digests or {})
    return 0 if all_ok else 1


def write_records(mode, spec_sha, counts, docker_ok, digests):
    """Append REPRO-family records (exact counts, METHODS.md section 5)."""
    store = os.path.join(HERE, "results", "conformance.jsonl")
    corpus_tag = "conformance-%s" % mode
    config = {"mode": mode, "spec_sha256": spec_sha,
              "n_blocks": counts.get("n_blocks"),
              "W": "per-block", "K": "per-block"}
    seed = metrics.derive_seed(EXPERIMENT_ID, 0, tag="run")
    n = counts.get("n_blocks", 0)
    legs_native = ["arm64-native-hashseed0", "arm64-native-hashseed0-rerun",
                   "arm64-native-hashseed1", "arm64-native-hashseed-unset"]
    common = dict(family="REPRO", experiment=EXPERIMENT_ID,
                  docker=digests or None,
                  statistics_note=("REPRO family: binary conformance, exact "
                                   "counts, no p-values/CI per METHODS.md "
                                   "section 5"))

    def rec(hypothesis, measured, legs):
        results_store.append(store, results_store.make_record(
            corpus_tag, "rnr1-ref-1.0.0", config, seed, measured,
            hypothesis=hypothesis, legs=legs, **common))

    cross = counts.get("cross_isa_identical_archives")
    rec("H8", {
        "total_blocks": n,
        "native_process_identity": counts.get("native_process_identity"),
        "hashseed_invariance": counts.get("hashseed_invariance"),
        "roundtrip_all": counts.get("roundtrip"),
        "cross_isa_identical_archives": cross,
        "cross_decode_x86_of_native": counts.get("cross_decode_d"),
        "cross_decode_native_of_x86": counts.get("cross_decode_e"),
        "byte_identity_rate": (None if cross is None else cross / max(1, n)),
    }, legs_native + (["amd64-linux-docker-emulated"] if cross is not None else []))
    rec("H9", {
        "total_blocks": n,
        "ladder_rungs": counts.get("n_rungs"),
        "cross_isa_identical_ladders": counts.get("cross_isa_identical_ladders"),
        "encoder_decoder_state_agreement": counts.get("roundtrip"),
    }, legs_native + (["amd64-linux-docker-emulated"]
                      if counts.get("cross_isa_identical_ladders") is not None
                      else []))
    rec("H14", {
        "total_blocks": n,
        "trace_positions": counts.get("n_positions"),
        "cross_isa_identical_traces": counts.get("cross_isa_identical_traces"),
        "identity_rate": (None if counts.get("cross_isa_identical_traces") is None
                          else counts["cross_isa_identical_traces"] / max(1, n)),
    }, legs_native + (["amd64-linux-docker-emulated"]
                      if counts.get("cross_isa_identical_traces") is not None
                      else []))
    rec("R-3.2-overflow", {
        "stress_max_total": counts.get("stress_max_total"),
        "bound": rnr1.PROB_TOTAL_BOUND,
        "headroom_x": (None if not counts.get("stress_max_total")
                       else rnr1.PROB_TOTAL_BOUND / counts["stress_max_total"]),
        "guard_trips": counts.get("overflow_guard_trips"),
    }, ["arm64-native-hashseed0"])
    rec("R-3.4-reduction-order", {
        "probe_ok": counts.get("reduction_probe"),
        "covered_by_cross_isa_trace": docker_ok,
    }, ["arm64-native-hashseed0"])
    print("results appended: %s" % store, flush=True)


if __name__ == "__main__":
    sys.exit(main())

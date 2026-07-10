#!/usr/bin/env python3
"""run_checks.py -- PASS/FAIL gate for the impl_cm deliverable (Task: CM
predictor + H1 first measurement).

Gates
-----
  G1  integer-table sanity           squash/stretch monotone, in-range,
                                     inverse-consistent
  G2  round-trip bit-exactness       CM pack -> unpack == input, on text,
                                     mixed binary, and incompressible
                                     (store-mode) slices; random access
  G3  two-process determinism        the CLI packs the same input to the
                                     byte-identical archive in two fresh
                                     interpreter processes (and matches
                                     the in-process archive)
  G4  CM beats n-gram                per-corpus model cross-entropy of
                                     the CM predictor strictly below the
                                     order-W n-gram's on EVERY text/code
                                     corpus of the H1 suite
  G5  coder redundancy < 0.01 bpb    measured repair bits minus model
                                     cross-entropy, per corpus and coder

Default scale: reduced (check slices; < 10 min total).
--full: gates G4/G5 evaluate on the exp-design-sized 'full' campaign.
Exit code 0 iff every gate PASSes.

Usage:  python -u run_checks.py [--full]
"""

import argparse
import hashlib
import subprocess
import sys
import tempfile
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "data"))
sys.path.insert(0, str(_HERE.parent / "impl"))

import loader  # noqa: E402
import rnr1  # noqa: E402
import cm_predictor  # noqa: E402
from cm_coder import pack_cm, read_cm, unpack_cm  # noqa: E402
import h1_measure  # noqa: E402

PY = sys.executable
RESULTS = []


def gate(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", name,
                         (" -- " + detail) if detail else ""), flush=True)
    return bool(ok)


def g1_tables():
    try:
        ok = cm_predictor.table_self_test()
        return gate("G1 table sanity", ok)
    except AssertionError as exc:
        return gate("G1 table sanity", False, str(exc))


def g2_roundtrip():
    text = loader.read_range("enwik8", 1 << 20, 96 * 1024)
    mixed = loader.read_range("calgary", 100 * 1024, 64 * 1024)
    rand = loader.read_range("ctrl_urandom", 0, 32 * 1024)
    ok = True
    for label, data, K in (("text", text, 16384), ("mixed", mixed, 16384),
                           ("incompressible", rand, 16384)):
        raw = pack_cm(data, W=3, K=K)
        back = unpack_cm(raw)
        good = back == data
        ok &= good
        print("    round-trip %-15s %6d bytes -> %6d bytes  %s"
              % (label, len(data), len(raw), "ok" if good else "MISMATCH"),
              flush=True)
    # store-mode fail-safe must engage on incompressible input
    raw = pack_cm(rand, W=3, K=16384)
    ok &= len(raw) <= len(rand) + 2048  # header+index margin only
    # random access
    raw = pack_cm(text, W=3, K=16384)
    piece, stats = read_cm(raw, 40000, 1234)
    ra_ok = piece == text[40000:41234] and stats["decoded_bytes"] <= \
        16384 + 1234 + 16384
    ok &= ra_ok
    print("    random access ok=%s warmup=%d decoded=%d"
          % (ra_ok, stats["warmup_bytes"], stats["decoded_bytes"]),
          flush=True)
    return gate("G2 round-trip bit-exactness (+RA, +store-mode)", ok)


def g3_two_process():
    data = loader.read_range("enwik8", 1 << 20, 128 * 1024)
    with tempfile.TemporaryDirectory() as td:
        src = Path(td) / "in.bin"
        src.write_bytes(data)
        digests = []
        for i in (1, 2):
            out = Path(td) / ("out%d.rnr1" % i)
            proc = subprocess.run(
                [PY, "-u", str(_HERE / "cm_coder.py"), "pack", str(src),
                 str(out), "--W", "3", "--K", "16384"],
                capture_output=True, text=True, timeout=1200)
            if proc.returncode != 0:
                return gate("G3 two-process determinism", False,
                            proc.stderr.strip()[:200])
            digests.append(hashlib.sha256(out.read_bytes()).hexdigest())
        inproc = hashlib.sha256(
            pack_cm(data, W=3, K=16384)).hexdigest()
        ok = digests[0] == digests[1] == inproc
        return gate("G3 two-process determinism",
                    ok, "sha256 %s" % digests[0][:16])


def g4_g5_measurement(scale):
    summary = h1_measure.run_campaign(
        scale, jsonl_path=_HERE / "results" / ("h1_%s.jsonl" % scale))
    ok4, ok5 = True, True
    for corpus, centry in summary["corpora"].items():
        cm = centry["coders"]["rnr1-cm"]["measured"]
        ng = centry["coders"]["rnr1-ngram"]["measured"]
        beats = cm["model_ce_bpb"] < ng["model_ce_bpb"]
        ok4 &= beats
        print("    %-15s CM CE %.4f  vs ngram CE %.4f  %s"
              % (corpus, cm["model_ce_bpb"], ng["model_ce_bpb"],
                 "ok" if beats else "NOT BETTER"), flush=True)
        for coder in ("rnr1-cm", "rnr1-ngram"):
            red = centry["coders"][coder]["measured"]["coder_redundancy_bpb"]
            good = 0.0 <= red < 0.01
            ok5 &= good
            if not good:
                print("    %-15s %s redundancy %.5f bpb OUT OF BOUND"
                      % (corpus, coder, red), flush=True)
        ok4 &= cm["roundtrip_bit_exact"] and ng["roundtrip_bit_exact"]
    gate("G4 CM strictly beats n-gram cross-entropy on every corpus", ok4)
    gate("G5 coder redundancy < 0.01 bpb (all rnr runs)", ok5)
    return ok4 and ok5


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--full", action="store_true",
                    help="evaluate G4/G5 on the exp-design-sized campaign")
    args = ap.parse_args(argv)
    t0 = time.perf_counter()
    g1_tables()
    g2_roundtrip()
    g3_two_process()
    g4_g5_measurement("full" if args.full else "check")
    ok = all(r for _, r in RESULTS)
    print("=" * 60, flush=True)
    for name, r in RESULTS:
        print("  %-58s %s" % (name, "PASS" if r else "FAIL"))
    print("OVERALL: %s (%.1f s)" % ("PASS" if ok else "FAIL",
                                    time.perf_counter() - t0), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

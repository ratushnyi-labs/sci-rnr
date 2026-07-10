#!/usr/bin/env python3
"""Gate for the Part II deviation-hierarchy simulation studies.

Runs all five studies (a)-(e) against Theorems 7.31-7.33, 7.35, 7.36 of the
random-access paper at reduced scale (default, < 10 min) or full scale
(--full).  Prints one PASS/FAIL line per check and OVERALL at the end, in the
style of scripts/verify/.

Usage:
    python run_checks.py [--full]

Each study is also runnable standalone:
    python study_a_overflow_ldp.py [--full]     # etc.
"""
import argparse
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import study_a_overflow_ldp
import study_b_universality
import study_c_moderate_fade
import study_d_query_evt
import study_e_buffer_fclt

STUDIES = [
    ("A (overflow LDP, Thm 7.31)", study_a_overflow_ldp),
    ("B (universality dividend, Thm 7.32)", study_b_universality),
    ("C (moderate-deviation RA fade, Thm 7.33)", study_c_moderate_fade),
    ("D (worst-case query EVT, Thm 7.35)", study_d_query_evt),
    ("E (buffer functional CLT, Thm 7.36)", study_e_buffer_fclt),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="full-scale runs (larger ladders / sample sizes)")
    args = ap.parse_args()

    t0 = time.time()
    all_results = []
    for name, mod in STUDIES:
        t1 = time.time()
        try:
            res = mod.run(full=args.full)
        except Exception:
            traceback.print_exc()
            res = [(f"{name} crashed", False, "unhandled exception")]
        all_results.extend(res)
        print(f"-- study {name}: {time.time() - t1:.1f}s\n", flush=True)

    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    n_pass = 0
    for tag, ok, _msg in all_results:
        print(f"[{'PASS' if ok else 'FAIL'}] {tag}")
        n_pass += bool(ok)
    n = len(all_results)
    ok_all = n_pass == n and n > 0
    print("=" * 72)
    print(f"{n_pass}/{n} checks passed in {time.time() - t0:.1f}s "
          f"({'full' if args.full else 'reduced'} scale)")
    print("OVERALL:", "PASS" if ok_all else "FAIL")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()

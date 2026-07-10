#!/usr/bin/env python3
"""
run_checks.py
============================================================================
Gate for the Part III dispersion simulation studies (experiments/
part3_dispersion/).  Runs studies (a)-(d) at reduced scale (< 10 min
total; the individual studies accept --full for the full-scale configs)
and prints PASS/FAIL lines per check plus an OVERALL verdict, in the
style of scripts/verify/.

  (a) study_a_finite_n_dispersion   -- Var(j_n)/n -> V_lossless on Gray,
                                       exact + Monte Carlo + rate fit
                                       (Remark 7.34b).
  (b) study_b_third_order_recentring-- j_n = i_n - nc recentring and the
                                       +1/2 log n third-order drift
                                       (Remark 7.34n').
  (c) study_c_beyond_gray_law       -- D_c^(n) - D_c = K(p)^2/n^2 with the
                                       closed-form K(p) (Remark 7.34o/o').
  (d) study_d_replica_margins       -- R_A(s;D) - 3/2 > 0 on Gray, A=2..5,
                                       corner law (A-2)/(A-1) p
                                       (Theorem 7.34m / Remark 7.34m').

Usage:
  /Users/para/.venvs/rnr/bin/python -u run_checks.py [--full]
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import study_a_finite_n_dispersion as study_a
import study_b_third_order_recentring as study_b
import study_c_beyond_gray_law as study_c
import study_d_replica_margins as study_d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="full-scale configs (longer ladders, more samples)")
    args = ap.parse_args()

    t0 = time.time()
    checkers = []
    for name, mod in (("a", study_a), ("b", study_b),
                      ("c", study_c), ("d", study_d)):
        ts = time.time()
        ck = mod.run(full=args.full)
        checkers.append((name, ck))
        print(f"  [study ({name}) wall time: {time.time() - ts:.1f} s]")

    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    total = ok_count = 0
    for name, ck in checkers:
        for check_name, ok in ck.results:
            total += 1
            ok_count += ok
        verdict = "PASS" if ck.ok else "FAIL"
        print(f"  study ({name}): "
              f"{sum(o for _, o in ck.results)}/{len(ck.results)} checks "
              f"{verdict}")
    overall = all(ck.ok for _, ck in checkers)
    print(f"  {ok_count}/{total} checks passed; wall time "
          f"{time.time() - t0:.1f} s")
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if overall else 'FAIL'}")
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())

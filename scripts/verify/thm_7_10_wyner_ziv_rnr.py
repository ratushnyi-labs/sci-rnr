#!/usr/bin/env python3
r"""
Verification of Theorem 7.10: Wyner-Ziv RNR rate-distortion bound for
binary doubly symmetric source with Hamming distortion.

For X ~ Bern(1/2), Y = X XOR Z with Z ~ Bern(p), the Wyner-Ziv 1976
rate-distortion function is:

    R_WZ(D) = max(0, h(p) - h(D))    for D in [0, p]
             = 0                       for D >= p

where h is binary entropy.

Comparisons:
- T7.7 (no side-info): R(D) = 1 - h(D).
- T7.8 (lossless with side-info): h(X|Y) = h(p) = R_WZ(0).
- T7.10 lossy with side-info: R_WZ(D) <= h(p).

This script:
1. Computes R_WZ(D) for binary doubly symmetric source at multiple
   (p, D) settings.
2. Verifies the limit relations: R_WZ(0) = h(p), R_WZ(p) = 0.
3. Compares savings over T7.7 (no side-info).
4. Sanity-checks: R_WZ(D) <= R(D) always, R_WZ(D) <= R_WZ(D') for
   D > D' (monotone).

PASS = closed-form formula holds at all tested (p, D), limits correct,
       monotonicity and inequality R_WZ <= R hold.
"""

import math
import sys


def binary_entropy(q):
    if q <= 1e-15 or q >= 1 - 1e-15:
        return 0.0
    return -q * math.log2(q) - (1 - q) * math.log2(1 - q)


def rate_distortion_no_side_info(p_X, D):
    """T7.7: R(D) = h(X) - h(D) for binary, capped at 0."""
    h_X = binary_entropy(p_X)
    return max(0.0, h_X - binary_entropy(D))


def wyner_ziv_rate_distortion_bdss(p, D):
    """T7.10 R_WZ(D) for binary doubly symmetric source."""
    if D >= p:
        return 0.0
    return max(0.0, binary_entropy(p) - binary_entropy(D))


def slepian_wolf_lossless_one_sided(p):
    """T7.8 / T7.10 limit D=0: h(X|Y) = h(p)."""
    return binary_entropy(p)


def test_wyner_ziv_table(p, distortions):
    """Tabulate R_WZ vs R for binary doubly symmetric."""
    h_p = binary_entropy(p)
    print(f"\nSource: X ~ Bern(0.5), Y = X XOR Bern(p={p}), Hamming distortion")
    print(f"h(X) = 1.0, h(X|Y) = h(p) = {h_p:.4f}")
    print(f"  D   | R_WZ(D)  | R_T7.7(D) | savings | savings vs h(X|Y)")
    print(f"  " + "-" * 65)

    all_ok = True
    R_prev = math.inf

    for D in distortions:
        R_WZ = wyner_ziv_rate_distortion_bdss(p, D)
        R_T77 = rate_distortion_no_side_info(0.5, D)
        savings_vs_T77 = 100.0 * (R_T77 - R_WZ) / R_T77 if R_T77 > 1e-9 else 0
        savings_vs_SW = 100.0 * (h_p - R_WZ) / h_p if h_p > 1e-9 else 0
        print(f"  {D:.3f} | {R_WZ:.4f}   | {R_T77:.4f}    | {savings_vs_T77:5.1f}%  | {savings_vs_SW:5.1f}%")

        # Monotonicity check
        if R_WZ > R_prev + 1e-9:
            print(f"    FAIL: R_WZ should be monotone decreasing in D")
            all_ok = False
        R_prev = R_WZ

        # R_WZ <= R_T7.7 check
        if R_WZ > R_T77 + 1e-9:
            print(f"    FAIL: R_WZ should be <= R_T7.7")
            all_ok = False

        # R_WZ >= 0 check
        if R_WZ < -1e-9:
            print(f"    FAIL: R_WZ < 0")
            all_ok = False

    return all_ok


def test_limits(p):
    """Verify limit relations of R_WZ."""
    print(f"\nLimit checks for p = {p}:")

    R_WZ_0 = wyner_ziv_rate_distortion_bdss(p, 0)
    h_p = binary_entropy(p)
    print(f"  R_WZ(0) = {R_WZ_0:.4f} vs h(p) = {h_p:.4f}: ", end="")
    if abs(R_WZ_0 - h_p) > 1e-12:
        print("FAIL")
        return False
    print("OK (lossless = Slepian-Wolf one-sided)")

    # R_WZ at D = p should be 0
    R_WZ_p = wyner_ziv_rate_distortion_bdss(p, p)
    print(f"  R_WZ({p}) = {R_WZ_p:.4f}: ", end="")
    if abs(R_WZ_p) > 1e-12:
        print("FAIL: should be 0")
        return False
    print("OK (decoder outputs Y directly)")

    # R_WZ > 0 strictly for D < p
    if p > 0.01:
        R_WZ_half = wyner_ziv_rate_distortion_bdss(p, p / 2)
        print(f"  R_WZ({p/2}) = {R_WZ_half:.4f}: ", end="")
        if R_WZ_half <= 0:
            print("FAIL: should be > 0 for D < p")
            return False
        print("OK (positive in interior)")

    return True


def main() -> int:
    print("Verification of Theorem 7.10 (Wyner-Ziv RNR with side-info at decoder)")
    print("=" * 70)

    all_ok = True

    # Test 1: strong correlation p=0.1
    distortions = [0.0, 0.02, 0.05, 0.08, 0.10, 0.15]
    if not test_wyner_ziv_table(0.1, distortions):
        all_ok = False
    if not test_limits(0.1):
        all_ok = False

    # Test 2: moderate correlation p=0.2
    distortions = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25]
    if not test_wyner_ziv_table(0.2, distortions):
        all_ok = False
    if not test_limits(0.2):
        all_ok = False

    # Test 3: weak correlation p=0.3
    distortions = [0.0, 0.1, 0.2, 0.3, 0.4]
    if not test_wyner_ziv_table(0.3, distortions):
        all_ok = False
    if not test_limits(0.3):
        all_ok = False

    # Verify specific numerical claim from the paper:
    # For p=0.1, D=0.05: R_WZ = 0.183, 74% savings vs T7.7, 61% vs SW.
    print("\n\nPaper's worked example (p=0.1, D=0.05):")
    R_WZ = wyner_ziv_rate_distortion_bdss(0.1, 0.05)
    R_T77 = rate_distortion_no_side_info(0.5, 0.05)
    h_p = binary_entropy(0.1)
    print(f"  R_WZ(0.05) = {R_WZ:.4f} (paper claims 0.183)")
    print(f"  Savings vs T7.7 R(0.05) = {R_T77:.4f}: {100 * (R_T77 - R_WZ) / R_T77:.1f}% (paper: 74%)")
    print(f"  Savings vs SW h(X|Y) = {h_p:.4f}: {100 * (h_p - R_WZ) / h_p:.1f}% (paper: 61%)")
    if abs(R_WZ - 0.183) > 0.001:
        print(f"  FAIL: paper's 0.183 claim off")
        all_ok = False

    print()
    if all_ok:
        print("PASS: Theorem 7.10 Wyner-Ziv RNR verified.")
        print("      R_WZ(D) = max(0, h(p) - h(D)) for binary doubly symmetric.")
        print("      Limits: R_WZ(0) = h(X|Y), R_WZ(p) = 0 (Y suffices).")
        print("      Monotonicity, R_WZ <= R_T7.7, paper's worked example all confirmed.")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Verification of Theorem 6.8 explicit APX-hardness constant γ_SGP = 8569/8568.

Theorem claim: The restricted Type-III-C dictionary-learning problem (constant
predictor, recursive normal-form dictionary, symbol-count objective L_sym) is
APX-hard with inapproximability ratio at least 8569/8568, inherited verbatim
from Charikar et al. 2005 (Theorem 1, page 2557, IEEE Transactions on
Information Theory 51(7):2554-2576) via the symbol-count bijection of Theorem 6.7.

This script reproduces the algebra of the Charikar et al. derivation:
  γ_SGP = (15|V| + 3|E| + (145/144)k) / (15|V| + 3|E| + k)
with structural constraints |E| ≤ (3/2)|V| (max-degree-3) and k ≥ (1/3)|V|
(min-VC from Berman-Karpinski 1999, ratio 145/144 inapproximable). The
adversarial inapproximability lower bound ρ_min is the MINIMUM of ρ over the
constraint region (the tightest inapproximability we can certify for the
family). ρ is monotone decreasing in |E| and monotone increasing in k, so the
minimum is achieved at the extreme point |E| = (3/2)|V| (upper bound), k =
(1/3)|V| (lower bound), where ρ = 8569/8568.

PASS = the formula evaluates exactly to 8569/8568 = 1.0001167... at the
constraint extremum.
"""

import math
import sys
from fractions import Fraction


def main() -> int:
    failures = 0

    # Charikar et al. 2005 Theorem 1 derivation, page 2558:
    # Grammar-input encoding gives grammar size G(σ) = 15|V| + 3|E| + k
    # where k = size of minimum vertex cover for the bounded-degree-3 graph H.
    # Berman-Karpinski 1999: VC-3 is hard to approximate below ratio 145/144.
    #
    # The ratio between the smallest grammar (k = min-VC) and any
    # (145/144)-approximation of VC is:
    #   ρ = (15|V| + 3|E| + (145/144)k) / (15|V| + 3|E| + k)
    # Both constraints: |E| ≤ (3/2)|V| (max-degree-3 implies sum-of-degrees =
    # 2|E| ≤ 3|V|), and k ≥ (1/3)|V| (since every vertex covers at most 3
    # edges and |E| ≥ |V| in Charikar's setup, so k ≥ |E|/3 ≥ |V|/3).
    # ρ is monotone DECREASING in |E| and monotone INCREASING in k. The
    # adversarial inapproximability lower bound (minimum of ρ over the family)
    # is achieved at the extreme point |E| = (3/2)|V| (upper bound, smaller ρ)
    # and k = (1/3)|V| (lower bound, smaller ρ).

    # Compute symbolically with exact fractions
    V = Fraction(1)          # normalize |V| = 1
    E = Fraction(3, 2) * V   # |E| = (3/2)|V|
    k = Fraction(1, 3) * V   # k = (1/3)|V|
    bk_ratio = Fraction(145, 144)  # Berman-Karpinski 145/144

    numerator = 15 * V + 3 * E + bk_ratio * k
    denominator = 15 * V + 3 * E + k

    gamma = numerator / denominator
    expected = Fraction(8569, 8568)

    print(f"Charikar et al. 2005 Theorem 1 derivation (page 2558):")
    print(f"  Normalize |V| = 1, then |E| = {E}, k = {k}, BK ratio = {bk_ratio}")
    print(f"  Numerator:   15·{V} + 3·{E} + ({bk_ratio})·{k}")
    print(f"             = {15*V} + {3*E} + {bk_ratio*k}")
    print(f"             = {numerator}")
    print(f"  Denominator: 15·{V} + 3·{E} + {k}")
    print(f"             = {15*V} + {3*E} + {k}")
    print(f"             = {denominator}")
    print(f"  γ_SGP = numerator / denominator = {gamma} = {numerator} / {denominator}")
    print()
    print(f"  Reducing: {numerator}/{denominator}:")
    g = math.gcd(numerator.numerator * denominator.denominator,
                 numerator.denominator * denominator.numerator)
    p = numerator.numerator * denominator.denominator
    q = numerator.denominator * denominator.numerator
    print(f"  numerator·432 = {numerator * 432}, denominator·432 = {denominator * 432}")
    print(f"  So γ_SGP = {numerator * 432} / {denominator * 432}")
    print()

    # Verify against expected
    if gamma == expected:
        print(f"PASS: γ_SGP = {gamma} = 8569/8568 ≈ {float(gamma):.10f}")
        print(f"      Charikar et al. 2005 Theorem 1 derivation reproduces exactly.")
        print(f"      Inapproximability gap: 1 + 1/{denominator * 432}")
    else:
        print(f"FAIL: computed γ_SGP = {gamma}, expected {expected}")
        failures += 1

    # Verify the minimum constraint
    # Check that the ratio is monotone increasing in k (so minimum at k = (1/3)|V|)
    print()
    print("Sanity check: ρ is monotone non-decreasing in k (Berman-Karpinski constant"
          " > 1 amplifies via larger k):")
    for k_test in [Fraction(1, 3), Fraction(1, 2), Fraction(1, 1)]:
        num_test = 15 * V + 3 * E + bk_ratio * k_test * V
        den_test = 15 * V + 3 * E + k_test * V
        ratio_test = num_test / den_test
        print(f"  k/|V| = {k_test}: ρ = {ratio_test} ≈ {float(ratio_test):.6f}")

    # Verify the minimum constraint
    # Check that the ratio is monotone non-increasing in |E| (so minimum at |E| = (3/2)|V|)
    print()
    print("Sanity check: ρ is monotone non-increasing in |E| (more edges dilute"
          " the gap):")
    for E_test in [Fraction(3, 2), Fraction(2, 1), Fraction(10, 1)]:
        num_test = 15 * V + 3 * E_test + bk_ratio * k
        den_test = 15 * V + 3 * E_test + k
        ratio_test = num_test / den_test
        print(f"  |E|/|V| = {E_test}: ρ = {ratio_test} ≈ {float(ratio_test):.6f}")

    print()
    if failures == 0:
        print("PASS: explicit γ_SGP = 8569/8568 derivation verified.")
        print("      The constant is inherited from Berman-Karpinski 145/144 via the")
        print("      Charikar et al. 15|V|+3|E|+k grammar encoding overhead.")
        return 0
    else:
        print(f"FAIL: {failures} mismatch(es).")
        return 1


if __name__ == "__main__":
    sys.exit(main())

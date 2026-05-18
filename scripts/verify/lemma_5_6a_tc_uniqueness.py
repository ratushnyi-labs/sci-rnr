#!/usr/bin/env python3
"""
Verification of Lemma 5.6a: TC(D_N) = Θ(N) for the uniform-over-distinct-PK
distribution.

Setting: D_N is uniform over the K!/(K-N)! injections from [N] (rows) to [K]
(possible primary-key values), where K = α·N for some constant α > 1.

Claim: TC(D_N) = N · c(α) + O(log N), where
  c(α) = log_2 e - (α-1) · log_2(α/(α-1)).

The coefficient c(α) is monotone DECREASING on (1, ∞):
  lim_{α → 1+} c(α) = log_2 e ≈ 1.4427 (tightest, K = N+1)
  c(1.5) ≈ 0.650, c(2) ≈ 0.4427, c(3) ≈ 0.273, c(10) ≈ 0.075
  c(α) → 0 as α → ∞ (UNIQUE constraint becomes vacuous when K >> N)

Verification: compute the exact joint H(D_N) = log_2(K!/(K-N)!) and exact sum
of marginal entropies N · log_2(K), then compare the residual to the formula.

PASS = the asymptotic formula matches within O(1) at moderate N.
"""

import math
import sys


def log_factorial(n: int) -> float:
    """log_2(n!) via math.lgamma for large n."""
    if n <= 1:
        return 0.0
    return math.lgamma(n + 1) / math.log(2)


def joint_H(N: int, K: int) -> float:
    """H(D_N) = log_2(K!/(K-N)!) for uniform-over-distinct-PK distribution."""
    return log_factorial(K) - log_factorial(K - N)


def sum_marginal_H(N: int, K: int) -> float:
    """Sum of per-row marginal entropies. Each row is uniform over K by symmetry."""
    return N * math.log2(K)


def TC(N: int, K: int) -> float:
    """TC(D_N) = sum of marginals - joint H."""
    return sum_marginal_H(N, K) - joint_H(N, K)


def predicted_c(alpha: float) -> float:
    """c(α) from the formula."""
    if alpha == 1.0:
        return float("inf")  # boundary
    return math.log2(math.e) - (alpha - 1) * math.log2(alpha / (alpha - 1))


def main() -> int:
    failures = 0

    # Test multiple α values and N values
    test_alphas = [2, 3, 4, 5, 10]
    test_Ns = [50, 100, 500, 1000]

    print(f"{'α':>3} {'N':>5} {'TC(D_N)':>12} {'N·c(α)':>12} "
          f"{'residual':>10} {'predicted c(α)':>14}")
    print("-" * 65)

    for alpha in test_alphas:
        c_alpha = predicted_c(float(alpha))
        for N in test_Ns:
            K = alpha * N
            tc = TC(N, K)
            predicted = N * c_alpha
            residual = tc - predicted
            print(f"{alpha:>3} {N:>5} {tc:>12.4f} {predicted:>12.4f} "
                  f"{residual:>10.4f} {c_alpha:>14.6f}")

            # Tolerance: residual should be O(1) (i.e., small relative to N · c(α))
            # As N → ∞, residual stays bounded
            if abs(residual) > 5.0:  # tolerance for finite-N effects
                print(f"  WARNING: residual {residual:.4f} exceeds 5.0 at α={alpha}, N={N}")
                failures += 1

    print()
    print("Special cases:")
    # α = 2: c ≈ 0.4427
    expected_c_2 = math.log2(math.e) - math.log2(2.0)
    print(f"  α = 2: c(2) = log_2 e - log_2(2/1) = log_2(e/2) = {expected_c_2:.6f}")
    if abs(predicted_c(2.0) - expected_c_2) > 1e-9:
        print(f"  FAIL: formula gives c(2) = {predicted_c(2.0)}")
        failures += 1

    # α → 1+: c → log_2 e (tightest constraint at K = N+1)
    # (α-1) log(α/(α-1)) → 0 as α → 1+ (L'Hopital)
    print(f"  α → 1+ (tightest): c(α) → log_2 e ≈ {math.log2(math.e):.6f}")
    print(f"    verified at α = 1.01: c = {predicted_c(1.01):.6f}")
    print(f"    verified at α = 1.001: c = {predicted_c(1.001):.6f}")

    # α → ∞: c → 0 (use L'Hopital or Taylor)
    # (α-1) log(α/(α-1)) = (α-1) log(1 + 1/(α-1)) ≈ (α-1) · 1/((α-1) ln 2) = 1/ln 2 = log_2 e
    # So c(α) → log_2 e - log_2 e = 0
    print(f"  α → ∞ (vacuous): c(α) → 0 (verified at α = 1000: c = {predicted_c(1000.0):.6f})")

    # Monotonicity check
    print(f"\nMonotonicity of c(α) on (1, ∞):")
    prev_c = None
    for a in [1.01, 1.5, 2.0, 3.0, 10.0, 100.0]:
        c = predicted_c(a)
        if prev_c is not None and c >= prev_c:
            print(f"  FAIL: c({a}) = {c:.6f} not less than c(prev) = {prev_c:.6f}")
            failures += 1
        else:
            mark = "" if prev_c is None else " (decreasing ✓)"
            print(f"  c({a}) = {c:.6f}{mark}")
        prev_c = c

    # α = 2 at large N: residual should converge to a small constant (-1/2 - 0.5)
    print()
    print("Asymptotic residual at α = 2:")
    for N in [100, 1000, 10000]:
        K = 2 * N
        tc = TC(N, K)
        predicted = N * predicted_c(2.0)
        residual = tc - predicted
        print(f"  N = {N:>6}: TC = {tc:.4f}, predicted N·c(2) = {predicted:.4f},"
              f" residual = {residual:.6f}")

    print()
    if failures == 0:
        print("PASS: Lemma 5.6a verified — TC(D_N) for the uniform-over-distinct-PK")
        print(f"      distribution scales linearly in N with c(α) = log_2 e -")
        print(f"      (α-1) log_2(α/(α-1)). Specifically:")
        for alpha in test_alphas:
            print(f"        α = {alpha}: c({alpha}) ≈ {predicted_c(float(alpha)):.6f}")
        return 0
    else:
        print(f"FAIL: {failures} mismatch(es).")
        return 1


if __name__ == "__main__":
    sys.exit(main())

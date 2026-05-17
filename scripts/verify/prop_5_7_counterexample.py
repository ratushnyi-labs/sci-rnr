#!/usr/bin/env python3
"""
Counterexample to the original Proposition 5.7 (Optimal support is a Q-prefix).

The original proposition (and its weakening to "monotonicity-of-A condition")
claimed: among supports of size s, the top-s positions sorted by Q_v(i)
minimize the expected Theorem 5.3 residual A_s(alpha) = log2 C(s, k-alpha)
+ log2 C(N-s, alpha), under the swap-argument that "stochastically decreases
alpha cannot increase expected cost".

Codex's counterexample (verified here): with N=7, k=2, s=2, A non-decreasing
on alpha in {0,1,2}, a specific atom distribution makes top-Q support {4,2}
strictly worse than support {4,0}. The marginal-only swap argument fails
because joint structure of M_v matters.

This file demonstrates the failure constructively. PASS means counterexample
reproduced; FAIL means counterexample no longer holds (which would mean codex's
attack was wrong).
"""

import math
import sys
import itertools


def log_binom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.log2(math.comb(n, k))


def main() -> int:
    N, k, s = 7, 2, 2

    # A_s(alpha) for alpha in feasible range
    feasible_alpha = [a for a in range(k + 1) if a <= N - s and k - a <= s]
    A = {a: log_binom(s, k - a) + log_binom(N - s, a) for a in feasible_alpha}

    print(f"N={N}, k={k}, s={s}")
    print(f"feasible alpha: {feasible_alpha}")
    print(f"A_s(alpha): {[(a, round(A[a], 4)) for a in feasible_alpha]}")

    # Monotonicity of A: check non-decreasing
    monotone = all(A[feasible_alpha[i + 1]] >= A[feasible_alpha[i]] - 1e-9
                   for i in range(len(feasible_alpha) - 1))
    print(f"A non-decreasing on feasible alpha: {monotone}")
    assert monotone, "Counterexample requires A non-decreasing, broken"

    # Codex's atom distribution: atoms are 2-element subsets of [0,N) with weights
    atoms = {
        frozenset({4, 5}): 0.12,
        frozenset({3, 4}): 0.04,
        frozenset({0, 4}): 0.36,
        frozenset({2, 4}): 0.20,
        frozenset({1, 2}): 0.28,
    }
    assert abs(sum(atoms.values()) - 1.0) < 1e-9

    # Marginal Q_v(i) = Pr(i in M_v) = sum of weights of atoms containing i
    q = [0.0] * N
    for atom, w in atoms.items():
        for i in atom:
            q[i] += w
    print(f"Marginals q[i]: {[(i, round(q[i], 4)) for i in range(N)]}")

    # Sort positions by decreasing Q
    order = sorted(range(N), key=lambda i: (-q[i], i))
    print(f"Order by Q desc: {order}")

    # Top-s prefix (codex's prediction: top-2 by Q is {4, 2})
    top_support = tuple(sorted(order[:s]))
    print(f"Top-s support (claim: optimal): {top_support}")

    # Expected A under top-s support
    def expected_cost(support: tuple) -> float:
        S = set(support)
        total = 0.0
        for atom, w in atoms.items():
            alpha = k - len(S.intersection(atom))
            total += w * A[alpha]
        return total

    top_cost = expected_cost(top_support)
    print(f"Expected cost of top-s support: {round(top_cost, 6)}")

    # Brute-force best over all C(N, s) supports
    all_supports = list(itertools.combinations(range(N), s))
    best_cost, best_support = min((expected_cost(S), S) for S in all_supports)
    print(f"Optimal support (brute-force): {best_support}, cost {round(best_cost, 6)}")

    # Verify counterexample: top is NOT optimal
    if top_cost > best_cost + 1e-6:
        print()
        print(f"PASS: counterexample reproduced.")
        print(f"  Top-s {top_support} costs {round(top_cost, 6)},")
        print(f"  but support {best_support} costs only {round(best_cost, 6)}.")
        print(f"  Proposition 5.7 (top-Q prefix optimal under monotone A) is FALSE.")
        return 0
    else:
        print()
        print(f"FAIL: top-s is optimal here, counterexample no longer holds.")
        print(f"  Either the atoms changed or the marginal-swap argument is recoverable.")
        return 1


if __name__ == "__main__":
    sys.exit(main())

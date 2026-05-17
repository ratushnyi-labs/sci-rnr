#!/usr/bin/env python3
"""
Verification of Theorem 4.3 (Varying-parameter NP-hardness of learned-code
labeling) via the Wagner-Corneil 1990 reduction.

Construction (from the proof in rnr_coding.tex):
  Given a tree T on n vertices with max degree Δ, build:
    F_{uv} = 1/Δ   if u != v and {u,v} ∈ E(T)
    F_{uu} = 1 - deg(u)/Δ
    F_{uv} = 0     otherwise
  Set threshold B = 2(n-1)/Δ.

  Then T embeds as a subgraph of Q_m iff there exists an injective
  labeling E : [n] → {0,1}^m with C_F(E) := sum_{i,j} F_{ij} d_H(E(i), E(j)) <= B.

This script verifies the reduction on small trees:
- F is symmetric and row-stochastic (theorem precondition)
- Hand-built embeddings of trees into hypercubes have C_F(E) = B (saturation)
- Non-embeddable instances (T not a subgraph of Q_m) have all
  C_F(E) > B (the YES iff embedding direction of Lemma 4.3c)

PASS = all sanity checks pass.

Reference: Wagner, A. and Corneil, D.G. (1990), "Embedding Trees in a
Hypercube is NP-Complete", SIAM J. Comput. 19(3), 570-590.
"""

import itertools
import math
import sys


def hamming(a: tuple, b: tuple) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def build_F(tree_edges: list, n: int) -> tuple:
    """Build symmetric row-stochastic F per the Thm 4.3 construction."""
    deg = [0] * n
    for u, v in tree_edges:
        deg[u] += 1
        deg[v] += 1
    Delta = max(deg)
    F = [[0.0] * n for _ in range(n)]
    for u, v in tree_edges:
        F[u][v] = 1.0 / Delta
        F[v][u] = 1.0 / Delta
    for u in range(n):
        F[u][u] = 1.0 - deg[u] / Delta
    return F, Delta


def C_F(F: list, E: tuple) -> float:
    n = len(F)
    total = 0.0
    for i in range(n):
        for j in range(n):
            if F[i][j] != 0:
                total += F[i][j] * hamming(E[i], E[j])
    return total


def best_C_F(F: list, n: int, m: int) -> tuple:
    """Brute-force minimum C_F over all injective labelings [n] -> {0,1}^m."""
    cube = list(itertools.product((0, 1), repeat=m))
    if n > len(cube):
        return float("inf"), None
    best = float("inf")
    best_E = None
    for E in itertools.permutations(cube, n):
        c = C_F(F, E)
        if c < best:
            best = c
            best_E = E
    return best, best_E


def main() -> int:
    failures = 0

    # Test 1: A path graph T = 0-1-2 (3 vertices, 2 edges) embeds in Q_2.
    # Embedding: 0 -> (0,0), 1 -> (0,1), 2 -> (1,1).
    print("Test 1: path P_3 embeds in Q_2")
    tree = [(0, 1), (1, 2)]
    n, m = 3, 2
    F, Delta = build_F(tree, n)
    print(f"  Delta={Delta}, F=\n   {F[0]}\n   {F[1]}\n   {F[2]}")
    assert all(abs(sum(row) - 1.0) < 1e-9 for row in F), "F not row-stochastic"
    assert all(F[i][j] == F[j][i] for i in range(n) for j in range(n)), "F not symmetric"
    B = 2 * (n - 1) / Delta
    print(f"  B = 2(n-1)/Delta = {B}")
    best, best_E = best_C_F(F, n, m)
    print(f"  best C_F over injective labelings = {best}, witness = {best_E}")
    if abs(best - B) > 1e-9:
        print(f"  FAIL: expected best == B = {B}, got {best}")
        failures += 1
    else:
        print(f"  PASS: best == B, embedding achievable")

    # Test 2: A 4-cycle is NOT a tree, skip. But K_{1,3} (star with 3 leaves)
    # embeds in Q_3 (center at 000, leaves at 001, 010, 100).
    print()
    print("Test 2: star K_{1,3} embeds in Q_3")
    tree = [(0, 1), (0, 2), (0, 3)]
    n, m = 4, 3
    F, Delta = build_F(tree, n)
    print(f"  Delta={Delta}")
    B = 2 * (n - 1) / Delta
    best, best_E = best_C_F(F, n, m)
    print(f"  B = {B}, best C_F = {best}, witness = {best_E}")
    if abs(best - B) > 1e-9:
        print(f"  FAIL: expected best == B")
        failures += 1
    else:
        print(f"  PASS: star embeds in Q_3 (Delta=3 leaves fit in 3-cube)")

    # Test 3: K_{1,4} (star with 4 leaves) does NOT embed in Q_3.
    # Q_3 has degree 3; the center would need 4 cube-neighbors, impossible.
    print()
    print("Test 3: star K_{1,4} does NOT embed in Q_3")
    tree = [(0, 1), (0, 2), (0, 3), (0, 4)]
    n, m = 5, 3
    F, Delta = build_F(tree, n)
    B = 2 * (n - 1) / Delta
    best, best_E = best_C_F(F, n, m)
    print(f"  Delta={Delta}, B={B}, best C_F = {best}, witness = {best_E}")
    if best <= B + 1e-9:
        print(f"  FAIL: expected best > B (no embedding), got best={best}")
        failures += 1
    else:
        print(f"  PASS: best > B, no embedding exists (Q_3 max degree 3 < {Delta})")

    # Test 4: K_{1,4} DOES embed in Q_4 (4-cube has degree 4).
    print()
    print("Test 4: star K_{1,4} embeds in Q_4")
    n, m = 5, 4
    F, Delta = build_F([(0, 1), (0, 2), (0, 3), (0, 4)], n)
    B = 2 * (n - 1) / Delta
    best, best_E = best_C_F(F, n, m)
    print(f"  B={B}, best C_F = {best}, witness = {best_E}")
    if abs(best - B) > 1e-9:
        print(f"  FAIL: expected best == B")
        failures += 1
    else:
        print(f"  PASS: K_{{1,4}} fits in Q_4")

    # Test 5: 6-vertex path P_6 embeds in Q_3 via a Gray-code-like walk.
    # Q_3 has 8 vertices, P_6 needs 6.
    print()
    print("Test 5: path P_6 embeds in Q_3")
    tree = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]
    n, m = 6, 3
    F, Delta = build_F(tree, n)
    B = 2 * (n - 1) / Delta
    best, best_E = best_C_F(F, n, m)
    print(f"  Delta={Delta}, B={B}, best C_F = {best}")
    if abs(best - B) > 1e-9:
        print(f"  FAIL: expected best == B")
        failures += 1
    else:
        print(f"  PASS: P_6 embeds in Q_3 via Gray-code path")

    # Test 6: Y-shaped tree with leaves at distance 3 - check non-trivial structure
    print()
    print("Test 6: tree with branching point - degree 3 internal node")
    tree = [(0, 1), (0, 2), (0, 3), (1, 4), (2, 5)]
    n, m = 6, 3
    F, Delta = build_F(tree, n)
    B = 2 * (n - 1) / Delta
    best, best_E = best_C_F(F, n, m)
    print(f"  Delta={Delta}, B={B}, best C_F = {best}")
    if abs(best - B) > 1e-9:
        print(f"  FAIL: expected best == B")
        failures += 1
    else:
        print(f"  PASS: tree embeds")

    print()
    if failures == 0:
        print(f"PASS: all reduction sanity checks passed (Wagner-Corneil precondition verified for our F construction).")
        return 0
    else:
        print(f"FAIL: {failures} sanity checks failed; reduction construction is wrong.")
        return 1


if __name__ == "__main__":
    sys.exit(main())

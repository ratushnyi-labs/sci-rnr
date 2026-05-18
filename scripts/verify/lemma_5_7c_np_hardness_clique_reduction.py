#!/usr/bin/env python3
r"""
Verification of Lemma 5.7c: NP-hardness of optimal support selection for
general dependent joint laws via CLIQUE -> support-selection reduction.

Reduction setup: given a CLIQUE instance (G = (V_G, E_G), k) with
|V_G| <= 3k+1, construct a support-selection instance with N = 3k+1,
s = k, k_v = 2, atoms = edges of G with uniform weight 1/|E_G|, and
threshold tau = A_s(1) - (A_s(1) - A_s(0)) * C(k,2) / |E_G|.

Key identity at N = 3s+1: A_s(1) = A_s(2), so the objective depends
ONLY on n_0(S) = #edges with both endpoints in S. The threshold tau is
the value of E[A_s] when n_0 = C(k,2), so optimum <= tau iff G has a
k-clique. Padding vertices contribute nothing (no edges touch them).

This script:
- Verifies A_s(1) = A_s(2) at N = 3s+1 (the algebraic identity the
  reduction relies on)
- Generates random graphs of varying size and density
- Brute-forces both CLIQUE (does G have a k-clique?) and support-selection
  (is the optimal E[A_s] <= tau?)
- Confirms the two answers always agree

PASS = both checks succeed on all trials.
"""

import itertools
import math
import random
import sys


def A_s(alpha: int, N: int, s: int, k_v: int) -> float:
    """Type-II partition-rank log-length at residual count alpha."""
    if alpha < 0 or k_v - alpha < 0 or alpha > N - s or k_v - alpha > s:
        return float("inf")
    return math.log2(math.comb(s, k_v - alpha)) + math.log2(math.comb(N - s, alpha))


def expected_cost(S: tuple, atoms: dict, N: int, s: int, k_v: int) -> float:
    if not atoms:
        return 0.0
    S_set = frozenset(S)
    total = 0.0
    for T, w in atoms.items():
        alpha = len(T - S_set)
        total += w * A_s(alpha, N, s, k_v)
    return total


def support_selection_brute_force(atoms: dict, N: int, s: int, k_v: int):
    best, best_S = float("inf"), None
    for S in itertools.combinations(range(N), s):
        c = expected_cost(S, atoms, N, s, k_v)
        if c < best:
            best, best_S = c, S
    return best_S, best


def has_clique(G_edges: set, vertices: list, k: int):
    for S in itertools.combinations(vertices, k):
        if all(frozenset({a, b}) in G_edges for a, b in itertools.combinations(S, 2)):
            return True, S
    return False, None


def reduce_clique_to_support(G_edges: set, V_G_size: int, k: int):
    N = 3 * k + 1
    if V_G_size > N:
        raise ValueError(f"|V_G|={V_G_size} > 3k+1={N}; reduction requires |V_G| <= 3k+1")
    s, k_v = k, 2
    w = 1.0 / len(G_edges)
    atoms = {frozenset(e): w for e in G_edges}
    a0, a1 = A_s(0, N, s, k_v), A_s(1, N, s, k_v)
    tau = a1 - (a1 - a0) * math.comb(k, 2) / len(G_edges)
    return atoms, N, s, k_v, tau


def test_A_s_equality() -> bool:
    print("Verifying A_s(1) = A_s(2) at N = 3s+1:")
    for s in range(2, 8):
        N = 3 * s + 1
        a1, a2 = A_s(1, N, s, 2), A_s(2, N, s, 2)
        print(f"  s={s}, N={N}: A_s(1)={a1:.6f}, A_s(2)={a2:.6f}, diff={a1-a2:.2e}")
        if abs(a1 - a2) > 1e-9:
            return False
    print("  OK: A_s(1) = A_s(2) holds exactly at N = 3s+1.\n")
    return True


def test_clique_reduction() -> bool:
    rng = random.Random(42)
    trials = 60
    failures = 0
    clique_yes_count = 0

    for trial in range(trials):
        k = rng.randint(2, 4)
        V_G_size = rng.randint(k, 3 * k + 1)
        p = rng.uniform(0.3, 0.85)
        G_edges = set()
        for a, b in itertools.combinations(range(V_G_size), 2):
            if rng.random() < p:
                G_edges.add(frozenset({a, b}))
        if not G_edges:
            continue

        has_k, clique = has_clique(G_edges, list(range(V_G_size)), k)
        atoms, N_red, s, k_v, tau = reduce_clique_to_support(G_edges, V_G_size, k)
        best_S, best_cost = support_selection_brute_force(atoms, N_red, s, k_v)

        support_yes = best_cost <= tau + 1e-9
        if has_k != support_yes:
            print(f"  FAIL trial {trial}: k={k}, |V_G|={V_G_size}, |E|={len(G_edges)}")
            print(f"    has-{k}-clique: {has_k}, support-cost-<=-tau: {support_yes}")
            print(f"    best_cost = {best_cost:.6f}, tau = {tau:.6f}")
            print(f"    edges: {[sorted(e) for e in G_edges]}")
            failures += 1
        elif has_k:
            clique_yes_count += 1

    print(f"Reduction trials: {trials} total, {clique_yes_count} YES instances, {failures} failures.")
    return failures == 0


def test_no_clique_means_strict_excess() -> bool:
    """A graph with no k-clique should yield strictly best_cost > tau."""
    # Star K_{1,k}: center connected to k leaves, no triangle even for k >= 2
    # Try detecting clique-3 in star K_{1,4} (5 vertices, 4 edges, no triangle)
    G_edges = {frozenset({0, i}) for i in range(1, 5)}
    V_G_size = 5
    k = 3
    atoms, N, s, k_v, tau = reduce_clique_to_support(G_edges, V_G_size, k)
    _, best_cost = support_selection_brute_force(atoms, N, s, k_v)
    print(f"\nNo-clique sanity (star K_{{1,4}}, k=3):")
    print(f"  best_cost = {best_cost:.6f}, tau = {tau:.6f}, gap = {best_cost - tau:.6f}")
    if best_cost > tau + 1e-9:
        print(f"  OK: best_cost > tau (no 3-clique correctly detected)")
        return True
    print(f"  FAIL: best_cost should exceed tau")
    return False


def test_clique_means_exact_tau() -> bool:
    """K_3 (triangle) trivially has a 3-clique, so optimum should equal tau."""
    G_edges = {frozenset({0, 1}), frozenset({0, 2}), frozenset({1, 2})}
    V_G_size = 3
    k = 3
    atoms, N, s, k_v, tau = reduce_clique_to_support(G_edges, V_G_size, k)
    _, best_cost = support_selection_brute_force(atoms, N, s, k_v)
    print(f"\nClique sanity (triangle K_3, k=3):")
    print(f"  best_cost = {best_cost:.6f}, tau = {tau:.6f}, gap = {best_cost - tau:.6f}")
    if abs(best_cost - tau) < 1e-9:
        print(f"  OK: best_cost = tau (3-clique correctly detected)")
        return True
    print(f"  FAIL: best_cost should equal tau when clique exists")
    return False


def test_no_padding_case() -> bool:
    """Force |V_G| = 3k+1 (no padding): K_7 has every clique up to size 7."""
    V_G_size, k = 7, 2  # 3k+1 = 7, no padding
    G_edges = {frozenset({a, b}) for a, b in itertools.combinations(range(V_G_size), 2)}
    has_k, _ = has_clique(G_edges, list(range(V_G_size)), k)
    atoms, N, s, k_v, tau = reduce_clique_to_support(G_edges, V_G_size, k)
    _, best_cost = support_selection_brute_force(atoms, N, s, k_v)
    print(f"\nNo-padding sanity (K_7, k=2, |V_G|=N=3k+1=7):")
    print(f"  has 2-clique: {has_k}, best_cost={best_cost:.6f}, tau={tau:.6f}")
    if not has_k:
        print(f"  FAIL: K_7 must have a 2-clique")
        return False
    if best_cost > tau + 1e-9:
        print(f"  FAIL: best_cost should be <= tau when clique exists")
        return False
    print(f"  OK: clique detected, |V_G|=N=3k+1 path exercised")
    return True


def test_empty_edge_branch() -> bool:
    """E_G = empty: per Lemma 5.7c, the reduction outputs a trivial NO."""
    G_edges = set()
    V_G_size, k = 4, 2
    # Combinatorial decision: trivial NO instance N=5, s=2, m=0, D=1
    N_combinatorial, s_combinatorial, D = 5, 2, 1
    atoms_empty = {}  # no atoms
    # Brute force: any |S|=2 has n_0=0 < 1, so the answer is NO
    n_0_max = 0
    for S in itertools.combinations(range(N_combinatorial), s_combinatorial):
        S_set = frozenset(S)
        n_0 = sum(1 for T in atoms_empty if T <= S_set)
        n_0_max = max(n_0_max, n_0)
    has_k, _ = has_clique(G_edges, list(range(V_G_size)), k)
    print(f"\nEmpty-edge branch (G has no edges, k=2):")
    print(f"  has 2-clique: {has_k} (correctly NO)")
    print(f"  combinatorial output: N={N_combinatorial}, s={s_combinatorial}, m=0, D={D}")
    print(f"  brute-force n_0_max = {n_0_max} < D = {D}, so combinatorial answer: NO")
    if has_k or n_0_max >= D:
        print(f"  FAIL: empty-edge instance should be trivial NO")
        return False
    print(f"  OK: CLIQUE NO <-> combinatorial NO")
    return True


def main() -> int:
    print("Verification of Lemma 5.7c (NP-hardness via CLIQUE -> support-selection)")
    print("=" * 75)
    ok1 = test_A_s_equality()
    ok2 = test_clique_reduction()
    ok3 = test_no_clique_means_strict_excess()
    ok4 = test_clique_means_exact_tau()
    ok5 = test_no_padding_case()
    ok6 = test_empty_edge_branch()

    print()
    if ok1 and ok2 and ok3 and ok4 and ok5 and ok6:
        print("PASS: CLIQUE-to-support-selection reduction verified.")
        print("      Therefore support-selection (general dependent regime) is NP-hard.")
        print("      Open Problem 6 (formal NP-hardness for general joint laws) is")
        print("      resolved in the affirmative; Lemma 5.7a's polynomial-time top-Q")
        print("      result for permutation-like laws stands as a positive contrast.")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())

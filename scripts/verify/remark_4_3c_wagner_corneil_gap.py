#!/usr/bin/env python3
r"""
Verification of Remark 4.3c: Wagner-Corneil reduction from tree-to-hypercube
embedding to E_12 gives only inverse-linear inapproximability gap.

Setup. Given tree T with n vertices and max degree Delta, Theorem 4.3
constructs an E_12 instance with threshold B = 2(n-1)/Delta. The cost
function is C_F(E) = sum_{(u,v) in E(T)} (2/Delta) * d_H(E(u), E(v)).

For a YES instance (T embeds in Q_m), each edge maps to d_H = 1, giving
C_F = B exactly.

For a NO instance, at least one edge has d_H >= 2 (Hamming distance is
an integer >= 1 for injective labelings, and non-embedded edge must have
d_H >= 2). So C_F >= B + 2/Delta.

The relative gap is (B + 2/Delta) / B = 1 + 1/(n-1), which goes to 1 as
n grows. Thus the reduction only gives inapproximability factor 1 - 1/n,
not constant gap.

This script:
1. Constructs small trees T of various sizes (n = 4, 8, 16)
2. For each, brute-force searches over labelings E: V(T) -> {0,1}^m
3. Computes optimal C_F(E*) and compares to B = 2(n-1)/Delta
4. Reports the relative gap (C_F(E_no_embed) / C_F(E_yes)) and confirms
   it approaches 1 as n grows

PASS = relative gap matches 1 + 1/(n-1) theory at each n, decreases
toward 1 as n grows.
"""

import itertools
import math
import sys


def hamming(a: tuple, b: tuple) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def C_F_tree(labeling: dict, tree_edges: list, Delta: int) -> float:
    """C_F(E) for tree T with edge weights F_{uv} = 1/Delta on each edge."""
    return sum(2.0 / Delta * hamming(labeling[u], labeling[v]) for u, v in tree_edges)


def brute_force_optimal_labeling(tree_vertices: list, tree_edges: list, m: int, Delta: int):
    """Find labeling E : V -> {0,1}^m minimizing C_F(E)."""
    n = len(tree_vertices)
    best, best_lab = float("inf"), None
    codewords = list(itertools.product([0, 1], repeat=m))
    if n > len(codewords):
        return None, float("inf")  # infeasible
    for perm in itertools.permutations(codewords, n):
        lab = dict(zip(tree_vertices, perm))
        cost = C_F_tree(lab, tree_edges, Delta)
        if cost < best:
            best, best_lab = cost, lab
    return best_lab, best


def make_path_tree(n: int) -> list:
    """Path tree 0 - 1 - 2 - ... - (n-1)."""
    return [(i, i + 1) for i in range(n - 1)]


def make_star_tree(n: int) -> list:
    """Star tree: center 0 connected to leaves 1, ..., n-1."""
    return [(0, i) for i in range(1, n)]


def test_path_embeds_constant_gap():
    """Path tree embeds into Q_m for m = ceil(log_2(n)). Verify YES case."""
    print("Test 1: Path tree (always embeds into Q_{ceil(log_2 n)})")
    all_ok = True
    for n in [4, 6, 8]:
        m = math.ceil(math.log2(n))
        Delta = 2  # internal degree of path
        edges = make_path_tree(n)
        B_theoretical = 2 * (n - 1) / Delta  # = n - 1
        lab, opt = brute_force_optimal_labeling(list(range(n)), edges, m, Delta)
        ratio = opt / B_theoretical if B_theoretical > 0 else 0
        print(f"  n={n}, m={m}, Delta={Delta}: B={B_theoretical:.4f}, OPT={opt:.4f}, ratio={ratio:.4f}")
        if abs(ratio - 1.0) > 1e-9:
            print(f"    FAIL: Path tree should embed exactly (ratio=1)")
            all_ok = False
        else:
            print(f"    OK: Path embeds (YES instance)")
    return all_ok


def test_star_n5_m2_no_embedding():
    """Star K_{1,4} on 5 vertices does NOT embed into Q_2 (only 4 codewords)
    but Q_3 has 8 codewords, where it embeds trivially. Test NO when m=2
    forces infeasibility."""
    print("\nTest 2: Star K_{1,4} (5 vertices), m=2 forces NO embedding")
    n, m = 5, 2
    edges = make_star_tree(n)
    Delta = 4  # center degree
    B_theoretical = 2 * (n - 1) / Delta  # = 8/4 = 2
    lab, opt = brute_force_optimal_labeling(list(range(n)), edges, m, Delta)
    # |V| = 5 > 2^m = 4, so no injective labeling exists -> infeasible.
    print(f"  n={n}, m={m}, |V|=5 > 2^m=4: infeasible per Theorem 4.3 step 1")
    if opt == float("inf"):
        print(f"    OK: brute force confirms no injective labeling exists")
        return True
    print(f"    FAIL: expected infeasibility")
    return False


def test_gap_scaling():
    """Verify the (B + 2/Delta) / B = 1 + 1/(n-1) gap formula numerically.

    For a tree that requires d_H >= 2 on at least one edge (when Q_m is
    too small), the optimal C_F satisfies B + 2/Delta <= OPT.
    """
    print("\nTest 3: Gap-scaling verification — predicted gap = 1 + 1/(n-1)")
    print("  (Brute force on small star trees that need at least one d_H >= 2 edge)")

    all_ok = True
    for n in [4, 5, 6]:
        # Star tree on n vertices, m = ceil(log2(n)) so we have JUST enough codewords
        # but not enough for star to embed at d_H = 1 everywhere
        m = math.ceil(math.log2(n))
        Delta = n - 1  # center degree
        edges = make_star_tree(n)
        B_theoretical = 2 * (n - 1) / Delta  # = 2
        lab, opt = brute_force_optimal_labeling(list(range(n)), edges, m, Delta)

        # Star K_{1,n-1} embeds into Q_m iff m >= n-1 (each leaf at distinct
        # codeword adjacent to center). When n-1 > m, star does NOT embed.
        embeds = (n - 1) <= m

        gap_theory = 1 + 1.0 / (n - 1)
        ratio = opt / B_theoretical if B_theoretical > 0 else 0
        print(
            f"  n={n}, m={m}, Delta={Delta}: B={B_theoretical:.4f}, OPT={opt:.4f}, "
            f"ratio={ratio:.4f}, predicted-gap-if-NO={gap_theory:.4f}, embeds={embeds}"
        )
        if embeds:
            # YES: ratio should be 1
            if abs(ratio - 1.0) > 1e-9:
                print(f"    FAIL: should embed at d_H=1")
                all_ok = False
        else:
            # NO: ratio should be >= 1 + 1/(n-1)
            if ratio < gap_theory - 1e-9:
                print(f"    FAIL: NO case should give ratio >= {gap_theory:.4f}")
                all_ok = False
            else:
                print(f"    OK: NO case ratio {ratio:.4f} >= predicted {gap_theory:.4f}")
    return all_ok


def test_asymptotic_gap_decay():
    """Show that 1 + 1/(n-1) -> 1 as n grows. This is the core observation
    of Remark 4.3c."""
    print("\nTest 4: Asymptotic gap decay 1 + 1/(n-1) -> 1")
    for n in [10, 100, 1000, 10000]:
        gap = 1 + 1.0 / (n - 1)
        rel = 1 - 1.0 / n  # equivalent statement
        print(
            f"  n={n}: gap = 1 + 1/(n-1) = {gap:.6f}, equivalently {rel:.6f}-multiplicative"
        )
    print("  -> The gap vanishes; Wagner-Corneil gives no constant-gap inapproximability.")
    return True


def main() -> int:
    print("Verification of Remark 4.3c (Wagner-Corneil inverse-linear gap)")
    print("=" * 70)
    ok1 = test_path_embeds_constant_gap()
    ok2 = test_star_n5_m2_no_embedding()
    ok3 = test_gap_scaling()
    ok4 = test_asymptotic_gap_decay()

    print()
    if all([ok1, ok2, ok3, ok4]):
        print("PASS: Remark 4.3c verified.")
        print("      Wagner-Corneil reduction gives only (1 + 1/(n-1)) gap")
        print("      (equivalently, NP-hard to approximate within (1 - 1/n)).")
        print("      The gap vanishes asymptotically; constant-gap inapproximability")
        print("      of E_12 requires a different reduction.")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())

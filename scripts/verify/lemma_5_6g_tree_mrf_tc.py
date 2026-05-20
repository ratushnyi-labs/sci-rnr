#!/usr/bin/env python3
r"""
Verification of Lemma 5.6g: tree-structured MRF total correlation
equals exactly the sum of edge pairwise mutual informations.

Claim (Pearl decomposition):
    TC(D_N) = sum_{(u,v) in E} I(X_u; X_v)

for any tree-structured Markov random field on N vertices with strictly
positive pairwise factors over a finite alphabet.

This script:
1. Generates several random tree-structured MRF instances of varying
   sizes and topologies (path / star / balanced binary / random tree).
2. Computes joint distribution by belief propagation / chain rule.
3. Computes TC(D_N) = sum H(X_v) - H(D_N).
4. Computes sum of edge MIs.
5. Verifies the identity TC = sum edge MIs at machine precision.

PASS = identity holds at machine precision (< 1e-9) for every test.
"""

import itertools
import math
import random
import sys
from collections import defaultdict


def entropy_of_dist(probs):
    return -sum(p * math.log2(p) for p in probs if p > 1e-15)


def random_factor(alphabet_size, eps=0.5):
    """Random strictly-positive pairwise factor psi(x, y) on alphabet."""
    sigma = alphabet_size
    f = {}
    for x in range(sigma):
        for y in range(sigma):
            # log-linear with random correlation strength
            f[(x, y)] = math.exp(random.gauss(0, eps))
    return f


def joint_distribution_tree(N, edges, factors, alphabet_size, root=0):
    """
    Compute exact joint distribution of tree MRF by enumeration.
    Only feasible for small N · alphabet_size combinations.
    """
    sigma = alphabet_size
    joint = {}
    total = 0.0
    for assignment in itertools.product(range(sigma), repeat=N):
        weight = 1.0
        for (u, v) in edges:
            weight *= factors[(u, v)][(assignment[u], assignment[v])]
        joint[assignment] = weight
        total += weight
    # Normalize
    for key in joint:
        joint[key] /= total
    return joint


def marginal_single(joint, v, N):
    """Marginal distribution of X_v."""
    marg = defaultdict(float)
    for assignment, p in joint.items():
        marg[assignment[v]] += p
    return dict(marg)


def marginal_pair(joint, u, v, N):
    """Joint marginal of (X_u, X_v)."""
    marg = defaultdict(float)
    for assignment, p in joint.items():
        marg[(assignment[u], assignment[v])] += p
    return dict(marg)


def mutual_information(joint, u, v, N):
    """I(X_u; X_v) from joint distribution."""
    pair_marg = marginal_pair(joint, u, v, N)
    marg_u = marginal_single(joint, u, N)
    marg_v = marginal_single(joint, v, N)
    mi = 0.0
    for (xu, xv), p_uv in pair_marg.items():
        if p_uv > 1e-15:
            mi += p_uv * math.log2(p_uv / (marg_u[xu] * marg_v[xv]))
    return mi


def total_correlation(joint, N):
    """TC(D_N) = sum_v H(X_v) - H(D_N)."""
    H_sum = 0.0
    for v in range(N):
        marg = marginal_single(joint, v, N)
        H_sum += entropy_of_dist(marg.values())
    H_joint = entropy_of_dist(joint.values())
    return H_sum - H_joint


def test_path_tree(N, sigma):
    """Path graph: 0-1-2-...-(N-1)."""
    edges = [(i, i + 1) for i in range(N - 1)]
    factors = {e: random_factor(sigma) for e in edges}
    joint = joint_distribution_tree(N, edges, factors, sigma)
    tc = total_correlation(joint, N)
    sum_mi = sum(mutual_information(joint, u, v, N) for (u, v) in edges)
    diff = abs(tc - sum_mi)
    return tc, sum_mi, diff


def test_star_tree(N, sigma):
    """Star graph: 0 connected to 1, 2, ..., N-1."""
    edges = [(0, i) for i in range(1, N)]
    factors = {e: random_factor(sigma) for e in edges}
    joint = joint_distribution_tree(N, edges, factors, sigma)
    tc = total_correlation(joint, N)
    sum_mi = sum(mutual_information(joint, u, v, N) for (u, v) in edges)
    diff = abs(tc - sum_mi)
    return tc, sum_mi, diff


def test_balanced_binary_tree(N, sigma):
    """Balanced binary tree with N = 2^h - 1 vertices."""
    edges = []
    for i in range((N - 1) // 2 + 1):
        left = 2 * i + 1
        right = 2 * i + 2
        if left < N:
            edges.append((i, left))
        if right < N:
            edges.append((i, right))
    factors = {e: random_factor(sigma) for e in edges}
    joint = joint_distribution_tree(N, edges, factors, sigma)
    tc = total_correlation(joint, N)
    sum_mi = sum(mutual_information(joint, u, v, N) for (u, v) in edges)
    diff = abs(tc - sum_mi)
    return tc, sum_mi, diff


def test_random_tree(N, sigma, seed=42):
    """Random tree via uniform random spanning tree (Prufer-like)."""
    random.seed(seed)
    parents = [random.randint(0, i - 1) for i in range(1, N)]
    edges = [(parents[i - 1], i) for i in range(1, N)]
    factors = {e: random_factor(sigma) for e in edges}
    joint = joint_distribution_tree(N, edges, factors, sigma)
    tc = total_correlation(joint, N)
    sum_mi = sum(mutual_information(joint, u, v, N) for (u, v) in edges)
    diff = abs(tc - sum_mi)
    return tc, sum_mi, diff


def main() -> int:
    print("Verification of Lemma 5.6g (tree-structured MRF: TC = sum edge MIs)")
    print("=" * 70)

    random.seed(2026)
    all_ok = True

    print("\nTest 1: Path graph, N=4, sigma=3")
    tc, sum_mi, diff = test_path_tree(4, 3)
    print(f"  TC = {tc:.10f}, sum_MI = {sum_mi:.10f}, diff = {diff:.2e}")
    if diff > 1e-9:
        all_ok = False
        print("  FAIL")

    print("\nTest 2: Path graph, N=6, sigma=2")
    tc, sum_mi, diff = test_path_tree(6, 2)
    print(f"  TC = {tc:.10f}, sum_MI = {sum_mi:.10f}, diff = {diff:.2e}")
    if diff > 1e-9:
        all_ok = False
        print("  FAIL")

    print("\nTest 3: Star graph, N=5, sigma=3")
    tc, sum_mi, diff = test_star_tree(5, 3)
    print(f"  TC = {tc:.10f}, sum_MI = {sum_mi:.10f}, diff = {diff:.2e}")
    if diff > 1e-9:
        all_ok = False
        print("  FAIL")

    print("\nTest 4: Balanced binary tree, N=7, sigma=2")
    tc, sum_mi, diff = test_balanced_binary_tree(7, 2)
    print(f"  TC = {tc:.10f}, sum_MI = {sum_mi:.10f}, diff = {diff:.2e}")
    if diff > 1e-9:
        all_ok = False
        print("  FAIL")

    print("\nTest 5: Random tree, N=5, sigma=3")
    tc, sum_mi, diff = test_random_tree(5, 3, seed=42)
    print(f"  TC = {tc:.10f}, sum_MI = {sum_mi:.10f}, diff = {diff:.2e}")
    if diff > 1e-9:
        all_ok = False
        print("  FAIL")

    print("\nTest 6: Random tree, N=6, sigma=2, multiple seeds")
    for s in [1, 7, 13, 100, 999]:
        tc, sum_mi, diff = test_random_tree(6, 2, seed=s)
        print(f"  seed={s}: diff={diff:.2e}", end=" ")
        if diff > 1e-9:
            all_ok = False
            print("FAIL")
        else:
            print("OK")

    print()
    if all_ok:
        print("PASS: Lemma 5.6g verified across path/star/binary/random topologies.")
        print("      TC(tree-MRF) = sum of edge pairwise mutual informations")
        print("      holds at machine precision for all tested instances.")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())

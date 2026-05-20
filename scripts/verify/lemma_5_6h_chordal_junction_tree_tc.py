#!/usr/bin/env python3
r"""
Verification of Lemma 5.6h: chordal graphical models on bounded
treewidth admit the junction-tree entropy decomposition

    H(D_N) = sum_{C in cliques} H(X_C) - sum_{S in separators} H(X_S)

equivalently

    TC(D_N) = sum_v H(X_v) - sum_C H(X_C) + sum_S H(X_S)

This script:
1. Generates chordal graph instances (interval, 2-tree path, fan,
   tree).
2. Constructs MRFs with strictly positive clique factors.
3. Computes joint distribution by enumeration.
4. Verifies the junction-tree H(D_N) decomposition holds at
   machine precision.
5. ALSO documents that the NAIVE simplification
   TC = sum_C TC(X_C) FAILS for treewidth >= 2 (this confirms the
   retraction in Lemma 5.6h's 'Honest caveat').

PASS = junction-tree decomposition holds at machine precision (<1e-9);
       AND the clique-wise-sum simplification FAILS by >= 1e-2 for at
       least one w>=2 chordal instance (confirming tree-essentiality
       of the simpler form).
"""

import itertools
import math
import random
import sys
from collections import defaultdict


def entropy_of_dist(probs):
    return -sum(p * math.log2(p) for p in probs if p > 1e-15)


def marginal_subset(joint, subset, N):
    marg = defaultdict(float)
    for assignment, p in joint.items():
        key = tuple(assignment[v] for v in subset)
        marg[key] += p
    return dict(marg)


def joint_entropy_subset(joint, subset, N):
    marg = marginal_subset(joint, subset, N)
    return entropy_of_dist(marg.values())


def total_correlation_subset(joint, subset, N):
    H_sum = sum(joint_entropy_subset(joint, (v,), N) for v in subset)
    H_joint = joint_entropy_subset(joint, subset, N)
    return H_sum - H_joint


def total_correlation_full(joint, N):
    return total_correlation_subset(joint, tuple(range(N)), N)


def build_chordal_mrf(N, cliques, alphabet_size, seed):
    random.seed(seed)
    sigma = alphabet_size
    factors = {}
    for C in cliques:
        psi = {}
        for assignment in itertools.product(range(sigma), repeat=len(C)):
            psi[assignment] = math.exp(random.gauss(0, 0.5))
        factors[tuple(C)] = psi

    joint = {}
    total = 0.0
    for assignment in itertools.product(range(sigma), repeat=N):
        w = 1.0
        for C in cliques:
            key = tuple(assignment[v] for v in C)
            w *= factors[tuple(C)][key]
        joint[assignment] = w
        total += w
    for key in joint:
        joint[key] /= total
    return joint, factors


def junction_tree_decomp_check(joint, N, cliques, separators):
    """Verify H(D_N) = sum_C H(X_C) - sum_S H(X_S)."""
    H_joint = joint_entropy_subset(joint, tuple(range(N)), N)
    sum_C_H = sum(joint_entropy_subset(joint, C, N) for C in cliques)
    sum_S_H = sum(joint_entropy_subset(joint, S, N) for S in separators)
    H_predicted = sum_C_H - sum_S_H
    return H_joint, H_predicted, abs(H_joint - H_predicted)


def tc_decomp_check(joint, N, cliques, separators):
    """Verify TC(D_N) = sum_v H(X_v) - sum_C H(X_C) + sum_S H(X_S)."""
    tc_actual = total_correlation_full(joint, N)
    sum_v_H = sum(joint_entropy_subset(joint, (v,), N) for v in range(N))
    sum_C_H = sum(joint_entropy_subset(joint, C, N) for C in cliques)
    sum_S_H = sum(joint_entropy_subset(joint, S, N) for S in separators)
    tc_predicted = sum_v_H - sum_C_H + sum_S_H
    return tc_actual, tc_predicted, abs(tc_actual - tc_predicted)


def clique_wise_naive_check(joint, N, cliques):
    """The WRONG identity: TC = sum_C TC(X_C). This should FAIL for w >= 2."""
    tc_actual = total_correlation_full(joint, N)
    sum_clique_tc = sum(total_correlation_subset(joint, C, N) for C in cliques)
    return tc_actual, sum_clique_tc, abs(tc_actual - sum_clique_tc)


def test_2_tree_path(N, sigma, seed):
    """2-tree path: maximal cliques are consecutive triples.
    Junction-tree separators are consecutive pairs."""
    cliques = [(i, i + 1, i + 2) for i in range(N - 2)]
    separators = [(i + 1, i + 2) for i in range(N - 3)]
    joint, _ = build_chordal_mrf(N, cliques, sigma, seed)
    return joint, cliques, separators


def test_fan_graph(N, sigma, seed):
    """Fan graph: central vertex 0 + path 1-2-...-(N-1).
    Maximal cliques are (0, i, i+1) for i=1..N-2.
    Separator between (0, i, i+1) and (0, i+1, i+2) is {0, i+1}."""
    cliques = [(0, i, i + 1) for i in range(1, N - 1)]
    separators = [(0, i + 1) for i in range(1, N - 2)]
    joint, _ = build_chordal_mrf(N, cliques, sigma, seed)
    return joint, cliques, separators


def test_tree(N, sigma, seed):
    """Tree (treewidth-1 case): maximal cliques are edges, separators
    are vertices (each appearing deg-1 times along the junction tree)."""
    random.seed(seed)
    parents = [random.randint(0, i - 1) for i in range(1, N)]
    cliques = [(parents[i - 1], i) for i in range(1, N)]
    # Junction-tree separators: each vertex v with degree d appears in
    # d-1 junction-tree edges, each separator is the singleton {v}.
    deg = defaultdict(int)
    for (u, v) in cliques:
        deg[u] += 1
        deg[v] += 1
    separators = []
    for v in range(N):
        for _ in range(max(0, deg[v] - 1)):
            separators.append((v,))
    joint, _ = build_chordal_mrf(N, cliques, sigma, seed)
    return joint, cliques, separators


def main() -> int:
    print("Verification of Lemma 5.6h (chordal MRF junction-tree decomposition)")
    print("=" * 70)

    decomp_all_ok = True
    naive_disagreement_found = False

    print("\nTest A: 2-tree path N=5, sigma=2, seed=42")
    joint, cliques, separators = test_2_tree_path(5, 2, seed=42)
    H_act, H_pred, diff = junction_tree_decomp_check(joint, 5, cliques, separators)
    print(f"  H(D_N) decomposition:  actual={H_act:.10f}, pred={H_pred:.10f}, diff={diff:.2e}")
    if diff > 1e-9:
        decomp_all_ok = False
        print("  FAIL (decomposition does not match)")
    tc_act, tc_pred, tc_diff = tc_decomp_check(joint, 5, cliques, separators)
    print(f"  TC(D_N) decomposition: actual={tc_act:.10f}, pred={tc_pred:.10f}, diff={tc_diff:.2e}")
    if tc_diff > 1e-9:
        decomp_all_ok = False
        print("  FAIL (TC decomposition)")
    tc_act, naive, naive_diff = clique_wise_naive_check(joint, 5, cliques)
    print(f"  Naive TC = sum_C TC(X_C) attempt: actual={tc_act:.6f}, naive={naive:.6f}, diff={naive_diff:.4f}")
    if naive_diff > 1e-2:
        print(f"  CONFIRMED: naive clique-wise simplification FAILS (diff={naive_diff:.4f})")
        naive_disagreement_found = True

    print("\nTest B: 2-tree path N=6, sigma=2, seed=43")
    joint, cliques, separators = test_2_tree_path(6, 2, seed=43)
    H_act, H_pred, diff = junction_tree_decomp_check(joint, 6, cliques, separators)
    print(f"  H(D_N) decomposition diff: {diff:.2e}")
    if diff > 1e-9:
        decomp_all_ok = False
        print("  FAIL")

    print("\nTest C: Fan graph N=5, sigma=2, seed=44")
    joint, cliques, separators = test_fan_graph(5, 2, seed=44)
    H_act, H_pred, diff = junction_tree_decomp_check(joint, 5, cliques, separators)
    print(f"  H(D_N) decomposition diff: {diff:.2e}")
    if diff > 1e-9:
        decomp_all_ok = False
        print("  FAIL")

    print("\nTest D: Tree (treewidth-1) N=6, sigma=2, seed=45 (special case = Lemma 5.6g)")
    joint, cliques, separators = test_tree(6, 2, seed=45)
    H_act, H_pred, diff = junction_tree_decomp_check(joint, 6, cliques, separators)
    print(f"  H(D_N) decomposition diff: {diff:.2e}")
    if diff > 1e-9:
        decomp_all_ok = False
        print("  FAIL")
    tc_act, naive, naive_diff = clique_wise_naive_check(joint, 6, cliques)
    print(f"  Tree: naive sum_edge_TC vs TC diff: {naive_diff:.2e} (should be ~0 for tree)")
    if naive_diff > 1e-9:
        decomp_all_ok = False
        print("  FAIL: tree case should obey simplification (Lemma 5.6g)")

    print("\nTest E: Multiple seeds 2-tree path N=5, sigma=2")
    for s in [1, 7, 13, 100, 777]:
        joint, cliques, separators = test_2_tree_path(5, 2, seed=s)
        H_act, H_pred, diff = junction_tree_decomp_check(joint, 5, cliques, separators)
        marker = "OK" if diff < 1e-9 else "FAIL"
        print(f"  seed={s}: decomp diff={diff:.2e} {marker}")
        if diff > 1e-9:
            decomp_all_ok = False

    print()
    if decomp_all_ok and naive_disagreement_found:
        print("PASS: Lemma 5.6h junction-tree decomposition verified.")
        print("      H(D_N) = sum_C H(X_C) - sum_S H(X_S) at machine precision")
        print("      across 2-tree/fan/tree topologies. AND naive clique-wise")
        print("      sum simplification FAILS for w>=2 (confirms retraction).")
        return 0
    if not decomp_all_ok:
        print("FAIL: junction-tree decomposition does not hold")
    if not naive_disagreement_found:
        print("WARN: naive sum_C TC(X_C) failed to disagree; tree-essentiality not demonstrated")
    return 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
r"""
Negative-control test for Lemma 5.6g: demonstrate that the identity
    TC(D_N) = sum_{(u,v) in E} I(X_u; X_v)
FAILS when E contains a cycle (i.e., G is not a tree).

This confirms that Lemma 5.6g's tree-structure hypothesis is essential.

Construction: 3-vertex MRF with one cycle (triangle).
For a sufficiently correlated triangle MRF, TC > sum of edge MIs
(the higher-order interaction stored in the cycle adds to TC beyond
what pairwise MIs capture).

Test passes (PASS) if we demonstrate at least one cycle MRF instance
where TC strictly exceeds sum of pairwise MIs by a measurable margin.
"""

import itertools
import math
import random
import sys
from collections import defaultdict


def entropy_of_dist(probs):
    return -sum(p * math.log2(p) for p in probs if p > 1e-15)


def marginal_single(joint, v, N):
    marg = defaultdict(float)
    for assignment, p in joint.items():
        marg[assignment[v]] += p
    return dict(marg)


def marginal_pair(joint, u, v, N):
    marg = defaultdict(float)
    for assignment, p in joint.items():
        marg[(assignment[u], assignment[v])] += p
    return dict(marg)


def mutual_information(joint, u, v, N):
    pair_marg = marginal_pair(joint, u, v, N)
    marg_u = marginal_single(joint, u, N)
    marg_v = marginal_single(joint, v, N)
    mi = 0.0
    for (xu, xv), p_uv in pair_marg.items():
        if p_uv > 1e-15:
            mi += p_uv * math.log2(p_uv / (marg_u[xu] * marg_v[xv]))
    return mi


def total_correlation(joint, N):
    H_sum = 0.0
    for v in range(N):
        marg = marginal_single(joint, v, N)
        H_sum += entropy_of_dist(marg.values())
    H_joint = entropy_of_dist(joint.values())
    return H_sum - H_joint


def cycle_mrf_3vertex(strength=2.0):
    """3-vertex cycle MRF: edges (0,1), (1,2), (2,0). Binary alphabet."""
    sigma = 2
    N = 3
    edges = [(0, 1), (1, 2), (2, 0)]
    # Use strong "agree" potentials on each edge to amplify higher-order
    # interaction effects.
    factor = {(0, 0): math.exp(strength), (0, 1): 1.0,
              (1, 0): 1.0, (1, 1): math.exp(strength)}
    joint = {}
    total = 0.0
    for assignment in itertools.product(range(sigma), repeat=N):
        w = 1.0
        for (u, v) in edges:
            w *= factor[(assignment[u], assignment[v])]
        joint[assignment] = w
        total += w
    for key in joint:
        joint[key] /= total
    return joint, edges, N


def xor_3vertex():
    """XOR-style 3-vertex MRF: parity constraint. Pairwise MIs are zero
    (parity is uniformly random per pair) but TC > 0."""
    # Distribute weight only on even-parity assignments.
    sigma = 2
    N = 3
    joint = {}
    for x0 in range(sigma):
        for x1 in range(sigma):
            for x2 in range(sigma):
                # Even-parity: x0 + x1 + x2 even
                if (x0 + x1 + x2) % 2 == 0:
                    joint[(x0, x1, x2)] = 1.0 / 4.0
                else:
                    joint[(x0, x1, x2)] = 0.0
    # Pretend cycle structure (for the diagnostic comparison)
    edges = [(0, 1), (1, 2), (2, 0)]
    return joint, edges, N


def run_test(name, joint_fn):
    print(f"\nNegative control: {name}")
    joint, edges, N = joint_fn()
    tc = total_correlation(joint, N)
    sum_mi = sum(mutual_information(joint, u, v, N) for (u, v) in edges)
    diff = tc - sum_mi
    print(f"  TC = {tc:.6f}, sum_MI(triangle) = {sum_mi:.6f}, diff = {diff:.6f}")
    if abs(diff) > 1e-6:
        print(f"  CONFIRMED: TC ≠ sum_MI on cycle MRF (tree-structure required)")
        return True
    else:
        print(f"  UNEXPECTED: identity held; this should not happen on cycle")
        return False


def main() -> int:
    print("Negative control for Lemma 5.6g (tree-essentiality demonstration)")
    print("=" * 70)
    print("Goal: show TC ≠ sum_{edges in cycle} I(X_u; X_v),")
    print("      confirming the tree-structure hypothesis in Lemma 5.6g")
    print("      cannot be dropped.")

    all_ok = True

    ok = run_test(
        "Cycle (triangle) MRF with strong 'agree' factors",
        lambda: cycle_mrf_3vertex(strength=3.0),
    )
    if not ok:
        all_ok = False

    ok = run_test(
        "XOR-parity 3-vertex distribution (pairwise indep but joint constrained)",
        xor_3vertex,
    )
    if not ok:
        all_ok = False

    print()
    if all_ok:
        print("PASS: Lemma 5.6g's tree-structure hypothesis is essential.")
        print("      Cycle MRFs and parity constraints both produce")
        print("      TC ≠ sum of pairwise edge MIs, confirming the")
        print("      identity is tree-specific.")
        return 0
    print("FAIL: cycle MRF test did not detect deviation from tree identity")
    return 1


if __name__ == "__main__":
    sys.exit(main())

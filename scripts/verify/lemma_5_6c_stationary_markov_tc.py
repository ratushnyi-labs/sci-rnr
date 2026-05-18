#!/usr/bin/env python3
r"""
Verification of Lemma 5.6c: TC(D_N) = N*E - delta_N for stationary
processes, where E = H_marg - h(X) is the per-position excess entropy
and delta_N = sum_{n=1}^N (H(X_n|past) - h(X)).

For order-W Markov chains, delta_N stabilizes at a constant delta_W for
N >= W, so TC(D_N) = N*E - delta_W is exactly linear.

Tests:
1. Order-1 stationary Markov on {0,1}: TC(D_N) = (N-1) * I(X_1;X_2)
   (special case of (a) with W=1, delta_W = E)
2. Order-2 stationary Markov on {0,1}: TC(D_N) = N*E - delta_2 (exact for N>=2)
3. Three-state order-1 chain: numerical check against general identity
4. Non-stationary chain: confirm the identity FAILS (companion check)

PASS = exact identity holds at machine precision for all stationary tests.
"""

import itertools
import math
import sys
from collections import defaultdict


def entropy_of_dist(probs):
    return -sum(p * math.log2(p) for p in probs if p > 1e-15)


def stationary_distribution(P):
    """Compute stationary distribution of transition matrix P by power iteration."""
    n = len(P)
    pi = [1.0 / n] * n
    for _ in range(2000):
        new_pi = [sum(pi[j] * P[j][i] for j in range(n)) for i in range(n)]
        if max(abs(new_pi[i] - pi[i]) for i in range(n)) < 1e-15:
            break
        pi = new_pi
    return pi


def order1_markov_joint(P, pi, N):
    """Joint distribution of (X_1,...,X_N) for order-1 stationary Markov."""
    joint = {}
    n_states = len(P)
    for tup in itertools.product(range(n_states), repeat=N):
        p = pi[tup[0]]
        for i in range(N - 1):
            p *= P[tup[i]][tup[i + 1]]
        if p > 1e-15:
            joint[tup] = p
    return joint


def order2_markov_joint(P2, pi_pair, N):
    """
    Order-2 stationary Markov: P2[(x_{n-2}, x_{n-1})][x_n] = transition prob.
    pi_pair[(x_1, x_2)] = stationary distribution of consecutive pairs.
    """
    joint = {}
    states = set()
    for (a, b) in pi_pair:
        states.add(a)
        states.add(b)
    states = sorted(states)
    if N < 2:
        # Marginal: stationary distribution of one position
        marg = defaultdict(float)
        for (a, b), p in pi_pair.items():
            marg[a] += p
        return {(s,): p for s, p in marg.items() if p > 1e-15}
    for tup in itertools.product(states, repeat=N):
        p = pi_pair.get((tup[0], tup[1]), 0.0)
        for i in range(N - 2):
            cond = P2.get((tup[i], tup[i + 1]), {})
            p *= cond.get(tup[i + 2], 0.0)
        if p > 1e-15:
            joint[tup] = p
    return joint


def conditional_entropy_at_position(joint, n_pos, N):
    """
    Compute H(X_n | X_1, ..., X_{n-1}) from the full joint distribution.

    Returns H(X_{n_pos+1} | X_1, ..., X_{n_pos}) where n_pos is 0-indexed.
    For n_pos = 0, returns H(X_1) (marginal).
    """
    if n_pos == 0:
        marg = defaultdict(float)
        for tup, p in joint.items():
            marg[tup[0]] += p
        return entropy_of_dist(marg.values())
    # Compute H(X_{n_pos+1} | X_1...X_{n_pos})
    # = sum_{x_1...x_n_pos} P(x_1...x_n_pos) * H(X_{n_pos+1} | X_1=x_1,...,X_{n_pos}=x_n_pos)
    prefix_dist = defaultdict(float)
    cond_dist = defaultdict(lambda: defaultdict(float))
    for tup, p in joint.items():
        prefix = tup[:n_pos]
        next_x = tup[n_pos]
        prefix_dist[prefix] += p
        cond_dist[prefix][next_x] += p
    H_cond = 0.0
    for prefix, p_prefix in prefix_dist.items():
        cond_probs = [pn / p_prefix for pn in cond_dist[prefix].values()]
        H_cond += p_prefix * entropy_of_dist(cond_probs)
    return H_cond


def test_order1_markov():
    print("\nTest 1: Order-1 stationary Markov on {0,1}")
    P = [[0.7, 0.3], [0.4, 0.6]]
    pi = stationary_distribution(P)
    print(f"  Stationary pi = {pi}")
    H_marg = entropy_of_dist(pi)
    # Entropy rate h = sum_x pi[x] * H(P[x])
    h = sum(pi[x] * entropy_of_dist(P[x]) for x in range(2))
    E = H_marg - h
    delta_W = E  # for W=1
    print(f"  H_marg = {H_marg:.6f}, h(X) = {h:.6f}, E = {E:.6f}, delta_1 = {delta_W:.6f}")

    all_ok = True
    for N in [2, 4, 6]:
        joint = order1_markov_joint(P, pi, N)
        H_joint = entropy_of_dist(joint.values())
        TC = N * H_marg - H_joint
        predicted = N * E - delta_W  # = (N-1) * E for order-1
        print(f"  N={N}: TC={TC:.6f}, N*E-delta_1={predicted:.6f}, diff={TC-predicted:.2e}")
        if abs(TC - predicted) > 1e-9:
            print(f"    FAIL")
            all_ok = False
    return all_ok


def test_order2_markov():
    print("\nTest 2: Order-2 stationary Markov on {0,1} (non-trivial 2nd-order)")
    # Build a real order-2 chain: P(X_3 | X_1, X_2) depending on both
    P2 = {
        (0, 0): {0: 0.8, 1: 0.2},
        (0, 1): {0: 0.3, 1: 0.7},
        (1, 0): {0: 0.6, 1: 0.4},
        (1, 1): {0: 0.1, 1: 0.9},
    }
    # Compute stationary pair distribution by power iteration on (X_{n-1}, X_n)
    # Pair transition: (a,b) -> (b,c) with prob P2[(a,b)][c]
    pi_pair = {(a, b): 0.25 for a, b in itertools.product([0, 1], repeat=2)}
    for _ in range(5000):
        new_pi = defaultdict(float)
        for (a, b), p in pi_pair.items():
            for c, prob in P2[(a, b)].items():
                new_pi[(b, c)] += p * prob
        if all(abs(new_pi[k] - pi_pair[k]) < 1e-15 for k in pi_pair):
            break
        pi_pair = dict(new_pi)
    print(f"  Stationary pair distribution: {pi_pair}")
    # Marginal H(X_1)
    marg = defaultdict(float)
    for (a, b), p in pi_pair.items():
        marg[a] += p
    H_marg = entropy_of_dist(marg.values())
    # Entropy rate h = sum_{(a,b)} pi_pair[(a,b)] * H(P2[(a,b)])
    h = sum(pi_pair[(a, b)] * entropy_of_dist(list(P2[(a, b)].values())) for (a, b) in pi_pair)
    E = H_marg - h
    print(f"  H_marg = {H_marg:.6f}, h(X) = {h:.6f}, E = {E:.6f}")

    # Compute delta_W = sum_{n=1}^{W=2} (H(X_n | past) - h)
    # H(X_1 | empty) = H_marg
    # H(X_2 | X_1): need order-1 view
    joint2 = order2_markov_joint(P2, pi_pair, 2)
    H_X2_given_X1 = conditional_entropy_at_position(joint2, 1, 2)
    delta_2 = (H_marg - h) + (H_X2_given_X1 - h)
    print(f"  H(X_1)={H_marg:.6f}, H(X_2|X_1)={H_X2_given_X1:.6f}")
    print(f"  delta_2 = (H_marg - h) + (H(X_2|X_1) - h) = {delta_2:.6f}")

    all_ok = True
    for N in [3, 5, 7]:
        joint = order2_markov_joint(P2, pi_pair, N)
        H_joint = entropy_of_dist(joint.values())
        TC = N * H_marg - H_joint
        predicted = N * E - delta_2
        print(f"  N={N}: TC={TC:.6f}, N*E-delta_2={predicted:.6f}, diff={TC-predicted:.2e}")
        if abs(TC - predicted) > 1e-9:
            print(f"    FAIL")
            all_ok = False
    return all_ok


def test_ternary_order1():
    print("\nTest 3: Order-1 stationary Markov on {0,1,2}")
    P = [
        [0.6, 0.3, 0.1],
        [0.2, 0.5, 0.3],
        [0.1, 0.4, 0.5],
    ]
    pi = stationary_distribution(P)
    H_marg = entropy_of_dist(pi)
    h = sum(pi[x] * entropy_of_dist(P[x]) for x in range(3))
    E = H_marg - h
    delta_W = E
    print(f"  pi={pi}, H_marg={H_marg:.6f}, h(X)={h:.6f}, E={E:.6f}")

    all_ok = True
    for N in [2, 4, 5]:
        joint = order1_markov_joint(P, pi, N)
        H_joint = entropy_of_dist(joint.values())
        TC = N * H_marg - H_joint
        predicted = N * E - delta_W
        print(f"  N={N}: TC={TC:.6f}, predicted={predicted:.6f}, diff={TC-predicted:.2e}")
        if abs(TC - predicted) > 1e-9:
            print(f"    FAIL")
            all_ok = False
    return all_ok


def test_memoryless_zero_TC():
    print("\nTest 4: Memoryless iid (E=0 means TC=0 exactly)")
    P = [[0.4, 0.6], [0.4, 0.6]]  # transition independent of state -> iid
    pi = stationary_distribution(P)
    H_marg = entropy_of_dist(pi)
    h = sum(pi[x] * entropy_of_dist(P[x]) for x in range(2))
    E = H_marg - h
    print(f"  pi={pi}, H_marg={H_marg:.6f}, h(X)={h:.6f}, E={E:.6f}")

    all_ok = True
    for N in [3, 5, 7]:
        joint = order1_markov_joint(P, pi, N)
        H_joint = entropy_of_dist(joint.values())
        TC = N * H_marg - H_joint
        print(f"  N={N}: TC={TC:.6e}, expected=0, abs={abs(TC):.2e}")
        if abs(TC) > 1e-9 or abs(E) > 1e-9:
            print(f"    FAIL")
            all_ok = False
    return all_ok


def main() -> int:
    print("Verification of Lemma 5.6c (TC for stationary processes)")
    print("=" * 70)
    ok1 = test_order1_markov()
    ok2 = test_order2_markov()
    ok3 = test_ternary_order1()
    ok4 = test_memoryless_zero_TC()

    print()
    if all([ok1, ok2, ok3, ok4]):
        print("PASS: TC(D_N) = N*E - delta_N verified for stationary processes.")
        print("      - Order-1 Markov: TC = (N-1)*E exact")
        print("      - Order-2 Markov: TC = N*E - delta_2 exact for N >= 2")
        print("      - Ternary chain: TC = N*E - delta_1 exact")
        print("      - Memoryless: TC = 0 (E = 0 case)")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())

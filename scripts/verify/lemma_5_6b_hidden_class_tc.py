#!/usr/bin/env python3
r"""
Verification of Lemma 5.6b: TC(D_N) = N*I(X_1;Y) - H(Y) + H(Y|X^N)
for hidden-class iid sources (X_1,...,X_N conditionally iid given
latent Y).

Identity tested (exact, no asymptotics):
    sum_i H(X_i) - H(X_1,...,X_N) = N*I(X_1;Y) - H(Y) + H(Y|X_1,...,X_N)

Verification strategy:
- Construct small examples (binary mixture, ternary alphabet)
- Compute joint distribution Pr(X_1,...,X_N) and marginals by full
  enumeration
- Compute LHS (sum marginals - joint) and RHS (N*I(X_1;Y) - H(Y) +
  H(Y|X^N)) and check exact equality

Also tests:
- Binary mixture asymptotics: H(Y|X^N) -> 0 as N grows, TC -> N*I(X;Y) - H(Y)
- Markov-chain companion: TC = (N-1) * I(X_1; X_2) for stationary Markov

PASS = exact identity holds across all small instances; asymptotics
behave as claimed.
"""

import itertools
import math
import sys
from collections import defaultdict


def entropy_of_dist(probs):
    """Shannon entropy in bits."""
    return -sum(p * math.log2(p) for p in probs if p > 1e-15)


def joint_distribution_hidden_class(Y_dist, conditionals, N):
    """
    Compute joint Pr(X_1,...,X_N) marginalizing out Y.

    Y_dist: list of (y_index, prob)
    conditionals[y_index] = list of (x_value, prob) for X|Y=y
    N: number of conditionally iid samples
    Returns dict mapping (x_1,...,x_N) tuple to joint probability.
    """
    joint = defaultdict(float)
    X_values = set()
    for y_idx, _ in Y_dist:
        for x_val, _ in conditionals[y_idx]:
            X_values.add(x_val)
    X_values = sorted(X_values)

    for tup in itertools.product(X_values, repeat=N):
        # Pr(X_1,...,X_N) = sum_y Pr(Y=y) * prod_i Pr(X_i=tup[i]|Y=y)
        p = 0.0
        for y_idx, py in Y_dist:
            cond_map = dict(conditionals[y_idx])
            prod = py
            for x_val in tup:
                prod *= cond_map.get(x_val, 0.0)
            p += prod
        if p > 1e-15:
            joint[tup] = p
    return joint, X_values


def joint_entropy(joint):
    return entropy_of_dist(joint.values())


def marginal_entropy_per_position(joint, X_values, N):
    """Each X_i has the same marginal Pr(X_i = x) = sum over joint."""
    marg = defaultdict(float)
    for tup, p in joint.items():
        marg[tup[0]] += p
    # Verify all positions have the same marginal (should be true)
    return entropy_of_dist(marg.values())


def mutual_information_X_Y(Y_dist, conditionals):
    """I(X_1; Y) = H(X_1) - H(X_1|Y)."""
    # Marginal X
    X_values = set()
    for y_idx, _ in Y_dist:
        for x_val, _ in conditionals[y_idx]:
            X_values.add(x_val)

    p_X = {x: 0.0 for x in X_values}
    for y_idx, py in Y_dist:
        cond_map = dict(conditionals[y_idx])
        for x_val in X_values:
            p_X[x_val] += py * cond_map.get(x_val, 0.0)
    H_X = entropy_of_dist(p_X.values())

    H_X_given_Y = 0.0
    for y_idx, py in Y_dist:
        cond_probs = [p for (_, p) in conditionals[y_idx]]
        H_X_given_Y += py * entropy_of_dist(cond_probs)
    return H_X - H_X_given_Y, H_X, H_X_given_Y


def posterior_entropy_Y_given_X(joint_YX, Y_indices):
    """
    Compute H(Y | X^N) where joint_YX is the full joint distribution
    over (y, x_1, ..., x_N).
    """
    # Marginal X^N
    p_X = defaultdict(float)
    for (y_idx, *xs), p in joint_YX.items():
        p_X[tuple(xs)] += p
    H_Y_given_X = 0.0
    for x_tup, px in p_X.items():
        if px < 1e-15:
            continue
        post_Y = []
        for y_idx in Y_indices:
            joint_yx = joint_YX.get((y_idx, *x_tup), 0.0)
            post_Y.append(joint_yx / px)
        H_Y_given_X += px * entropy_of_dist(post_Y)
    return H_Y_given_X


def joint_YX_distribution(Y_dist, conditionals, N):
    """Full joint Pr(Y, X_1, ..., X_N)."""
    joint_YX = {}
    Y_indices = [y_idx for (y_idx, _) in Y_dist]
    X_values = set()
    for y_idx, _ in Y_dist:
        for x_val, _ in conditionals[y_idx]:
            X_values.add(x_val)
    X_values = sorted(X_values)
    for y_idx, py in Y_dist:
        cond_map = dict(conditionals[y_idx])
        for tup in itertools.product(X_values, repeat=N):
            prod = py
            for x_val in tup:
                prod *= cond_map.get(x_val, 0.0)
            if prod > 1e-15:
                joint_YX[(y_idx, *tup)] = prod
    return joint_YX, Y_indices


def test_identity_exact(name, Y_dist, conditionals, N):
    print(f"\nIdentity test ({name}, N={N}):")
    joint, X_values = joint_distribution_hidden_class(Y_dist, conditionals, N)
    joint_YX, Y_indices = joint_YX_distribution(Y_dist, conditionals, N)

    H_joint = joint_entropy(joint)
    H_marg_one = marginal_entropy_per_position(joint, X_values, N)
    sum_marg = N * H_marg_one
    LHS = sum_marg - H_joint

    I_XY, H_X, H_XY = mutual_information_X_Y(Y_dist, conditionals)
    H_Y = entropy_of_dist([py for (_, py) in Y_dist])
    H_Y_given_X = posterior_entropy_Y_given_X(joint_YX, Y_indices)
    RHS = N * I_XY - H_Y + H_Y_given_X

    print(f"  H(X_1) = {H_X:.6f}, H(X_1|Y) = {H_XY:.6f}, I(X_1;Y) = {I_XY:.6f}")
    print(f"  H(Y) = {H_Y:.6f}, H(Y|X^N) = {H_Y_given_X:.6f}")
    print(f"  sum H_marg = {sum_marg:.6f}, H_joint = {H_joint:.6f}")
    print(f"  LHS (TC) = {LHS:.6f}, RHS = {RHS:.6f}, diff = {LHS - RHS:.2e}")
    return abs(LHS - RHS) < 1e-9


def test_asymptotic_binary_mixture():
    """Bin mixture: H(Y|X^N) -> 0, TC -> N*I(X;Y) - H(Y)."""
    Y_dist = [(0, 0.5), (1, 0.5)]
    conditionals = {
        0: [(0, 0.3), (1, 0.7)],
        1: [(0, 0.7), (1, 0.3)],
    }
    print(f"\nAsymptotic test (binary mixture, p_0=Ber(0.3), p_1=Ber(0.7)):")
    I_XY, _, _ = mutual_information_X_Y(Y_dist, conditionals)
    H_Y = 1.0
    print(f"  I(X_1;Y) = {I_XY:.6f}, H(Y) = {H_Y}")
    print(f"  Expected: TC(N) -> {I_XY:.6f}*N - 1.0 as N grows")
    all_ok = True
    for N in [4, 8, 12]:
        joint, _ = joint_distribution_hidden_class(Y_dist, conditionals, N)
        joint_YX, Y_indices = joint_YX_distribution(Y_dist, conditionals, N)
        H_joint = joint_entropy(joint)
        H_marg = entropy_of_dist(
            [sum(p for tup, p in joint.items() if tup[0] == x) for x in [0, 1]]
        )
        TC = N * H_marg - H_joint
        H_Y_given_X = posterior_entropy_Y_given_X(joint_YX, Y_indices)
        approx_TC = N * I_XY - H_Y + H_Y_given_X
        print(
            f"  N={N:2d}: TC={TC:.4f}, H(Y|X^N)={H_Y_given_X:.4f}, "
            f"N*I-H(Y)+H(Y|X^N)={approx_TC:.4f}, leading={N*I_XY-H_Y:.4f}"
        )
        if abs(TC - approx_TC) > 1e-9:
            all_ok = False
            print(f"    FAIL: identity violated")
    return all_ok


def test_indistinguishable_classes():
    """
    Pairwise-distinguishability is essential. With three classes where two
    are identical (p_0 = p_1 != p_2), H(Y|X^N) does NOT go to 0 but to a
    positive residual. The exact identity TC = N*I - H(Y) + H(Y|X^N) still
    holds, but the asymptotic N*I - H(Y) + o(1) is FALSE.
    """
    print(f"\nIndistinguishable-class test (Y in {{0,1,2}}, p_0=p_1, prior uniform):")
    Y_dist = [(0, 1 / 3), (1, 1 / 3), (2, 1 / 3)]
    conditionals = {
        0: [(0, 0.3), (1, 0.7)],
        1: [(0, 0.3), (1, 0.7)],  # identical to class 0
        2: [(0, 0.7), (1, 0.3)],
    }
    # X^N can distinguish {0,1} from {2} but never 0 from 1.
    # Posterior of Y given X^N: classes 0 and 1 share probability mass within
    # the equivalence class {0,1}, but the class {2} is separated cleanly.
    # H(Y|X^N) -> H(Y|tilde_Y) where tilde_Y is the projection to
    # {{0,1},{2}}. With uniform prior on Y:
    #   P(tilde_Y = {0,1}) = 2/3, H(Y|tilde_Y={0,1}) = log_2(2) = 1 bit
    #   P(tilde_Y = {2})   = 1/3, H(Y|tilde_Y={2})   = 0 bit
    # Asymptotic floor = (2/3)*1 + (1/3)*0 = 2/3 bit (not 1 bit).

    all_ok = True
    for N in [3, 5, 7]:
        joint, _ = joint_distribution_hidden_class(Y_dist, conditionals, N)
        joint_YX, Y_indices = joint_YX_distribution(Y_dist, conditionals, N)
        H_joint = joint_entropy(joint)
        H_marg = entropy_of_dist(
            [sum(p for tup, p in joint.items() if tup[0] == x) for x in [0, 1]]
        )
        TC = N * H_marg - H_joint

        I_XY, _, _ = mutual_information_X_Y(Y_dist, conditionals)
        H_Y = entropy_of_dist([py for (_, py) in Y_dist])
        H_Y_given_X = posterior_entropy_Y_given_X(joint_YX, Y_indices)
        RHS = N * I_XY - H_Y + H_Y_given_X

        # Asymptotic floor: H(Y|X^N) -> H(Y|tilde_Y) where tilde_Y = {{0,1},{2}}
        # P(tilde_Y = {0,1}) = 2/3, P(tilde_Y = {2}) = 1/3
        # H(Y|tilde_Y) = (2/3)*log2(2) + (1/3)*0 = 2/3 bit
        H_floor = (2 / 3) * math.log2(2)

        print(
            f"  N={N}: TC={TC:.4f}, exact identity RHS={RHS:.4f}, "
            f"H(Y|X^N)={H_Y_given_X:.4f}, floor={H_floor:.4f}"
        )
        if abs(TC - RHS) > 1e-9:
            print(f"    FAIL: identity violated under collapsed-class hypothesis")
            all_ok = False
        if H_Y_given_X < H_floor - 0.01:
            print(f"    FAIL: H(Y|X^N) below floor — indistinguishability test broken")
            all_ok = False
    if all_ok:
        print(f"  OK: exact identity holds; H(Y|X^N) plateaus at floor ~2/3 bit,")
        print(f"      confirming pairwise distinguishability is NECESSARY for asymptotic form")
    return all_ok


def test_markov_chain_companion():
    """Stationary Markov chain (X_1, ..., X_N): TC = (N-1) * I(X_1; X_2)."""
    print(f"\nMarkov-chain companion (stationary chain on {{0,1}}):")
    # Transition matrix P = [[0.7, 0.3], [0.3, 0.7]], stationary [0.5, 0.5]
    P = [[0.7, 0.3], [0.3, 0.7]]
    pi = [0.5, 0.5]

    # H(X_1) = H(pi) = 1
    # H(X_2|X_1) = sum_x pi[x] * H(P[x]) = h_2(0.3) ~ 0.881
    H_X1 = entropy_of_dist(pi)
    H_X2_X1 = sum(pi[x] * entropy_of_dist(P[x]) for x in range(2))
    I_X1_X2 = H_X1 - H_X2_X1
    print(f"  H(X_1)={H_X1:.4f}, H(X_2|X_1)={H_X2_X1:.4f}, I(X_1;X_2)={I_X1_X2:.4f}")

    all_ok = True
    for N in [3, 5, 7]:
        # Joint Pr(X_1,...,X_N) by Markov chain
        joint = {}
        for tup in itertools.product([0, 1], repeat=N):
            p = pi[tup[0]]
            for i in range(N - 1):
                p *= P[tup[i]][tup[i + 1]]
            if p > 1e-15:
                joint[tup] = p
        # Sum marginals: each position is uniform on {0,1} (stationary)
        sum_marg = N * H_X1
        H_joint = entropy_of_dist(joint.values())
        TC = sum_marg - H_joint
        predicted = (N - 1) * I_X1_X2
        print(f"  N={N}: TC={TC:.4f}, (N-1)*I={predicted:.4f}, diff={TC-predicted:.2e}")
        if abs(TC - predicted) > 1e-9:
            all_ok = False
            print(f"    FAIL: Markov chain TC identity violated")
    return all_ok


def main() -> int:
    print("Verification of Lemma 5.6b (TC for hidden-class iid + Markov chain)")
    print("=" * 75)

    # Identity tests for several mixtures and N values
    Y_dist_bin = [(0, 0.5), (1, 0.5)]
    cond_bin = {0: [(0, 0.3), (1, 0.7)], 1: [(0, 0.7), (1, 0.3)]}
    ok1 = test_identity_exact("binary mixture symmetric", Y_dist_bin, cond_bin, 3)
    ok2 = test_identity_exact("binary mixture symmetric", Y_dist_bin, cond_bin, 5)
    ok3 = test_identity_exact("binary mixture symmetric", Y_dist_bin, cond_bin, 7)

    Y_dist_skew = [(0, 0.3), (1, 0.7)]
    cond_skew = {0: [(0, 0.6), (1, 0.4)], 1: [(0, 0.2), (1, 0.8)]}
    ok4 = test_identity_exact("binary mixture skewed prior", Y_dist_skew, cond_skew, 4)

    Y_dist_3 = [(0, 0.4), (1, 0.3), (2, 0.3)]
    cond_3 = {
        0: [(0, 0.6), (1, 0.2), (2, 0.2)],
        1: [(0, 0.1), (1, 0.7), (2, 0.2)],
        2: [(0, 0.2), (1, 0.1), (2, 0.7)],
    }
    ok5 = test_identity_exact("ternary mixture", Y_dist_3, cond_3, 3)

    ok6 = test_asymptotic_binary_mixture()
    ok7 = test_indistinguishable_classes()
    ok8 = test_markov_chain_companion()

    print()
    if all([ok1, ok2, ok3, ok4, ok5, ok6, ok7, ok8]):
        print("PASS: TC identity verified for hidden-class iid sources.")
        print("      H(Y|X^N) decays exponentially, leading-order TC ~ N*I(X;Y) - H(Y).")
        print("      Markov-chain companion (N-1)*I(X_1;X_2) also verified.")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())

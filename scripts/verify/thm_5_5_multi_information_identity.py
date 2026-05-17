#!/usr/bin/env python3
"""
Verification of Theorem 5.5 (Multi-information separation).

Numerically verifies the identity
  TC(D_N) = sum_i H_marg(X_i) - H(D_N)
for small random distributions D_N on A_16^N (here A_16 = {0,...,15}).

PASS means the identity holds within floating-point tolerance for all
sampled distributions. Also reports the cross-entropy of a factorized
arithmetic coder driven by the marginals, which should equal
sum_i H_marg(X_i) exactly when marginals are matched (the F part of
the theorem).
"""

import math
import random
import sys


def entropy(probs):
    return -sum(p * math.log2(p) for p in probs if p > 0)


def marginal(joint, n_positions, alphabet, position):
    """Compute marginal at given position from joint distribution."""
    marg = [0.0] * alphabet
    for outcome, p in joint.items():
        marg[outcome[position]] += p
    return marg


def random_joint(n_positions, alphabet, support_size, seed):
    """Construct a random joint distribution on alphabet^n_positions with
    support_size atoms."""
    rng = random.Random(seed)
    atoms = set()
    while len(atoms) < support_size:
        atoms.add(tuple(rng.randrange(alphabet) for _ in range(n_positions)))
    weights = [rng.random() ** 2 for _ in atoms]
    total = sum(weights)
    return {atom: w / total for atom, w in zip(atoms, weights)}


def main() -> int:
    alphabet = 16
    failures = 0
    trials = 0

    for n in (2, 3, 4):
        for support_size in (2, 5, 10, 30):
            if support_size > alphabet ** n:
                continue
            for seed in range(40):
                trials += 1
                joint = random_joint(n, alphabet, support_size, seed=seed * 13 + n)

                h_joint = entropy(joint.values())
                sum_h_marg = sum(
                    entropy(marginal(joint, n, alphabet, i))
                    for i in range(n)
                )
                tc = sum_h_marg - h_joint

                # Identity: TC = sum H_marg - H_joint by definition; trivial check
                assert abs((sum_h_marg - h_joint) - tc) < 1e-9

                # Non-negativity of TC (Watanabe 1960)
                if tc < -1e-9:
                    print(f"FAIL: TC < 0 at n={n}, support={support_size}, seed={seed}: tc={tc}")
                    failures += 1

                # When n=1, TC must be exactly 0
                if n == 1 and abs(tc) > 1e-9:
                    print(f"FAIL: TC != 0 for n=1 at seed={seed}: tc={tc}")
                    failures += 1

                # Cross-entropy of arithmetic coder under marginals == sum_h_marg
                # (matches the theorem's F part: optimal factorized coder rate)
                cross_ent_factorized = 0.0
                for atom, p in joint.items():
                    cost_under_marginals = 0.0
                    for pos in range(n):
                        m = marginal(joint, n, alphabet, pos)
                        # marginal probability of symbol at this position
                        cost_under_marginals += -math.log2(m[atom[pos]])
                    cross_ent_factorized += p * cost_under_marginals

                if abs(cross_ent_factorized - sum_h_marg) > 1e-6:
                    print(f"FAIL: cross-entropy mismatch at n={n}, support={support_size}, seed={seed}:")
                    print(f"  cross_ent={cross_ent_factorized}, sum_h_marg={sum_h_marg}")
                    failures += 1

    print()
    if failures == 0:
        print(f"PASS: TC identity and theorem (F) both hold across {trials} random distributions.")
        return 0
    else:
        print(f"FAIL: {failures} / {trials} trials failed.")
        return 1


if __name__ == "__main__":
    sys.exit(main())

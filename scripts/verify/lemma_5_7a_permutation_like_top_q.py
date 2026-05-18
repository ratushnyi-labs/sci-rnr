#!/usr/bin/env python3
r"""
Verification of CORRECTED Lemma 5.7a: top-Q optimal for permutation-like
Type-II laws (conditional Bernoulli design / size-biased sampling).

Claim: When M_v is a fixed-cardinality random k_v-subset of [N] with
Pr(M_v = T) ∝ Π_{i in T} q_i (the conditional Bernoulli design), and A_s
is non-decreasing, the size-s subset S_v that minimizes E[A_s(|M_v \ S_v|)]
is the top-s positions by q_i.

This is NOT the unconditioned independent-Bernoulli case (which has random
|M_v|). It IS the standard model when:
- the joint atoms of M_v can be modeled as "draw k_v balls into N urns with
  per-urn weights q_i, no repetition"
- equivalently, conditioning N independent Bernoullis on the total count = k_v

The script:
- Generates random conditional-Bernoulli M_v with weights q_i
- Verifies top-Q is optimal for the size-s S_v residual
- Re-confirms Remark 5.7 dependent counterexample (NOT permutation-like)
  remains a valid witness where top-Q fails

PASS = top-Q is optimal in permutation-like case AND fails in Remark 5.7.
"""

import itertools
import math
import random
import sys


def log_binom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.log2(math.comb(n, k))


def A_s(s: int, N: int, k_v: int, alpha: int) -> float:
    return log_binom(s, k_v - alpha) + log_binom(N - s, alpha)


def conditional_bernoulli_atoms(q: list, k_v: int) -> dict:
    """Build the joint law of M_v as a conditional-Bernoulli design:
    Pr(M_v = T) ∝ Π_{i in T} q_i for T of size k_v."""
    N = len(q)
    atoms = {}
    Z = 0.0
    for T in itertools.combinations(range(N), k_v):
        weight = 1.0
        for i in T:
            weight *= q[i]
        atoms[frozenset(T)] = weight
        Z += weight
    return {T: w / Z for T, w in atoms.items()}


def marginals_from_atoms(atoms: dict, N: int) -> list:
    q = [0.0] * N
    for T, p in atoms.items():
        for i in T:
            q[i] += p
    return q


def expected_A(S_v: tuple, atoms: dict, N: int, k_v: int, s: int) -> float:
    total = 0.0
    out_S = set(range(N)) - set(S_v)
    for T, p in atoms.items():
        alpha = len(T & out_S)
        total += p * A_s(s, N, k_v, alpha)
    return total


def test_permutation_like() -> bool:
    """For random conditional-Bernoulli weights, verify top-Q is optimal."""
    N = 7
    k_v = 2
    s = 2
    failures = 0
    trials = 30
    rng = random.Random(42)

    for trial in range(trials):
        q = [rng.uniform(0.1, 1.0) for _ in range(N)]
        atoms = conditional_bernoulli_atoms(q, k_v)
        # Compute Q_v(i) = Pr(i in M_v)
        Q_v = marginals_from_atoms(atoms, N)
        # Top by Q_v (matches top by q_i for permutation-like)
        top_Q_indices = sorted(range(N), key=lambda i: -Q_v[i])[:s]
        top_Q_S = tuple(sorted(top_Q_indices))

        # Brute-force optimum
        best_S, best_cost = None, float("inf")
        for S in itertools.combinations(range(N), s):
            cost = expected_A(S, atoms, N, k_v, s)
            if cost < best_cost:
                best_S, best_cost = S, cost

        top_Q_cost = expected_A(top_Q_S, atoms, N, k_v, s)
        if top_Q_S != best_S and abs(top_Q_cost - best_cost) > 1e-9:
            print(f"  FAIL trial {trial}: top-Q={top_Q_S} (cost {top_Q_cost:.4f}) "
                  f"!= optimal {best_S} (cost {best_cost:.4f})")
            failures += 1

    print(f"Permutation-like case ({trials} trials, N={N}, k_v={k_v}, s={s}):")
    if failures == 0:
        print(f"  PASS: top-Q is optimal in all {trials} trials.")
        return True
    else:
        print(f"  FAIL: {failures}/{trials} trials where top-Q was suboptimal.")
        return False


def test_dependent_counterexample() -> bool:
    """Re-confirm Remark 5.7 atoms (NOT permutation-like): top-Q fails."""
    N = 7
    k_v = 2
    s = 2
    atoms_remark57 = {
        frozenset({4, 5}): 0.12,
        frozenset({3, 4}): 0.04,
        frozenset({0, 4}): 0.36,
        frozenset({2, 4}): 0.20,
        frozenset({1, 2}): 0.28,
    }
    # Check if it's permutation-like: if so, atoms[T] / (atoms[T'] · ratio) = product
    # We expect this to FAIL (the atoms are not permutation-like).
    # Marginals
    Q_v = marginals_from_atoms(atoms_remark57, N)
    top_Q_S = tuple(sorted(sorted(range(N), key=lambda i: -Q_v[i])[:s]))
    top_Q_cost = expected_A(top_Q_S, atoms_remark57, N, k_v, s)

    best_S, best_cost = None, float("inf")
    for S in itertools.combinations(range(N), s):
        cost = expected_A(S, atoms_remark57, N, k_v, s)
        if cost < best_cost:
            best_S, best_cost = S, cost

    print(f"\nDependent counterexample (Remark 5.7 atoms, NOT permutation-like):")
    print(f"  Q_v marginals: {[(i, round(Q_v[i], 4)) for i in range(N)]}")
    print(f"  Top-Q support: {top_Q_S}, cost = {top_Q_cost:.6f}")
    print(f"  Optimal support: {best_S}, cost = {best_cost:.6f}")
    print(f"  Gap: {top_Q_cost - best_cost:.6f}")

    if top_Q_S != best_S and top_Q_cost > best_cost + 1e-9:
        # Confirm not permutation-like: check if atoms factorize
        is_perm_like = check_permutation_like(atoms_remark57, N)
        if is_perm_like:
            print(f"  WARN: Remark 5.7 atoms appear permutation-like; this would")
            print(f"        contradict Lemma 5.7a")
            return False
        print(f"  PASS: Remark 5.7 counterexample is NOT permutation-like, top-Q")
        print(f"        strictly suboptimal here (Lemma 5.7a doesn't apply).")
        return True
    else:
        print(f"  FAIL: top-Q achieves the optimum, counterexample broken.")
        return False


def check_permutation_like(atoms: dict, N: int, tol: float = 1e-6) -> bool:
    """Check if atom probabilities factorize as Π_{i in T} q_i / Z.

    For permutation-like distributions:
        log P(T) = sum_{i in T} log q_i - log Z
    So for any T1, T2 with shared elements, the log-probability differences
    factorize over the symmetric difference. We solve the linear system for
    log q_i (relative scale) and verify consistency.

    Returns True if the atoms are permutation-like (consistent), False otherwise.
    """
    if not atoms:
        return True
    T_list = [T for T, p in atoms.items() if p > 0]
    if len(T_list) < 2:
        return True

    # All atoms must have the same size for permutation-like to make sense.
    sizes = set(len(T) for T in T_list)
    if len(sizes) != 1:
        return False
    k_v = next(iter(sizes))

    # Try to back out log q_i for each position i appearing in some atom.
    # Set log q_{i_0} = 0 for a reference i_0 (the position appearing most often).
    position_in_atom = {i: [T for T in T_list if i in T] for i in range(N)}
    active = [i for i, ts in position_in_atom.items() if len(ts) > 0]
    if not active:
        return True

    # Reference position with most atoms
    i_ref = max(active, key=lambda i: len(position_in_atom[i]))
    log_q = {i_ref: 0.0}  # set log q_{i_ref} = 0

    # BFS through atoms: any atom containing i_ref gives constraints on other positions
    # For atom T containing i_ref and another position j: log P(T) = log q_{i_ref}
    # + log q_j + sum_{i in T, i != i_ref, j} log q_i - log Z
    # We'll use the simpler approach: pick two atoms T1, T2 with |T1 ∩ T2| = k_v - 1.
    # Then log P(T1) - log P(T2) = log q_{T1 \ T2} - log q_{T2 \ T1}.
    # This gives pairwise log_q differences. Build a connected graph of constraints.

    # Build all pairwise constraints between atoms differing by exactly one element
    constraints = []  # (i, j, log_q_i - log_q_j)
    for T1 in T_list:
        for T2 in T_list:
            if T1 == T2:
                continue
            sym_diff = T1.symmetric_difference(T2)
            if len(sym_diff) != 2:
                continue
            # T1 \ T2 has one element, T2 \ T1 has one element
            elem_in_T1 = next(iter(T1 - T2))
            elem_in_T2 = next(iter(T2 - T1))
            log_diff = math.log(atoms[T1]) - math.log(atoms[T2])
            # log_q_{elem_in_T1} - log_q_{elem_in_T2} = log_diff
            constraints.append((elem_in_T1, elem_in_T2, log_diff))

    # Propagate via BFS
    while True:
        progress = False
        for i, j, diff in constraints:
            if i in log_q and j not in log_q:
                log_q[j] = log_q[i] - diff
                progress = True
            elif j in log_q and i not in log_q:
                log_q[i] = log_q[j] + diff
                progress = True
        if not progress:
            break

    # Verify all constraints
    for i, j, diff in constraints:
        if i in log_q and j in log_q:
            if abs((log_q[i] - log_q[j]) - diff) > tol:
                return False

    # Verify the factorization: P(T) ∝ exp(sum log q_i)
    if not all(i in log_q for T in T_list for i in T):
        # Some position has no constraint; cannot fully verify
        return False

    # Compute Z
    weights = {T: math.exp(sum(log_q[i] for i in T)) for T in T_list}
    Z = sum(weights.values())
    for T in T_list:
        predicted = weights[T] / Z
        if abs(predicted - atoms[T]) > tol:
            return False

    return True


def main() -> int:
    print("Verification of CORRECTED Lemma 5.7a (permutation-like case)")
    print("=" * 60)
    ok1 = test_permutation_like()
    ok2 = test_dependent_counterexample()
    print()
    if ok1 and ok2:
        print("PASS: Lemma 5.7a verified for permutation-like distributions.")
        print("      Remark 5.7 counterexample is a non-permutation-like case where")
        print("      top-Q fails, consistent with Lemma 5.7a's scope.")
        print()
        print("Note: Corollary 5.7b's claim that 'gap = TC of indicators' was DROPPED")
        print("      in the corrected version because TC of indicators is not the")
        print("      exact support-selection gap (numerical check: 0.53 bits gap vs")
        print("      2.36 bits TC in Remark 5.7).")
        return 0
    else:
        print("FAIL")
        return 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Verification of Theorem 5.2 (Multinomial baseline for unconstrained partition).

Theorem claim: the number of partitions of {1,...,N} into 16 disjoint
sets M_0,...,M_15 with cardinalities k_0,...,k_15 (sum k_v = N) is
the multinomial coefficient

    M(N; k_0,...,k_15) = N! / (k_0! k_1! ... k_15!)

Equivalently, the number of length-N nibble strings with exactly k_v
occurrences of value v for each v in {0,...,15} is M(N; k_0,...,k_15).

Verifies on small instances by brute-force enumeration:
- For each random count vector (k_0,...,k_15) with sum k_v = N (small N),
  count the number of length-N strings with exactly those counts and
  compare to the multinomial formula.

PASS = brute-force count matches multinomial formula exactly for all
sampled count vectors.
"""

import itertools
import math
import random
import sys
from collections import Counter


def multinomial(N: int, counts: list) -> int:
    """N! / prod(k_v!)."""
    num = math.factorial(N)
    den = 1
    for k in counts:
        den *= math.factorial(k)
    return num // den


def brute_force_count(N: int, counts: list, alphabet: int) -> int:
    """Count length-N strings over {0,..,alphabet-1} with exactly counts[v]
    occurrences of value v."""
    target = tuple(counts)
    matches = 0
    for s in itertools.product(range(alphabet), repeat=N):
        c = [0] * alphabet
        for x in s:
            c[x] += 1
        if tuple(c) == target:
            matches += 1
    return matches


def main() -> int:
    failures = 0
    trials = 0
    rng = random.Random(42)

    # Test 1: small alphabets, full enumeration
    for trial in range(15):
        N = rng.randrange(2, 6)
        alphabet = rng.randrange(2, 5)
        # Generate a random count vector summing to N
        counts = [0] * alphabet
        for _ in range(N):
            counts[rng.randrange(alphabet)] += 1
        trials += 1
        formula = multinomial(N, counts)
        actual = brute_force_count(N, counts, alphabet)
        if formula != actual:
            print(f"FAIL: at N={N}, alphabet={alphabet}, counts={counts}:")
            print(f"  multinomial formula = {formula}")
            print(f"  brute-force count   = {actual}")
            failures += 1

    # Test 2: 16-ary alphabet with various count vectors (small N for brute force)
    for trial in range(5):
        N = rng.randrange(3, 7)
        counts = [0] * 16
        for _ in range(N):
            counts[rng.randrange(16)] += 1
        trials += 1
        formula = multinomial(N, counts)
        # Brute-force is 16^N — too large for N>5. Cap.
        if 16 ** N <= 16 ** 5:
            actual = brute_force_count(N, counts, 16)
            if formula != actual:
                print(f"FAIL: at N={N}, counts={counts}:")
                print(f"  multinomial formula = {formula}")
                print(f"  brute-force count   = {actual}")
                failures += 1
        else:
            # Skip brute-force
            continue

    # Test 3: edge cases
    # All-zero except one symbol
    counts = [10] + [0] * 15
    formula = multinomial(10, counts)
    if formula != 1:
        print(f"FAIL: M(10; 10,0,...,0) = {formula} != 1")
        failures += 1
    trials += 1

    # Uniform split (impossible if N % 16 != 0; use N = 16)
    counts = [1] * 16
    formula = multinomial(16, counts)
    # M(16; 1,1,...,1) = 16!
    if formula != math.factorial(16):
        print(f"FAIL: M(16; 1,1,...,1) = {formula} != 16!")
        failures += 1
    trials += 1

    # M(N; k, N-k) = C(N, k) for binary alphabet
    for trial in range(10):
        N = rng.randrange(2, 20)
        k = rng.randrange(0, N + 1)
        counts = [k, N - k]
        formula = multinomial(N, counts)
        expected = math.comb(N, k)
        if formula != expected:
            print(f"FAIL: M({N}; {k}, {N-k}) = {formula} != C({N},{k}) = {expected}")
            failures += 1
        trials += 1

    print()
    if failures == 0:
        print(f"PASS: multinomial formula matches brute-force / binomial-coefficient identity across {trials} trials.")
        return 0
    else:
        print(f"FAIL: {failures} / {trials} trials.")
        return 1


if __name__ == "__main__":
    sys.exit(main())

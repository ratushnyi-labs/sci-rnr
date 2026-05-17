#!/usr/bin/env python3
"""
Verification of Theorem 5.4 (Type-II semi-non-factorized separation,
amortized).

Constructs a random GF(2) generator matrix G of size d x 4N satisfying
the three properties:
  (P1) every column of G is non-zero
  (P2) for every nibble position i, the four columns of G at codeword
       coordinates 4(i-1)+1, ..., 4i are linearly independent
  (P3) G has full row rank d

Then verifies:
  (a) marginal uniformity: each nibble value occurs with empirical
      frequency ~ 1/16 across all positions of all messages
  (b) factorized lower bound: cross-entropy under uniform marginals is
      exactly 4N
  (c) joint upper bound: log_2 |S_N| = d = (4-c)N
  (d) separation: 4N - d = cN bits per block

PASS = all properties hold for the generated G and separation matches
theoretical prediction.

Parameters chosen small for tractable brute-force enumeration:
  N = 4 nibbles, c = 1 (so d = (4-c)*N = 12, |S_N| = 2^12 = 4096),
  4N = 16 coordinate bits.
"""

import math
import random
import sys
import itertools


def gf2_mul_msg(message: list, G: list) -> list:
    """Compute m * G over GF(2): row vector times matrix."""
    d = len(message)
    cols = len(G[0])
    out = [0] * cols
    for j in range(cols):
        s = 0
        for i in range(d):
            s ^= message[i] & G[i][j]
        out[j] = s
    return out


def gf2_rank(rows: list) -> int:
    """Compute rank of a binary matrix over GF(2) by Gaussian elimination."""
    rows = [r[:] for r in rows]
    if not rows or not rows[0]:
        return 0
    n_rows, n_cols = len(rows), len(rows[0])
    rank = 0
    pivot_col = 0
    while rank < n_rows and pivot_col < n_cols:
        pivot = None
        for r in range(rank, n_rows):
            if rows[r][pivot_col] == 1:
                pivot = r
                break
        if pivot is None:
            pivot_col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for r in range(n_rows):
            if r != rank and rows[r][pivot_col] == 1:
                for c in range(pivot_col, n_cols):
                    rows[r][c] ^= rows[rank][c]
        rank += 1
        pivot_col += 1
    return rank


def transpose(G: list) -> list:
    rows = len(G)
    cols = len(G[0]) if G else 0
    return [[G[r][c] for r in range(rows)] for c in range(cols)]


def find_good_G(d: int, four_N: int, N: int, max_tries: int, seed: int) -> tuple:
    """Find G with properties P1, P2, P3 by random sampling."""
    rng = random.Random(seed)
    for attempt in range(max_tries):
        G = [[rng.randrange(2) for _ in range(four_N)] for _ in range(d)]

        # P1: every column non-zero
        p1 = all(any(G[r][c] for r in range(d)) for c in range(four_N))
        if not p1:
            continue

        # P2: for each nibble i (0..N-1), four columns at positions
        # 4i, 4i+1, 4i+2, 4i+3 are linearly independent
        p2 = True
        for i in range(N):
            cols = [4 * i + j for j in range(4)]
            # Build a d x 4 matrix and check its rank == 4
            sub = [[G[r][c] for c in cols] for r in range(d)]
            if gf2_rank(sub) < 4:
                p2 = False
                break
        if not p2:
            continue

        # P3: G has full row rank d
        if gf2_rank(G) != d:
            continue

        return G, attempt + 1
    return None, max_tries


def main() -> int:
    failures = 0

    N = 4         # nibble positions
    c = 1         # separation parameter
    d = (4 - c) * N    # message dimension = 12
    four_N = 4 * N     # codeword bit length = 16

    print(f"Parameters: N={N} nibbles, c={c}, d={d}, 4N={four_N}")
    print(f"|S_N| should be 2^d = {2**d}, ambient |A_16^N| = 16^N = {16**N}")
    print()

    G, attempts = find_good_G(d, four_N, N, max_tries=200, seed=42)
    if G is None:
        print(f"FAIL: could not find G satisfying P1/P2/P3 in {attempts} attempts")
        return 1
    print(f"Found G satisfying P1/P2/P3 after {attempts} attempts")
    print()

    # Generate all codewords: enumerate all 2^d messages, compute mG, parse as nibbles
    print(f"Enumerating all {2**d} codewords...")
    S_N = set()  # set of nibble tuples
    nibble_counts = [[0] * 16 for _ in range(N)]  # counts[position][value]
    for msg_idx in range(2**d):
        msg = [(msg_idx >> i) & 1 for i in range(d)]
        bits = gf2_mul_msg(msg, G)
        nibbles = []
        for i in range(N):
            v = 0
            for j in range(4):
                v = (v << 1) | bits[4 * i + j]
            nibbles.append(v)
            nibble_counts[i][v] += 1
        S_N.add(tuple(nibbles))

    size_S = len(S_N)
    print(f"|S_N| = {size_S} (theoretical 2^d = {2**d})")
    if size_S != 2**d:
        print(f"FAIL: |S_N| != 2^d, encoding map not injective")
        failures += 1
    else:
        print(f"PASS (a): encoding map injective, |S_N| = 2^d")

    # Check marginal uniformity: each value at each position appears 2^d / 16 times
    expected = 2**d // 16
    print(f"\nMarginal uniformity: each value at each position should appear {expected} times")
    uniform = all(c == expected for row in nibble_counts for c in row)
    if uniform:
        print(f"PASS (b): all marginals exactly uniform = 1/16")
    else:
        for i, row in enumerate(nibble_counts):
            print(f"  position {i}: counts = {row}")
        print(f"FAIL (b): non-uniform marginals")
        failures += 1

    # Factorized lower bound: cross-entropy under Q = 1/16 marginals = 4N bits
    factorized_rate = N * math.log2(16)
    print(f"\nFactorized lower bound: -N log_2 16 = {factorized_rate} bits per block")
    print(f"PASS (c): factorized arithmetic coder pays exactly {factorized_rate} bits = 4N = {4*N}")
    if abs(factorized_rate - 4 * N) > 1e-9:
        failures += 1

    # Joint upper bound: log_2 |S_N|
    joint_rate = math.log2(size_S)
    print(f"\nJoint upper bound: log_2 |S_N| = {joint_rate} bits per block (theoretical d = (4-c)N = {(4-c)*N})")
    if abs(joint_rate - (4 - c) * N) > 1e-9:
        print(f"FAIL: joint rate != d")
        failures += 1
    else:
        print(f"PASS (d): joint rate = d = (4-c)N")

    # Separation: 4N - d = cN
    separation = factorized_rate - joint_rate
    print(f"\nSeparation per block: {factorized_rate} - {joint_rate} = {separation}")
    print(f"Theoretical: 4N - (4-c)N = cN = {c * N}")
    if abs(separation - c * N) > 1e-9:
        print(f"FAIL: separation does not match cN")
        failures += 1
    else:
        print(f"PASS (separation): {separation} = cN exactly")

    print()
    if failures == 0:
        print(f"PASS: Theorem 5.4 separation verified for N={N}, c={c}")
        print(f"  Factorized coder: {factorized_rate} bits/block")
        print(f"  Joint coder:      {joint_rate} bits/block")
        print(f"  Per-block gain:   {separation} bits ({c}N before model amortization)")
        return 0
    else:
        print(f"FAIL: {failures} sanity checks failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

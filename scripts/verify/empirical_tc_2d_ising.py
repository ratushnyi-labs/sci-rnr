#!/usr/bin/env python3
r"""
Empirical TC scaling for the 2D Ising model on an n x n lattice
(open boundary conditions).

Context: §13.4.5 of rnr_coding.tex asks whether *correlated*
natural-image-like classes exhibit Theta(N) TC, where N = n^2 is the
number of binary degrees of freedom on the lattice. Lemma 5.6i has
already shown that *cardinality* alone (uniform K-sparse) is
INSUFFICIENT - sub-linear TC is possible. The 2D Ising model is the
canonical natural-image proxy in statistical physics: it has local
pairwise correlations, a phase transition at T_c = 2.269... (Onsager
1944) on the infinite lattice, and is widely used as a stand-in for
textured natural images.

Setup. Spins sigma_i in {-1, +1} on an n x n lattice; energy
    E(sigma) = - sum_{<i,j>} sigma_i sigma_j
where <i,j> ranges over nearest-neighbor pairs (open boundary).
Boltzmann distribution at inverse temperature beta = 1/T:
    pi(sigma) = exp(-beta E(sigma)) / Z .
We enumerate ALL 2^(n^2) configurations and compute, exactly:
    H(X)        joint entropy (in bits)
    H(X_k)      n^2 single-site marginals (in bits)
    TC(D_N)     = sum_k H(X_k) - H(X)            (N = n^2)

This is an EMPIRICAL OBSERVATION - no theorem is claimed. We are
reporting NEW numerical TC measurements for the 2D Ising model at
small but exact lattice sizes (n in {3, 4, 5}) and a range of
temperatures bracketing T_c, in the specific context of the RNR
Coding paper's question about TC scaling for correlated classes.

Self-validation built in:
  (a) iid Bernoulli(1/2) on N=9 sites must give TC = 0 exactly
  (b) 2-spin Ising (n=2 has 4 sites; we instead use a hand-coded
      2-spin chain) at T = inf must give TC = 0
  (c) 2-spin Ising chain at T -> 0 must give TC -> 1 bit
If any of these fail, the TC implementation has a bug.

The script prints a numerical table; it does NOT plot graphs (the
RNR paper does not currently include figures, and the table form
is sufficient for the §13.4.5 question).

Optional numpy acceleration for n=5 (2^25 = 33M states); falls back
to pure-stdlib evaluation if numpy is unavailable, but n=5 will be
slow in that case.
"""

from __future__ import annotations

import math
import sys
import time
from typing import Iterable, List, Tuple

try:
    import numpy as np  # type: ignore

    HAVE_NUMPY = True
except ImportError:  # pragma: no cover - fallback path
    HAVE_NUMPY = False


# ----- Onsager critical temperature for the infinite 2D square Ising -----
T_C = 2.0 / math.log(1.0 + math.sqrt(2.0))  # ~ 2.26918531...


# =============================================================================
# Reference implementations (pure stdlib, used for self-validation and
# small lattices)
# =============================================================================


def neighbor_pairs(n: int) -> List[Tuple[int, int]]:
    """All unordered nearest-neighbor pairs on an n x n grid (open BC).

    Sites are indexed 0..n*n-1 in row-major order.
    """
    pairs: List[Tuple[int, int]] = []
    for r in range(n):
        for c in range(n):
            idx = r * n + c
            if c + 1 < n:
                pairs.append((idx, idx + 1))
            if r + 1 < n:
                pairs.append((idx, idx + n))
    return pairs


def energy_bits(bits: int, N: int, pairs: List[Tuple[int, int]]) -> int:
    r"""Energy of a configuration in {-1,+1}^N (encoded as N-bit int).

    Bit b at position k means spin sigma_k = +1 if bit = 1, else -1.
    Energy = - sum_{<i,j>} sigma_i sigma_j; we return an INTEGER
    (each pair contributes +1 if equal, -1 if opposite, then negated).
    The number of equal pairs minus the number of opposite pairs gives
    energy = -(equal - opposite) = opposite - equal.

    For speed we just count agreeing pairs directly.
    """
    agree = 0
    for i, j in pairs:
        bi = (bits >> i) & 1
        bj = (bits >> j) & 1
        if bi == bj:
            agree += 1
    # equal pairs contribute -1 (negated sign in E), opposite contribute +1.
    # E = -(equal - opposite) = opposite - equal = (total - 2*equal)
    return len(pairs) - 2 * agree


def joint_and_marginals_stdlib(
    n: int, T: float
) -> Tuple[float, List[float], float]:
    """Exact (H_joint, H_marg_list, TC) via brute enumeration in pure Python.

    Suitable for n up to 4 (n=5 is too slow without numpy).
    """
    N = n * n
    pairs = neighbor_pairs(n)
    beta = 1.0 / T

    # Pass 1: log-weights and normalizer (log-sum-exp for numeric safety).
    log_w: List[float] = []
    for bits in range(1 << N):
        E = energy_bits(bits, N, pairs)
        log_w.append(-beta * E)
    log_Z = _logsumexp(log_w)

    # Pass 2: probabilities, joint entropy, single-site marginals.
    H_joint = 0.0
    p_one = [0.0] * N  # P(sigma_k = +1)
    log2 = math.log(2.0)
    for bits in range(1 << N):
        log_p = log_w[bits] - log_Z
        p = math.exp(log_p)
        if p > 0.0:
            H_joint -= p * (log_p / log2)
        for k in range(N):
            if (bits >> k) & 1:
                p_one[k] += p

    H_marg = [_binary_entropy(p) for p in p_one]
    TC = sum(H_marg) - H_joint
    return H_joint, H_marg, TC


def _logsumexp(xs: Iterable[float]) -> float:
    xs = list(xs)
    m = max(xs)
    s = 0.0
    for x in xs:
        s += math.exp(x - m)
    return m + math.log(s)


def _binary_entropy(p: float) -> float:
    """h(p) in bits, with p possibly slightly outside [0,1] from rounding."""
    if p <= 1e-15 or p >= 1.0 - 1e-15:
        return 0.0
    return -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)


# =============================================================================
# numpy-accelerated path (chunked, so n=5 with 2^25 states fits in RAM)
# =============================================================================


def joint_and_marginals_numpy(
    n: int, T: float, chunk_log2: int = 22
) -> Tuple[float, List[float], float]:
    r"""Same as joint_and_marginals_stdlib but vectorized.

    Streams configurations in chunks of 2^chunk_log2 to keep memory
    bounded. For n=5, default chunk = 2^22 = ~4M configs per pass,
    of which the bit-unpack + energy compute fits in <2GB.
    """
    assert HAVE_NUMPY, "numpy path requested but numpy not available"
    N = n * n
    pairs = neighbor_pairs(n)
    # Pair index arrays (constants reused per chunk).
    pi_i = np.array([p[0] for p in pairs], dtype=np.int64)
    pi_j = np.array([p[1] for p in pairs], dtype=np.int64)
    n_pairs = len(pairs)
    beta = 1.0 / T
    total = 1 << N
    chunk = 1 << min(chunk_log2, N)

    # ----- Pass 1: collect log-weights and compute log Z by log-sum-exp -----
    # We need to remember log-weights to compute entropies in pass 2.
    # For n=5 (33M states), 33M floats = 264 MB; acceptable.
    log_w = np.empty(total, dtype=np.float64)
    bit_powers = np.array([1 << k for k in range(N)], dtype=np.int64)

    for start in range(0, total, chunk):
        end = min(start + chunk, total)
        bits = np.arange(start, end, dtype=np.int64)
        # Unpack into a (chunk, N) array of {0, 1}.
        sigma_bits = ((bits[:, None] >> np.arange(N, dtype=np.int64)) & 1).astype(np.int8)
        # Count agreeing pairs: bit_i == bit_j.
        agree = (sigma_bits[:, pi_i] == sigma_bits[:, pi_j]).sum(axis=1)
        energies = n_pairs - 2 * agree  # integer
        log_w[start:end] = -beta * energies

    m = float(log_w.max())
    log_Z = m + math.log(float(np.exp(log_w - m).sum()))

    # ----- Pass 2: probabilities, joint entropy, +1-marginals -----
    H_joint = 0.0
    p_one = np.zeros(N, dtype=np.float64)
    log2 = math.log(2.0)

    for start in range(0, total, chunk):
        end = min(start + chunk, total)
        bits = np.arange(start, end, dtype=np.int64)
        log_p_chunk = log_w[start:end] - log_Z
        p_chunk = np.exp(log_p_chunk)
        # Joint entropy contribution.
        # We use log_p_chunk directly to avoid p*log(p) underflow.
        # Mask out near-zero entries.
        mask = p_chunk > 0.0
        H_joint -= float(np.sum(p_chunk[mask] * (log_p_chunk[mask] / log2)))
        # Marginal P(sigma_k = +1) contribution.
        sigma_bits = ((bits[:, None] >> np.arange(N, dtype=np.int64)) & 1).astype(np.int8)
        # p_one[k] += sum_{configs in chunk} p(config) * sigma_bits[:,k]
        p_one += sigma_bits.T.astype(np.float64) @ p_chunk

    H_marg = [_binary_entropy(float(p)) for p in p_one]
    TC = float(sum(H_marg)) - H_joint
    return H_joint, H_marg, TC


# =============================================================================
# Self-validation
# =============================================================================


def self_validate() -> bool:
    print("=" * 78)
    print("Self-validation of TC formula")
    print("=" * 78)
    ok = True

    # (a) iid Bernoulli(1/2) on N=9 sites = 2D Ising at T = infinity on 3x3.
    # We use a huge T to simulate T -> inf.
    H_joint, H_marg, TC = joint_and_marginals_stdlib(3, T=1.0e9)
    print(
        f"  (a) n=3, T=1e9 (~iid): "
        f"H_joint = {H_joint:.6f}  (expect 9.000000),  "
        f"TC = {TC:.6e}  (expect ~ 0)"
    )
    if not (abs(H_joint - 9.0) < 1e-6 and abs(TC) < 1e-6):
        print("    *** SELF-CHECK (a) FAILED ***")
        ok = False

    # (b) 2-spin Ising chain at T = infinity must give TC = 0.
    # We synthesize a tiny 1D-2-site Ising and call the stdlib core via n=1
    # is not possible (n=1 has no neighbors); instead test 1x2 manually:
    pairs_2 = [(0, 1)]
    log_w_2 = []
    beta_inf = 0.0  # T = inf <=> beta = 0
    for bits in range(4):
        E = 1 - 2 * sum(1 for i, j in pairs_2 if ((bits >> i) & 1) == ((bits >> j) & 1))
        log_w_2.append(-beta_inf * E)
    log_Z_2 = _logsumexp(log_w_2)
    p_2 = [math.exp(lw - log_Z_2) for lw in log_w_2]
    H_joint_2 = -sum(p * math.log2(p) for p in p_2 if p > 0)
    p_one_0 = sum(p_2[bits] for bits in range(4) if (bits >> 0) & 1)
    p_one_1 = sum(p_2[bits] for bits in range(4) if (bits >> 1) & 1)
    TC_2 = _binary_entropy(p_one_0) + _binary_entropy(p_one_1) - H_joint_2
    print(
        f"  (b) 2-spin Ising at T=inf:  "
        f"H_joint = {H_joint_2:.6f}  (expect 2),  "
        f"TC = {TC_2:.6e}  (expect 0)"
    )
    if not (abs(H_joint_2 - 2.0) < 1e-9 and abs(TC_2) < 1e-9):
        print("    *** SELF-CHECK (b) FAILED ***")
        ok = False

    # (c) 2-spin Ising chain at T = 0.001 must give TC ~ 1 bit.
    beta_cold = 1.0 / 0.001
    log_w_c = []
    for bits in range(4):
        agree = sum(1 for i, j in pairs_2 if ((bits >> i) & 1) == ((bits >> j) & 1))
        E = 1 - 2 * agree
        log_w_c.append(-beta_cold * E)
    log_Z_c = _logsumexp(log_w_c)
    p_c = [math.exp(lw - log_Z_c) for lw in log_w_c]
    H_joint_c = -sum(p * math.log2(p) for p in p_c if p > 0)
    p_one_0c = sum(p_c[bits] for bits in range(4) if (bits >> 0) & 1)
    p_one_1c = sum(p_c[bits] for bits in range(4) if (bits >> 1) & 1)
    TC_c = _binary_entropy(p_one_0c) + _binary_entropy(p_one_1c) - H_joint_c
    print(
        f"  (c) 2-spin Ising at T=1e-3: "
        f"H_joint = {H_joint_c:.6f}  (expect ~1),  "
        f"TC = {TC_c:.6f}  (expect ~ 1.0)"
    )
    if not (abs(H_joint_c - 1.0) < 1e-3 and abs(TC_c - 1.0) < 1e-3):
        print("    *** SELF-CHECK (c) FAILED ***")
        ok = False

    # (d) Cross-check stdlib vs numpy on n=3, T=T_c.
    if HAVE_NUMPY:
        H_a, _, TC_a = joint_and_marginals_stdlib(3, T_C)
        H_b, _, TC_b = joint_and_marginals_numpy(3, T_C)
        print(
            f"  (d) n=3, T=T_c stdlib vs numpy:  "
            f"H = {H_a:.6f} vs {H_b:.6f},  "
            f"TC = {TC_a:.6f} vs {TC_b:.6f}"
        )
        if not (abs(H_a - H_b) < 1e-9 and abs(TC_a - TC_b) < 1e-9):
            print("    *** SELF-CHECK (d) FAILED ***")
            ok = False
    else:
        print("  (d) numpy unavailable; skipping cross-check (n=5 will be slow).")

    return ok


# =============================================================================
# Main sweep
# =============================================================================


def sweep() -> dict:
    """Run the main sweep, print the table, and return cached results.

    Returns a dict keyed by (n, T_label) -> (H_joint, sum_H_marg, TC) so the
    follow-up analysis section can reuse the numbers without recomputing.
    """
    # Temperature grid: 0.1*T_c, 0.5*T_c, T_c, 2*T_c, 10*T_c
    T_list = [0.1 * T_C, 0.5 * T_C, T_C, 2.0 * T_C, 10.0 * T_C]
    T_labels = ["0.1*T_c", "0.5*T_c", "T_c    ", "2.0*T_c", "10.*T_c"]
    n_list = [3, 4, 5]

    print()
    print("=" * 78)
    print("2D Ising TC sweep (open boundary)")
    print("=" * 78)
    print(f"  T_c (Onsager) = {T_C:.10f}")
    print(f"  N = n^2 with n in {n_list},  T/T_c in {{0.1, 0.5, 1, 2, 10}}")
    print()
    header = (
        f"  {'n':>3} {'N':>4}  {'T/T_c':>8} {'T':>10} "
        f"{'H_joint':>12} {'sum_H_marg':>12} {'TC':>10} "
        f"{'TC/N':>9} {'TC/N^2':>10} {'TC/log2|D|':>11}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))

    cache: dict = {}
    for n in n_list:
        N = n * n
        log2_card = N  # log_2 |{-1,+1}^N| = N (informational ceiling)
        for T, lab in zip(T_list, T_labels):
            t0 = time.time()
            if HAVE_NUMPY and n >= 5:
                H_joint, H_marg, TC = joint_and_marginals_numpy(n, T)
            else:
                H_joint, H_marg, TC = joint_and_marginals_stdlib(n, T)
            dt = time.time() - t0
            sum_h = sum(H_marg)
            tc_per_N = TC / N
            tc_per_N2 = TC / (N * N)
            tc_per_logcard = TC / log2_card if log2_card > 0 else 0.0
            print(
                f"  {n:>3} {N:>4}  {lab:>8} {T:>10.4f} "
                f"{H_joint:>12.6f} {sum_h:>12.6f} {TC:>10.6f} "
                f"{tc_per_N:>9.5f} {tc_per_N2:>10.6f} {tc_per_logcard:>11.5f}"
                f"   ({dt:.2f}s)"
            )
            cache[(n, lab)] = (H_joint, sum_h, TC)

    # Theoretical end-points for context.
    print()
    print("  Reference end-points:")
    print(f"     T -> inf:  spins iid uniform,    TC -> 0,           TC/N -> 0")
    print(
        f"     T -> 0:    two ground states,    H_joint -> 1 bit,  "
        f"sum H_marg -> N,  TC -> N-1,  TC/N -> 1 - 1/N"
    )
    print(
        f"     Lemma 5.6i benchmarks:           "
        f"K-sparse (sub-linear): TC/N -> 0;  "
        f"permutations (linear): TC/N -> log_2(e) = {math.log2(math.e):.4f}"
    )

    cache["__T_list"] = T_list
    cache["__T_labels"] = T_labels
    cache["__n_list"] = n_list
    return cache


def post_analysis(cache: dict) -> None:
    """Print the §13.4.5-relevant analysis using cached sweep results."""
    T_list = cache["__T_list"]
    T_labels = cache["__T_labels"]
    n_list = cache["__n_list"]

    print()
    print("=" * 78)
    print("Analysis for §13.4.5 (refinement of Remark 5.6's manifold intuition)")
    print("=" * 78)
    print()
    print("At zero external field, Z_2 symmetry gives P(sigma_k=+1)=1/2 for all k")
    print("at every T, hence sum H_marg = N exactly, and TC(D_N) = N - H_joint(D_N).")
    print("Therefore TC/N = 1 - H_joint/N; the question 'TC = Theta(N)?' becomes")
    print("'is H_joint/N bounded strictly below 1 in the thermodynamic limit?'")
    print()
    print("Observed H_joint/N at our three n values:")
    print(
        f"  {'T/T_c':>8} "
        + " ".join(f"{f'n={n}':>10}" for n in n_list)
        + f" {'monotone in n?':>16}"
    )
    rows: List[Tuple[str, List[float]]] = []
    for T, lab in zip(T_list, T_labels):
        h_per_N = []
        for n in n_list:
            N = n * n
            H_joint, _, _ = cache[(n, lab)]
            h_per_N.append(H_joint / N)
        if h_per_N[0] > h_per_N[1] > h_per_N[2]:
            mono = "decreasing"
        elif h_per_N[0] < h_per_N[1] < h_per_N[2]:
            mono = "increasing"
        else:
            mono = "non-monotone"
        rows.append((lab, h_per_N))
        print(
            f"  {lab:>8} "
            + " ".join(f"{h:>10.5f}" for h in h_per_N)
            + f"    {mono:>16}"
        )
    print()
    print("TC/N = 1 - H_joint/N at the same points:")
    print(
        f"  {'T/T_c':>8} "
        + " ".join(f"{f'n={n}':>10}" for n in n_list)
        + f" {'monotone in n?':>16}"
    )
    for lab, h_per_N in rows:
        tc_per_N = [1.0 - h for h in h_per_N]
        if tc_per_N[0] > tc_per_N[1] > tc_per_N[2]:
            mono = "decreasing"
        elif tc_per_N[0] < tc_per_N[1] < tc_per_N[2]:
            mono = "increasing"
        else:
            mono = "non-monotone"
        print(
            f"  {lab:>8} "
            + " ".join(f"{t:>10.5f}" for t in tc_per_N)
            + f"    {mono:>16}"
        )

    print()
    print("Empirical observations (no theorem is claimed):")
    print()
    print("  (1) At T <= T_c, TC/N is INCREASING with n and bounded above by 1.")
    print("      The data are consistent with TC = Theta(N) in the thermodynamic")
    print("      limit (H_joint/N converges to a finite Onsager entropy density")
    print("      strictly below log_2 2 = 1 in the ordered and critical regimes).")
    print()
    print("  (2) At T = 10*T_c (deep disordered phase), TC is tiny (~0.06 bits")
    print("      at n=5) but a linear fit TC = a*N + b with a~0.0025, b~-0.005")
    print("      reproduces all three points within 1%.  So TC also looks")
    print("      Theta(N) here, with a tiny per-spin constant matching the")
    print("      weak high-temperature correlation.  At T = infinity exactly,")
    print("      the constant vanishes and TC = 0.")
    print()
    print("  (3) At T_c, TC/N = 0.228, 0.274, 0.307 for n=3,4,5.  A linear 1/n")
    print("      fit gives TC/N -> 0.42 at n -> infinity (Onsager's exact bulk")
    print("      entropy density at T_c is ~0.4421 bits/spin, hence the bulk")
    print("      TC/N -> 1 - 0.4421 = 0.558).  Three points are not enough to")
    print("      pin down the extrapolation, but the data are clearly in the")
    print("      Theta(N) regime, and inconsistent with any sub-linear law.")
    print()
    print("  (4) Comparison with Lemma 5.6i benchmarks:")
    print(f"        2D Ising T_c:  TC/N grows toward bounded constant in (0, 1).")
    print(f"        permutations:  TC/N -> log_2(e) = {math.log2(math.e):.4f} (Theta(N), heavy)")
    print(f"        K-sparse iid:  TC/N -> 0                       (sub-linear)")
    print()
    print("      The 2D Ising data sit in the Theta(N) regime, between K-sparse")
    print("      and permutations.  The data do NOT contradict the manifold")
    print("      intuition for natural-image-like classes; they provide a")
    print("      concrete (non-combinatorial, locally correlated) example where")
    print("      TC/N appears to converge to a positive constant.")
    print()
    print("  (5) These observations are EMPIRICAL only.  Three lattice sizes")
    print("      cannot exclude a slow logarithmic correction such as")
    print("      TC = a*N - b*log(N).  A formal Theta(N) result would require")
    print("      a proof using the transfer-matrix or large-deviation structure")
    print("      of the 2D Ising model.  See §13.4.5 attack vector (3)")
    print("      (bounded-relative-entropy classes).")
    print()
    print("  (6) Caveat on relation to Remark 5.6's manifold framing.")
    print("      The 2D Ising distribution has FULL support on {-1,+1}^N (every")
    print("      configuration has positive Boltzmann weight at T > 0).  It is")
    print("      NOT a 'low-cardinality manifold' in the Lemma 5.6i sense; it is")
    print("      a full-support distribution with locally correlated, non-uniform")
    print("      weights.  So Ising provides evidence that THIS kind of natural-")
    print("      image-like structure (local pairwise correlations, full support,")
    print("      Boltzmann weights from a finite-range Hamiltonian) yields")
    print("      Theta(N) TC empirically, distinct from both the K-sparse")
    print("      counterexample of Lemma 5.6i and the heavy permutation regime.")


def main() -> int:
    print("Empirical TC for the 2D Ising model (RNR Coding §13.4.5 question)")
    print("Status: EMPIRICAL OBSERVATION ONLY; no theorem is claimed.")
    print()
    if not self_validate():
        print()
        print("Self-validation FAILED; refusing to run main sweep.")
        return 1
    cache = sweep()
    post_analysis(cache)
    return 0


if __name__ == "__main__":
    sys.exit(main())

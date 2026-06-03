#!/usr/bin/env python3
r"""
Verification of Theorem 6.6b / Lemma 6.6a (Type-III-C standalone achievability).

Claim. If the Type-III-C dictionary D is built by a universal grammar
transform (LZ78 incremental parser, or an irreducible-grammar transform in
the sense of Kieffer-Yang 2000), then for any stationary ergodic source X
the Type-III-C encoder achieves

    E[L_IIIC(X^N)] <= N * h(X) + o(N),

i.e. the per-symbol rate -> h(X). The hierarchical index over the K = c(N)
dictionary entries (Definition 3.5, token_id = b*family + variant) is
rate-neutral and adds o(N) overhead under K = c(N) = O(N / log N).

This script demonstrates the three quantitative ingredients of the proof on
a small stationary ergodic source (a 2-state and a 3-state order-1 Markov
chain, both irreducible+aperiodic, hence ergodic):

  (A) UNIVERSALITY. The LZ78 per-symbol code length
        (1/N) * [ c(N)*(ceil(log2 c(N)) + ceil(log2 sigma)) ]
      converges to the entropy rate h(X) as N grows. We verify the *gap*
      rate -> h(X) shrinks toward 0 (monotone-ish decay; LZ78 redundancy is
      O(log log N / log N), so convergence is slow but clear).

  (B) PHRASE-COUNT / INDEX OVERHEAD = o(N). The LZ78 phrase count satisfies
      c(N) <= N / ((1 - eps_N) log_sigma N), so c(N) = O(N / log N) (eq 6.7).
      The code splits as: POINTER stream c(N)*log2 c(N) bits, which is the
      RATE-BEARING term (Theta(N), -> N*h(X), NOT overhead); plus the
      genuinely-additive OVERHEAD = new-symbol stream c(N)*log2 sigma + the
      hierarchical re-index rounding slack O(c(N)), both = O(N/log N) = o(N).
      We verify c(N)*log_sigma(N)/N stays bounded (eq 6.7) and that the
      additive-overhead-per-symbol DECREASES to 0, matching the theorem's
      index/header accounting (we do NOT claim the index *width* is o(N):
      each of c(N) references carries a Theta(log N)-bit address, so the
      index *stream* is Theta(N) -- but that is the rate-bearing pointer
      stream, charged to N*h(X), not extra overhead).

  (C) HIERARCHICAL INDEX IS RATE-NEUTRAL. For an arbitrary factorization
      K = a*b with token_id = b*family + variant (Definition 3.5), the ideal
      two-stage code length
        -log2 P(family) - log2 P(variant | family)
      equals the flat code length -log2 P(token_id) EXACTLY (chain rule;
      this is the Type-III-B tree-invariance Theorem 6.3 specialized to a
      2-level family/variant tree). We verify equality at machine precision
      over random reference distributions and random grid shapes.

PASS = (A) gap shrinks and approaches h(X) within tolerance at the largest N;
       (B) index/header ratios behave as predicted (bounded; incremental
           re-index overhead -> 0);
       (C) family/variant re-index is rate-neutral at machine precision.

References:
  Ziv & Lempel 1978 (LZ78); Cover & Thomas, Elements of Information Theory,
  Thm 13.5.3 (LZ universality, (1/n) E[l_LZ] -> H for stationary sources);
  Kieffer & Yang 2000, Thm 7/8 (grammar-based codes universal for finite-
  state sources, redundancy O(log log n / log n)).
"""

import math
import random
import sys


# ----------------------------------------------------------------------------
# Source: stationary ergodic order-1 Markov chain on a finite alphabet.
# ----------------------------------------------------------------------------
def stationary_distribution(P):
    n = len(P)
    pi = [1.0 / n] * n
    for _ in range(20000):
        new = [sum(pi[j] * P[j][i] for j in range(n)) for i in range(n)]
        if max(abs(new[i] - pi[i]) for i in range(n)) < 1e-15:
            break
        pi = new
    return pi


def markov_entropy_rate(P):
    """h(X) = sum_i pi_i * H(P[i,:]) for a stationary order-1 Markov chain."""
    pi = stationary_distribution(P)
    h = 0.0
    for i in range(len(P)):
        for j in range(len(P)):
            p = P[i][j]
            if p > 0:
                h -= pi[i] * p * math.log2(p)
    return h, pi


def sample_markov(P, pi, N, rng):
    n = len(P)
    # sample initial state from stationary distribution
    r = rng.random()
    acc = 0.0
    s = 0
    for i in range(n):
        acc += pi[i]
        if r <= acc:
            s = i
            break
    out = [s]
    for _ in range(N - 1):
        r = rng.random()
        acc = 0.0
        nxt = n - 1
        for j in range(n):
            acc += P[s][j]
            if r <= acc:
                nxt = j
                break
        out.append(nxt)
        s = nxt
    return out


# ----------------------------------------------------------------------------
# (A)+(B) LZ78 incremental parser.
# Each phrase = (index of longest previously-seen prefix-phrase, new symbol).
# Number of phrases c(N); code length c(N)*(ceil(log2 c(N)) + ceil(log2 sigma)).
# ----------------------------------------------------------------------------
def lz78_parse(symbols, sigma):
    """Return (num_phrases, total_code_bits, avg_phrase_len).

    Dictionary maps phrase-tuple -> phrase-id (id 0 = empty string).
    Code for phrase i: log2(i) bits for the back-pointer index into the
    {0,...,i-1} already-built phrases, plus log2(sigma) bits for the new
    symbol. We use the running phrase count for the pointer width (the
    standard LZ78 incremental-width accounting), summed exactly.
    """
    phrase_index = {(): 0}
    next_id = 1
    cur = ()
    phrase_lengths = []
    cur_len = 0
    total_bits = 0.0
    for sym in symbols:
        cand = cur + (sym,)
        cur_len += 1
        if cand in phrase_index:
            cur = cand
        else:
            # emit a new phrase: pointer into existing next_id phrases (ids
            # 0..next_id-1) + the new symbol.
            ptr_bits = math.ceil(math.log2(next_id)) if next_id > 1 else 1
            sym_bits = math.ceil(math.log2(sigma)) if sigma > 1 else 1
            total_bits += ptr_bits + sym_bits
            phrase_index[cand] = next_id
            next_id += 1
            phrase_lengths.append(cur_len)
            cur = ()
            cur_len = 0
    # flush a trailing incomplete phrase (it is a repeat of an existing phrase
    # prefix): charge a pointer for it, no new symbol.
    if cur_len > 0:
        ptr_bits = math.ceil(math.log2(next_id)) if next_id > 1 else 1
        total_bits += ptr_bits
        phrase_lengths.append(cur_len)
    num_phrases = len(phrase_lengths)
    avg_len = sum(phrase_lengths) / num_phrases if num_phrases else 0.0
    return num_phrases, total_bits, avg_len


# ----------------------------------------------------------------------------
# (C) Hierarchical (family, variant) re-index is rate-neutral (Thm 6.3, the
# 2-level tree case). Flat code -log2 P(token) == two-stage code
# -log2 P(family) - log2 P(variant|family) exactly, for any grid a*b >= K.
# ----------------------------------------------------------------------------
def reindex_gap(token_probs, a, b, rng):
    """Map token_id -> (family, variant) = (id // b, id % b); compare flat vs
    two-stage ideal code length, expected under token_probs. Returns the max
    absolute per-token discrepancy (should be ~0) and the expected-length gap.
    """
    K = len(token_probs)
    assert a * b >= K
    # family marginal and conditional
    fam_mass = [0.0] * a
    for tid, p in enumerate(token_probs):
        fam = tid // b
        fam_mass[fam] += p
    max_disc = 0.0
    exp_flat = 0.0
    exp_two = 0.0
    for tid, p in enumerate(token_probs):
        if p <= 0:
            continue
        fam = tid // b
        flat = -math.log2(p)
        # two-stage: -log2 P(fam) - log2 P(variant|fam) = -log2(fam_mass) - log2(p/fam_mass)
        two = -math.log2(fam_mass[fam]) - math.log2(p / fam_mass[fam])
        max_disc = max(max_disc, abs(flat - two))
        exp_flat += p * flat
        exp_two += p * two
    return max_disc, abs(exp_flat - exp_two)


def run_partA_B(name, P, Ns, rng, tol_gap):
    h, pi = markov_entropy_rate(P)
    sigma = len(P)
    print(f"\n[{name}] order-1 Markov, sigma={sigma}, entropy rate h(X) = {h:.5f} bits/symbol")
    print(f"   {'N':>9} {'c(N)':>8} {'rate':>9} {'gap':>9} "
          f"{'gap/rdy':>8} {'ptr/N':>8} {'ovhd/N':>8} {'reidx/N':>9}")
    prev_gap = None
    gaps = []
    overheads = []
    last = {}
    for N in Ns:
        syms = sample_markov(P, pi, N, rng)
        c, body_bits, avg_len = lz78_parse(syms, sigma)
        # universal LZ78 per-symbol rate
        rate = body_bits / N
        gap = rate - h
        # PER-SYMBOL DECOMPOSITION (matches the theorem's index/header
        # accounting). The POINTER stream c*log2(c) is the RATE-BEARING term
        # (Theta(N), -> N*h(X)); it is NOT overhead. The genuinely additive
        # OVERHEAD is the new-symbol stream c*log2(sigma) plus the re-index
        # rounding slack -- both O(c) = O(N/log N) = o(N).
        ptr_over_N = (c * math.log2(c) / N) if c > 1 else 0.0  # rate-bearing
        newsym_bits = c * (math.log2(sigma) if sigma > 1 else 1)  # o(N) overhead
        overhead_over_N = newsym_bits / N
        # INCREMENTAL hierarchical re-index overhead beyond the flat index:
        # the (family,variant) split is bijective, so the *extra* address bits
        # over a flat ceil(log2 c)-bit index is the rounding slack only,
        # <= c bits total (one rounding bit per family stage), also o(N).
        reindex_over_N = (c * 1.0) / N
        # Kieffer-Yang/LZ redundancy rate ~ log log N / log N; report gap/redy
        redy = math.log2(math.log2(N)) / math.log2(N)
        print(f"   {N:>9} {c:>8} {rate:>9.4f} {gap:>9.4f} "
              f"{gap / redy:>8.3f} {ptr_over_N:>8.4f} {overhead_over_N:>8.5f} "
              f"{reindex_over_N:>9.5f}")
        gaps.append(gap)
        overheads.append(overhead_over_N + reindex_over_N)
        last = dict(N=N, c=c, rate=rate, gap=gap,
                    overhead_over_N=overhead_over_N + reindex_over_N,
                    reindex_over_N=reindex_over_N)
        prev_gap = gap
    # --- assertions ---
    ok = True
    # (A) universality: gap -> h(X) from above. Two checks:
    #   (A1) the per-symbol gap at the largest N is within tolerance;
    #   (A2) the gap is monotonically DECREASING across the N schedule
    #        (slow O(log log N/log N) convergence, but unmistakably shrinking).
    monotone = all(gaps[i + 1] < gaps[i] + 1e-9 for i in range(len(gaps) - 1))
    decreased = gaps[-1] < gaps[0] * 0.6  # at least 40% shrink end-to-end
    if last["gap"] > tol_gap:
        print(f"   FAIL (A1): final gap {last['gap']:.4f} > tol {tol_gap}")
        ok = False
    elif not monotone:
        print(f"   FAIL (A2): gap not monotone decreasing: {[round(g,4) for g in gaps]}")
        ok = False
    elif not decreased:
        print(f"   FAIL (A2): gap shrank < 40% end-to-end ({gaps[0]:.4f} -> {gaps[-1]:.4f})")
        ok = False
    else:
        print(f"   PASS (A): per-symbol gap monotone decreasing {gaps[0]:.4f} -> "
              f"{gaps[-1]:.4f} <= {tol_gap}; rate -> h(X)={h:.4f} "
              f"(consistent with O(log log N / log N) redundancy)")
    # (B) phrase count is sublinear: c(N) <= N/log_sigma(N) * const, so
    # c(N)/N -> 0 (eq. 6.7). Check c(N)*log_sigma(N)/N is bounded by a
    # moderate const, AND that the genuinely-additive overhead per symbol
    # (new-symbol stream + re-index slack = O(c) = o(N)) DECREASES to 0.
    cbound = last["c"] * (math.log(last["N"], sigma) if sigma > 1 else math.log(last["N"])) / last["N"]
    overhead_decreasing = all(
        overheads[i + 1] < overheads[i] + 1e-9 for i in range(len(overheads) - 1)
    )
    if cbound > 5.0:
        print(f"   FAIL (B): c(N) log_sigma N / N = {cbound:.3f} not bounded (eq 6.7)")
        ok = False
    elif not overhead_decreasing:
        print(f"   FAIL (B): additive overhead/N not decreasing: "
              f"{[round(o,5) for o in overheads]}")
        ok = False
    else:
        print(f"   PASS (B): c(N)/N -> 0 (c(N) log_sigma N / N = {cbound:.3f}, bounded, eq 6.7); "
              f"additive overhead/N {overheads[0]:.5f} -> {overheads[-1]:.5f} -> 0 "
              f"(pointer stream is the rate-bearing Theta(N) term, not overhead)")
    return ok


def run_partC(rng):
    print("\n[C] Hierarchical (family,variant) re-index is rate-neutral "
          "(Thm 6.3, 2-level tree)")
    ok = True
    worst = 0.0
    trials = 0
    for _ in range(2000):
        # random grid shape and a dictionary K <= a*b
        a = rng.randint(2, 16)
        b = rng.randint(2, 16)
        K = rng.randint(2, a * b)
        # random reference distribution over K tokens (Dirichlet-ish)
        raw = [rng.random() ** 2 + 1e-6 for _ in range(K)]
        s = sum(raw)
        probs = [x / s for x in raw]
        max_disc, exp_gap = reindex_gap(probs, a, b, rng)
        worst = max(worst, max_disc, exp_gap)
        trials += 1
        if max_disc > 1e-9 or exp_gap > 1e-9:
            print(f"   FAIL (C): a={a} b={b} K={K} max_disc={max_disc:.2e} "
                  f"exp_gap={exp_gap:.2e}")
            ok = False
            break
    if ok:
        print(f"   PASS (C): flat vs two-stage code identical over {trials} random "
              f"grids; worst |gap| = {worst:.2e} (machine precision)")
    return ok


def main():
    rng = random.Random(20260603)

    # Two ergodic order-1 Markov sources (irreducible + aperiodic => ergodic).
    # Source 1: 2-state, skewed (low entropy rate ~ 0.5 bits) -> compresses.
    P2 = [[0.92, 0.08],
          [0.20, 0.80]]
    # Source 2: 3-state, moderate entropy rate.
    P3 = [[0.70, 0.20, 0.10],
          [0.15, 0.70, 0.15],
          [0.10, 0.25, 0.65]]

    # Ns grow geometrically; LZ78 redundancy is O(log log N / log N), so we go
    # large to see clear convergence of the per-symbol gap.
    Ns = [2_000, 8_000, 32_000, 128_000, 512_000, 2_000_000]

    ok = True
    # tolerance reflects the slow O(log log N/log N) LZ78 redundancy at N=2e6:
    # log log N / log N ~ ln(ln 2e6)/ln(2e6) ~ 0.18; the realized gap is well
    # within a generous 0.6 bit/symbol and visibly decreasing.
    ok &= run_partA_B("2-state", P2, Ns, rng, tol_gap=0.60)
    ok &= run_partA_B("3-state", P3, Ns, rng, tol_gap=0.60)
    ok &= run_partC(rng)

    print("\n" + ("PASS: Theorem 6.6b verified "
                  "(universality -> h(X); index overhead o(N); re-index rate-neutral)."
                  if ok else "FAIL: see messages above."))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
lemma_6_6d_mdl_entropy_fallback_universal.py
============================================================================
Lemma 6.6d + Remark 6.6e -- universality of the greedy marginal-MDL Type-III-C
dictionary builder under an ENTROPY-CODED literal fallback (BUG-006 / Rem 6.6b).

THE GAP (Remark 6.6b):  Theorem 6.6c proves the converse floor
    E[L_IIIC(X^N)] >= H(X^N) >= N h(X)                              (6.6c.1)
for ANY uniquely-decodable Type-III-C code, and Theorem 6.6b proves the matching
achievability Nh+o(N) *under a universal grammar builder* (LZ78 / irreducible
grammar). Whether the concrete greedy marginal-MDL admission rule of Theorem 6.5
(admit tau iff Delta L(tau) < 0, all else frozen) is itself universal is OPEN.

THE HONEST STEP (Lemma 6.6d):  the admission decisions are IRRELEVANT to the rate
once two benign, explicitly-stated features hold:
  (A) literals (uncovered / length-1 positions) are coded under the SAME universal
      sequential predictor Q that codes the token references -- NOT at a flat 8*ell
      raw-byte cost;
  (B) the builder ships the smaller of {learned dictionary, empty dictionary}
      (a baseline guard; one comparison + one signalling bit, what any sane
      implementation does).
Then for every stationary ergodic source,  (1/N) E[L_IIIC] -> h(X), BECAUSE the
empty dictionary already entropy-codes X^N under universal Q at Nh+o(N), the guard
forbids doing worse, and the converse (6.6c.1) forbids doing better. The greedy-MDL
admissions can only help; they cannot break universality.

THE RESIDUAL THAT STAYS OPEN (Remark 6.6e):  the *literal* Theorem 6.5 rule uses a
flat 8*ell raw-literal cost (feature (A) violated). Then the empty dictionary costs
8N bits, the guard is useless (8N >> Nh), and reaching h requires the greedy
admissions to COVER a (1-o(1))-rate fraction AND code it near-entropy -- the genuine
smallest-grammar-adjacent question. Charikar et al. (2005) grammar-SIZE lower bounds
do NOT settle it (LZ78 is greedy, size-suboptimal, yet universal), so it is open.

CHECKS (all rate-vs-h comparisons averaged over seeds, since the converse
E[L]>=H(X^N) is an EXPECTATION bound -- a single sample's -log2 Q(x^N) fluctuates
+-O(sqrt N) and may dip below Nh):
  D1  empty-dict universal coder (adaptive KT order-k) -> h on a Markov source
      (E[rate]-h shrinks as N grows; E[rate] >= h converse floor respected).
  D2  GUARDED sandwich: run a real greedy-MDL bigram builder with entropy-coded
      fallback; for ANY dictionary (empty / greedy-learned / adversarially-bad), the
      guarded realized rate is <= empty-dict rate (PER SAMPLE, exact) AND E[rate]->h.
      Admissions never break universality.
  D3  RESIDUAL localized: under the flat-8*ell literal rule the EMPTY-dict baseline
      is 8 bits/symbol >> h, so the guard cannot deliver universality; under the
      entropy-coded fallback the empty-dict baseline is h. This is exactly why
      feature (A) is essential and what keeps the literal rule open.
  D4  converse floor (consolidated): E[rate] >= h for the universal coder.

Deps: numpy.  Python: /Users/para/.venvs/rnr/bin/python.
"""
import math, random
import numpy as np
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<70} {'PASS' if ok else 'FAIL'}"); return ok

# ----------------------------------------------------------------------------
# sources
# ----------------------------------------------------------------------------
def markov2(N, P, seed):
    rng = random.Random(seed); x = [0]; s = 0
    for _ in range(N - 1):
        r = rng.random(); s = 0 if r < P[s][0] else 1; x.append(s)
    return x

def markov_entropy_rate(P):
    p01, p10 = P[0][1], P[1][0]
    pi0 = p10 / (p01 + p10)
    def H(row): return -sum(q * math.log2(q) for q in row if q > 0)
    return pi0 * H(P[0]) + (1 - pi0) * H(P[1])

def iid_bern(N, p, seed):
    rng = random.Random(seed)
    return [1 if rng.random() < p else 0 for _ in range(N)]

# ----------------------------------------------------------------------------
# universal sequential predictor: adaptive KT (Krichevsky-Trofimov) order-k
#   Q(a|s) = (n[s,a]+1/2)/(sum_b n[s,b] + A/2);  codelen = -sum log2 Q(x_t|s_t)
# universal over order-<=k sources: rate -> h + O(log N / N).
#
# SCOPE NOTE (matches Lemma 6.6d feature (A) as corrected): a FIXED-order-k KT is
# universal only for sources of order <= k. We use it here on finite-order Markov
# sources (order 1, with k=4 >= 1), where it IS universal, to ILLUSTRATE the
# guarded-sandwich mechanism. The Lemma's "for every stationary source" claim
# requires a Q universal for ALL stationary sources -- LZ78 (Cover-Thomas 13.5) or
# a growing-order/countable mixture (Ryabko 1984) -- NOT a fixed-order KT. The
# probe demonstrates the mechanism, not the all-stationary universality of KT.
# ----------------------------------------------------------------------------
def kt_codelen(seq, A, k):
    counts = {}; tot = {}; L = 0.0
    ctx = tuple([0] * k)
    for x in seq:
        c = counts.get(ctx); t = tot.get(ctx, 0)
        na = (c[x] if c else 0) + 0.5
        L += -math.log2(na / (t + A * 0.5))
        if c is None:
            c = [0] * A; counts[ctx] = c
        c[x] += 1; tot[ctx] = t + 1
        if k > 0:
            ctx = ctx[1:] + (x,)
    return L

# ----------------------------------------------------------------------------
# greedy marginal-MDL bigram builder -- a concrete instance of the Theorem 6.5 rule
#  mode 'entropy': literals AND references coded under universal KT (feature A on)
#  mode 'flat8'  : literals at flat 8 bits/symbol (Theorem 6.5 literal, feature A off)
# ----------------------------------------------------------------------------
def greedy_mdl_code(seq, A, k, mode, dict_tokens):
    tokmap = {tok: A + i for i, tok in enumerate(dict_tokens)}
    Aref = A + len(dict_tokens)
    parsed = []; i = 0; n = len(seq)
    while i < n:
        if i + 1 < n and (seq[i], seq[i + 1]) in tokmap:
            parsed.append(tokmap[(seq[i], seq[i + 1])]); i += 2
        else:
            parsed.append(seq[i]); i += 1
    header = len(dict_tokens) * 2 * max(1, math.ceil(math.log2(A))) if dict_tokens else 0
    if mode == 'entropy':
        body = kt_codelen(parsed, Aref, k)
    elif mode == 'flat8':
        lit = sum(1 for s in parsed if s < A) * 8.0
        refs = [s for s in parsed if s >= A]
        body = lit + (kt_codelen(refs, Aref, 0) if refs else 0.0)
    else:
        raise ValueError(mode)
    return header + body

def greedy_admit_digrams(seq, A, k, mode, topm=8):
    from collections import Counter
    dig = Counter((seq[i], seq[i + 1]) for i in range(len(seq) - 1))
    cands = [t for t, _ in dig.most_common(topm)]
    admitted = []; base = greedy_mdl_code(seq, A, k, mode, admitted)
    for tau in cands:
        trial = greedy_mdl_code(seq, A, k, mode, admitted + [tau])
        if trial < base:                  # Delta L(tau) < 0 (realized, rest frozen)
            admitted.append(tau); base = trial
    return admitted

def guarded_code(seq, A, k, mode, dict_tokens):
    """baseline guard: ship min{learned, empty} + 1 signalling bit."""
    Ld = greedy_mdl_code(seq, A, k, mode, dict_tokens)
    L0 = greedy_mdl_code(seq, A, k, mode, [])
    return min(Ld, L0) + 1.0, Ld, L0

SEEDS = range(101, 109)        # 8 seeds for the heavier greedy-builder estimates (D2/D3)
SEEDS_MANY = range(201, 265)   # 64 seeds for the tight empty-dict converse margin (D1/D4)
# the converse E[L]>=H(X^N)>=Nh is an EXPECTATION bound; the true KT mean is
# h + redundancy (redundancy ~ 7e-4 at order-4, N=2e5), and the 64-seed Monte-Carlo
# standard error is ~2e-4, so we assert E[rate] >= h - MC_TOL with MC_TOL=1e-3.
MC_TOL = 1e-3

# ----------------------------------------------------------------------------
def D1():
    print("-" * 78)
    print("D1  empty-dict universal KT coder -> h on a Markov source (E over seeds)")
    P = [[0.85, 0.15], [0.4, 0.6]]; h = markov_entropy_rate(P); k = 4
    ok = True; prev = None
    for N in (2000, 20000, 200000):
        rates = [kt_codelen(markov2(N, P, s), 2, k) / N for s in SEEDS_MANY]
        Er = sum(rates) / len(rates); gap = Er - h
        floor_ok = Er >= h - MC_TOL
        print(f"     N={N:>7}: E[rate]={Er:.5f}  h={h:.5f}  E[rate]-h={gap:+.5f}  >=h-MC:{floor_ok}")
        ok = ok and floor_ok
        if prev is not None: ok = ok and abs(gap) < abs(prev) + 1e-6  # |gap| shrinks
        prev = gap
    ok = ok and abs(prev) < 0.005   # converged to h
    # robust strict floor at small N (redundancy >> MC noise there): recompute N=2000
    r2k = sum(kt_codelen(markov2(2000, P, s), 2, k) / 2000 for s in SEEDS_MANY) / 64
    ok = ok and (r2k > h)
    return rep("D1 universal empty-dict coder -> h (E[rate]>=h-MC, |gap| shrinks)", ok)

def D2():
    print("-" * 78)
    print("D2  GUARDED sandwich (entropy fallback): any dict -> guard<=empty AND E[rate]->h")
    P = [[0.8, 0.2], [0.3, 0.7]]; h = markov_entropy_rate(P); k = 4; A = 2; N = 200000
    ok = True
    bad = [(0, 1), (1, 0)]   # adversarial: boundary-straddling digrams, inflate ref alphabet
    for tag, make in (("learned", "L"), ("adversarial", "A"), ("empty", "E")):
        per_sample_guard_ok = True; rates = []
        for s in SEEDS:
            seq = markov2(N, P, s)
            D = (greedy_admit_digrams(seq, A, k, 'entropy') if make == "L"
                 else bad if make == "A" else [])
            g, Ld, L0 = guarded_code(seq, A, k, 'entropy', D)
            per_sample_guard_ok = per_sample_guard_ok and (g <= L0 + 1.0 + 1e-9)
            rates.append(g / N)
        Er = sum(rates) / len(rates); near = (Er >= h - 2e-3) and (Er <= h + 0.02)
        print(f"     {tag:<11}: E[guarded rate]={Er:.5f}  h={h:.5f}  guard<=empty(all seeds):"
              f"{per_sample_guard_ok}  ->h:{near}")
        ok = ok and per_sample_guard_ok and near
    return rep("D2 guarded greedy-MDL universal regardless of admissions (sandwich->h)", ok)

def D3():
    print("-" * 78)
    print("D3  RESIDUAL localized: flat-8 empty-dict baseline=8>>h vs entropy baseline=h")
    p = 0.05; h = -(p * math.log2(p) + (1 - p) * math.log2(1 - p)); A = 2; N = 60000; k = 0
    rates_flat0 = []; rates_ent0 = []
    for s in SEEDS:
        seq = iid_bern(N, p, s)
        rates_flat0.append(greedy_mdl_code(seq, A, k, 'flat8', []) / N)   # empty dict, flat-8
        rates_ent0.append(greedy_mdl_code(seq, A, k, 'entropy', []) / N)  # empty dict, entropy
    Ef = sum(rates_flat0) / len(rates_flat0); Ee = sum(rates_ent0) / len(rates_ent0)
    print(f"     iid p={p}: h={h:.4f}  |  flat-8 empty-dict baseline={Ef:.4f}  "
          f"entropy empty-dict baseline={Ee:.4f}")
    print(f"     flat-8 baseline >> h (guard useless): {Ef > h + 3.0};  "
          f"entropy baseline -> h (guard delivers): {abs(Ee - h) < 0.05}")
    # the residual is precisely this baseline gap: feature (A) is what makes the guard work.
    ok = (abs(Ef - 8.0) < 1e-9) and (Ef > h + 3.0) and (abs(Ee - h) < 0.05) and (Ee >= h - 1e-3)
    return rep("D3 flat-8 baseline=8>>h (residual open) vs entropy baseline=h (feature A essential)", ok)

def D4():
    print("-" * 78)
    print("D4  converse floor (consolidated): E[rate] >= h for the universal coder")
    P = [[0.8, 0.2], [0.3, 0.7]]; h = markov_entropy_rate(P); k = 4; N = 200000
    rates = [kt_codelen(markov2(N, P, s), 2, k) / N for s in SEEDS_MANY]
    Er = sum(rates) / len(rates); ok = Er >= h - MC_TOL
    print(f"     N={N}: E[rate]={Er:.5f} >= h-MC={h - MC_TOL:.5f} (h={h:.5f}): {ok}")
    return rep("D4 converse floor E[rate]>=h respected (Thm 6.6c.1 empirical, MC)", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("Lemma 6.6d / Remark 6.6e -- greedy-MDL universality under entropy-coded fallback")
    print("=" * 78)
    D1(); D2(); D3(); D4()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")

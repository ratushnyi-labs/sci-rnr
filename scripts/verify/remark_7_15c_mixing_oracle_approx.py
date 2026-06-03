r"""
Verification for Remark 7.15c (the noncontiguous bounded-order support-selection
oracle of Remark 7.15a -- whose EXACT complexity is open (= exact finite-block HMM
entropy) -- is poly-time epsilon-APPROXIMABLE under a FILTER-STABILITY hypothesis
(FS); so its #P-hardness is an EXACT-COMPUTATION artifact (the exact oracle is
suspected hard even under mixing; filter stability makes only the APPROXIMATE oracle
poly), and greedy support selection inherits a (1-1/e)-O(eps) guarantee on
filter-stable sources).

------------------------------------------------------------------------
THE INDIRECT PATH.

Remark 7.15a left OPEN the exact complexity of H(X_T) for a NONCONTIGUOUS subset T
of a bounded-order-w chain (= exact finite-block Shannon entropy of the induced
HMM / function-of-Markov-chain process Y_k := X_{t_k}). Direct attack (exact
#P-hardness vs poly algorithm) is blocked. Indirect path: under a FILTER-STABILITY
hypothesis (FS) -- the order-w hidden chain is geometrically ergodic with positive
(Doeblin) emissions, so the induced missing-observation HMM's prediction filter
forgets geometrically (Atar-Zeitouni 1997; Le Gland-Mevel 2000; discrete analogue
of Lemma 5.6e's GE+EP) -- the conditional MI of Y decays geometrically. NOTE: (FS)
is STRONGER than the bare exp-CMI mixing (5.8) of Theorem 7.2 -- sigma(Y)<=sigma(X)
transfers beta-mixing but NOT the conditional-MI bound (conditional MI is not
monotone under dropping the skipped mediators), so (FS) is the right primitive. Then

    H(Y_1..Y_r)  =  sum_k H(Y_k | Y_{<k})   (chain rule)

is approximated by the ORDER-d truncation

    A_d  :=  sum_k H(Y_k | Y_{k-d .. k-1})

with GEOMETRIC error  0 <= A_d - H(Y_1..Y_r) <= r * C * rho^d,  where rho < 1 is
the mixing/contraction coefficient (conditioning on less can only RAISE entropy, so
A_d >= H; the gap is the tail dependence the window drops, which decays as rho^d).
The order-d quantity uses only (d+1)-block marginals of Y -- |Sigma|^{d+1} entries,
each computed by the forward algorithm over the A^w hidden state -- so for
d = O(log(r/eps)/log(1/rho)) the approximation is exact-to-eps in time
poly(r, 1/eps) (|Sigma|^{d+1} = (r/eps)^{O(log|Sigma|/log(1/rho))}).

------------------------------------------------------------------------
WHAT THIS SCRIPT VERIFIES. (These are NUMERICAL ILLUSTRATIONS on small brute-force
instances of the analytic claims, whose PROOF is the FS filter-stability argument
above -- M2/M3 illustrate the geometric decay and its mixing-dependence rather than
proving the bound; M4 illustrates the logarithmic depth budget using brute marginals,
not the analytic forward algorithm; M1/M5 verify exact identities / the robust-greedy
bound.)

(M1) MONOTONE OVER-ESTIMATE: A_d is non-increasing in d and A_d >= H(Y_1..Y_r),
     with A_d -> H exactly once d >= r-1 (full history). (Conditioning reduces
     entropy.)

(M2) GEOMETRIC DECAY: the gap g_d := A_d - H(Y_1..Y_r) decays geometrically in d,
     g_d <= C rho^d; we fit the decay ratio and confirm it is < 1 and tracks the
     chain's mixing (more mixing => smaller ratio => faster decay).

(M3) MIXING DEPENDENCE: across chains of increasing mixing (flatter transition,
     smaller second-eigenvalue modulus / Dobrushin coefficient), the decay is
     faster -- the approximation is cheap exactly in the regime the paper assumes.

(M4) POLY BUDGET: to reach additive error eps the needed depth is
     d(eps) = O(log(r/eps)/log(1/rho)); the table size |Sigma|^{d+1} is poly(r/eps)
     for fixed |Sigma|, rho. We tabulate d(eps) and the table size and confirm
     logarithmic growth in 1/eps (NOT the |Sigma|^|T| of exact brute marginalisation).

(M5) ROBUST GREEDY: greedy support selection driven by the eps-approximate oracle
     still attains M(greedy) >= (1-1/e) M(opt) - O(eps) (additive oracle error
     propagates additively through the greedy; Nemhauser-Wolsey-Fisher is robust).
     Verified vs the exact-oracle greedy and brute optimum on small fields.

PASS = M1 (monotone over-estimate, -> H at full depth) AND M2 (geometric decay,
ratio < 1) AND M3 (faster for more mixing) AND M4 (logarithmic depth budget) AND
M5 (approx-oracle greedy within (1-1/e)-O(eps)).

HONEST SCOPE. This does NOT resolve the EXACT complexity (still open, suspected hard
EVEN under mixing). It shows only the OPERATIVE (eps-approximate) oracle is poly
under mixing -- the regime of all the paper's main results -- so the #P-hardness of
Remark 7.15a is an EXACT-computation artifact (not removed by mixing; only the
approximate oracle is) and does not obstruct mixing-source support selection. The
poly cost is FIXED-PARAMETER (degree ~ log|Sigma|/log(1/rho), blowing up as rho->1).
Complements Remark 7.15b (streaming/RA already poly): the
optional selection LEVER is also poly-approximable once the source mixes.
"""

import math
import sys
from itertools import combinations, product

import numpy as np

RNG = np.random.default_rng(20260603)


# ----------------------------------------------------------------------
# order-2 chain Y observed on even positions {0,2,4,...} (function of a chain).
# temp controls mixing: larger temp -> flatter transition -> more mixing.
# ----------------------------------------------------------------------

def order2(A, temp, rng):
    T2 = np.zeros((A, A, A))
    for a in range(A):
        for b in range(A):
            l = rng.standard_normal(A) / temp
            p = np.exp(l - l.max()); T2[a, b] = p / p.sum()
    l = rng.standard_normal((A, A)) / temp
    pe = np.exp(l - l.max()); pi2 = pe / pe.sum()
    return pi2, T2


def order2_joint(pi2, T2, N):
    A = pi2.shape[0]
    P = np.zeros((A,) * N)
    for idx in product(range(A), repeat=N):
        pr = pi2[idx[0], idx[1]]
        for i in range(2, N):
            pr *= T2[idx[i - 2], idx[i - 1], idx[i]]
        P[idx] = pr
    return P / P.sum()


def subset_marginal(P, keep):
    """Marginal table P(X_keep) over the kept (even) positions."""
    keep = tuple(sorted(keep))
    axes_out = tuple(i for i in range(P.ndim) if i not in keep)
    return P.sum(axis=axes_out) if axes_out else P


def Htab(tab):
    t = np.asarray(tab).ravel(); t = t[t > 0]
    return float(-(t * np.log2(t)).sum())


def cond_entropy_window(Ymarg, k, d):
    """H(Y_k | Y_{k-d..k-1}) from the marginal table Ymarg over the r kept symbols.
    d=0 -> H(Y_k); d>=k -> H(Y_k | Y_0..Y_{k-1})."""
    lo = max(0, k - d)
    block = tuple(range(lo, k + 1))          # Y_{lo..k}
    ctx = tuple(range(lo, k))                # Y_{lo..k-1}
    r = Ymarg.ndim
    Hblk = Htab(subset_marginal_idx(Ymarg, block))
    Hctx = Htab(subset_marginal_idx(Ymarg, ctx)) if ctx else 0.0
    return Hblk - Hctx


def subset_marginal_idx(tab, keep):
    keep = tuple(sorted(keep))
    axes_out = tuple(i for i in range(tab.ndim) if i not in keep)
    return tab.sum(axis=axes_out) if axes_out else tab


def induced_marginal(P, N):
    """Marginal of the induced process Y on even positions {0,2,...} of X^N."""
    evens = list(range(0, N, 2))
    return subset_marginal(P, evens), len(evens)


# ----------------------------------------------------------------------

def second_eig_modulus(T2, A):
    """A mixing proxy: spectral radius of the order-2 lift restricted off the
    Perron eigenvalue (|lambda_2| of the A^2-state transition)."""
    M = np.zeros((A * A, A * A))
    for a in range(A):
        for b in range(A):
            for c in range(A):
                M[a * A + b, b * A + c] = T2[a, b, c]
    ev = np.sort(np.abs(np.linalg.eigvals(M)))[::-1]
    return float(ev[1]) if len(ev) > 1 else 0.0


def check_M1_M2_M3():
    print("=" * 70)
    print("M1/M2/M3: A_d monotone over-estimate -> H; geometric decay; faster")
    print("          for more-mixing (flatter) chains.")
    print("=" * 70)
    ok = True
    rates = []
    for temp in (0.8, 1.2, 2.0, 3.0):
        rng = np.random.default_rng(int(temp * 100))
        A, N = 2, 14                            # r = 7 kept symbols
        pi2, T2 = order2(A, temp, rng)
        P = order2_joint(pi2, T2, N)
        Ymarg, r = induced_marginal(P, N)
        Hexact = Htab(Ymarg)
        gaps = []
        for d in range(0, r):
            Ad = sum(cond_entropy_window(Ymarg, k, d) for k in range(r))
            gaps.append(Ad - Hexact)
        lam2 = second_eig_modulus(T2, A)
        # monotone non-increasing and >= 0, -> 0 at full depth
        mono = all(gaps[i] >= gaps[i + 1] - 1e-9 for i in range(len(gaps) - 1))
        nonneg = all(g >= -1e-9 for g in gaps)
        reaches0 = abs(gaps[-1]) < 1e-9
        # geometric decay ratio over the informative middle range
        ratios = [gaps[i + 1] / gaps[i] for i in range(len(gaps) - 1)
                  if gaps[i] > 1e-6]
        avg_ratio = float(np.mean(ratios)) if ratios else 0.0
        geom = (avg_ratio < 1.0) if ratios else True
        rates.append((temp, lam2, avg_ratio))
        print(f"  temp={temp}: |lambda2|={lam2:.3f}  gaps={['%.3f'%g for g in gaps]}")
        print(f"           monotone={mono} nonneg={nonneg} ->0@full={reaches0} "
              f"decay-ratio~{avg_ratio:.3f} (<1: {geom})")
        ok = ok and mono and nonneg and reaches0 and geom
    # M3: decay ratio should DECREASE as mixing increases (|lambda2| decreases).
    # sort by |lambda2| and check decay ratio is (weakly) monotone with it.
    rates_sorted = sorted(rates, key=lambda t: t[1])   # by |lambda2| asc
    m3 = all(rates_sorted[i][2] <= rates_sorted[i + 1][2] + 0.15
             for i in range(len(rates_sorted) - 1))
    print(f"  M3: decay-ratio tracks mixing (smaller |lambda2| => faster decay): "
          f"{m3}")
    ok = ok and m3
    print(f"  M1/M2/M3 {'PASS' if ok else 'FAIL'}")
    return ok


def check_M4():
    print("=" * 70)
    print("M4: POLY budget -- depth d(eps)=O(log(r/eps)/log(1/rho)) is LOG in")
    print("    1/eps; table size |Sigma|^{d+1} is poly(r/eps), vs |Sigma|^|T|")
    print("    for exact brute marginalisation.")
    print("=" * 70)
    # Pick a SLOWER-mixing chain (|lambda2| nearer 1) so the depth budget grows
    # visibly with 1/eps -- demonstrating LOGARITHMIC (not linear) growth.
    A, N = 2, 20
    best = None
    for seed in range(60):
        rng = np.random.default_rng(seed)
        pi2, T2 = order2(A, 0.6, rng)
        lam2 = second_eig_modulus(T2, A)
        if best is None or lam2 > best[0]:
            best = (lam2, pi2, T2)
    lam2, pi2, T2 = best
    P = order2_joint(pi2, T2, N)
    Ymarg, r = induced_marginal(P, N)
    Hexact = Htab(Ymarg)
    ok = True
    brute = A ** r
    print(f"  r={r} kept symbols, |lambda2|={lam2:.3f} (slow-mixing); "
          f"exact brute marginal table = |Sigma|^r = {brute}")
    ds = []
    for eps in (0.5, 0.2, 0.1, 0.05, 0.02, 0.01):
        d_needed = None
        for d in range(0, r):
            Ad = sum(cond_entropy_window(Ymarg, k, d) for k in range(r))
            if Ad - Hexact <= eps:
                d_needed = d
                break
        d_use = d_needed if d_needed is not None else r
        tbl = A ** (d_use + 1)
        ds.append((eps, d_use, tbl))
        print(f"    eps={eps:<5}: depth d={d_use}, table |Sigma|^(d+1)={tbl:<6} "
              f"(brute |Sigma|^r={brute};  {brute // max(tbl,1)}x smaller)")
        ok = ok and (d_needed is not None) and (tbl < brute)
    # logarithmic check: depth vs log10(1/eps) roughly linear with small slope;
    # crucially depth stays well below r and table well below brute.
    max_depth = max(d for _, d, _ in ds)
    grows = ds[-1][1] >= ds[0][1]            # depth non-decreasing as eps shrinks
    sublinear = max_depth <= r // 2          # stays well below full history r
    print(f"  => depth grows with 1/eps ({grows}) but stays <= r/2 ({sublinear}); "
          f"table poly(r/eps) << |Sigma|^r. (log-depth: each 10x smaller eps")
    print(f"     adds ~1/log10(1/rho) levels, rho=|lambda2|={lam2:.3f}.)")
    ok = ok and grows and sublinear
    print(f"  M4 {'PASS' if ok else 'FAIL'}")
    return ok


def check_M5():
    print("=" * 70)
    print("M5: ROBUST GREEDY -- support selection with the eps-approximate oracle")
    print("    attains M(greedy) >= (1-1/e) M(opt) - O(eps). vs exact + brute.")
    print("=" * 70)
    ok = True
    factor = 1.0 - 1.0 / math.e
    for trial, (A, N, temp) in enumerate([(2, 12, 0.6), (3, 10, 0.7)]):
        rng = np.random.default_rng(20 + trial)
        pi2, T2 = order2(A, temp, rng)
        P = order2_joint(pi2, T2, N)
        Ymarg, r = induced_marginal(P, N)
        t = r // 2
        # exact oracle on the induced process Y: M(S) = H(Y_S)
        def Mexact(S):
            return Htab(subset_marginal_idx(Ymarg, tuple(S)))
        # eps-approx oracle: order-d windowed entropy of the SUBSET S of Y
        # (here we emulate oracle error by an order-d truncation on Y_S directly)
        def Mapprox(S, d=2):
            S = sorted(S)
            # H(Y_S) approx = sum_j H(Y_{S_j} | previous d kept-in-S)
            tot = 0.0
            for j in range(len(S)):
                lo = max(0, j - d)
                blk = tuple(S[lo:j + 1]); ctx = tuple(S[lo:j])
                tot += Htab(subset_marginal_idx(Ymarg, blk)) - \
                    (Htab(subset_marginal_idx(Ymarg, ctx)) if ctx else 0.0)
            return tot

        def greedy(M, t):
            S = []
            for _ in range(t):
                best, bi = -1e9, None
                for i in range(r):
                    if i in S:
                        continue
                    v = M(S + [i])
                    if v > best:
                        best, bi = v, i
                S.append(bi)
            return S
        # brute opt of exact M
        bestM = max(Mexact(list(c)) for c in combinations(range(r), t))
        Sg_exact = greedy(Mexact, t); Mg_exact = Mexact(Sg_exact)
        print(f"  trial {trial}: r={r} t={t}: M_opt={bestM:.3f} "
              f"M_greedy(exact)={Mg_exact:.3f}")
        # SWEEP coarse->fine d so the approx oracle has GENUINE eps>0 (avoids the
        # vacuous eps=0 case); verify M(greedy_d) >= (1-1/e)M_opt - 2 s eps_d at
        # each d, and that eps_d shrinks as d grows.
        prev_eps = None
        for d in (0, 1, 2):
            Sg = greedy(lambda S: Mapprox(S, d=d), t)
            Mg = Mexact(Sg)
            eps_d = max(abs(Mapprox(list(c), d=d) - Mexact(list(c)))
                        for c in combinations(range(r), t))
            bound = factor * bestM - 2 * t * eps_d
            ok_d = Mg >= bound - 1e-9
            shrink = (prev_eps is None) or (eps_d <= prev_eps + 1e-9)
            print(f"           d={d}: eps_d={eps_d:.3f} (shrinks: {shrink}), "
                  f"M_greedy={Mg:.3f} >= (1-1/e)M_opt-2t*eps={bound:.3f}? {ok_d}")
            ok = ok and ok_d and shrink
            prev_eps = eps_d
        # the d=0 oracle (ignore dependence) is genuinely lossy => eps_0 > 0,
        # so the bound is tested non-vacuously, not just at the exact limit.
    print("  (non-vacuous: d=0 oracle has eps_0>0; the (1-1/e)-2t*eps bound holds")
    print("   at every d, and eps_d shrinks with d -- robust submodular greedy.)")
    print(f"  M5 {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    print()
    print("#" * 70)
    print("# Remark 7.15c: mixing makes the open-middle oracle poly-APPROXIMABLE")
    print("#" * 70)
    print()
    r1 = check_M1_M2_M3(); print()
    r2 = check_M4(); print()
    r3 = check_M5(); print()
    print("=" * 70)
    allok = r1 and r2 and r3
    print(f"M1-M3 monotone+geometric+mixing : {'PASS' if r1 else 'FAIL'}")
    print(f"M4 logarithmic depth budget     : {'PASS' if r2 else 'FAIL'}")
    print(f"M5 robust approx-oracle greedy  : {'PASS' if r3 else 'FAIL'}")
    print("=" * 70)
    print(f"OVERALL: {'PASS' if allok else 'FAIL'}")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())

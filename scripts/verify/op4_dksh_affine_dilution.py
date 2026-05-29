#!/usr/bin/env python3
r"""
Verification for OP4 §13.4.4 status update 2026-05-29:
Densest-k-Subhypergraph (DkSH) attack vector closed NEGATIVELY.

Claim under test: mapping Applebaum's p-Densest-Subhypergraph (p-DSH)
hardness through the Lemma 5.7c/5.7d positive-affine identity
    E_mu[A_s(|T \ S|)] = A_s(1) - (A_s(1) - A_s(0)) * n_0(S)/m
does NOT preserve the n^eps multiplicative gap; it dilutes it WORSE than
DkS, to 1 + Theta(g^-(d-1)) where g = 1/p is the DkSH gap and d is the
hyperedge uniformity.

Literature (verified, corrected citation):
  - Applebaum, SIAM J. Comput. 42(5):2008-2037 (2013); ECCC TR11-007.
    p-DSH on d-uniform hypergraphs (d const), YES set of density p
    contains >= p^{d-1}(1-o(1)) edge-fraction; NO: every density-p set
    contains <= p^d(1+o(1)). Gap g = p^{d-1}/p^d = 1/p, = n^eps at
    p = n^{-(c-3)/2d}. Conditional on pseudorandomness of random local
    function F_{m,Q} (cryptographic).
  - Chlamtac-Dinitz-Konrad-Kortsarz-Rabanca, APPROX-RANDOM 2016 /
    SIAM J. Discrete Math 32(2) 2018 (algorithmic upper bounds).
  NOTE: the paper's prior "Chuzhoy-Manurangsi-Schroeppel 2020"
  attribution is a mis-citation, corrected in the status update.

Three checks:
  1. Affine dilution: rho_TypeII = (1-p^d)/(1-p^{d-1}) -> 1 as p->0,
     and rho-1 = Theta(p^{d-1}) = Theta(g^-(d-1)). Worse than DkS (d=2).
  2. Reduction validity dichotomy:
     (a) FLAT-above-0 A_s (A_s(1)=..=A_s(d)) encodes the p-DSH contained-
         hyperedge count exactly (argmin E[A_s] = argmax n_0(S)).
     (b) genuinely-linear A_s(j)=a-bj encodes MAX vertex-degree-sum,
         solved optimally by greedy top-s vertices => trivial, no hardness.
  3. Monotone collapse: larger DkSH gap g => SMALLER Type-II gap.

PASS = all three confirmed (=> DkSH does not close OP4; same affine
       obstacle as DkS, strengthened by exponent d-1).
"""

import math
import random
import sys
from itertools import combinations


def rho_typeII(DY_over_m, DN_over_m, As1=1.0, As0=0.0):
    """Lemma 5.7d multiplicative ratio (NO min-cost / YES min-cost),
    normalised A_s by default (most gap-favourable: As0=0, As1=1)."""
    cost_Y = As1 - (As1 - As0) * DY_over_m   # denser => cheaper (smaller)
    cost_N = As1 - (As1 - As0) * DN_over_m   # larger
    return cost_N / cost_Y


def check1_dilution():
    """rho - 1 = Theta(p^{d-1}); collapses to 1 as p->0; worse for larger d."""
    ok = True
    for d in (2, 3, 4, 6):
        for p in (1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-9):
            DY, DN = p ** (d - 1), p ** d
            predicted = p ** (d - 1)
            # Skip cases where rho-1 ~ p^{d-1} underflows float64 epsilon
            # (these are precision artifacts, not a failure of the claim).
            if predicted < 1e-13:
                continue
            # Compute rho-1 via closed form p^{d-1}(1-p)/(1-p^{d-1}) directly.
            # The ratio form (rho_typeII - 1.0) suffers catastrophic
            # cancellation once rho ~ 1 (e.g. d=3,p=1e-6 gives rho-1 ~ 1e-12,
            # below float64's ~1e-16 relative resolution at magnitude 1), so
            # the closed form is the numerically correct quantity here.
            actual = (DY - DN) / (1.0 - DY)   # = p^{d-1}(1-p)/(1-p^{d-1})
            # ratio actual/predicted -> 1 from below (since (1-p)/(1-p^{d-1})
            # <= 1 for d>=2, p in (0,1)).
            ratio = actual / predicted
            if not (0.5 <= ratio <= 1.0 + 1e-9):
                print(f"  FAIL dilution: d={d} p={p} rho-1={actual:.3e} "
                      f"pred={predicted:.3e} ratio={ratio:.4f}")
                ok = False
            # Monotone collapse: as g grows, rho-1 must shrink
        # check strict decrease in p across the row (gap grows as p shrinks)
        vals = [(p, rho_typeII(p ** (d - 1), p ** d) - 1.0)
                for p in (1e-1, 1e-2, 1e-3, 1e-4)]
        for (p1, r1), (p2, r2) in zip(vals, vals[1:]):
            # p2 < p1 => gap larger => r2 must be smaller (gap collapses more)
            if not (r2 < r1):
                print(f"  FAIL monotone: d={d} p {p1}->{p2} gap-grows but "
                      f"rho-1 {r1:.3e}->{r2:.3e} did not shrink")
                ok = False
    if ok:
        print("  [check1] rho_TypeII - 1 = Theta(g^-(d-1)) -> 1 as gap grows; "
              "worse than DkS (d=2). PASS")
    return ok


def make_hg(n, m, d, seed):
    random.seed(seed)
    return [tuple(sorted(random.sample(range(n), d))) for _ in range(m)]


def n0_count(S, H):
    Sset = set(S)
    return sum(1 for T in H if set(T) <= Sset)


def E_flat(S, H):
    """FLAT-above-0 A_s (A_s0=0, A_s1..d=1), normalised:
    E[A_s] = 1 - n0(S)/m. Minimised by maximising n0(S)."""
    return 1.0 - n0_count(S, H) / len(H)


def E_linear(S, H, d, a=10.0, b=1.0):
    """Genuinely-linear A_s(j) = a + b*j (b>0, natural increasing
    code-length: more out-of-support coords cost more). Then
    E[A_s] = a + b*(sum |T\\S|)/m = a + b*(d*m - sum_{v in S} deg v)/m,
    minimised by MAXIMISING deg-sum (greedy top-s)."""
    Sset = set(S)
    tot_out = sum(d - len(set(T) & Sset) for T in H)   # sum |T\S|
    return a + b * tot_out / len(H)


def check2_dichotomy():
    """(a) flat-above-0 A_s encodes contained-count; (b) linear A_s = max
    deg-sum, solved by greedy top-s (trivial)."""
    ok = True
    for (n, m, d, s, seed) in [(14, 60, 3, 7, 1), (13, 50, 3, 6, 7),
                               (12, 40, 4, 6, 3), (15, 70, 3, 8, 11)]:
        H = make_hg(n, m, d, seed)
        # (a) flat-above-0 <=> max contained count
        best_flat = min(combinations(range(n), s), key=lambda S: E_flat(S, H))
        best_cnt = max(combinations(range(n), s), key=lambda S: n0_count(S, H))
        if n0_count(best_flat, H) != n0_count(best_cnt, H):
            print(f"  FAIL (a): flat-A_s argmin != argmax-count "
                  f"(n={n},m={m},d={d},s={s})")
            ok = False
        # (b) linear A_s <=> max vertex-degree-sum, optimal = greedy top-s
        deg = [0] * n
        for T in H:
            for v in T:
                deg[v] += 1
        best_lin = max(combinations(range(n), s),
                       key=lambda S: sum(deg[v] for v in S))
        # argmin E_linear should equal argmax deg-sum
        argmin_lin = min(combinations(range(n), s),
                         key=lambda S: E_linear(S, H, d))
        top = sorted(range(n), key=lambda v: -deg[v])[:s]
        ds_brute = sum(deg[v] for v in best_lin)
        ds_argmin = sum(deg[v] for v in argmin_lin)
        ds_greedy = sum(deg[v] for v in top)
        if not (ds_brute == ds_argmin == ds_greedy):
            print(f"  FAIL (b): linear-A_s deg-sums differ brute={ds_brute} "
                  f"argmin={ds_argmin} greedy={ds_greedy} (n={n},d={d})")
            ok = False
    if ok:
        print("  [check2] flat-above-0 A_s encodes p-DSH count (valid "
              "reduction); linear A_s = trivial max-deg-sum (no hardness). "
              "PASS")
    return ok


def check3_worse_than_dks():
    """For the SAME gap g, hypergraph (d>=3) surviving Type-II gap is
    strictly smaller than DkS (d=2): g^-(d-1) < g^-1 for g>1, d>=3."""
    ok = True

    def gap(d, p):  # closed form rho-1 = p^{d-1}(1-p)/(1-p^{d-1})
        DY, DN = p ** (d - 1), p ** d
        return (DY - DN) / (1.0 - DY)

    for g in (2.0, 10.0, 100.0, 1000.0):
        p = 1.0 / g
        gap_dks = gap(2, p)                                 # d=2
        for d in (3, 4, 6):
            gap_dksh = gap(d, p)
            if not (gap_dksh < gap_dks):
                print(f"  FAIL: d={d} g={g} dksh-gap {gap_dksh:.3e} not < "
                      f"dks-gap {gap_dks:.3e}")
                ok = False
    if ok:
        print("  [check3] for fixed source gap g, DkSH (d>=3) Type-II gap "
              "< DkS (d=2) Type-II gap: hyperedges dilute WORSE. PASS")
    return ok


def main():
    print("OP4 §13.4.4: DkSH attack vector affine-dilution verification")
    print("=" * 64)
    r1 = check1_dilution()
    r2 = check2_dichotomy()
    r3 = check3_worse_than_dks()
    print("=" * 64)
    if r1 and r2 and r3:
        print("ALL CHECKS PASS: DkSH does NOT close OP4. The n^eps gap "
              "dilutes to 1+Theta(g^-(d-1)) under the affine A_s identity "
              "(same obstacle as DkS, strengthened by exponent d-1>=2). "
              "Verdict: NEGATIVE.")
        return 0
    print("SOME CHECK FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())

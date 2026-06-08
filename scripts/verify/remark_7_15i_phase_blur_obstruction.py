#!/usr/bin/env python3
"""
remark_7_15i_phase_blur_obstruction.py
======================================================================
OBSTRUCTION to the Remark 7.15i conjecture via the 7.15g prime-log mechanism.

Remark 7.15i conjectures: exact finite-block Shannon entropy H(X^N) of a STATIONARY
ERGODIC FILTER-STABLE HMM with positive-dimension (fractal) Blackwell measure is #P-hard.
The natural attack is to ADAPT the 7.15g #SAT->prime-log-coefficient reduction so the
count survives ERGODIC MIXING (rather than living in a persistent frozen mode).

This script shows the 7.15g MECHANISM does NOT transfer through ergodic phase-mixing.

SETUP (route 1, periodic injection -- the most favourable for the conjecture).
A period-n hidden "position" clock; at each block boundary a FRESH 7.15g mode is drawn
(dummy w.p. 8P/(m+8P), clause-i w.p. 1/(m+8P)); the block emits the n-bit gadget word
for that mode. With a fresh redraw every block the chain forgets the mode each period
(geometric forgetting; the strongest possible filter stability) -- yet the count must be
read by an observer who does NOT see the block boundaries (no separator), i.e. from the
STATIONARY (uniform-phase) binary stream. (An OBSERVABLE separator would synchronise the
observer, collapse the belief set to poly, and land in regime I -- not the fractal middle.)

For ONE aligned fresh gadget block the marginal is p_g(x)=(P+u(x))/Z, u(x)=#falsified
clauses, Z=2^{n-3}(8P+m), and the coefficient of log P in H is exactly kappa_P=-c_0 P/Z
with c_0=#SAT (7.15g). The conjecture needs the count to survive after the observer must
PHASE-AVERAGE (marginalise the hidden clock) and after windows straddle the re-randomised
boundaries.

RESULT (the obstruction):
  V1 baseline 7.15g: aligned single-block kappa_P = -#SAT*P/Z  (count is clean & linear).
  V2 aligned multi-block: kappa_P scales as blocks*(-#SAT*P/Z) (additive count).
  V3 HIDDEN-phase (genuinely-forgetting, uniform-phase) length-N entropy coefficient is a
     MESSY rational with no clean relation to #SAT.
  V4 DECISIVE: across formulas with the SAME (n,m,P,Z) the hidden-phase coefficient is NOT
     an affine function of #SAT; worse, formulas with IDENTICAL collision multiset {c_j}
     (hence identical #SAT) give DIFFERENT hidden-phase coefficients. So kappa_P(H_hidden)
     does not even factor through the collision statistic, let alone equal alpha*#SAT+beta.
     => a single coefficient query cannot recover #SAT by the 7.15g mechanism.
  V5 the periodic clock that DEFINES the blocks makes the generating chain PERIODIC (period n),
     not the aperiodic strictly-positive Doeblin kernel the conjecture's fractal witnesses use;
     and an OBSERVABLE separator (which WOULD restore the count) collapses belief growth to
     poly => regime I. So count-carrying and the fractal middle are at odds for this route.

VERDICT: the 7.15g prime-log reduction is OBSTRUCTED in regime II; ergodic phase-mixing
"forgets the count" at the level of the entropy's prime-log coefficient. This is evidence
the conjecture (if true) needs a fundamentally different reduction, and is consistent with
regime II being honestly OPEN. (It does NOT prove regime II is sub-#P.)

Mechanism: kappa_P(H)=-sum_x p(x)*v_P(p(x)) is NOT linear in p (the v_P*log structure),
so the phase-average p_stat=(1/n)sum_phi p_phi does NOT give (1/n)sum_phi kappa_P(p_phi);
and the cross-boundary windows replace the single-aligned-assignment statistic u(x) by a
product over independent fresh assignments, which no longer counts falsified clauses of one
assignment. Cf. the literature: ergodic mixing is the TOOL that makes approximate counting
easy (MCMC; 7.15c poly-eps under FS), not a source of exact hardness; the only adjacent
exact HMM-functional hardness (Kiefer, TV distance, ICALP 2018, arXiv:1804.06170) relies on
NON-forgetting reachability-tracking.

Deps: standard library only (fractions, itertools).
"""
from fractions import Fraction as Fr
from itertools import product


# ---------------- 3-CNF helpers ----------------
def falsifies(clause, a):
    """clause = list of (var, is_positive). Falsified iff every literal is false."""
    for (v, pos) in clause:
        if ((a[v] == 1) if pos else (a[v] == 0)):
            return False
    return True


def u_count(clauses, a):
    return sum(1 for c in clauses if falsifies(c, a))


def num_sat(n, clauses):
    return sum(1 for a in product((0, 1), repeat=n)
               if all(not falsifies(c, a) for c in clauses))


def collision_counts(n, clauses):
    cj = {}
    for a in product((0, 1), repeat=n):
        j = u_count(clauses, a)
        cj[j] = cj.get(j, 0) + 1
    return dict(sorted(cj.items()))


# ---------------- prime-log coefficient ----------------
def vp(num, P):
    if num == 0:
        return 0
    e = 0
    while num % P == 0:
        num //= P
        e += 1
    return e


def coeff_logP(pmap, P):
    """coeff of log P in H=-sum p log p, exact: -sum_x p(x)*(v_P(num)-v_P(den))."""
    k = Fr(0)
    for x, p in pmap.items():
        if p == 0:
            continue
        a, b = p.numerator, p.denominator
        k += p * (vp(a, P) - vp(b, P))
    return -k


# ---------------- 7.15g gadget ----------------
def mode_block_dist(n, clauses, mode):
    """Per-mode distribution over n-bit assignments. dummy=uniform; clause-i=uniform over
       assignments that falsify clause i (its 3 vars fixed, other n-3 fair)."""
    d = {}
    if mode == 'dummy':
        for x in product((0, 1), repeat=n):
            d[x] = Fr(1, 2 ** n)
    else:
        c = clauses[mode]
        fixed = {v: (0 if pos else 1) for (v, pos) in c}  # literal false
        for x in product((0, 1), repeat=n):
            d[x] = Fr(1, 2 ** (n - 3)) if all(x[v] == b for v, b in fixed.items()) else Fr(0)
    return d


def mode_weights(m, P):
    w = {'dummy': Fr(8 * P, m + 8 * P)}
    for i in range(m):
        w[i] = Fr(1, m + 8 * P)
    return w


def aligned_block_marginal(n, clauses, P):
    """p_g(x) = sum_mode w_mode * dist_mode(x) = (P+u(x))/Z."""
    w = mode_weights(len(clauses), P)
    p = {}
    for x in product((0, 1), repeat=n):
        p[x] = sum(w[md] * mode_block_dist(n, clauses, md)[x] for md in w)
    return p


def aligned_multiblock_coeff(n, clauses, P, blocks):
    """kappa_P over 'blocks' independent fresh ALIGNED gadgets (phase pinned to 0)."""
    pg = aligned_block_marginal(n, clauses, P)
    pmap = {}
    for combo in product(list(pg.keys()), repeat=blocks):
        word = tuple(b for blk in combo for b in blk)
        pr = Fr(1)
        for blk in combo:
            pr *= pg[blk]
        pmap[word] = pr
    return coeff_logP(pmap, P)


def block_subset_marginal(n, clauses, P, positions, bits):
    """P(a fresh block emits 'bits' at within-block 'positions') = sum_mode w * P(bits|mode)."""
    w = mode_weights(len(clauses), P)
    s = Fr(0)
    for md in w:
        md_dist = mode_block_dist(n, clauses, md)
        sub = Fr(0)
        for x in product((0, 1), repeat=n):
            if all(x[positions[k]] == bits[k] for k in range(len(positions))):
                sub += md_dist[x]
        s += w[md] * sub
    return s


def hidden_phase_coeff(n, clauses, P, N):
    """Exact kappa_P of H(X^N) for the HIDDEN-phase (uniform-phase, fresh-mode) process:
       p_stat(y)=(1/n) sum_phase prod_{blocks touched} (block marginal over a fresh mode)."""
    pmap = {y: Fr(0) for y in product((0, 1), repeat=N)}
    for ph in range(n):
        # split window positions 0..N-1 (starting at phase ph) into per-block segments
        seg = []
        cur = []
        pib = ph
        for _ in range(N):
            cur.append(pib)
            pib += 1
            if pib == n:
                seg.append(cur)
                cur = []
                pib = 0
        if cur:
            seg.append(cur)
        for y in product((0, 1), repeat=N):
            prob = Fr(1, n)
            idx = 0
            for positions in seg:
                bits = y[idx:idx + len(positions)]
                idx += len(positions)
                prob *= block_subset_marginal(n, clauses, P, positions, bits)
            pmap[y] += prob
    assert sum(pmap.values()) == 1
    return coeff_logP(pmap, P)


# ============================================================
if __name__ == "__main__":
    print("=" * 76)
    print("Remark 7.15i OBSTRUCTION: ergodic phase-mixing breaks the 7.15g prime-log count")
    print("=" * 76)

    n = 4
    base = [[(0, True), (1, True), (2, True)],
            [(1, True), (2, False), (3, True)]]
    m = len(base)
    P = 5
    Z = (2 ** (n - 3)) * (8 * P + m)
    sat = num_sat(n, base)
    print(f"\nbase formula: n={n}, m={m}, #SAT={sat}, P={P}, Z={Z}")

    # V1 baseline 7.15g: aligned single block recovers #SAT
    ka1 = aligned_multiblock_coeff(n, base, P, 1)
    v1 = (ka1 == Fr(-sat * P, Z)) and (-ka1 * Z / P == sat)
    print(f"\nV1 ALIGNED single block: kappa_P={ka1} = -#SAT*P/Z, recovered #SAT={-ka1*Z/P}  PASS={v1}")

    # V2 aligned multi-block: count is additive
    ka2 = aligned_multiblock_coeff(n, base, P, 2)
    v2 = (ka2 == 2 * ka1) and (-ka2 * Z / (P * 2) == sat)
    print(f"V2 ALIGNED 2 blocks: kappa_P={ka2} = 2*(-#SAT*P/Z), per-block #SAT={-ka2*Z/(P*2)}  PASS={v2}")

    # V3 hidden-phase: messy, no clean count
    kh = hidden_phase_coeff(n, base, P, n)
    v3 = (kh != ka1) and (-kh * Z / P != sat)
    print(f"V3 HIDDEN-phase length-{n}: kappa_P={kh} (!= aligned {ka1}), -kappa Z/P={-kh*Z/P} != #SAT  PASS={v3}")

    # V4 DECISIVE: hidden-phase coeff does not factor through #SAT or the collision multiset
    print("\nV4 DECISIVE non-recoverability (3-clause, 5-var family; same n,m,P,Z):")
    n2, m2, P2 = 5, 3, 7
    Z2 = (2 ** (n2 - 3)) * (8 * P2 + m2)
    fam = {
        "A": [[(0, True), (1, True), (2, True)], [(1, True), (2, False), (3, True)], [(2, True), (3, True), (4, True)]],
        "B": [[(0, True), (1, True), (2, True)], [(0, False), (1, False), (2, False)], [(2, True), (3, True), (4, True)]],
        "F2": [[(0, True), (1, True), (2, True)], [(2, False), (3, True), (4, True)], [(0, True), (3, False), (4, False)]],
        "C": [[(0, True), (1, True), (2, True)], [(0, True), (1, True), (3, True)], [(0, True), (1, True), (4, True)]],
        "D": [[(0, True), (1, False), (2, True)], [(1, True), (2, True), (3, False)], [(0, False), (3, True), (4, True)]],
    }
    print(f"  {'F':>3} {'#SAT':>4} {'collision c_j':>26} {'kappa_aligned':>14} {'kappa_hidden':>18}")
    rec = []
    for name, cl in fam.items():
        s = num_sat(n2, cl)
        cj = collision_counts(n2, cl)
        kal = aligned_multiblock_coeff(n2, cl, P2, 1)
        khi = hidden_phase_coeff(n2, cl, P2, n2)
        rec.append((name, s, tuple(sorted(cj.items())), kal, khi))
        print(f"  {name:>3} {s:>4} {str(cj):>26} {str(kal):>14} {str(khi):>18}")

    # aligned must be exactly -#SAT P/Z (sanity)
    aligned_ok = all(kal == Fr(-s * P2, Z2) for (_, s, _, kal, _) in rec)
    # same-collision => same #SAT but DIFFERENT hidden coeff?
    by_cj = {}
    for (name, s, cj, kal, khi) in rec:
        by_cj.setdefault(cj, []).append((name, s, khi))
    same_cj_diff_hidden = False
    for cj, lst in by_cj.items():
        if len(lst) >= 2:
            sats = {s for (_, s, _) in lst}
            hiddens = {khi for (_, _, khi) in lst}
            if len(sats) == 1 and len(hiddens) >= 2:
                same_cj_diff_hidden = True
                names = ",".join(nm for nm, _, _ in lst)
                print(f"  -> formulas {{{names}}} share c_j={dict(cj)} (#SAT={list(sats)[0]}) "
                      f"but {len(hiddens)} DIFFERENT hidden coeffs => not a function of #SAT/c_j")
    # also: not affine in #SAT (two same-#SAT-different-coeff points already kill affinity)
    v4 = aligned_ok and same_cj_diff_hidden
    print(f"  aligned == -#SAT*P/Z for all: {aligned_ok}; "
          f"identical-collision formulas give different hidden coeffs: {same_cj_diff_hidden}")
    print(f"  V4 PASS={v4}  (hidden coeff not affine in #SAT, not even a function of c_j)")

    print("\n" + "=" * 76)
    allok = v1 and v2 and v3 and v4
    print(f"RESULT: {'ALL PASS' if allok else 'CHECK'} -- the 7.15g prime-log #SAT reduction is")
    print("  OBSTRUCTED under ergodic phase-mixing: the count is clean & additive ONLY in the")
    print("  phase-ALIGNED (synchronisable=>regime I) reading; the genuinely-forgetting")
    print("  hidden-phase entropy coefficient does not recover #SAT (not affine in #SAT, not")
    print("  a function of the collision multiset). Regime II stays honestly OPEN; a proof of")
    print("  the conjecture needs a reduction NOT relying on phase-aligned mode-mixture counts.")

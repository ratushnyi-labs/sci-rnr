#!/usr/bin/env python3
r"""
lemma_7_34_realness_pi_w_closure.py
============================================================================
BUG-009-D / Route B (A=3): REALNESS OF THE REPLICA SPECTRUM AT s = pi,
proven on the whole open superdomain -- fully self-contained.

THEOREM (this script certifies every step, exact rational arithmetic):
  For A = 3, all 0 < p < 2/3 and 0 < D <= Dbar(p), the 5x5
  pattern-quotient replica matrix Q(p, D, w=-1) has 5 REAL eigenvalues
  (Dbar(p) = (2/3)(1-sg)^2/(1+sg^2), sg = sqrt(1-(3/2)p), the proven
  upper bound of the Gray threshold, lemma_7_34_dbar_upper_bound.py --
  so this covers the whole Gray region at the binding angle s = pi).

PROOF STRUCTURE (all derived in-script from the replica matrix itself):
  (a) SPLIT: with c = p(p-1)(2D^2-3D+2)^2, M = cQ is a polynomial
      matrix; chiM = charpoly(M) splits exactly as f1 * f4 with f1
      linear in X (explicit eigenvalue) and f4 quartic.  Q's remaining
      spectrum is real iff f4's four roots are real (f1's is rational).
  (b) DISC: disc_X(f4), recomputed here, equals
        (2D^2-3D+2)^2 (D-1)^8 D^8 (3p-4)^2 (3p-2)^12 * W(p, D)
      by exact division (W = 572 terms, deg (20, 28)).  Every cofactor
      is an even power of a polynomial that does not vanish on the open
      strip (2D^2-3D+2 has negative discriminant; D=0, p=2/3, p=4/3 are
      outside), so  sign(disc f4) = sign(W)  on the whole open strip.
  (c) W-CLOSURE (the analytic core, regions + blow-ups below):
        W(p, D) > 0   for all  0 < p < 2/3,  0 < D <= Dbar(p).
  (d) ANCHOR + CONTINUITY: at (p,D) = (1/2, 1/10), f4 has 4 distinct
      real roots (exact Sturm).  disc f4 > 0 on the connected open strip
      by (b)+(c), so no root collision occurs anywhere on it; since
      complex roots of a real polynomial can only appear through a
      collision (disc = 0) along a path, the root count 4-real is
      constant on the strip.  With f1's rational eigenvalue: 5 real.

METHOD (everything exact rational; Bernstein = closed-box strict
positivity: min coefficient > 0 AND zero exact-zero coefficients)

  (1) Rationalize the domain:  p = (2/3)(1 - sg^2),  D = th * Dbar(sg),
      Dbar(sg) = (2/3)(1-sg)^2/(1+sg^2);  (sg, th) in (0,1) x (0,1]
      is a bijective parametrization of the open strip.  Exactly:
          W(p(sg), D(sg,th)) = (1-sg)^4 * core(sg,th)
                               / (3^28 * (1+sg^2)^26),
      core in Z[sg,th], 2511 terms, deg (88, 28).  All cofactors are > 0
      on sg in (0,1), so W > 0 on the strip  <=>  core > 0 on
      (0,1) x (0,1].  core(0,th) = c0(th) vanishes in [0,1] exactly at
      th = 3/4 and th = 1 (both order 4) -- excluded boundary touches
      (sg = 0 is p = 2/3, not in the open strip).

  (2) Decompose [0,1] x [0,1] (sg0 = 1/64, window half-width h = 1/16):
        R1  = [1/64, 1] x [0, 1]          Bernstein, subdivision
        R2a = [0, 1/64] x [0, 11/16]      Bernstein, root box
        R2b = [0, 1/64] x [13/16, 15/16]  Bernstein, root box
        K1  = [0, 1/64] x [11/16, 13/16]  blow-up at (sg,th) = (0, 3/4)
        K2  = [0, 1/64] x [15/16, 1]      blow-up at (sg,th) = (0, 1)
      R1 u R2a u R2b covers everything outside the corner windows;
      K1/K2 cover the windows minus the two touch points themselves.

  (3) Corner blow-up (the key step).  For th0 in {3/4, 1} write
      C(sg,mu) = core(sg, th0 + mu) = sum c_ij sg^i mu^j.  VERIFIED
      exactly: every monomial has total degree i + j >= 4, and the
      degree-4 form F4 is positive definite:
        corner 1:  F4 = 3^28 [ (1/2) sg^4 + (1/81) (sg - (2/3) mu)^4 ]
        corner 2:  F4 = 2^34 [ (mu - 2 sg)^4 + 128 sg^4 ]
      (definiteness is not assumed anywhere -- it is what makes the
      chart certificates below succeed on their closed boxes).
      Define the order-4 projective blow-up charts (POLYNOMIAL identities,
      monomial bijections (i,j) -> (i+j-4, j) resp. (i+j-4, i)):
        HA (sg, m) := sum c_ij sg^(i+j-4) m^j
                       ==>  C(sg, sg*m)   = sg^4 * HA(sg, m)
        HB+-(nu, s) := sum c_ij (+-1)^j s^i nu^(i+j-4)
                       ==>  C(nu*s, +-nu) = nu^4 * HB+-(nu, s)
      SOUNDNESS: if HA > 0 on the CLOSED box [0,sg0] x [-1,1] and
      HB+- > 0 on the CLOSED box [0,h] x [0,1], then for every
      (sg, mu) != (0,0) with 0 <= sg <= sg0, |mu| <= h:
        |mu| <= sg  (then sg > 0):  core = sg^4 HA(sg, mu/sg)      > 0,
        |mu| >  sg:                 core = |mu|^4 HB_sgn(|mu|, sg/|mu|) > 0.
      Hence core > 0 on the whole window except the touch point (0, th0),
      which is outside the open strip.  Corner 2 needs only mu <= 0
      (th <= 1): charts HAm (m in [-1,0], stored as m' = -m in [0,1])
      and HB-.  All chart edges at sg = 0 / nu = 0 equal restrictions of
      the positive definite F4 -- this is why closed-strict certification
      terminates (and it does so at the ROOT box for every chart).

  (4) Every Bernstein certificate uses the FIXED tensor engine
      (bernstein_tensor.py): a box is certified only if min > 0 AND
      nzeros == 0 (exact-zero coefficients are counted; a zero vertex
      coefficient equals a zero of the polynomial at that corner).

CHECKS (each prints PASS/FAIL; OVERALL at the end):
  E0   engine self-test on knowable examples: exact Bernstein coefficients
       of x^2 on the lo != 0 box [1/2,1] (= {1/4, 1/2, 1}: pins the
       repaired h**b affine map); exact-zero accounting for f = x on [0,1]
       (min > 0 but nzeros = 1 => NOT closed-strict, f(0) = 0); a strict
       certificate and a False witness (f = xy - 1/2)
  DRV  full in-script derivation from the replica matrix: build Q(p,D,-1)
       (pattern-quotient, eta real at w = -1); M = cQ polynomial;
       chiM = charpoly(M); split chiM == f1 * f4 by exact division with
       the explicit linear factor; disc_X(f4) via sympy discriminant
       (1398 terms, deg (34,48)); W := disc / (explicit even cofactors)
       by exact division with zero remainder; W = 572 terms, deg (20,28);
       cofactor nonvanishing on the strip (2D^2-3D+2 negative disc).
  ANCH exact Sturm anchor: f4 at (p,D) = (1/2, 1/10) has 4 distinct real
       roots; disc f4 != 0 there.
  W1   reduction identity: cleared denominator == 3^28 (1+sg^2)^26,
       numerator == (1-sg)^4 * core (exact reconstruction), core != 0
       at sg=1 (no further (1-sg) factor); parametrization sanity
       p(0)=2/3, p(1)=0, p strictly decreasing; Dbar(0)=2/3, Dbar(1)=0
  W2   c0 = core(0,.) factorization: 1048576 (th-1)^4 (4th-3)^4 Q4^2 Q12
       with Q4, Q12 having 0 roots in [0,1] (Sturm) => c0 zeros in [0,1]
       are exactly {3/4, 1}; also W(p,0) sextic has 0 real roots
  R1   certify core > 0 closed-strict on [1/64,1] x [0,1] (subdivision)
  R2   certify core > 0 closed-strict on [0,1/64] x [0,11/16] and
       [0,1/64] x [13/16,15/16]
  K1   corner (0,3/4): min total degree == 4; the degree-4 form equals
       3^28 [(1/2) sg^4 + (1/81)(sg - (2/3) mu)^4] exactly (positive
       definite); chart construction collision-free; EXACT POLYNOMIAL
       blow-up identities (expand(C(sg, sg*m) - sg^4 HA) == 0 etc., not
       just spot checks); certify HA on [0,1/64] x [-1,1],
       HB+ and HB- on [0,1/16] x [0,1]
  K2   corner (0,1):   same with degree-4 form == 2^34 [(mu - 2 sg)^4
       + 128 sg^4] exactly; charts HAm on [0,1/64] x [0,1] (m' = -m)
       and HB- on [0,1/16] x [0,1]
  COV  coverage arithmetic: 11/16 = 3/4 - 1/16, 13/16 = 3/4 + 1/16,
       15/16 = 1 - 1/16; shared sg-boundary 1/64; union = (0,1) x [0,1]
  NUM  independent sanity (not part of the proof): exact rational
       evaluation of the ORIGINAL W at 200 random strip points > 0

Runtime: ~4 min single-threaded (R1 needs only 9 subdivision nodes; every
other Bernstein certificate passes at its ROOT box).
Deps: sympy only.  NO external inputs -- everything is derived in-script
from the replica-matrix construction (same recipe as
lemma_7_34_midband_chi_certificates.py).
Python: /Users/para/.venvs/rnr/bin/python.
"""
import itertools
import random
import sys
import time
from math import comb

import sympy as sp

T0 = time.time()
PASS = True


def rep(name, ok):
    global PASS
    PASS = PASS and bool(ok)
    print(f"  [{time.time()-T0:7.1f}s] {name:<62} {'PASS' if ok else 'FAIL'}",
          flush=True)
    return bool(ok)


def log(msg):
    print(f"  [{time.time()-T0:7.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Exact tensor Bernstein engine (closed-box strict semantics: a box is
# certified only if min > 0 AND nzeros == 0 -- exact-zero coefficients are
# dropped from the sparse dict and recovered from the dense count; a zero
# VERTEX coefficient equals a zero of the polynomial at that box corner).
# ---------------------------------------------------------------------------
def to_tensor(polyE, vs):
    P = sp.Poly(polyE, *vs)
    degs = [P.degree(v) for v in vs]
    return {mon: c for mon, c in zip(P.monoms(), P.coeffs())}, degs


def axis_apply(T, M, axis, newdim):
    out = {}
    for mon, c in T.items():
        a = mon[axis]
        for b in range(newdim):
            f = M.get((b, a))
            if f:
                nm = list(mon); nm[axis] = b; nm = tuple(nm)
                out[nm] = out.get(nm, sp.Rational(0)) + f * c
    return {k: v for k, v in out.items() if v != 0}


def box_tensor(T, degs, box):
    for ax, (lo, hi) in enumerate(box):
        d = degs[ax]
        h = hi - lo
        # v -> lo + h v:  v^a = sum_b C(a,b) lo^(a-b) h^b v^b  (h power is b)
        M = {}
        for a in range(d + 1):
            for b in range(a + 1):
                M[(b, a)] = sp.binomial(a, b) * lo ** (a - b) * h ** b
        T = axis_apply(T, M, ax, d + 1)
    return T


def bern_coeff_min(T, degs):
    total = 1
    for d in degs:
        total *= (d + 1)
    for ax, d in enumerate(degs):
        M = {}
        for b in range(d + 1):
            for a in range(b + 1):
                M[(b, a)] = sp.Rational(comb(b, a), comb(d, a))
        T = axis_apply(T, M, ax, d + 1)
    nzeros = total - len(T)
    if not T:
        return sp.Rational(0), sp.Rational(0), nzeros
    vals = list(T.values())
    return min(vals), max(vals), nzeros


def certify_box(polyE, vs, box, maxdepth=24, nodecap=4000):
    T0_, degs = to_tensor(polyE, vs)
    n = len(vs)
    work = [(list(box), 0)]
    nodes = 0
    stuck = []
    while work:
        bx, dep = work.pop()
        nodes += 1
        if nodes > nodecap:
            return None, nodes, stuck, 'nodecap'
        Tb = box_tensor(dict(T0_), degs, bx)
        mn, mx, nz = bern_coeff_min(Tb, degs)
        if mn > 0 and nz == 0:
            continue
        center = {v: sp.Rational(lo + hi, 2) if isinstance(lo + hi, int)
                  else (lo + hi) / 2 for v, (lo, hi) in zip(vs, bx)}
        cv = polyE.subs(center)
        if cv <= 0:
            return False, nodes, [(bx, cv)], 'witness'
        if dep >= maxdepth:
            stuck.append((bx, mn))
            if len(stuck) >= 6:
                return None, nodes, stuck, 'stuck'
            continue
        i = dep % n
        lo, hi = bx[i]
        mid = sp.Rational(lo + hi, 2) if isinstance(lo + hi, int) else (lo + hi) / 2
        b1 = list(bx); b1[i] = (lo, mid)
        b2 = list(bx); b2[i] = (mid, hi)
        work.append((b1, dep + 1))
        work.append((b2, dep + 1))
    return (True if not stuck else None), nodes, stuck, 'done'


p, D = sp.symbols('p D')
sg, th, mu = sp.symbols('sigma theta mu')
m, nu, s2 = sp.symbols('m nu s2')

SG0 = sp.Rational(1, 64)     # corner sg-thickness = R1 left edge
H = sp.Rational(1, 16)       # corner window half-width in theta
W1LO, W1HI = sp.Rational(11, 16), sp.Rational(13, 16)
W2LO = sp.Rational(15, 16)

PSUB = sp.Rational(2, 3) * (1 - sg**2)
DBAR = sp.Rational(2, 3) * (1 - sg)**2 / (1 + sg**2)
DSUB = th * DBAR


def engine_selftest():
    """Pin the repaired engine semantics on examples with known exact answers."""
    x, y = sp.symbols('x y')
    # (a) lo != 0 affine map: Bernstein coeffs of x^2 on [1/2, 1] are exactly
    #     {b0, b1, b2} = {lo^2, lo*hi, hi^2} = {1/4, 1/2, 1}
    T, dg = to_tensor(sp.expand(x**2), (x,))
    Tb = box_tensor(dict(T), dg, [(sp.Rational(1, 2), sp.Integer(1))])
    coef = dict(Tb)
    for ax, d in enumerate(dg):
        M = {}
        for b in range(d + 1):
            for a in range(b + 1):
                M[(b, a)] = sp.Rational(comb(b, a), comb(d, a))
        coef = axis_apply(coef, M, ax, d + 1)
    ok_a = (coef.get((0,)) == sp.Rational(1, 4)
            and coef.get((1,)) == sp.Rational(1, 2)
            and coef.get((2,)) == 1)
    rep("E0a Bernstein coeffs of x^2 on [1/2,1] == {1/4, 1/2, 1} exactly", ok_a)
    # (b) exact-zero accounting: f = x on [0,1] has b0 = f(0) = 0 -> nzeros=1,
    #     min of stored coefficients > 0: NOT closed-strict (f(0) = 0)
    T, dg = to_tensor(x, (x,))
    mn, mx, nz = bern_coeff_min(box_tensor(dict(T), dg,
                                           [(sp.Integer(0), sp.Integer(1))]), dg)
    rep("E0b f = x on [0,1]: min > 0 but nzeros == 1 (vertex zero caught)",
        mn > 0 and nz == 1)
    # (c) strict certificate and a False witness
    res1, _, _, _ = certify_box((x - 1)**2 + sp.Rational(1, 9), (x,),
                                [(sp.Integer(0), sp.Integer(1))])
    res2, _, wit, why2 = certify_box(x * y - sp.Rational(1, 2), (x, y),
                                     [(sp.Integer(0), sp.Integer(1)),
                                      (sp.Integer(0), sp.Integer(1))])
    rep("E0c certify_box: (x-1)^2 + 1/9 True; x*y - 1/2 False witness",
        res1 is True and res2 is False and why2 == 'witness')


def main():
    print("=" * 78)
    print("W-closure: W(p,D) > 0 on {0 < p < 2/3, 0 < D <= Dbar(p)}  (A=3, w=-1)")
    print("exact Bernstein + order-4 corner blow-ups, rational arithmetic only")
    print("=" * 78)

    # ---------------- E0: engine self-test ----------------
    engine_selftest()

    # ---------------- DRV: derive everything from the replica matrix -------
    X = sp.Symbol('X')
    log("DRV building Q(p, D, w=-1) (pattern-quotient replica matrix) ...")
    Aval = 3

    def _pat(tr):
        x0, u_, up = tr
        if x0 == u_ == up: return 0
        if x0 == u_ and u_ != up: return 1
        if x0 == up and u_ != up: return 2
        if u_ == up and x0 != u_: return 3
        return 4

    # eta at w = -1 is real rational: E = -D/((A-1)(1-D))
    Ew = -D / ((Aval - 1) * (1 - D))
    etaw = sp.cancel(((Aval - 1) * D * Ew + D - (Aval - 1) * Ew)
                     / ((Aval - 1) * D * Ew + D - (Aval - 1)))
    Tm = [[(1 - p) if i == j else p / (Aval - 1) for j in range(Aval)]
          for i in range(Aval)]
    states = list(itertools.product(range(Aval), repeat=3))
    reps_ = [None] * 5
    for k_, tr_ in enumerate(states):
        if reps_[_pat(tr_)] is None:
            reps_[_pat(tr_)] = k_
    Qm = sp.zeros(5, 5)
    for a_ in range(5):
        x0, xp_, xq = states[reps_[a_]]
        for (y, yp, yq) in states:
            term = Tm[xp_][yp] * Tm[xq][yq] / Tm[x0][y]
            if yp != y: term *= etaw
            if yq != y: term *= etaw
            Qm[a_, _pat((y, yp, yq))] += term
    c_ = p * (p - 1) * (2 * D**2 - 3 * D + 2)**2
    Mm = sp.zeros(5, 5)
    ok_poly = True
    for i in range(5):
        for j in range(5):
            n2, d2 = sp.fraction(sp.cancel(c_ * Qm[i, j]))
            ok_poly &= d2.is_number
            Mm[i, j] = sp.expand(n2 / d2)
    rep("DRV M = p(p-1)(2D^2-3D+2)^2 Q is a polynomial matrix", ok_poly)
    chiM = sp.expand(Mm.charpoly(X).as_expr())
    log("DRV charpoly(M) computed")
    # explicit linear factor: the split-off eigenvalue of M
    lam1M = sp.expand((3 * p - 2)**2 * D * (D - 1) * (2 * D**2 - 3 * D + 2))
    f1p = X - lam1M
    q4, r4 = sp.div(chiM, f1p, X)
    if r4 != 0:      # sign convention fallback
        f1p = X + lam1M
        q4, r4 = sp.div(chiM, f1p, X)
    f4p = sp.expand(q4)
    rep("DRV chiM == f1 * f4 exactly (explicit linear factor divides)",
        r4 == 0 and sp.expand(chiM - f1p * f4p) == 0)
    rep("DRV f4 is quartic in X with polynomial coefficients",
        sp.degree(f4p, X) == 4)
    disc = sp.discriminant(sp.Poly(f4p, X))
    discE = sp.expand(disc.as_expr() if isinstance(disc, sp.Poly) else disc)
    Pd = sp.Poly(discE, p, D)
    rep("DRV disc_X(f4) recomputed: 1398 terms, deg (34, 48)",
        len(Pd.terms()) == 1398 and Pd.degree(p) == 34 and Pd.degree(D) == 48)
    # W by exact division with the EXPLICIT even cofactors
    evens = ((2 * D**2 - 3 * D + 2)**2 * (D - 1)**8 * D**8
             * (3 * p - 4)**2 * (3 * p - 2)**12)
    qW, rW = sp.div(discE, sp.expand(evens), p, D)
    W = sp.expand(qW)
    PW = sp.Poly(W, p, D)
    rep("DRV disc == evens * W exactly (zero remainder)",
        rW == 0 and sp.expand(discE - evens * W) == 0)
    rep("DRV W: 572 terms, deg_p 20, deg_D 28",
        len(PW.terms()) == 572 and PW.degree(p) == 20 and PW.degree(D) == 28)
    # cofactors never vanish on the open strip: 2D^2-3D+2 has disc 9-16 < 0;
    # D = 0, D = 1, p = 2/3, p = 4/3 are all outside {0<p<2/3, 0<D<=Dbar<2/3}
    rep("DRV cofactor 2D^2-3D+2 has no real roots (disc = -7 < 0)",
        sp.discriminant(sp.Poly(2 * D**2 - 3 * D + 2, D)) == -7)

    # ---------------- ANCH: Sturm anchor for the continuity step ----------
    f4_anchor = sp.Poly(f4p.subs({p: sp.Rational(1, 2), D: sp.Rational(1, 10)}), X)
    n_real = f4_anchor.count_roots()
    disc_anchor = discE.subs({p: sp.Rational(1, 2), D: sp.Rational(1, 10)})
    rep("ANCH f4 at (p,D)=(1/2,1/10): 4 distinct real roots, disc != 0",
        n_real == 4 and disc_anchor != 0)

    # ---------------- W1: reduction to core ----------------
    n_, d_ = sp.fraction(sp.cancel(sp.together(W.subs({p: PSUB, D: DSUB}))))
    ok_den = sp.expand(d_ - 3**28 * (1 + sg**2)**26) == 0
    rep("W1 cleared denominator == 3^28 (1+sigma^2)^26 (positive)", ok_den)

    Pn = sp.Poly(n_, sg, th)
    f1s = sp.Poly(1 - sg, sg, th)
    corePoly, mult = Pn, 0
    while mult < 4:
        q_, r_ = sp.div(corePoly, f1s, sg, th)
        if r_ == 0 and not q_.is_zero:
            corePoly, mult = sp.Poly(q_, sg, th), mult + 1
        else:
            break
    core = corePoly.as_expr()
    ok_strip = (mult == 4
                and sp.expand(n_ - (1 - sg)**4 * core) == 0
                and sp.expand(core.subs(sg, 1)) != 0)
    rep("W1 numerator == (1-sigma)^4 * core exactly; core(1,.) != 0", ok_strip)
    log(f"core: {len(corePoly.terms())} terms, deg_sg {corePoly.degree(sg)}, "
        f"deg_th {corePoly.degree(th)}")
    rep("W1 core size: 2511 terms, deg (88, 28)",
        len(corePoly.terms()) == 2511 and corePoly.degree(sg) == 88
        and corePoly.degree(th) == 28)
    ok_par = (PSUB.subs(sg, 0) == sp.Rational(2, 3) and PSUB.subs(sg, 1) == 0
              and sp.diff(PSUB, sg) == -sp.Rational(4, 3) * sg
              and DBAR.subs(sg, 0) == sp.Rational(2, 3)
              and DBAR.subs(sg, 1) == 0)
    rep("W1 parametrization: p(0)=2/3, p(1)=0, p'=-(4/3)sg; Dbar(0)=2/3", ok_par)

    # ---------------- W2: boundary-slice sanity ----------------
    c0 = sp.expand(core.subs(sg, 0))
    Q4 = 16 * th**4 - 64 * th**3 + 113 * th**2 - 102 * th + 45
    Q12 = (16384 * th**12 - 184320 * th**11 + 973056 * th**10
           - 3175424 * th**9 + 7104992 * th**8 - 11424528 * th**7
           + 13441257 * th**6 - 11537370 * th**5 + 7055343 * th**4
           - 2920860 * th**3 + 745767 * th**2 - 100602 * th + 6561)
    c0_claim = 1048576 * (th - 1)**4 * (4 * th - 3)**4 * Q4**2 * Q12
    rep("W2 c0 == 1048576 (th-1)^4 (4th-3)^4 Q4^2 Q12 exactly",
        sp.expand(c0 - c0_claim) == 0)
    rep("W2 Q4, Q12 have 0 roots in [0,1] (Sturm)",
        sp.Poly(Q4, th).count_roots(0, 1) == 0
        and sp.Poly(Q12, th).count_roots(0, 1) == 0)
    sext = sp.Poly(9 * p**6 - 24 * p**5 + 16 * p**4 + 4 * p**2 - 8 * p + 4, p)
    ok_D0 = (sp.expand(W.subs(D, 0)
                       - 16384 * p**4 * (p - 1)**4 * sext.as_expr()) == 0
             and sext.count_roots() == 0)
    rep("W2 W(p,0) = 16384 p^4 (p-1)^4 * sextic, sextic has no real roots",
        ok_D0)

    # ---------------- R1: bulk, subdivision ----------------
    log("R1 certify_box on [1/64,1] x [0,1] (subdivision; the slow step) ...")
    res, nodes, stuck, why = certify_box(
        core, (sg, th),
        [(SG0, sp.Integer(1)), (sp.Integer(0), sp.Integer(1))],
        maxdepth=44, nodecap=40000)
    log(f"R1: result={res}, nodes={nodes}, stuck={len(stuck)}, why={why}")
    rep("R1 core > 0 closed-strict on [1/64,1] x [0,1]", res is True)

    # ---------------- R2: thin strips outside windows, root boxes ----------
    Tc, dgc = to_tensor(core, (sg, th))
    for tag, box in (("R2a [0,1/64]x[0,11/16]",
                      [(sp.Integer(0), SG0), (sp.Integer(0), W1LO)]),
                     ("R2b [0,1/64]x[13/16,15/16]",
                      [(sp.Integer(0), SG0), (W1HI, W2LO)])):
        Tb = box_tensor(dict(Tc), dgc, box)
        mn, mx, nz = bern_coeff_min(Tb, dgc)
        log(f"{tag}: root-box Bernstein min ~ {float(mn):.4g}, zeros {nz}")
        rep(f"{tag} closed-strict at root box", mn > 0 and nz == 0)

    # ---------------- K1/K2: corner blow-ups ----------------
    random.seed(20260704)
    for th0, name, mrange in ((sp.Rational(3, 4), "K1 corner (0,3/4)", "full"),
                              (sp.Integer(1), "K2 corner (0,1)", "neg")):
        C = sp.Poly(sp.expand(core.subs(th, th0 + mu)), sg, mu)
        tab = list(zip(C.monoms(), C.coeffs()))
        d0 = min(i + j for (i, j), _ in tab)
        rep(f"{name}: min total degree of core(sg, th0+mu) == 4", d0 == 4)

        # exact closed form of the (positive definite) degree-4 leading form
        F4 = sum(c * sg**i * mu**j for (i, j), c in tab if i + j == 4)
        if th0 == sp.Rational(3, 4):
            F4ref = 3**28 * (sp.Rational(1, 2) * sg**4
                             + sp.Rational(1, 81)
                             * (sg - sp.Rational(2, 3) * mu)**4)
        else:
            F4ref = 2**34 * ((mu - 2 * sg)**4 + 128 * sg**4)
        rep(f"{name}: degree-4 form == positive definite closed form",
            sp.expand(F4 - F4ref) == 0)

        HA, HAm, HBp, HBm = {}, {}, {}, {}
        for (i, j), c in tab:
            HA[(i + j - 4, j)] = c
            HAm[(i + j - 4, j)] = c * (-1) ** j
            HBp[(i + j - 4, i)] = c
            HBm[(i + j - 4, i)] = c * (-1) ** j
        ok_coll = all(len(dd) == len(tab) for dd in (HA, HAm, HBp, HBm))
        rep(f"{name}: chart monomial maps are collision-free", ok_coll)
        HAe = sp.Poly.from_dict(HA, sg, m).as_expr()
        HAme = sp.Poly.from_dict(HAm, sg, m).as_expr()
        HBpe = sp.Poly.from_dict(HBp, nu, s2).as_expr()
        HBme = sp.Poly.from_dict(HBm, nu, s2).as_expr()

        # EXACT polynomial blow-up identities (monomial-map substitutions
        # keep these cheap; they eliminate any construction-bug risk)
        Ce = C.as_expr()
        if mrange == "full":
            id_jobs = [
                ("C(sg, sg m) == sg^4 HA", sp.expand(Ce.subs(mu, sg * m))
                 - sg**4 * HAe),
                ("C(nu s, +nu) == nu^4 HB+",
                 sp.expand(Ce.subs({sg: nu * s2, mu: nu})) - nu**4 * HBpe),
                ("C(nu s, -nu) == nu^4 HB-",
                 sp.expand(Ce.subs({sg: nu * s2, mu: -nu})) - nu**4 * HBme),
            ]
        else:
            id_jobs = [
                ("C(sg, -sg m') == sg^4 HAm",
                 sp.expand(Ce.subs(mu, -sg * m)) - sg**4 * HAme),
                ("C(nu s, -nu) == nu^4 HB-",
                 sp.expand(Ce.subs({sg: nu * s2, mu: -nu})) - nu**4 * HBme),
            ]
        for idtag, dexpr in id_jobs:
            rep(f"{name}: exact identity {idtag}", sp.expand(dexpr) == 0)

        if mrange == "full":
            chart_jobs = [
                ("HA  on [0,1/64] x [-1,1]", HAe, (sg, m),
                 [(sp.Integer(0), SG0), (sp.Integer(-1), sp.Integer(1))]),
                ("HB+ on [0,1/16] x [0,1]", HBpe, (nu, s2),
                 [(sp.Integer(0), H), (sp.Integer(0), sp.Integer(1))]),
                ("HB- on [0,1/16] x [0,1]", HBme, (nu, s2),
                 [(sp.Integer(0), H), (sp.Integer(0), sp.Integer(1))]),
            ]
        else:  # corner 2: only mu <= 0 is needed (th <= 1)
            chart_jobs = [
                ("HAm on [0,1/64] x [0,1]  (m in [-1,0])", HAme, (sg, m),
                 [(sp.Integer(0), SG0), (sp.Integer(0), sp.Integer(1))]),
                ("HB- on [0,1/16] x [0,1]", HBme, (nu, s2),
                 [(sp.Integer(0), H), (sp.Integer(0), sp.Integer(1))]),
            ]
        for ctag, ce, vs, box in chart_jobs:
            Tch, dch = to_tensor(ce, vs)
            Tb = box_tensor(dict(Tch), dch, box)
            mn, mx, nz = bern_coeff_min(Tb, dch)
            log(f"{name} {ctag}: min ~ {float(mn):.4g}, zeros {nz}")
            rep(f"{name} {ctag} closed-strict at root box", mn > 0 and nz == 0)

    # ---------------- COV: coverage arithmetic ----------------
    ok_cov = (W1LO == sp.Rational(3, 4) - H and W1HI == sp.Rational(3, 4) + H
              and W2LO == 1 - H and SG0 < H)
    rep("COV window arithmetic: [11/16,13/16] = 3/4 -+ 1/16; 15/16 = 1-1/16",
        ok_cov)
    log("COV union: sg in [1/64,1]: R1;  sg in (0,1/64]:")
    log("    th in [0,11/16] R2a | [11/16,13/16] K1 | [13/16,15/16] R2b |")
    log("    [15/16,1] K2  ==> core > 0 on (0,1) x [0,1] \\ {(0,3/4),(0,1)}")
    log("    touch points have sg = 0 (p = 2/3): OUTSIDE the open strip")

    # ---------------- NUM: independent numeric sanity ----------------
    okn, npts = True, 0
    for _ in range(200):
        sv = sp.Rational(random.randint(1, 511), 512)
        tv = sp.Rational(random.randint(1, 256), 256)
        pv = PSUB.subs(sg, sv)
        Dv = (DSUB.subs({sg: sv, th: tv}))
        val = W.subs({p: pv, D: Dv})  # exact rational
        okn &= val > 0
        npts += 1
    rep(f"NUM exact W > 0 at {npts} random strip points (sanity)", okn)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("REALNESS AT s = pi PROVEN, self-contained: the spectrum of the")
    print("replica matrix Q(p, D, w=-1) is real (5 real eigenvalues) for all")
    print("0 < p < 2/3, 0 < D <= Dbar(p) -- split f1*f4, disc f4 = evens*W,")
    print("W > 0 (regions + order-4 blow-ups), Sturm anchor + connectedness.")
    return PASS


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)

#!/usr/bin/env python3
r"""
lemma_7_34_a5_realness_pi.py
============================================================================
BUG-009-D / Route B (A=5): REALNESS OF THE REPLICA SPECTRUM AT s = pi,
proven on the TWO-PIECE CURVE-RESTRICTED Gray domain -- fully
self-contained.

THEOREM (this script certifies every step, exact rational arithmetic):
  For A = 5 and all (sg, D) with sg in (0,1) and
      0 < D <= curve*(sg) * Dbar(sg),
  the 5x5 pattern-quotient replica matrix Q(p, D, w = -1) has 5 REAL
  eigenvalues, where
      p = (4/5)(1 - sg^2)                  (all p in (0, 4/5)),
      Dbar(sg) = (4/5)(1-sg)^2/(1+sg^2),
      curve*(sg) = curveM(sg) := 1 - (3/4) sg^2      on (0, 7/10],
                   curveP(sg) := (9 - 5 sg^2)/10     on [7/10, 1),
  (curveP >= curveM on [7/10, 3/4]: curveP - curveM = (5sg^2-2)/20 > 0
  there -- gated).  The claim is proven PER PIECE:
      piece M:  sg in (0, 3/4],  0 < D <= curveM(sg) Dbar(sg),
      piece P:  sg in [7/10, 1), 0 < D <= curveP(sg) Dbar(sg),
  and the union of the two pieces covers the full curve* domain.

WHY TWO PIECES AND A CURVE (vs the A=3 full-strip lemma): full-strip
realness is FALSE at A=5 -- the complexification island is worse than at
A=4 (it reaches DOWN to th ~ 0.727 at sg ~ 0.64 at s = pi, scout map).
This script EXHIBITS the island exactly (check ISL): at (sg, th) =
(16/25, 87/100) -- strictly inside the open superdomain strip, strictly
ABOVE curveM (87/100 > 433/625 = curveM(16/25)) -- the sign object W
below is < 0 (exact rational) and f4 has exactly 2 real roots (exact
Sturm): two eigenvalues of Q are genuinely complex there.  A SINGLE
quadratic curve 1 - a sg^2 cannot both clear the island and thread the
(A)/(B) corridor at A=5 (scout G4: the corridor pinches to
a in (0.5988, 0.6036) near sg -> 1 while the island forces a >= ~0.72
at sg ~ 0.62): the two-piece splice is necessary, not an artifact.

PROOF STRUCTURE (same architecture as the committed A=4 lemma
lemma_7_34_a4_realness_pi.py; everything derived in-script from the
replica matrix):
  (a) SPLIT: at w = -1, eta = 2D(D-1)/(2D^2-5D+4) is real rational; with
      c = p(p-1)(2D^2-5D+4)^2, M = cQ is a polynomial matrix; chiM =
      charpoly(M) splits exactly as f1 * f4 with
        f1 = X - lam1M,  lam1M = (1/2) D(D-1)(5p-4)^2 (2D^2-5D+4)
      (the swap eigenvector e1 - e2 of the pattern basis) and f4 monic
      quartic.  Q's spectrum is real iff f4's four roots are real.
  (b) DISC: disc_X(f4), recomputed here (1398 terms, deg (34,48)),
      satisfies
        262144 disc = (2D^2-5D+4)^2 (D-1)^8 D^8 (15p-16)^2 (5p-4)^12 * W
      by exact division (W = 572 terms, deg (20,28), primitive).  Every
      cofactor is an even power of a polynomial that does not vanish on
      the open curve domain (2D^2-5D+4 has discriminant -7 < 0; D = 0,
      D = 1, p = 4/5, p = 16/15 are outside), so sign(disc f4) = sign(W)
      there.
  (c) W-CLOSURE ON THE CURVE, PER PIECE (the analytic core):
      piece M:  W(p(sg), u curveM(sg) Dbar(sg)) > 0 on (0,3/4] x (0,1],
      piece P:  W(p(sg), u curveP(sg) Dbar(sg)) > 0 on [7/10,1) x (0,1]
      -- regions + order-4 corner blow-ups below.
  (d) ANCHOR + CONTINUITY per piece: piece M at (sg,u) = (1/2,1/2),
      i.e. (p,D) = (3/5, 13/200); piece P at (sg,u) = (4/5,1/2), i.e.
      (p,D) = (36/125, 29/5125): f4 has 4 distinct real roots (exact
      Sturm) and disc != 0.  Each piece domain is connected and
      disc f4 != 0 on it by (b)+(c); complex roots of a real monic
      polynomial can only appear through a root collision (disc = 0)
      along a path, so the count 4-real is constant per piece.  With
      f1's rational eigenvalue: 5 real.                            QED

METHOD (everything exact rational; Bernstein = closed-box strict
positivity: min coefficient > 0 AND zero exact-zero coefficients):

  (1) Rationalize each piece:  p = (4/5)(1-sg^2), D = u curve(sg)
      Dbar(sg), u in (0,1].  Exactly (checks W1M/W1P):
        piece M:  W = 2^22 (1-sg)^4 coreM(sg,u) / (5^28 (1+sg^2)^26),
        piece P:  W = 2^36 (1-sg)^4 coreP(sg,u) / (5^56 (1+sg^2)^26),
      coreM, coreP in Z[sg,u] PRIMITIVE, 3321 terms each, deg (144,28).
      All cofactors > 0 on sg in (0,1), so W > 0 on a piece <=> core > 0
      on the piece's (sg,u) domain.
  (2) Boundary geography (checks W2M/W2P).
      piece M: coreM(0,u) = c0(u) vanishes in [0,1] exactly at u = 5/8
      and u = 1 (both order 4):
        c0 = 2^28 (u-1)^4 (8u-5)^4 Q4(u)^2 Q12(u),
        Q4 = 64u^4 - 256u^3 + 443u^2 - 410u + 175,
      Q4 and Q12 (deg 12) root-free on [0,1] (Sturm).  These are the
      A=5 avatars of the A=4 sg=0 touches (there u = 2/3, 1; here
      u = 5/8, 1 -- both map to D = 1/2 resp. D = Dbar(0) at sg = 0,
      p = pmax = 4/5: OUTSIDE the open domain).  Also coreM(sg,0)
      (deg 76) and coreM(3/4,u) (deg 28) are root-free on the piece
      (Sturm).
      piece P: coreP(1,u) == 2^44 5^52 (a positive constant -- the
      sg=1 edge is uniformly clear, the exact A=4 pattern);
      coreP(7/10,u) (deg 28) and coreP(sg,0) (deg 76) root-free
      (Sturm).  coreP has NO zeros on [7/10,1] x [0,1]: piece P is one
      closed-strict Bernstein certificate, no blow-ups needed.
  (3) Piece M decomposition (sg0 = 1/64, window half-width h = 1/16):
        R1  = [1/64, 3/4] x [0, 1]         Bernstein, subdivision
        R2a = [0, 1/64] x [0, 9/16]        Bernstein, subdivision
        R2b = [0, 1/64] x [11/16, 15/16]   Bernstein, subdivision
        K1  = [0, 1/64] x [9/16, 11/16]    blow-up at (sg,u) = (0, 5/8)
        K2  = [0, 1/64] x [15/16, 1]       blow-up at (sg,u) = (0, 1)
      R1 u R2a u R2b covers everything outside the corner windows;
      K1/K2 cover the windows minus the two touch points themselves.
  (4) Corner blow-ups (order 4, the committed A=3/A=4 technique).  For
      u0 in {5/8, 1} write C(sg,mu) = coreM(sg, u0 + mu) = sum c_ij
      sg^i mu^j.  VERIFIED exactly: every monomial has total degree
      i + j >= 4, and the degree-4 form F4 is positive definite with
      EXACT closed forms
        K1:  F4 = (5^20/4) [ (12 mu - 15 sg)^4 + (25 sg)^4 ],
        K2:  F4 = 2^46    [ (3 mu - 6 sg)^4  + (4 sg)^4  ]
      (plus exact Sturm definiteness gates).  Define the order-4
      projective blow-up charts (POLYNOMIAL identities, monomial
      bijections (i,j) -> (i+j-4, j) resp. (i+j-4, i)):
        HA (sg, m) := sum c_ij sg^(i+j-4) m^j
                       ==>  C(sg, sg*m)   = sg^4 * HA(sg, m)
        HB+-(nu, s) := sum c_ij (+-1)^j s^i nu^(i+j-4)
                       ==>  C(nu*s, +-nu) = nu^4 * HB+-(nu, s)
      SOUNDNESS: if HA > 0 on the CLOSED box [0,sg0] x [-1,1] and HB+-
      > 0 on the CLOSED box [0,h] x [0,1], then for every (sg,mu) !=
      (0,0) with 0 <= sg <= sg0, |mu| <= h:
        |mu| <= sg (then sg > 0):  core = sg^4 HA(sg, mu/sg)        > 0,
        |mu| >  sg:                core = |mu|^4 HB_sgn(|mu|, sg/|mu|) > 0.
      Hence coreM > 0 on the whole window except the touch point
      (0,u0), which is outside the open domain.  K2 needs only mu <= 0
      (u <= 1): charts HAm (m' = -m in [0,1]) and HB-.
  (5) Every Bernstein certificate uses the FIXED tensor engine of the
      committed A=3/A=4 lemmas (closed-box strict semantics; the
      repaired h**b affine map), run in Fraction arithmetic and
      cross-checked coefficient-exactly against the sympy engine (E0d)
      on a real certificate box.  Subdivision is sound: closed-strict
      sub-boxes glue.

CHECKS (each prints PASS/FAIL; OVERALL at the end):
  E0   engine self-tests (repaired affine map; exact-zero accounting;
       strict certificate + False witness; sympy == Fraction on a real
       box).
  DRV  full in-script derivation: Q(p,D,-1); M = cQ polynomial; chiM =
       f1 * f4 by exact division; disc_X(f4) (1398 terms, deg (34,48));
       262144 disc == evens * W by exact division with zero remainder;
       W = 572 terms, deg (20,28), primitive; cofactor nonvanishing
       (2D^2-5D+4 disc = -7).
  ANCH exact Sturm anchors: piece M (sg,u) = (1/2,1/2) <=> (p,D) =
       (3/5, 13/200); piece P (4/5,1/2) <=> (36/125, 29/5125): 4
       distinct real roots each, disc != 0.
  ISL  the island is real and is excluded by the curve: at (sg,th) =
       (16/25, 87/100): th < 1, th > curveM = 433/625, D = 7047/110125
       exactly, W < 0 exactly, f4 has exactly 2 real roots (Sturm)
       => full-strip realness FALSE at A=5; the claim is sharp-in-kind.
  W1M/W1P  reduction identities per piece (cleared denominator, exact
       numerator reconstruction, core primitive, sizes, edge constants,
       parametrization sanity).
  W2M/W2P  boundary-slice geography per piece (c0 factorization exact;
       Sturm root-freeness of Q4, Q12, u=0 and edge slices; coreP(1,.)
       == 2^44 5^52).
  R1/R2a/R2b   piece M region certificates (closed-strict, subdivision)
  K1/K2        piece M corner blow-ups (exact identities + closed-strict
       chart certificates)
  R1P  piece P single certificate on [7/10,1] x [0,1] (closed-strict)
  COV  coverage arithmetic per piece + the two-piece union (curveP >=
       curveM on the overlap; (0,3/4] u [7/10,1) = (0,1))
  NUM  independent sanity: exact rational evaluation of the ORIGINAL W
       at 100 random points per piece

HONEST SCOPE:
  * The claim is the two-piece curve domain.  Whether it contains the
    whole A=5 Gray region (D_c < curve* Dbar) is the (B)-side of the
    separate cert5-on-gray deliverable (the scout's piecewise Sturm
    inclusion passed for both pieces; that lemma is NOT this script).
  * On the full superdomain strip the statement is FALSE (check ISL);
    no full-strip claim is made.
  * Interior angles 0 < s < pi are the SEPARATE deliverable
    lemma_7_34_a5_realness_interior.py; this script is s = pi only.
  * All load-bearing steps are exact rational; no floats anywhere.

RESULT (numbers from the recorded run of this script): every Bernstein
certificate except R1 passes AT ITS ROOT BOX (R1 needs 13 subdivision
nodes -- the island-adjacent tight zone near (sg,u) ~ (0.62, 1), where
the island bottom clears the curve by only ~3.3% in u); piece P is ONE
root-box certificate.  The whole script runs in ~4 min (measured 238 s).
Adversarial cross-check (scratchpad dev6, not part of this script):
direct mpmath dps-60 eig of Q(w=-1) at 320 hostile curve-domain points
(island-adjacent u -> 1 on piece M, both sg -> 0 touch corners, piece P
sg -> 1) gives max |Im|/rho ~ 2e-55, while the island witness point
shows a genuine complex pair (|Im|/rho ~ 1.9e-2) -- matching ISL.

Runtime: ~4 min single-threaded.  Deps: sympy only.
NO external inputs -- everything is derived in-script from the
replica-matrix construction.  Python: /Users/para/.venvs/rnr/bin/python.
"""
import itertools
import random
import sys
import time
from fractions import Fraction
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
# Exact tensor Bernstein engine -- committed A=3/A=4 semantics (closed-box
# strict: a box is certified only if min > 0 AND nzeros == 0).  Two
# implementations of the SAME algorithm, cross-checked in E0d.
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
        # v -> lo + h v:  v^a = sum_b C(a,b) lo^(a-b) h^b v^b  (h power is b,
        # NOT a: the h**a version silently evaluates every lo != 0 sub-box on
        # the wrong region -- the historical A=3 cert5 bug class)
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


# ----------------------- Fraction fast path (same maps) --------------------
def _frx(x):
    if isinstance(x, Fraction):
        return x
    return Fraction(int(sp.numer(x)), int(sp.denom(x)))


def to_ftensor(polyE, vs):
    P = sp.Poly(polyE, *vs)
    degs = [P.degree(v) for v in vs]
    T = {}
    for mon, c in zip(P.monoms(), P.coeffs()):
        T[mon] = Fraction(int(sp.numer(c)), int(sp.denom(c)))
    return T, degs


def axis_apply_f(T, M, axis, newdim):
    out = {}
    for mon, c in T.items():
        a = mon[axis]
        for b in range(newdim):
            f = M.get((b, a))
            if f:
                nm = list(mon); nm[axis] = b; nm = tuple(nm)
                out[nm] = out.get(nm, Fraction(0)) + f * c
    return {k: v for k, v in out.items() if v != 0}


def box_tensor_f(T, degs, box):
    for ax, (lo, hi) in enumerate(box):
        d = degs[ax]
        lo = _frx(lo); hi = _frx(hi)
        h = hi - lo
        M = {}
        pl = [Fraction(1)] * (d + 1)
        ph = [Fraction(1)] * (d + 1)
        for i in range(1, d + 1):
            pl[i] = pl[i - 1] * lo
            ph[i] = ph[i - 1] * h
        for a in range(d + 1):
            for b in range(a + 1):
                M[(b, a)] = comb(a, b) * pl[a - b] * ph[b]
        T = axis_apply_f(T, M, ax, d + 1)
    return T


def bern_coeff_min_f(T, degs):
    total = 1
    for d in degs:
        total *= (d + 1)
    for ax, d in enumerate(degs):
        M = {}
        for b in range(d + 1):
            for a in range(b + 1):
                M[(b, a)] = Fraction(comb(b, a), comb(d, a))
        T = axis_apply_f(T, M, ax, d + 1)
    nzeros = total - len(T)
    if not T:
        return Fraction(0), Fraction(0), nzeros
    vals = list(T.values())
    return min(vals), max(vals), nzeros


def certify_box(polyE, vs, box, maxdepth=24, nodecap=4000, T0f=None,
                degs=None):
    """closed-strict certification with sound bisection (Fraction engine)."""
    if T0f is None:
        T0f, degs = to_ftensor(polyE, vs)
    n = len(vs)
    work = [([(_frx(lo), _frx(hi)) for lo, hi in box], 0)]
    nodes = 0
    stuck = []
    while work:
        bx, dep = work.pop()
        nodes += 1
        if nodes > nodecap:
            return None, nodes, stuck, 'nodecap'
        Tb = box_tensor_f(dict(T0f), degs, bx)
        mn, mx, nz = bern_coeff_min_f(Tb, degs)
        if mn > 0 and nz == 0:
            continue
        center = {v: sp.Rational((lo + hi).numerator,
                                 (lo + hi).denominator * 2)
                  for v, (lo, hi) in zip(vs, bx)}
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
        mid = (lo + hi) / 2
        b1 = list(bx); b1[i] = (lo, mid)
        b2 = list(bx); b2[i] = (mid, hi)
        work.append((b1, dep + 1))
        work.append((b2, dep + 1))
    return (True if not stuck else None), nodes, stuck, 'done'


p, D = sp.symbols('p D')
sg, u, mu = sp.symbols('sigma u mu')
m, nu, s2 = sp.symbols('m nu s2')

SG0 = sp.Rational(1, 64)     # corner sg-thickness = R1 left edge
H = sp.Rational(1, 16)       # corner window half-width in u
W1LO, W1HI = sp.Rational(9, 16), sp.Rational(11, 16)   # 5/8 -+ 1/16
W2LO = sp.Rational(15, 16)                             # 1 - 1/16

PSUB = sp.Rational(4, 5) * (1 - sg**2)
DBAR = sp.Rational(4, 5) * (1 - sg)**2 / (1 + sg**2)
CURVE_M = 1 - sp.Rational(3, 4) * sg**2
CURVE_P = (9 - 5 * sg**2) / sp.Integer(10)

# island witness (exact rationals): sg = 16/25, th = 87/100 (ABOVE curveM)
ISL_P = sp.Rational(1476, 3125)
ISL_D = sp.Rational(7047, 110125)


def engine_selftest(core=None):
    """Pin both engine implementations on examples with known exact answers."""
    x, y = sp.symbols('x y')
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
    Tf, dgf = to_ftensor(sp.expand(x**2), (x,))
    Tbf = box_tensor_f(dict(Tf), dgf, [(sp.Rational(1, 2), sp.Integer(1))])
    coeff = dict(Tbf)
    for ax, d in enumerate(dgf):
        M = {}
        for b in range(d + 1):
            for a in range(b + 1):
                M[(b, a)] = Fraction(comb(b, a), comb(d, a))
        coeff = axis_apply_f(coeff, M, ax, d + 1)
    ok_a &= (coeff.get((0,)) == Fraction(1, 4)
             and coeff.get((1,)) == Fraction(1, 2)
             and coeff.get((2,)) == 1)
    rep("E0a Bernstein coeffs of x^2 on [1/2,1] == {1/4,1/2,1} (both)", ok_a)
    T, dg = to_tensor(x, (x,))
    mn, mx, nz = bern_coeff_min(box_tensor(dict(T), dg,
                                           [(sp.Integer(0), sp.Integer(1))]),
                                dg)
    Tf, dgf = to_ftensor(x, (x,))
    mnf, mxf, nzf = bern_coeff_min_f(
        box_tensor_f(dict(Tf), dgf, [(sp.Integer(0), sp.Integer(1))]), dgf)
    rep("E0b f = x on [0,1]: min > 0 but nzeros == 1 (both engines)",
        mn > 0 and nz == 1 and mnf > 0 and nzf == 1)
    res1, _, _, _ = certify_box((x - 1)**2 + sp.Rational(1, 9), (x,),
                                [(sp.Integer(0), sp.Integer(1))])
    res2, _, wit, why2 = certify_box(x * y - sp.Rational(1, 2), (x, y),
                                     [(sp.Integer(0), sp.Integer(1)),
                                      (sp.Integer(0), sp.Integer(1))])
    rep("E0c certify_box: (x-1)^2 + 1/9 True; x*y - 1/2 False witness",
        res1 is True and res2 is False and why2 == 'witness')
    if core is not None:
        box = [(sp.Integer(0), SG0), (sp.Integer(0), W1LO)]
        Tc, dgc = to_tensor(core, (sg, u))
        mn_s, mx_s, nz_s = bern_coeff_min(box_tensor(dict(Tc), dgc, box), dgc)
        Tcf, dgcf = to_ftensor(core, (sg, u))
        mn_f, mx_f, nz_f = bern_coeff_min_f(
            box_tensor_f(dict(Tcf), dgcf, box), dgcf)
        okd = (Fraction(int(sp.numer(mn_s)), int(sp.denom(mn_s))) == mn_f
               and Fraction(int(sp.numer(mx_s)), int(sp.denom(mx_s))) == mx_f
               and nz_s == nz_f)
        rep("E0d sympy engine == Fraction engine on the R2a box (exact)", okd)


def reduce_on_curve(W, curve, denpow5, ncont2, tag):
    """Substitute p = PSUB, D = u*curve*DBAR into W; gate the exact shape
    W == ncont2 (1-sg)^4 core / (5^denpow5 (1+sg^2)^26); return core."""
    n_, d_ = sp.fraction(sp.cancel(sp.together(
        W.subs({p: PSUB, D: u * curve * DBAR}))))
    ok_den = sp.expand(d_ - 5**denpow5 * (1 + sg**2)**26) == 0
    rep(f"W1{tag} cleared denominator == 5^{denpow5} (1+sigma^2)^26", ok_den)
    Pn = sp.Poly(n_, sg, u)
    f1s = sp.Poly(1 - sg, sg, u)
    corePoly, mult = Pn, 0
    while mult < 4:
        q_, r_ = sp.div(corePoly, f1s, sg, u)
        if r_ == 0 and not q_.is_zero:
            corePoly, mult = sp.Poly(q_, sg, u), mult + 1
        else:
            break
    cont = sp.gcd_list(corePoly.coeffs())
    core = sp.expand(corePoly.as_expr() / cont)
    corePoly = sp.Poly(core, sg, u)
    ok_strip = (mult == 4 and cont == ncont2
                and sp.expand(n_ - ncont2 * (1 - sg)**4 * core) == 0
                and sp.gcd_list(corePoly.coeffs()) == 1)
    rep(f"W1{tag} numerator == {ncont2} (1-sigma)^4 * core; core primitive",
        ok_strip)
    rep(f"W1{tag} core size: 3321 terms, deg (144, 28)",
        len(corePoly.terms()) == 3321 and corePoly.degree(sg) == 144
        and corePoly.degree(u) == 28)
    return core


def main():
    print("=" * 78)
    print("A=5 realness at s = pi on the TWO-PIECE curve-restricted domain:")
    print("W(p,D) > 0 for p = (4/5)(1-sg^2), 0 < D <= curve*(sg) Dbar(sg),")
    print("curve* = curveM = 1-(3/4)sg^2 (sg <= 3/4) / curveP = (9-5sg^2)/10")
    print("(sg >= 7/10); exact Bernstein + order-4 corner blow-ups")
    print("=" * 78)

    # ---------------- DRV: derive everything from the replica matrix -------
    X = sp.Symbol('X')
    log("DRV building Q(p, D, w=-1) (pattern-quotient replica matrix) ...")
    Aval = 5

    def _pat(tr):
        x0, u_, up = tr
        if x0 == u_ == up: return 0
        if x0 == u_ and u_ != up: return 1
        if x0 == up and u_ != up: return 2
        if u_ == up and x0 != u_: return 3
        return 4

    # eta at w = -1 is real rational: E = -D/((A-1)(1-D))
    K5 = 2 * D**2 - 5 * D + 4
    Ew = -D / ((Aval - 1) * (1 - D))
    etaw = sp.cancel(((Aval - 1) * D * Ew + D - (Aval - 1) * Ew)
                     / ((Aval - 1) * D * Ew + D - (Aval - 1)))
    rep("DRV eta(w=-1) == 2D(D-1)/(2D^2-5D+4) (real rational)",
        sp.cancel(etaw - 2 * D * (D - 1) / K5) == 0)
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
    c_ = p * (p - 1) * K5**2
    Mm = sp.zeros(5, 5)
    ok_poly = True
    for i in range(5):
        for j in range(5):
            n2, d2 = sp.fraction(sp.cancel(c_ * Qm[i, j]))
            ok_poly &= bool(d2.is_number)
            Mm[i, j] = sp.expand(n2 / d2)
    rep("DRV M = p(p-1)(2D^2-5D+4)^2 Q is a polynomial matrix", ok_poly)
    chiM = sp.expand(Mm.charpoly(X).as_expr())
    log("DRV charpoly(M) computed")
    # explicit linear factor: the swap eigenvalue of M (eigenvector e1 - e2)
    lam1M = sp.expand(sp.Rational(1, 2) * D * (D - 1) * (5 * p - 4)**2 * K5)
    f1p = X - lam1M
    q4, r4 = sp.div(chiM, f1p, X)
    if r4 != 0:      # sign convention fallback
        f1p = X + lam1M
        q4, r4 = sp.div(chiM, f1p, X)
    f4p = sp.expand(q4)
    rep("DRV chiM == f1 * f4 exactly (explicit linear factor divides)",
        r4 == 0 and sp.expand(chiM - f1p * f4p) == 0)
    rep("DRV f4 is MONIC quartic in X with polynomial coefficients",
        sp.degree(f4p, X) == 4 and sp.Poly(f4p, X).LC() == 1)
    log("DRV disc_X(f4) (the slow discriminant) ...")
    disc = sp.discriminant(sp.Poly(f4p, X))
    discE = sp.expand(disc.as_expr() if isinstance(disc, sp.Poly) else disc)
    Pd = sp.Poly(discE, p, D)
    rep("DRV disc_X(f4) recomputed: 1398 terms, deg (34, 48)",
        len(Pd.terms()) == 1398 and Pd.degree(p) == 34 and Pd.degree(D) == 48)
    # W by exact division with the EXPLICIT even cofactors
    evens = sp.expand(K5**2 * (D - 1)**8 * D**8 * (15 * p - 16)**2
                      * (5 * p - 4)**12)
    qW, rW = sp.div(sp.expand(262144 * discE), evens, p, D)
    W = sp.expand(qW)
    PW = sp.Poly(W, p, D)
    rep("DRV 262144 disc == evens * W exactly (zero remainder)",
        rW == 0 and sp.expand(262144 * discE - evens * W) == 0)
    rep("DRV W: 572 terms, deg_p 20, deg_D 28, primitive",
        len(PW.terms()) == 572 and PW.degree(p) == 20
        and PW.degree(D) == 28 and sp.gcd_list(PW.coeffs()) == 1)
    # cofactors never vanish on the open curve domain: 2D^2-5D+4 has
    # disc = 25 - 32 = -7; D = 0, D = 1, p = 4/5, p = 16/15 are all outside
    # {0 < p < 4/5, 0 < D <= curve* Dbar < 4/5}
    rep("DRV cofactor 2D^2-5D+4 has no real roots (disc = -7 < 0)",
        sp.discriminant(sp.Poly(K5, D)) == -7)

    # ---------------- ANCH: Sturm anchors for the continuity step ----------
    ok_anch = True
    for tag, sgA, uA, curve, pA, DA in (
            ("M", sp.Rational(1, 2), sp.Rational(1, 2), CURVE_M,
             sp.Rational(3, 5), sp.Rational(13, 200)),
            ("P", sp.Rational(4, 5), sp.Rational(1, 2), CURVE_P,
             sp.Rational(36, 125), sp.Rational(29, 5125))):
        ok_on = (PSUB.subs(sg, sgA) == pA
                 and sp.cancel(uA * (curve * DBAR).subs(sg, sgA) - DA) == 0)
        f4_anchor = sp.Poly(f4p.subs({p: pA, D: DA}), X)
        n_real = f4_anchor.count_roots()
        disc_anchor = discE.subs({p: pA, D: DA})
        ok_anch &= (ok_on and n_real == 4 and disc_anchor != 0)
        rep(f"ANCH piece {tag} ({sgA},{uA}) -> ({pA},{DA}): 4 distinct "
            f"real roots", ok_on and n_real == 4 and disc_anchor != 0)

    # ---------------- ISL: the island is real and curve-excluded ----------
    sI = sp.Rational(16, 25)
    thI = sp.Rational(87, 100)
    ok_geo = (thI < 1 and thI > CURVE_M.subs(sg, sI)
              and CURVE_M.subs(sg, sI) == sp.Rational(433, 625)
              and PSUB.subs(sg, sI) == ISL_P
              and sp.cancel(thI * DBAR.subs(sg, sI) - ISL_D) == 0)
    rep("ISL witness (16/25, 87/100): inside strip, ABOVE curveM 433/625",
        ok_geo)
    WI = W.subs({p: ISL_P, D: ISL_D})
    rep("ISL exact W < 0 at the island point (full-strip realness FALSE)",
        WI < 0)
    f4_isl = sp.Poly(f4p.subs({p: ISL_P, D: ISL_D}), X)
    rep("ISL f4 at the island point has exactly 2 real roots (Sturm)",
        f4_isl.count_roots() == 2)

    # ---------------- W1: reductions to core, per piece -----------
    log("W1M substituting piece M into W (the slow expand) ...")
    coreM = reduce_on_curve(W, CURVE_M, 28, 2**22, "M")
    log("W1P substituting piece P into W ...")
    coreP = reduce_on_curve(W, CURVE_P, 56, 2**36, "P")
    ok_par = (PSUB.subs(sg, 0) == sp.Rational(4, 5) and PSUB.subs(sg, 1) == 0
              and DBAR.subs(sg, 0) == sp.Rational(4, 5)
              and DBAR.subs(sg, 1) == 0
              and CURVE_M.subs(sg, 0) == 1
              and CURVE_M.subs(sg, sp.Rational(3, 4)) == sp.Rational(37, 64)
              and CURVE_P.subs(sg, sp.Rational(7, 10)) == sp.Rational(131, 200)
              and CURVE_P.subs(sg, 1) == sp.Rational(2, 5))
    rep("W1 parametrization: p(0)=4/5, p(1)=0; Dbar(0)=4/5; curve values",
        ok_par)

    # ---------------- E0: engine self-tests (needs coreM for E0d) ----------
    engine_selftest(coreM)

    # ---------------- W2M: piece M boundary geography ----------------
    c0 = sp.expand(coreM.subs(sg, 0))
    Q4 = 64 * u**4 - 256 * u**3 + 443 * u**2 - 410 * u + 175
    pref = 2**28 * (u - 1)**4 * (8 * u - 5)**4 * Q4**2
    q12, r12 = sp.div(c0, sp.expand(pref), u)
    Q12 = sp.expand(q12)
    ok_c0 = (r12 == 0 and sp.expand(c0 - pref * Q12) == 0
             and sp.degree(Q12, u) == 12)
    rep("W2M c0 == 2^28 (u-1)^4 (8u-5)^4 Q4^2 * Q12 exactly", ok_c0)
    rep("W2M Q4 root-free on [0,1] (Sturm) and Q4(0) = 175 > 0",
        sp.Poly(Q4, u).count_roots(0, 1) == 0
        and sp.Poly(Q4, u).eval(0) == 175)
    rep("W2M Q12 root-free on [0,1] (Sturm) and Q12(0) > 0",
        sp.Poly(Q12, u).count_roots(0, 1) == 0
        and sp.Poly(Q12, u).eval(0) > 0)
    # u = 0 slice (D = 0): deg 76, root-free on the piece
    e0 = sp.Poly(sp.expand(coreM.subs(u, 0)), sg)
    rep("W2M coreM(sg,0) (deg 76) root-free on [0,3/4] (Sturm), >0 at 1/2",
        e0.degree() == 76
        and e0.count_roots(0, sp.Rational(3, 4)) == 0
        and e0.eval(sp.Rational(1, 2)) > 0)
    # sg = 3/4 edge: deg 28, root-free (also covered by R1's closed box)
    e34 = sp.Poly(sp.expand(coreM.subs(sg, sp.Rational(3, 4))), u)
    rep("W2M coreM(3/4,u) (deg 28) root-free on [0,1] (Sturm), >0 at 1/2",
        e34.count_roots(0, 1) == 0 and e34.eval(sp.Rational(1, 2)) > 0)
    # u = 1 edge: spot sanity only (deg-144 Sturm avoided as at A=4; the
    # closed edge is certified by R1 + K2)
    e1 = sp.Poly(sp.expand(coreM.subs(u, 1)), sg)
    rep("W2M coreM(sg,1) > 0 at sg = 1/64, 1/2, 47/64 (spot; edge closure "
        "is R1+K2)",
        e1.eval(sp.Rational(1, 64)) > 0 and e1.eval(sp.Rational(1, 2)) > 0
        and e1.eval(sp.Rational(47, 64)) > 0)

    # ---------------- W2P: piece P boundary geography ----------------
    e1P = sp.expand(coreP.subs(sg, 1))
    rep("W2P coreP(1,.) == 2^44 5^52 (positive constant: sg=1 edge clear)",
        sp.expand(e1P - 2**44 * 5**52) == 0)
    e7P = sp.Poly(sp.expand(coreP.subs(sg, sp.Rational(7, 10))), u)
    rep("W2P coreP(7/10,u) (deg 28) root-free on [0,1] (Sturm), >0 at 1/2",
        e7P.count_roots(0, 1) == 0 and e7P.eval(sp.Rational(1, 2)) > 0)
    e0P = sp.Poly(sp.expand(coreP.subs(u, 0)), sg)
    rep("W2P coreP(sg,0) (deg 76) root-free on [7/10,1] (Sturm), >0 at 4/5",
        e0P.count_roots(sp.Rational(7, 10), 1) == 0
        and e0P.eval(sp.Rational(4, 5)) > 0)

    # precompute the Fraction tensors once
    TcfM, dgcfM = to_ftensor(coreM, (sg, u))
    TcfP, dgcfP = to_ftensor(coreP, (sg, u))

    # ---------------- R1: piece M bulk, subdivision ----------------
    log("R1 certify_box on [1/64,3/4] x [0,1] (subdivision; slow) ...")
    res, nodes, stuck, why = certify_box(
        coreM, (sg, u),
        [(SG0, sp.Rational(3, 4)), (sp.Integer(0), sp.Integer(1))],
        maxdepth=44, nodecap=40000, T0f=TcfM, degs=dgcfM)
    log(f"R1: result={res}, nodes={nodes}, stuck={len(stuck)}, why={why}")
    rep("R1 coreM > 0 closed-strict on [1/64,3/4] x [0,1]", res is True)

    # ---------------- R2: piece M thin strips outside windows ----------
    for tag, box in (("R2a [0,1/64]x[0,9/16]",
                      [(sp.Integer(0), SG0), (sp.Integer(0), W1LO)]),
                     ("R2b [0,1/64]x[11/16,15/16]",
                      [(sp.Integer(0), SG0), (W1HI, W2LO)])):
        res, nodes, stuck, why = certify_box(
            coreM, (sg, u), box, maxdepth=44, nodecap=40000,
            T0f=TcfM, degs=dgcfM)
        log(f"{tag}: result={res}, nodes={nodes}, why={why}")
        rep(f"{tag} closed-strict", res is True)

    # ---------------- K1/K2: piece M corner blow-ups ----------------
    for u0, name, mrange in ((sp.Rational(5, 8), "K1 corner (0,5/8)", "full"),
                             (sp.Integer(1), "K2 corner (0,1)", "neg")):
        C = sp.Poly(sp.expand(coreM.subs(u, u0 + mu)), sg, mu)
        tab = list(zip(C.monoms(), C.coeffs()))
        d0 = min(i + j for (i, j), _ in tab)
        rep(f"{name}: min total degree of coreM(sg, u0+mu) == 4", d0 == 4)

        # exact closed form of the (positive definite) degree-4 leading form
        F4 = sp.expand(sum(c * sg**i * mu**j for (i, j), c in tab
                           if i + j == 4))
        if u0 == sp.Rational(5, 8):
            F4ref = sp.Rational(5**20, 4) * ((12 * mu - 15 * sg)**4
                                             + (25 * sg)**4)
        else:
            F4ref = 2**46 * ((3 * mu - 6 * sg)**4 + (4 * sg)**4)
        rep(f"{name}: degree-4 form == positive definite closed form",
            sp.expand(F4 - F4ref) == 0)
        # belt-and-suspenders definiteness by exact Sturm
        F41 = sp.Poly(F4.subs(sg, 1), mu)
        F404 = F4.subs({sg: 0, mu: 1})   # mu^4 coefficient
        rep(f"{name}: F4 positive definite (F4(1,m) rootfree, F4(1,0)>0, "
            f"F4(0,1)>0)",
            F41.count_roots() == 0 and F41.eval(0) > 0 and F404 > 0)

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

        # EXACT polynomial blow-up identities
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
        else:  # corner 2: only mu <= 0 is needed (u <= 1)
            chart_jobs = [
                ("HAm on [0,1/64] x [0,1]  (m in [-1,0])", HAme, (sg, m),
                 [(sp.Integer(0), SG0), (sp.Integer(0), sp.Integer(1))]),
                ("HB- on [0,1/16] x [0,1]", HBme, (nu, s2),
                 [(sp.Integer(0), H), (sp.Integer(0), sp.Integer(1))]),
            ]
        for ctag, ce, vs, box in chart_jobs:
            res, nodes, stuck, why = certify_box(ce, vs, box,
                                                 maxdepth=44, nodecap=40000)
            log(f"{name} {ctag}: result={res}, nodes={nodes}, why={why}")
            rep(f"{name} {ctag} closed-strict", res is True)

    # ---------------- R1P: piece P single certificate ----------------
    log("R1P certify_box on [7/10,1] x [0,1] (subdivision; slow) ...")
    res, nodes, stuck, why = certify_box(
        coreP, (sg, u),
        [(sp.Rational(7, 10), sp.Integer(1)), (sp.Integer(0), sp.Integer(1))],
        maxdepth=44, nodecap=40000, T0f=TcfP, degs=dgcfP)
    log(f"R1P: result={res}, nodes={nodes}, stuck={len(stuck)}, why={why}")
    rep("R1P coreP > 0 closed-strict on [7/10,1] x [0,1]", res is True)

    # ---------------- COV: coverage arithmetic ----------------
    ok_cov = (W1LO == sp.Rational(5, 8) - H and W1HI == sp.Rational(5, 8) + H
              and W2LO == 1 - H and SG0 < H)
    rep("COV window arithmetic: [9/16,11/16] = 5/8 -+ 1/16; 15/16 = 1-1/16",
        ok_cov)
    # two-piece union: curveP >= curveM on the overlap [7/10, 3/4]
    ok_uni = (sp.expand(CURVE_P - CURVE_M - (5 * sg**2 - 2) / sp.Integer(20))
              == 0
              and (5 * sp.Rational(7, 10)**2 - 2) == sp.Rational(9, 20)
              and sp.Rational(7, 10) < sp.Rational(3, 4))
    rep("COV union: curveP - curveM == (5sg^2-2)/20 >= 9/400 on overlap; "
        "(0,3/4] u [7/10,1) = (0,1)", ok_uni)
    log("COV piece M union: sg in [1/64,3/4]: R1;  sg in (0,1/64]:")
    log("    u in [0,9/16] R2a | [9/16,11/16] K1 | [11/16,15/16] R2b |")
    log("    [15/16,1] K2  ==> coreM > 0 on (0,3/4] x [0,1] minus the two")
    log("    sg=0 touch points (which have p = 4/5: OUTSIDE the domain)")
    log("COV piece P: single closed box [7/10,1] x [0,1]")

    # ---------------- NUM: independent numeric sanity ----------------
    TW, dgW = to_ftensor(W, (p, D))
    random.seed(20260704)
    okn, npts = True, 0
    for piece in ("M", "P"):
        for _ in range(100):
            if piece == "M":
                sv = Fraction(random.randint(1, 383), 512)  # (0, 3/4)
                cv = 1 - Fraction(3, 4) * sv**2
            else:
                sv = Fraction(7, 10) + Fraction(random.randint(1, 152), 512)
                cv = (9 - 5 * sv**2) / 10
            uv = Fraction(random.randint(1, 256), 256)
            pv = Fraction(4, 5) * (1 - sv**2)
            Dv = uv * cv * Fraction(4, 5) * (1 - sv)**2 / (1 + sv**2)
            pp = [Fraction(1)] * (dgW[0] + 1)
            for i in range(1, dgW[0] + 1):
                pp[i] = pp[i - 1] * pv
            dd = [Fraction(1)] * (dgW[1] + 1)
            for i in range(1, dgW[1] + 1):
                dd[i] = dd[i - 1] * Dv
            val = sum(cf * pp[mon[0]] * dd[mon[1]] for mon, cf in TW.items())
            okn &= val > 0
            npts += 1
    rep(f"NUM exact W > 0 at {npts} random curve-domain points (sanity)",
        okn)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("REALNESS AT s = pi (A=5) PROVEN on the two-piece curve domain:")
    print("the spectrum of Q(p, D, w=-1) is real (5 real eigenvalues) for")
    print("all 0 < p < 4/5, 0 < D <= curve*(sg) Dbar(sg), where curve* =")
    print("curveM = 1-(3/4)sg^2 on (0,3/4] and curveP = (9-5sg^2)/10 on")
    print("[7/10,1).  Full-strip realness is FALSE (exact island witness):")
    print("the curve restriction is necessary.")
    return PASS


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)

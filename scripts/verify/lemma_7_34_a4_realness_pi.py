#!/usr/bin/env python3
r"""
lemma_7_34_a4_realness_pi.py
============================================================================
BUG-009-D / Route B (A=4): REALNESS OF THE REPLICA SPECTRUM AT s = pi,
proven on the CURVE-RESTRICTED Gray domain -- fully self-contained.

THEOREM (this script certifies every step, exact rational arithmetic):
  For A = 4, all sg in (0,1) and all D with 0 < D <= thetatilde4(sg) *
  Dbar(sg), the 5x5 pattern-quotient replica matrix Q(p, D, w=-1) has
  5 REAL eigenvalues, where
      p = (3/4)(1 - sg^2)              (all p in (0, 3/4)),
      Dbar(sg)       = (3/4)(1-sg)^2/(1+sg^2),
      thetatilde4(sg) = 1 - sg^2/2     (the committed A=4 separating
                                        curve of cert5-on-Gray).
  Since D_c(p) < thetatilde4 * Dbar for all p (certificate (B) of the
  committed lemma_7_34_a4_cert5_on_gray.py), this covers the WHOLE A=4
  Gray region at the binding angle s = pi.

WHY A CURVE-RESTRICTED CLAIM (vs the A=3 full-strip lemma): full-strip
realness is FALSE at A=4.  A genuine complexification ISLAND lives inside
the superdomain at t = 2 (scout geography: sg in ~[0.614, 0.698] at
th = 1, floor th ~ 0.939 at sg ~ 0.649).  This script EXHIBITS the island
exactly (check ISL): at (sg, th) = (13/20, 97/100) -- strictly inside the
open superdomain strip, strictly ABOVE the curve (97/100 > 631/800 =
thetatilde4(13/20)) -- the sign object W below is < 0 (exact rational)
and f4 has exactly 2 real roots (exact Sturm): two eigenvalues of Q are
genuinely complex there.  The curve restriction is necessary, not an
artifact; the island clears the curve by ~0.15 in th.

PROOF STRUCTURE (all derived in-script from the replica matrix itself,
same architecture as the committed A=3 lemma
lemma_7_34_realness_pi_w_closure.py):
  (a) SPLIT: at w = -1, eta = 2D(D-1)/(2D^2-4D+3) is real rational; with
      c = p(p-1)(2D^2-4D+3)^2, M = cQ is a polynomial matrix; chiM =
      charpoly(M) splits exactly as f1 * f4 with
        f1 = X - lam1M,  lam1M = (2/3) D(D-1)(4p-3)^2 (2D^2-4D+3)
      (the swap eigenvector e1 - e2 of the pattern basis) and f4 monic
      quartic.  Q's spectrum is real iff f4's four roots are real.
  (b) DISC: disc_X(f4), recomputed here (1397 terms, deg (34,48)), equals
        (2D^2-4D+3)^2 (D-1)^8 D^8 (8p-9)^2 (4p-3)^12 * W(p, D)
      by exact division (W = 571 terms, deg (20, 28)).  Every cofactor is
      an even power of a polynomial that does not vanish on the open
      curve-domain (2D^2-4D+3 has discriminant -8 < 0; D = 0, D = 1,
      p = 3/4, p = 9/8 are outside), so sign(disc f4) = sign(W) there.
  (c) W-CLOSURE ON THE CURVE (the analytic core):
        W(p(sg), u * thetatilde4(sg) * Dbar(sg)) > 0
      for all (sg, u) in (0,1) x (0,1] -- regions + blow-ups below.
  (d) ANCHOR + CONTINUITY: at (sg,u) = (1/2,1/2), i.e. (p,D) =
      (9/16, 21/320), f4 has 4 distinct real roots (exact Sturm) and
      disc != 0.  The curve domain (0,1) x (0,1] is connected and
      disc f4 != 0 on it by (b)+(c); complex roots of a real monic
      polynomial can only appear through a root collision (disc = 0)
      along a path, so the count 4-real is constant.  With f1's rational
      eigenvalue: 5 real.                                          QED

METHOD (everything exact rational; Bernstein = closed-box strict
positivity: min coefficient > 0 AND zero exact-zero coefficients):

  (1) Rationalize the curve domain:  p = (3/4)(1-sg^2),
      D = u * (1 - sg^2/2) * Dbar(sg),  (sg,u) in (0,1) x (0,1].
      Exactly (check W1):
        W(p(sg), D(sg,u)) = 3^12 (1-sg)^4 * core(sg,u)
                            / (2^60 (1+sg^2)^22),
      core in Z[sg,u] PRIMITIVE (content 1), 3093 terms, deg (136, 28).
      All cofactors > 0 on sg in (0,1), so W > 0 on the curve domain
      <=> core > 0 on (0,1) x (0,1].
  (2) Boundary geography (check W2).  core(0,u) = c0(u) vanishes in
      [0,1] exactly at u = 2/3 and u = 1 (both order 4):
        c0 = 2^26 (u-1)^4 (3u-2)^4 (3u^2-8u+6)^2 (3u^2-4u+4)^2 Q12(u),
      the quadratics have negative discriminants (-8, -32) and Q12 has 0
      roots in [0,1] (Sturm).  These are the A=4 avatars of the A=3
      touches (there at th = 3/4, 1; here the sg=0 endpoint has
      theta_c -> 1 and the curve -> 1, and the touches sit at u = 2/3,
      1).  Both touch points have sg = 0, i.e. p = 3/4: OUTSIDE the open
      domain.  Also core(1, u) == 2^86 (a positive constant -- the sg=1
      edge is uniformly clear) and core(sg, 0) has 0 roots in [0,1]
      (u = 0 slice, exact factor shape + Sturm).
  (3) Decompose [0,1] x [0,1] (sg0 = 1/64, window half-width h = 1/16):
        R1  = [1/64, 1] x [0, 1]           Bernstein, subdivision
        R2a = [0, 1/64] x [0, 29/48]       Bernstein, subdivision
        R2b = [0, 1/64] x [35/48, 15/16]   Bernstein, subdivision
        K1  = [0, 1/64] x [29/48, 35/48]   blow-up at (sg,u) = (0, 2/3)
        K2  = [0, 1/64] x [15/16, 1]       blow-up at (sg,u) = (0, 1)
      R1 u R2a u R2b covers everything outside the corner windows;
      K1/K2 cover the windows minus the two touch points themselves.
  (4) Corner blow-ups (order 4, the committed A=3 technique).  For
      u0 in {2/3, 1} write C(sg,mu) = core(sg, u0 + mu) = sum c_ij
      sg^i mu^j.  VERIFIED exactly: every monomial has total degree
      i + j >= 4, and the degree-4 form F4 is positive definite
      (exact Sturm definiteness gates: F4(1,m) has 0 real roots,
      F4(1,0) > 0, F4(0,1) > 0 -- definiteness is what makes the chart
      certificates close at their boxes; the certificates themselves are
      the load-bearing step).  Define the order-4 projective blow-up
      charts (POLYNOMIAL identities, monomial bijections
      (i,j) -> (i+j-4, j) resp. (i+j-4, i)):
        HA (sg, m) := sum c_ij sg^(i+j-4) m^j
                       ==>  C(sg, sg*m)   = sg^4 * HA(sg, m)
        HB+-(nu, s) := sum c_ij (+-1)^j s^i nu^(i+j-4)
                       ==>  C(nu*s, +-nu) = nu^4 * HB+-(nu, s)
      SOUNDNESS: if HA > 0 on the CLOSED box [0,sg0] x [-1,1] and HB+-
      > 0 on the CLOSED box [0,h] x [0,1], then for every (sg,mu) !=
      (0,0) with 0 <= sg <= sg0, |mu| <= h:
        |mu| <= sg (then sg > 0):  core = sg^4 HA(sg, mu/sg)        > 0,
        |mu| >  sg:                core = |mu|^4 HB_sgn(|mu|, sg/|mu|) > 0.
      Hence core > 0 on the whole window except the touch point (0,u0),
      which is outside the open domain.  K2 needs only mu <= 0 (u <= 1):
      charts HAm (m' = -m in [0,1]) and HB-.
  (5) Every Bernstein certificate uses the FIXED tensor engine of the
      committed A=3 lemma (closed-box strict semantics; the repaired
      h**b affine map), run in Fraction arithmetic for speed and
      cross-checked coefficient-exactly against the sympy engine (E0d)
      on a real certificate box.  Subdivision is sound: closed-strict
      sub-boxes glue.

CHECKS (each prints PASS/FAIL; OVERALL at the end):
  E0   engine self-tests: exact Bernstein coefficients of x^2 on the
       lo != 0 box [1/2,1] == {1/4, 1/2, 1} (pins the repaired h**b
       affine map) in BOTH engines; exact-zero accounting for f = x on
       [0,1]; a strict certificate and a False witness (f = xy - 1/2);
       E0d sympy-vs-Fraction identical (min,max,nzeros) on the R2a box.
  DRV  full in-script derivation: Q(p,D,-1) (eta real rational);
       M = cQ polynomial; chiM = charpoly(M); chiM == f1 * f4 by exact
       division (f4 monic quartic); disc_X(f4) (1397 terms, deg
       (34,48)); W := disc / evens by exact division with zero
       remainder; W = 571 terms, deg (20,28); cofactor nonvanishing
       (2D^2-4D+3 disc = -8).
  ANCH exact Sturm anchor at (sg,u) = (1/2,1/2) <=> (p,D) = (9/16,
       21/320) (exact substitution identity): 4 distinct real roots,
       disc != 0.
  ISL  the island is real and is excluded by the curve: at (sg,th) =
       (13/20, 97/100): th < 1 (inside the open superdomain strip),
       th > thetatilde4 = 631/800 (above the curve), D = 14259/227600
       exactly, W < 0 exactly, and f4 has exactly 2 real roots (Sturm)
       => full-strip realness FALSE at A=4; the claim is sharp-in-kind.
  W1   reduction identity: cleared denominator == 2^60 (1+sg^2)^22,
       numerator == 3^12 (1-sg)^4 * core (exact reconstruction), core
       PRIMITIVE, 3093 terms, deg (136,28); core(1,.) == 2^86 != 0;
       parametrization sanity p(0)=3/4, p(1)=0, Dbar(0)=3/4, Dbar(1)=0,
       thetatilde4(0)=1, thetatilde4(1)=1/2, D(sg,1) == thetatilde4*Dbar.
  W2   c0 = core(0,.) factorization (exact division, quadratic
       discriminants < 0, Q12 Sturm-free on [0,1]) => c0 zeros in [0,1]
       are exactly {2/3, 1}, both order 4; core(sg,0) exact factor shape
       with 0 roots in [0,1]; core(sg,1) spot positivity (the closed
       u=1 edge itself is certified by R1 + K2, see COV; a deg-136
       Sturm sequence is a memory hazard and is deliberately avoided).
  R1   certify core > 0 closed-strict on [1/64,1] x [0,1] (subdivision)
  R2   certify core > 0 closed-strict on [0,1/64] x [0,29/48] and
       [0,1/64] x [35/48,15/16] (subdivision allowed)
  K1   corner (0,2/3): min total degree == 4; F4 positive definite
       (Sturm); chart construction collision-free; EXACT POLYNOMIAL
       blow-up identities (expand(C(sg, sg*m) - sg^4 HA) == 0 etc.);
       certify HA on [0,1/64] x [-1,1], HB+ and HB- on [0,1/16] x [0,1]
  K2   corner (0,1): same with charts HAm on [0,1/64] x [0,1] (m'=-m)
       and HB- on [0,1/16] x [0,1] (only mu <= 0 is needed: u <= 1)
  COV  coverage arithmetic: 29/48 = 2/3 - 1/16, 35/48 = 2/3 + 1/16,
       15/16 = 1 - 1/16; shared sg-boundary 1/64 < 1/16; union =
       (0,1) x [0,1] minus the two sg=0 touch points
  NUM  independent sanity (not part of the proof): exact rational
       evaluation of the ORIGINAL W at 200 random curve-domain points

HONEST SCOPE:
  * The claim is the curve-restricted domain th <= thetatilde4(sg) --
    which contains the whole Gray region D <= D_c by certificate (B) of
    the committed cert5-on-gray lemma.  On the full superdomain strip
    the statement is FALSE (check ISL); no full-strip claim is made.
  * Interior angles 0 < s < pi at A=4 are a SEPARATE deliverable (the
    A=3 analogue is lemma_7_34_realness_interior.py); this script is the
    s = pi endpoint only.
  * Numeric guards run at mpmath dps >= 40 where used (NUM is exact
    rational; float64 flips signs on these polynomials).

RESULT (numbers from the recorded run of this script): every Bernstein
certificate except R1 passes AT ITS ROOT BOX (R1 needs 9 subdivision
nodes, max depth 4 -- exactly the A=3 shape: the entire difficulty was
the domain geometry, not the polynomial); the whole script runs in
~3 min.  Adversarial cross-check (scratchpad explore5, not part of this
script): direct mpmath dps-60 eig of Q at 437 hostile curve-domain
points (u -> 1, island-adjacent sg at u = 1, degenerate corners) gives
max |Im|/rho ~ 1e-57, while the island point shows a genuine complex
pair (|Im|/rho ~ 8e-3) -- matching ISL.

Runtime: ~3 min single-threaded (measured 164 s; the A=4 core is deg
(136,28) vs A=3's (88,28)).
Deps: sympy only (all load-bearing steps exact rational; no floats).
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
# Exact tensor Bernstein engine -- committed A=3 semantics (closed-box
# strict: a box is certified only if min > 0 AND nzeros == 0; exact-zero
# coefficients are dropped from the sparse dict and recovered from the dense
# count; a zero VERTEX coefficient equals a zero of the polynomial at that
# box corner).  Two implementations of the SAME algorithm:
#   * sympy Rational (verbatim from lemma_7_34_realness_pi_w_closure.py),
#   * Fraction (fast path for the heavy runs),
# cross-checked coefficient-exactly in E0d on a real certificate box.
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
    """sympy Rational/Integer or Fraction -> Fraction (exact)."""
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
W1LO, W1HI = sp.Rational(29, 48), sp.Rational(35, 48)
W2LO = sp.Rational(15, 16)

PSUB = sp.Rational(3, 4) * (1 - sg**2)
DBAR = sp.Rational(3, 4) * (1 - sg)**2 / (1 + sg**2)
THETA4 = 1 - sg**2 / 2
DSUB = u * THETA4 * DBAR    # curve-restricted D(sg, u)

# island witness (exact rationals): sg = 13/20, th = 97/100 (ABOVE the curve)
ISL_P = sp.Rational(693, 1600)
ISL_D = sp.Rational(14259, 227600)


def engine_selftest(core=None):
    """Pin both engine implementations on examples with known exact answers."""
    x, y = sp.symbols('x y')
    # (a) lo != 0 affine map: Bernstein coeffs of x^2 on [1/2, 1] are exactly
    #     {b0, b1, b2} = {lo^2, lo*hi, hi^2} = {1/4, 1/2, 1} -- both engines
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
    # (b) exact-zero accounting: f = x on [0,1] has b0 = f(0) = 0 -> nzeros=1
    T, dg = to_tensor(x, (x,))
    mn, mx, nz = bern_coeff_min(box_tensor(dict(T), dg,
                                           [(sp.Integer(0), sp.Integer(1))]), dg)
    Tf, dgf = to_ftensor(x, (x,))
    mnf, mxf, nzf = bern_coeff_min_f(
        box_tensor_f(dict(Tf), dgf, [(sp.Integer(0), sp.Integer(1))]), dgf)
    rep("E0b f = x on [0,1]: min > 0 but nzeros == 1 (both engines)",
        mn > 0 and nz == 1 and mnf > 0 and nzf == 1)
    # (c) strict certificate and a False witness
    res1, _, _, _ = certify_box((x - 1)**2 + sp.Rational(1, 9), (x,),
                                [(sp.Integer(0), sp.Integer(1))])
    res2, _, wit, why2 = certify_box(x * y - sp.Rational(1, 2), (x, y),
                                     [(sp.Integer(0), sp.Integer(1)),
                                      (sp.Integer(0), sp.Integer(1))])
    rep("E0c certify_box: (x-1)^2 + 1/9 True; x*y - 1/2 False witness",
        res1 is True and res2 is False and why2 == 'witness')
    # (d) engine cross-check on a REAL certificate box (R2a) -- exact equality
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


def main():
    print("=" * 78)
    print("A=4 realness at s = pi on the CURVE-RESTRICTED Gray domain:")
    print("W(p,D) > 0 for p = (3/4)(1-sg^2), 0 < D <= (1-sg^2/2) Dbar(sg)")
    print("exact Bernstein + order-4 corner blow-ups, rational arithmetic only")
    print("=" * 78)

    # ---------------- DRV: derive everything from the replica matrix -------
    X = sp.Symbol('X')
    log("DRV building Q(p, D, w=-1) (pattern-quotient replica matrix) ...")
    Aval = 4

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
    rep("DRV eta(w=-1) == 2D(D-1)/(2D^2-4D+3) (real rational)",
        sp.cancel(etaw - 2 * D * (D - 1) / (2 * D**2 - 4 * D + 3)) == 0)
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
    K4 = 2 * D**2 - 4 * D + 3
    c_ = p * (p - 1) * K4**2
    Mm = sp.zeros(5, 5)
    ok_poly = True
    for i in range(5):
        for j in range(5):
            n2, d2 = sp.fraction(sp.cancel(c_ * Qm[i, j]))
            ok_poly &= bool(d2.is_number)
            Mm[i, j] = sp.expand(n2 / d2)
    rep("DRV M = p(p-1)(2D^2-4D+3)^2 Q is a polynomial matrix", ok_poly)
    chiM = sp.expand(Mm.charpoly(X).as_expr())
    log("DRV charpoly(M) computed")
    # explicit linear factor: the swap eigenvalue of M (eigenvector e1 - e2)
    lam1M = sp.expand(sp.Rational(2, 3) * D * (D - 1) * (4 * p - 3)**2 * K4)
    f1p = X - lam1M
    q4, r4 = sp.div(chiM, f1p, X)
    if r4 != 0:      # sign convention fallback
        f1p = X + lam1M
        q4, r4 = sp.div(chiM, f1p, X)
    f4p = sp.expand(q4)
    rep("DRV chiM == f1 * f4 exactly (explicit linear factor divides)",
        r4 == 0 and sp.expand(chiM - f1p * f4p) == 0)
    rep("DRV f4 is MONIC quartic in X with polynomial coefficients",
        sp.degree(f4p, X) == 4
        and sp.Poly(f4p, X).LC() == 1)
    disc = sp.discriminant(sp.Poly(f4p, X))
    discE = sp.expand(disc.as_expr() if isinstance(disc, sp.Poly) else disc)
    Pd = sp.Poly(discE, p, D)
    rep("DRV disc_X(f4) recomputed: 1397 terms, deg (34, 48)",
        len(Pd.terms()) == 1397 and Pd.degree(p) == 34 and Pd.degree(D) == 48)
    # W by exact division with the EXPLICIT even cofactors
    evens = (K4**2 * (D - 1)**8 * D**8 * (8 * p - 9)**2 * (4 * p - 3)**12)
    qW, rW = sp.div(discE, sp.expand(evens), p, D)
    W = sp.expand(qW)
    PW = sp.Poly(W, p, D)
    rep("DRV disc == evens * W exactly (zero remainder)",
        rW == 0 and sp.expand(discE - evens * W) == 0)
    rep("DRV W: 571 terms, deg_p 20, deg_D 28",
        len(PW.terms()) == 571 and PW.degree(p) == 20 and PW.degree(D) == 28)
    # cofactors never vanish on the open curve domain: 2D^2-4D+3 has
    # disc = 16 - 24 = -8; D = 0, D = 1, p = 3/4, p = 9/8 are all outside
    # {0 < p < 3/4, 0 < D <= thetatilde4*Dbar < 3/4}
    rep("DRV cofactor 2D^2-4D+3 has no real roots (disc = -8 < 0)",
        sp.discriminant(sp.Poly(K4, D)) == -8)

    # ---------------- ANCH: Sturm anchor for the continuity step ----------
    pA = sp.Rational(9, 16)
    DA = sp.Rational(21, 320)
    ok_on = (PSUB.subs(sg, sp.Rational(1, 2)) == pA
             and sp.Rational(1, 2) * (THETA4 * DBAR).subs(sg, sp.Rational(1, 2)) == DA)
    rep("ANCH (sg,u)=(1/2,1/2) maps exactly to (p,D)=(9/16,21/320)", ok_on)
    f4_anchor = sp.Poly(f4p.subs({p: pA, D: DA}), X)
    n_real = f4_anchor.count_roots()
    disc_anchor = discE.subs({p: pA, D: DA})
    rep("ANCH f4 at (9/16,21/320): 4 distinct real roots, disc != 0",
        n_real == 4 and disc_anchor != 0)

    # ---------------- ISL: the island is real and curve-excluded ----------
    sI = sp.Rational(13, 20)
    thI = sp.Rational(97, 100)
    ok_geo = (thI < 1 and thI > THETA4.subs(sg, sI)
              and THETA4.subs(sg, sI) == sp.Rational(631, 800)
              and PSUB.subs(sg, sI) == ISL_P
              and thI * DBAR.subs(sg, sI) == ISL_D)
    rep("ISL witness (13/20, 97/100): inside strip, ABOVE curve 631/800",
        ok_geo)
    WI = W.subs({p: ISL_P, D: ISL_D})
    rep("ISL exact W < 0 at the island point (full-strip realness FALSE)",
        WI < 0)
    f4_isl = sp.Poly(f4p.subs({p: ISL_P, D: ISL_D}), X)
    rep("ISL f4 at the island point has exactly 2 real roots (Sturm)",
        f4_isl.count_roots() == 2)

    # ---------------- W1: reduction to core on the curve domain -----------
    log("W1 substituting the curve domain into W (the slow expand) ...")
    n_, d_ = sp.fraction(sp.cancel(sp.together(W.subs({p: PSUB, D: DSUB}))))
    ok_den = sp.expand(d_ - 2**60 * (1 + sg**2)**22) == 0
    rep("W1 cleared denominator == 2^60 (1+sigma^2)^22 (positive)", ok_den)

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
    ok_strip = (mult == 4 and cont == 3**12
                and sp.expand(n_ - 3**12 * (1 - sg)**4 * core) == 0
                and sp.gcd_list(corePoly.coeffs()) == 1)
    rep("W1 numerator == 3^12 (1-sigma)^4 * core exactly; core primitive",
        ok_strip)
    log(f"core: {len(corePoly.terms())} terms, deg_sg {corePoly.degree(sg)}, "
        f"deg_u {corePoly.degree(u)}")
    rep("W1 core size: 3093 terms, deg (136, 28)",
        len(corePoly.terms()) == 3093 and corePoly.degree(sg) == 136
        and corePoly.degree(u) == 28)
    rep("W1 core(1,.) == 2^86 (positive constant: sg=1 edge clear)",
        sp.expand(core.subs(sg, 1) - 2**86) == 0)
    ok_par = (PSUB.subs(sg, 0) == sp.Rational(3, 4) and PSUB.subs(sg, 1) == 0
              and DBAR.subs(sg, 0) == sp.Rational(3, 4)
              and DBAR.subs(sg, 1) == 0
              and THETA4.subs(sg, 0) == 1
              and THETA4.subs(sg, 1) == sp.Rational(1, 2)
              and sp.cancel(DSUB.subs(u, 1) - THETA4 * DBAR) == 0)
    rep("W1 parametrization: p(0)=3/4, p(1)=0; Dbar(0)=3/4; "
        "D(sg,1)=thetatilde4*Dbar", ok_par)

    # ---------------- E0: engine self-tests (needs core for E0d) ----------
    engine_selftest(core)

    # ---------------- W2: boundary-slice geography ----------------
    c0 = sp.expand(core.subs(sg, 0))
    Q2a = 3 * u**2 - 8 * u + 6
    Q2b = 3 * u**2 - 4 * u + 4
    pref = 2**26 * (u - 1)**4 * (3 * u - 2)**4 * Q2a**2 * Q2b**2
    q12, r12 = sp.div(c0, sp.expand(pref), u)
    Q12 = sp.expand(q12)
    ok_c0 = (r12 == 0 and sp.expand(c0 - pref * Q12) == 0
             and sp.degree(Q12, u) == 12)
    rep("W2 c0 == 2^26 (u-1)^4 (3u-2)^4 Q2a^2 Q2b^2 * Q12 exactly", ok_c0)
    rep("W2 Q2a, Q2b have negative discriminants (-8, -32)",
        sp.discriminant(sp.Poly(Q2a, u)) == -8
        and sp.discriminant(sp.Poly(Q2b, u)) == -32)
    rep("W2 Q12 has 0 roots in [0,1] (Sturm) and Q12(0) > 0",
        sp.Poly(Q12, u).count_roots(0, 1) == 0
        and sp.Poly(Q12, u).eval(0) > 0)
    # u = 0 slice (D = 0): exact factor shape, no roots in [0,1]
    e0 = sp.expand(core.subs(u, 0))
    pref0 = 2**44 * (sg + 1)**4 * (3 * sg**2 + 1)**4 * (sg**2 + 1)**22
    q0, r0 = sp.div(e0, sp.expand(pref0), sg)
    S12 = sp.expand(q0)
    ok_u0 = (r0 == 0 and sp.expand(e0 - pref0 * S12) == 0
             and sp.degree(S12, sg) == 12
             and sp.Poly(S12, sg).count_roots(0, 1) == 0
             and sp.Poly(S12, sg).eval(sp.Rational(1, 2)) > 0)
    rep("W2 core(sg,0) == 2^44 (sg+1)^4 (3sg^2+1)^4 (sg^2+1)^22 * S12, "
        "S12 rootfree on [0,1]", ok_u0)
    # u = 1 edge (the curve itself): NOT Sturm-counted here -- the deg-136
    # Sturm sequence is a memory hazard and the closed edge is already
    # covered by the R1 box (sg >= 1/64, u = 1 face included) and the K2
    # charts (sg <= 1/64); spot values are sanity only
    e1 = sp.Poly(sp.expand(core.subs(u, 1)), sg)
    rep("W2 core(sg,1) > 0 at sg = 1/64, 1/2, 63/64 (spot sanity; edge "
        "closure is R1+K2)",
        e1.eval(sp.Rational(1, 64)) > 0 and e1.eval(sp.Rational(1, 2)) > 0
        and e1.eval(sp.Rational(63, 64)) > 0)

    # precompute the Fraction tensor of core once
    Tcf, dgcf = to_ftensor(core, (sg, u))

    # ---------------- R1: bulk, subdivision ----------------
    log("R1 certify_box on [1/64,1] x [0,1] (subdivision; the slow step) ...")
    res, nodes, stuck, why = certify_box(
        core, (sg, u),
        [(SG0, sp.Integer(1)), (sp.Integer(0), sp.Integer(1))],
        maxdepth=44, nodecap=40000, T0f=Tcf, degs=dgcf)
    log(f"R1: result={res}, nodes={nodes}, stuck={len(stuck)}, why={why}")
    rep("R1 core > 0 closed-strict on [1/64,1] x [0,1]", res is True)

    # ---------------- R2: thin strips outside windows ----------
    for tag, box in (("R2a [0,1/64]x[0,29/48]",
                      [(sp.Integer(0), SG0), (sp.Integer(0), W1LO)]),
                     ("R2b [0,1/64]x[35/48,15/16]",
                      [(sp.Integer(0), SG0), (W1HI, W2LO)])):
        res, nodes, stuck, why = certify_box(
            core, (sg, u), box, maxdepth=44, nodecap=40000,
            T0f=Tcf, degs=dgcf)
        log(f"{tag}: result={res}, nodes={nodes}, why={why}")
        rep(f"{tag} closed-strict", res is True)

    # ---------------- K1/K2: corner blow-ups ----------------
    for u0, name, mrange in ((sp.Rational(2, 3), "K1 corner (0,2/3)", "full"),
                             (sp.Integer(1), "K2 corner (0,1)", "neg")):
        C = sp.Poly(sp.expand(core.subs(u, u0 + mu)), sg, mu)
        tab = list(zip(C.monoms(), C.coeffs()))
        d0 = min(i + j for (i, j), _ in tab)
        rep(f"{name}: min total degree of core(sg, u0+mu) == 4", d0 == 4)

        # exact closed form of the (positive definite) degree-4 leading form
        F4 = sp.expand(sum(c * sg**i * mu**j for (i, j), c in tab
                           if i + j == 4))
        if u0 == sp.Rational(2, 3):
            F4ref = sp.Rational(2**46, 3**10) * ((3 * mu - 4 * sg)**4
                                                 + 3072 * sg**4)
        else:
            F4ref = 2**26 * 3**4 * ((mu - 2 * sg)**4 + 12 * sg**4)
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

    # ---------------- COV: coverage arithmetic ----------------
    ok_cov = (W1LO == sp.Rational(2, 3) - H and W1HI == sp.Rational(2, 3) + H
              and W2LO == 1 - H and SG0 < H)
    rep("COV window arithmetic: [29/48,35/48] = 2/3 -+ 1/16; 15/16 = 1-1/16",
        ok_cov)
    log("COV union: sg in [1/64,1]: R1;  sg in (0,1/64]:")
    log("    u in [0,29/48] R2a | [29/48,35/48] K1 | [35/48,15/16] R2b |")
    log("    [15/16,1] K2  ==> core > 0 on (0,1) x [0,1] \\ {(0,2/3),(0,1)}")
    log("    touch points have sg = 0 (p = 3/4): OUTSIDE the open domain")

    # ---------------- NUM: independent numeric sanity ----------------
    # exact rational evaluation of the ORIGINAL W (Fraction tensor; no
    # floats anywhere -- float64 flips signs on these polynomials)
    TW, dgW = to_ftensor(W, (p, D))
    random.seed(20260704)
    okn, npts = True, 0
    for _ in range(200):
        sv = Fraction(random.randint(1, 511), 512)
        uv = Fraction(random.randint(1, 256), 256)
        pv = Fraction(3, 4) * (1 - sv**2)
        Dv = uv * (1 - sv**2 / 2) * Fraction(3, 4) * (1 - sv)**2 / (1 + sv**2)
        pp = [Fraction(1)] * (dgW[0] + 1)
        for i in range(1, dgW[0] + 1):
            pp[i] = pp[i - 1] * pv
        dd = [Fraction(1)] * (dgW[1] + 1)
        for i in range(1, dgW[1] + 1):
            dd[i] = dd[i - 1] * Dv
        val = sum(cf * pp[mon[0]] * dd[mon[1]] for mon, cf in TW.items())
        okn &= val > 0
        npts += 1
    rep(f"NUM exact W > 0 at {npts} random curve-domain points (sanity)", okn)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("REALNESS AT s = pi (A=4) PROVEN on the curve-restricted domain:")
    print("the spectrum of Q(p, D, w=-1) is real (5 real eigenvalues) for all")
    print("0 < p < 3/4, 0 < D <= (1 - sg^2/2) Dbar -- covering ALL of Gray")
    print("(D_c < thetatilde4*Dbar, committed cert5-on-gray (B)).  Full-strip")
    print("realness is FALSE (exact island witness): the curve is necessary.")
    return PASS


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)

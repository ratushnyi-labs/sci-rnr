#!/usr/bin/env python3
r"""
lemma_7_34_a5_cert5_on_gray.py
============================================================================
CERT5-ON-GRAY (BUG-009-D / Route B, A=5): exact certificate that the k=0
charpoly sign  chi(Lambda) > 0  holds on a GRAY-CONTAINING subdomain of the
rationalized superdomain -- the certificate that is FALSE on the full
superdomain (real Perron pocket, exact witnesses chi(Lambda) < 0 at
sg = 1/2, th = 1, t = 2 and at sg = 4/5, th = 1, t = 2; check P5).

STATEMENT (all exact, A = 5).  Rationalized domain
    p = (4/5)(1 - sg^2),   Dbar = (4/5)(1-sg)^2/(1+sg^2),   D = th*Dbar,
    t = 1 - cos s in (0, 2],  sg in (0, 1),  th in (0, 1].
TWO-PIECE separating splice (rational, degree 2 each; a SINGLE quadratic
curve is DEAD at A=5 -- scout G4: the a-corridor for 1 - a sg^2 requires
a > 0.6467 at sg = 0.4 yet a < 0.6036 at sg = 0.995, an EMPTY
intersection; this is the A=3 on-gray + smallp splice pattern, not the
A=4 single-curve pattern):
    curveM(sg) = 1 - (3/4) sg^2          on  sg in (0, 3/4],
    curveP(sg) = (9 - 5 sg^2)/10         on  sg in [7/10, 1),
    overlap [7/10, 3/4] nonempty.
 (A-M) chi(Lambda) > 0  for all sg in (0, 3/4], th in (0, curveM(sg)],
       t in (0, 2].
 (A-P) chi(Lambda) > 0  for all sg in [7/10, 1), th in (0, curveP(sg)],
       t in (0, 2].
 (B-M) D_c(p) < curveM(sg) * Dbar(p)  for all sg in (0, 3/4],
 (B-P) D_c(p) < curveP(sg) * Dbar(p)  for all sg in [7/10, 1),
      where D_c is the A=5 Gray threshold = first sign change in D of the
      reduced-cubic discriminant disc(D,p) of the 3x3 symmetry-reduced
      alternating product (the committed threshold object of
      lemma_7_34_dbar_upper_bound.py, steps (i)-(iii)).
 =>   UNION LOGIC (exact splice, gate P14; stated exactly as in the A=3
      smallp splice): the sg-coverage (0, 3/4] u [7/10, 1) = (0, 1) is
      seamless because 7/10 < 3/4, and on the overlap [7/10, 3/4] the two
      certified slabs NEST: curveP - curveM = (5 sg^2 - 2)/20 > 0 there
      (5 sg^2 >= 49/20 > 2), i.e. min(curveM, curveP) = curveM, whose own
      (B-M) inclusion already places the Gray edge strictly below it.
      Hence the union region {th <= curveM(sg), sg <= 3/4} u
      {th <= curveP(sg), sg >= 7/10} contains the Gray region
      {0 < D <= D_c(p)} x {angles s in (0, pi]} for ALL p in (0, 4/5)
      (each piece by its own (B)-certificate), and chi(Lambda) > 0 holds
      on all of it (each piece by its own (A)-certificate).

SIGN LINKAGE (exact, no numerics needed).  With the t-domain assembly of
the committed pipeline (lemma_7_34_midband_chi_certificates.py at A=3,
lemma_7_34_a4_chi_certificates.py at A=4) run at A = 5:
    target := assemble(0) = (1+sg^2)^K * chat_t^5 * CrN_t^5 * chi(Lambda)
on the rationalized domain, where (check P3, exact expansions)
    chat_t = p(p-1) * B(D,t),      CrN_t = -B(D,t),
    B(D,t) = 2 D^2 t (1-D)(4-D) + (5D-4)^2  >  0   for 0 < D < 4/5
    (sum of nonnegative terms, second strictly positive),
hence  chat_t^5 CrN_t^5 = [p(1-p)]^5 B^10 > 0  on the claim region and
    sign(target) = sign(chi(Lambda))            (P4 = 140-pt numeric guard).
The identity is the charpoly expansion chi = sum (-1)^k e_k Lam^{5-k} with
e_k = Etld_k / chat^k, Lam = L*CrD_t/CrN_t, L = 1 - (3/2)D(1-D)t (target
rate R_5 >= 3/2, same L as A=3/A=4), valid once the grading gates P1 hold
(E_k symmetric at (-w)^k, chat at (-w)^1, mN = mD).

CERTIFICATE (A) STRUCTURE.  target = (1-sg)^7 * th * t * core  (exact
division, P6; all three stripped factors > 0 on the claim region).
PIECE M:  substitute  sg = (3/4) x,  th = u (64 - 27 x^2)/64
          [= u * curveM(sg), identity P6c],  t = 2 tau;
PIECE P:  substitute  sg = (7 + 3x)/10,  th = u (131 - 42x - 9x^2)/200
          [= u * curveP(sg), identity P6c],  t = 2 tau;
clear denominators (integer tensors, exact cross-check P7 + independent
Schwartz-Zippel identity check P7b, per piece).  Then on [0,1]^3 in
(x, u, tau):
  PIECE M (claim CLOSED at sg = 3/4, OPEN at sg = 0 <=> p = 4/5, the
  degenerate pmax corner -- the A=3 on-gray face pattern):
    P8   ALL scaled Bernstein coefficients >= 0, not all 0
         =>  core_u > 0 on the OPEN box (0,1)^3;
    P9   faces x=1, u=1, tau=1 (collapse): 2-var no-negatives + nonzero;
    P10  edges (x=1,u=1), (x=1,tau=1), (u=1,tau=1): 1-var no-negatives +
         nonzero; corner (1,1,1): exact positive value.
    Union = (0,1]^3 exactly = the claim region of (A-M).
  PIECE P (claim CLOSED at sg = 7/10, OPEN at sg = 1 <=> p = 0, the
  degenerate p -> 0 corner -- the A=3 smallp face pattern):
    P8   root-box no-negatives  =>  core_u > 0 on the OPEN box;
    P9   faces x=0 (restriction), u=1, tau=1: no-negatives + nonzero;
    P10  edges (x=0,u=1), (x=0,tau=1), (u=1,tau=1): 1-var no-negatives +
         nonzero; corner (0,1,1): exact positive value.
    Union = [0,1) x (0,1] x (0,1] exactly = the claim region of (A-P).
  Zero-coefficient counts are reported: zeros are allowed by the
  open-box+closure logic; the pieceM edge (u=1,tau=1) vanishes EXACTLY at
  x=0 (sg=0 <=> p=pmax) and the pieceP edge (u=1,tau=1) vanishes EXACTLY
  at x=1 (sg=1 <=> p=0) -- the two degenerate corners, which is why each
  piece's claim is open at its bad end and no closed-box certificate can
  exist there.

CERTIFICATE (B) (exact Sturm, per piece).  disc(D,p) = N(D,p) /
(2^18 (5D-4)^12)  (identity P11a; the denominator is > 0 for D < 4/5, so
sign(disc) = sign(N) there).  N = D^4 * C(D,p) with
    C(0,p) = 2^18 p^4 (1-p)^2 (5p-4)^4 (p^2-p+1)      (identity P11b),
and p^2-p+1 > 0 always (discriminant -3 < 0), so disc > 0 as D -> 0+
for every p in (0, 4/5).  On each curve D = curve*Dbar, p = (4/5)(1-sg^2):
after clearing the positive denominator const*(1+sg^2)^6 (P11c;
const = 5^12 on M, 5^22 on P) and stripping the positive-on-the-piece
factors (P11d; sg^13 (1-sg)^12 (1+sg)^2 on M, sg^8 (1-sg)^13 (1+sg)^2 on
P), the residual integer polynomial (deg 29 on M, deg 33 on P) has NO
roots in the piece interval ([0, 3/4] resp. [7/10, 1], Sturm) and is
negative there (P11e).  A continuous function of D positive as D -> 0+
and negative at D = curve*Dbar has its FIRST sign change strictly below
the curve; that first sign change IS the A=5 threshold object D_c
(committed convention: lemma_7_34_dbar_upper_bound.py steps (i)-(iii),
threshold mechanism probe_7_34_beyond_gray_structure B1).  Only the SAFE
direction is used: D_c <= (any D where disc has already turned negative
with disc > 0 before it), hence D_c < curve * Dbar on each piece.  The
identification of the first disc sign change with the all-word Gray
threshold inherits the grade of the committed threshold mechanism,
exactly as in the A=3/A=4 cert5-on-gray usage.

HONEST NOTES.
 *  The Perron pocket is REAL at A=5 (P5 exact witnesses INSIDE each
    piece's sg-range: sg=1/2 (pieceM), sg=4/5 (pieceP), both at th=1
    above the respective curve) -- the restriction to th <= curve is
    necessary, not an artifact.
 *  Near sg -> 1 the corridor pinches EXACTLY as at A=3 smallp: in
    om = 1 - sg^2 units, th_c ~ 2/5 + 0.2403 om and th_x ~ 2/5 + 0.7204
    om (scout 8), and curveP = 2/5 + om/2 rides through the pinch; P15
    checks th_c < curveP < th_x numerically at sg = 0.995.  The margins
    vanish linearly in om -- forced by the sharpness of R_5 >= 3/2 at
    p -> 0, not an artifact.
 *  The s = pi complexification island (W < 0, scouts 2b/4b: sg in
    ~[0.58, 0.77], bottom theta ~ 0.7275 at sg = 0.64, t = 2) sits INSIDE
    the superdomain strip at A=5 -- a new phenomenon vs A=3/A=4.  It is
    IRRELEVANT to this k=0 certificate (the target-sign linkage is
    algebraic, not spectral), but it kills any full-superdomain realness
    lemma: the eventual interior-angle realness certificate must be
    curve-restricted (the island clears curveM by ~+0.024 and curveP by
    ~+0.028 -- comfortable at A=4 grade, thin near sg -> 1 on curveP,
    ~0.0025 in a-units at sg = 0.995).
 *  chi(Lambda) > 0 is the k=0 member of the Budan-Fourier tower; the
    j >= 1 positive-side and all negative-side psi^(j) certificates hold
    on the FULL superdomain box at A=5 (companion scripts, root-box
    grade); curve-restricted interior-angle realness remains a separate
    open item (unchanged by this script).
 *  Numeric guards run at mpmath dps 40-80: double precision FLIPS signs
    of these polynomials near the degenerate corners (A=3 lesson,
    reproduced at A=4/A=5 in development).
 *  sp.nsimplify must NOT touch exact Rationals in the anchors P12 (A=4
    lesson: at ambient dps 40 it returns an only-approximately-equal
    radical-power form; P12a gates against recurrence).

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.
Runtime: ~2-3 min (measured 144 s; two towers + two Sturm inclusions).
Machinery:
committed A=3/A=4 pipeline (build_Q/build_Cr/Chebyshev t-domain assembly;
integer Bernstein tower of lemma_7_34_midband_cert5_on_gray.py /
lemma_7_34_a4_cert5_on_gray.py; disc_pieces of
lemma_7_34_dbar_upper_bound.py) run at AVAL = 5, TWO pieces.
"""
import itertools
import math
import random
import time
from fractions import Fraction as Fr
from math import comb

import sympy as sp
import mpmath as mp

mp.mp.dps = 40
T0 = time.time()
PASS = True


def rep(name, ok):
    global PASS
    PASS = PASS and bool(ok)
    print(f"  {name:<68} {'PASS' if ok else 'FAIL'}  [{time.time()-T0:.0f}s]",
          flush=True)
    return ok


# ---------------------------------------------------------------------------
# symbols and the rationalized domain (A = 5)
# ---------------------------------------------------------------------------
p, D, w, t = sp.symbols('p D w t')
sg, th, U, X = sp.symbols('sigma theta u x')
AVAL = 5
A = sp.Integer(AVAL)
PSUB = sp.Rational(4, 5) * (1 - sg ** 2)
DBAR = sp.Rational(4, 5) * (1 - sg) ** 2 / (1 + sg ** 2)
DSUB = th * DBAR
VS = (sg, th, t)
CURVM = 1 - sp.Rational(3, 4) * sg ** 2       # piece M curve, sg in (0, 3/4]
CURVP = (9 - 5 * sg ** 2) / 10                # piece P curve, sg in [7/10, 1)

# ---------------------------------------------------------------------------
# replica pattern-quotient machinery (committed pipeline, at A=5)
# ---------------------------------------------------------------------------
STATES = list(itertools.product(range(AVAL), repeat=3))


def _pat(tr):
    x_, u_, up = tr
    if x_ == u_ == up: return 0
    if x_ == u_ and u_ != up: return 1
    if x_ == up and u_ != up: return 2
    if u_ == up and x_ != u_: return 3
    return 4


REPS = [None] * 5
for _k, _tr in enumerate(STATES):
    if REPS[_pat(_tr)] is None:
        REPS[_pat(_tr)] = _k

eth = D / ((A - 1) * (1 - D))


def _etaf(W):
    E = eth * W
    return ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))


def build_Q():
    eta, etb = _etaf(w), _etaf(1 / w)
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(AVAL)]
         for i in range(AVAL)]
    Q = sp.zeros(5, 5)
    for a_ in range(5):
        x_, xp, xq = STATES[REPS[a_]]
        for (y, yp, yq) in STATES:
            term = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x_][y])
            if yp != y: term *= eta
            if yq != y: term *= etb
            Q[a_, _pat((y, yp, yq))] += term
    return Q


def build_Cr():
    def Cf(W):
        E = eth * W
        return ((A - 1) * D * E + D - (A - 1)) / (A * D - (A - 1))
    return sp.cancel(Cf(w) * Cf(1 / w) / (Cf(1) ** 2))


UU = sp.symbols('uu')


def laurent_sym_to_t(expr, K=60):
    """w-symmetric Laurent polynomial -> t-polynomial via Chebyshev in
    uu = w + 1/w = 2 - 2t on |w| = 1; None if asymmetric."""
    ex = sp.expand(expr)
    npoly = sp.Poly(sp.expand(ex * w ** K), w)
    co = {}
    for (md,), c in npoly.terms():
        co[md - K] = co.get(md - K, 0) + c
    ks = sorted(co.keys())
    maxk = max(abs(k) for k in ks) if ks else 0
    C = {0: sp.Integer(2), 1: UU}
    for k2 in range(2, maxk + 1):
        C[k2] = sp.expand(UU * C[k2 - 1] - C[k2 - 2])
    out = sp.Integer(0)
    done = set()
    for k2 in ks:
        if k2 in done: continue
        if k2 == 0:
            out += co[0]; done.add(0)
        else:
            kk = abs(k2)
            if sp.simplify(co.get(kk, 0) - co.get(-kk, 0)) != 0:
                return None
            out += co.get(kk, 0) * C[kk]
            done.add(kk); done.add(-kk)
    return sp.expand(out.subs(UU, 2 - 2 * t))


def sym_to_t_mid(expr, K=60):
    ex = sp.expand(expr)
    poly = sp.Poly(sp.expand(ex * w ** K), w)
    degs = [md - K for (md,), cc in poly.terms()]
    mid = sp.Rational(min(degs) + max(degs), 2)
    if int(mid) != mid:
        return None, None
    r = laurent_sym_to_t(sp.cancel(ex / (-w) ** int(mid)), K)
    return r, (int(mid) if r is not None else None)


def dom_piece(ex):
    """(Poly on (sg,th,t), (1+sg^2)-exponent, numeric const) with denominator
    exactly const*(1+sg^2)^k after the domain substitution, else None."""
    n_, d_ = sp.fraction(sp.cancel(sp.together(ex.subs({p: PSUB, D: DSUB}))))
    f = sp.factor_list(d_)
    const, expo = f[0], 0
    for b, e in f[1]:
        if sp.simplify(b - (1 + sg ** 2)) == 0:
            expo += e
        elif b.is_number:
            const *= b ** e
        else:
            return None
    return sp.Poly(sp.expand(n_), *VS), expo, const


def Q_numeric(pv, Dv, sv):
    thv = mp.log(Dv / ((AVAL - 1) * (1 - Dv)))
    E = mp.e ** (thv + 1j * sv)
    et = ((AVAL - 1) * Dv * E + Dv - (AVAL - 1) * E) / \
        ((AVAL - 1) * Dv * E + Dv - (AVAL - 1))
    Cn = ((AVAL - 1) * Dv * E + Dv - (AVAL - 1)) / (AVAL * Dv - (AVAL - 1))
    C0 = ((AVAL - 1) * Dv * mp.e ** thv + Dv - (AVAL - 1)) / \
        (AVAL * Dv - (AVAL - 1))
    Crn = abs(Cn / C0) ** 2
    Tn = [[(1 - pv) if i == j else pv / (AVAL - 1) for j in range(AVAL)]
          for i in range(AVAL)]
    etb = mp.conj(et)
    Qn = mp.zeros(5, 5)
    for a2 in range(5):
        x_, xp, xq = STATES[REPS[a2]]
        for (y, yp, yq) in STATES:
            Qn[a2, _pat((y, yp, yq))] += (Tn[xp][yp] * Tn[xq][yq] / Tn[x_][y]
                                          * (et if yp != y else 1)
                                          * (etb if yq != y else 1))
    return Qn, Crn


def chi_numeric(pv, Dv, sv, tvf):
    """chi(Lambda) from the numeric eigenvalues of Q (A=5)."""
    Qn, Crn = Q_numeric(pv, Dv, sv)
    ev = mp.eig(Qn)[0]
    Lm = (1 - 1.5 * Dv * (1 - Dv) * tvf) / Crn
    out = mp.mpc(1)
    for lam in ev:
        out *= (Lm - lam)
    return float(mp.re(out))


# ---------------------------------------------------------------------------
# the A=5 threshold object: reduced-cubic discriminant of the alternating
# product (committed convention, lemma_7_34_dbar_upper_bound.py)
# ---------------------------------------------------------------------------
def disc_pieces(Aval):
    Aq = sp.Integer(Aval)
    lam2 = 1 - D * Aq / (Aq - 1)
    al = sp.Rational(1, Aval) + (1 - sp.Rational(1, Aval)) / lam2
    be = sp.Rational(1, Aval) - sp.Rational(1, Aval) / lam2
    td = 1 - p
    to = p / (Aq - 1)

    def Bred(y):
        Ki = {0: (al if y == 0 else be), 1: (be if y == 0 else al), 2: be}
        Trow = {0: {0: td, 1: to, 2: (Aq - 2) * to},
                1: {0: to, 1: td, 2: (Aq - 2) * to},
                2: {0: to, 1: to, 2: td + (Aq - 3) * to}}
        M = sp.zeros(3, 3)
        for i in range(3):
            for j in range(3):
                M[i, j] = Ki[i] * Trow[i][j]
        return M
    M = (Bred(1) * Bred(0)).applyfunc(lambda e: sp.cancel(sp.together(e)))
    c2m = sp.cancel(M.trace())
    c1m = sp.cancel(sum(M[[i for i in r], [j for j in r]].det()
                        for r in ([0, 1], [0, 2], [1, 2])))
    c0m = sp.cancel(M.det())
    a = -c2m
    b = c1m
    c = -c0m
    disc = sp.cancel(sp.together(18 * a * b * c - 4 * a ** 3 * c
                                 + a ** 2 * b ** 2 - 4 * b ** 3 - 27 * c ** 2))
    N, Dn = sp.fraction(disc)
    return sp.expand(N), Dn


def Dc5_numeric(pv, prec_bits=60):
    """A=5 Gray threshold via alternating-product collision bisection
    (numeric consistency checks only; the exact work is Certificate B)."""
    Aq = AVAL

    def imflag(Dv):
        lam2 = 1 - Dv * Aq / (Aq - 1)
        al = mp.mpf(1) / Aq + (1 - mp.mpf(1) / Aq) / lam2
        be = mp.mpf(1) / Aq - (mp.mpf(1) / Aq) / lam2
        td = 1 - pv
        to = pv / (Aq - 1)

        def Bred(y):
            Ki = {0: (al if y == 0 else be), 1: (be if y == 0 else al), 2: be}
            Trow = {0: {0: td, 1: to, 2: (Aq - 2) * to},
                    1: {0: to, 1: td, 2: (Aq - 2) * to},
                    2: {0: to, 1: to, 2: td + (Aq - 3) * to}}
            Mx = mp.zeros(3, 3)
            for i in range(3):
                for jj in range(3):
                    Mx[i, jj] = Ki[i] * Trow[i][jj]
            return Mx
        Mx = Bred(1) * Bred(0)
        ev = mp.eig(Mx)[0]
        return max(abs(mp.im(e)) for e in ev)

    lo, hi = mp.mpf('1e-8'), mp.mpf(AVAL - 1) / AVAL - mp.mpf('1e-9')
    for _ in range(prec_bits):
        mid = (lo + hi) / 2
        if imflag(mid) < mp.mpf('1e-30'):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ---------------------------------------------------------------------------
# shared tower helpers (integer tensors)
# ---------------------------------------------------------------------------
def bern_axis(T, degs, ax):
    d_ = degs[ax]
    out = {}
    for mon, cval in T.items():
        a2 = mon[ax]
        lm = list(mon)
        for b2 in range(a2, d_ + 1):
            f = comb(d_ - a2, b2 - a2)
            lm[ax] = b2
            key = tuple(lm)
            out[key] = out.get(key, 0) + f * cval
    return {k: v for k, v in out.items() if v != 0}


def collapse_at1(T, ax, degs):
    out = {}
    for mon, cf in T.items():
        key = tuple(m for i, m in enumerate(mon) if i != ax)
        out[key] = out.get(key, 0) + cf         # value at var=1: sum coeffs
    return {k: v for k, v in out.items() if v != 0}, \
        [d_ for i, d_ in enumerate(degs) if i != ax]


def restrict_at0(T, ax, degs):
    out = {}
    for mon, cf in T.items():
        if mon[ax] != 0:
            continue                            # value at var=0: only a=0
        key = tuple(m for i, m in enumerate(mon) if i != ax)
        out[key] = out.get(key, 0) + cf
    return {k: v for k, v in out.items() if v != 0}, \
        [d_ for i, d_ in enumerate(degs) if i != ax]


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("CERT5-ON-GRAY (A=5): chi(Lambda) > 0 below the TWO-PIECE splice")
    print("curveM = 1-(3/4)sg^2 on sg in (0,3/4]  u  curveP = (9-5sg^2)/10 on")
    print("sg in [7/10,1)  (ALL p in (0,4/5)), all angles; + D_c < curve*Dbar")
    print("per piece (exact Sturm).  A = 5.  Exact integer/rational arithmetic.")
    print("=" * 78)

    # ---- build the t-domain pieces (committed assembly, A=5) --------------
    print("building Q, Cr, charpoly pieces ...", flush=True)
    Q = build_Q()
    Cr = build_Cr()
    CrN, CrD = sp.fraction(Cr)
    L = 1 - sp.Rational(3, 2) * D * (1 - D) * t
    G = D ** 2 * w + (D - (AVAL - 1)) * (1 - D)
    Gt = D ** 2 + (D - (AVAL - 1)) * (1 - D) * w
    c = sp.expand(p * (1 - p) * G * Gt)
    XX = sp.symbols('XX')
    M = sp.zeros(5, 5)
    for i in range(5):
        for j in range(5):
            n_, d_ = sp.fraction(sp.cancel(c * Q[i, j]))
            M[i, j] = sp.expand(n_ / d_)
    E = [sp.expand(co * (-1) ** k)
         for k, co in enumerate(M.charpoly(XX).all_coeffs())]
    Etld, mks = [sp.Integer(1)], [0]
    for k in range(1, 6):
        ek, mk = sym_to_t_mid(E[k])
        Etld.append(ek); mks.append(mk)
    chat, mc = sym_to_t_mid(c)
    CrNt, mN = sym_to_t_mid(CrN)
    CrDt, mD = sym_to_t_mid(CrD)

    print("-" * 78)
    print("P1  charpoly grading gates (exact)")
    rep("P1 E_k symmetric at (-w)^k, chat at (-w)^1, mN == mD",
        (None not in mks and mc == 1 and mN == mD
         and all(mks[k] == k * mc for k in range(6))))

    PC = {'L': dom_piece(L), 'CrN': dom_piece(CrNt), 'CrD': dom_piece(CrDt),
          'chat': dom_piece(chat)}
    for k in range(1, 6):
        PC[f'E{k}'] = dom_piece(Etld[k])
    print("-" * 78)
    print("P2  domain pieces")
    rep("P2 all pieces have const*(1+sigma^2)^k denominators",
        all(v is not None for v in PC.values()))

    # ---- P3: exact multiplier sign structure ------------------------------
    print("-" * 78)
    print("P3  multiplier identities (exact expansions)")
    Bpoly = 2 * D ** 2 * t * (1 - D) * (4 - D) + (5 * D - 4) ** 2
    ok3 = (sp.expand(chat - p * (p - 1) * Bpoly) == 0)
    ok3 &= (sp.expand(CrNt + Bpoly) == 0)
    ok3 &= (sp.expand(CrDt + (5 * D - 4) ** 2) == 0)
    rep("P3 chat_t = p(p-1)B, CrN_t = -B, B = 2D^2t(1-D)(4-D)+(5D-4)^2", ok3)
    rep("P3b => chat_t^5 CrN_t^5 = [p(1-p)]^5 B^10 > 0 for 0<p<1, D<4/5",
        ok3)   # algebraic consequence; recorded as a gate

    # ---- assemble target = (1+sg^2)^K * chat^5 CrN^5 * chi ---------------
    print("-" * 78)
    print("P4  target assembly + numeric sign guard")

    def assemble(j):
        terms = []
        for k in range(0, 6):
            if 5 - k - j < 0:
                continue
            fall = sp.Integer(1)
            for q_ in range(j):
                fall *= (5 - k - q_)
            pw_ch, pw_LCrD, pw_CrN = 5 - j - k, 5 - k - j, k
            Pk = PC[f'E{k}'][0] if k >= 1 else sp.Poly(sp.Integer(1), *VS)
            dE = PC[f'E{k}'][1] if k >= 1 else 0
            cE = PC[f'E{k}'][2] if k >= 1 else sp.Integer(1)
            dens = (dE + pw_ch * PC['chat'][1]
                    + pw_LCrD * (PC['L'][1] + PC['CrD'][1])
                    + pw_CrN * PC['CrN'][1])
            consts = (cE * PC['chat'][2] ** pw_ch
                      * (PC['L'][2] * PC['CrD'][2]) ** pw_LCrD
                      * PC['CrN'][2] ** pw_CrN)
            terms.append((k, fall, Pk, pw_ch, pw_LCrD, pw_CrN, dens, consts))
        maxden = max(tt[6] for tt in terms)
        onep = sp.Poly(1 + sg ** 2, *VS)
        num = sp.Poly(sp.Integer(0), *VS)
        for (k, fall, Pk, pw_ch, pw_LCrD, pw_CrN, dens, consts) in terms:
            term = Pk
            for _ in range(pw_ch):
                term = term * PC['chat'][0]
            for _ in range(pw_LCrD):
                term = term * PC['L'][0] * PC['CrD'][0]
            for _ in range(pw_CrN):
                term = term * PC['CrN'][0]
            for _ in range(maxden - dens):
                term = term * onep
            num = num + term * sp.Rational((-1) ** k, 1) * fall / consts
        return num

    Ptg = assemble(0)
    print(f"     target: {len(Ptg.terms())} terms, degs "
          f"{[Ptg.degree(v) for v in VS]}", flush=True)

    # Fraction tensor of the target (exact; reused by P4 eval and P5 witness)
    Tg = {}
    for mon, cc in zip(Ptg.monoms(), Ptg.coeffs()):
        Tg[mon] = Fr(int(sp.numer(cc)), int(sp.denom(cc)))

    def tg_eval_mp(ptf, dps=60):
        """mpmath evaluation of the target tensor at (sg,th,t) Fractions."""
        with mp.workdps(dps):
            vals = [mp.mpf(q.numerator) / q.denominator for q in ptf]
            mx = [max(m[i] for m in Tg) for i in range(3)]
            pows = [[vals[i] ** k2 for k2 in range(mx[i] + 1)]
                    for i in range(3)]
            s_ = mp.mpf(0)
            for m, cf in Tg.items():
                s_ += (mp.mpf(cf.numerator) / cf.denominator
                       * pows[0][m[0]] * pows[1][m[1]] * pows[2][m[2]])
            return s_

    # NOTE: the target has huge inter-term cancellation near the degenerate
    # corners; double precision FLIPS signs there.  Evaluate at dps 60 and
    # compare against chi from the numeric eigenvalues (ambient dps 40).
    random.seed(23)
    agree = disagree = npts = 0
    for _ in range(140):
        sv = Fr(random.randint(1, 31), 32)
        tv = Fr(random.randint(1, 32), 16)
        thv = Fr(random.randint(1, 16), 16)
        pv = 0.8 * (1.0 - float(sv) ** 2)
        Dv = float(thv) * 0.8 * (1.0 - float(sv)) ** 2 / (1.0 + float(sv) ** 2)
        if Dv <= 1e-9 or Dv >= 0.7999:
            continue
        tvf = float(tv)
        cn = chi_numeric(pv, Dv, math.acos(1 - tvf), tvf)
        tgv = tg_eval_mp((sv, thv, tv))
        if tgv == 0 or abs(cn) < 1e-24:
            continue
        npts += 1
        if (cn > 0) == (tgv > 0):
            agree += 1
        else:
            disagree += 1
    rep(f"P4 numeric sign chain sign(target)=sign(chi), {npts} pts, "
        f"{disagree} flips", npts >= 80 and disagree == 0)

    # ---- P5: pocket witnesses (why the restriction is necessary) ----------
    print("-" * 78)
    print("P5  pocket witnesses: target < 0 at (sg=1/2, th=1, t=2) [pieceM")
    print("    range] and (sg=4/5, th=1, t=2) [pieceP range]  (exact)")
    mxa = max(m[0] for m in Tg)
    mxc = max(m[2] for m in Tg)
    pow2 = [Fr(2) ** c2 for c2 in range(mxc + 1)]
    ok5 = True
    for sgw, tag in ((Fr(1, 2), "sg=1/2 (p=3/5,  D=Dbar=4/25)"),
                     (Fr(4, 5), "sg=4/5 (p=36/125, D=Dbar=4/205)")):
        pws = [Fr(1)] * (mxa + 1)
        for a2 in range(1, mxa + 1):
            pws[a2] = pws[a2 - 1] * sgw
        wit = sum(cf * pws[m[0]] * pow2[m[2]] for m, cf in Tg.items())
        print(f"     target = {float(wit):.6e}  (exact rational; {tag})")
        ok5 &= (wit < 0)
    rep("P5 exact target < 0 at both pocket points (th=1 above each curve)",
        ok5)

    # ---- P6: strip trivial factors, integerize ----------------------------
    print("-" * 78)
    print("P6  strip (1-sg)^7 * th * t and integerize (exact divisions)")
    corePoly = Ptg
    mults = {}
    for f in ((1 - sg), th, t, sg, (1 + sg), (1 + sg ** 2)):
        fP = sp.Poly(f, *VS)
        while True:
            q_, r_ = sp.div(corePoly, fP)
            if r_.is_zero and not q_.is_zero:
                corePoly = q_
                mults[f] = mults.get(f, 0) + 1
            else:
                break
    ok6 = (mults.get(1 - sg, 0) == 7 and mults.get(th, 0) == 1
           and mults.get(t, 0) == 1 and len(mults) == 3)
    print(f"     multiplicities: {[(str(k), v) for k, v in mults.items()]}; "
          f"core {len(corePoly.terms())} terms, degs "
          f"{[corePoly.degree(v) for v in VS]}")
    rep("P6 target = (1-sg)^7 * th * t * core exactly (all factors > 0 on "
        "claim region)", ok6)
    dens = [sp.fraction(cc)[1] for cc in corePoly.coeffs()]
    lc = sp.ilcm(*[int(d_) for d_ in dens])
    Tcore = {}
    intOK = True
    for mon, cc in zip(corePoly.monoms(), corePoly.coeffs()):
        val = cc * lc
        if not (val.is_rational and sp.Integer(val) == val):
            intOK = False
            break
        Tcore[mon] = int(val)
    rep(f"P6b core * lcm({lc}) has integer coefficients", intOK)

    # ---- P6c: x-domain identities of the piece substitutions --------------
    print("-" * 78)
    print("P6c exact identities: curveM(3x/4) == (64-27x^2)/64;")
    print("    curveP((7+3x)/10) == (131-42x-9x^2)/200;")
    print("    1-sg^2 |P == 3(1-x)(17+3x)/100")
    idM = sp.expand(64 * CURVM.subs(sg, sp.Rational(3, 4) * X)
                    - (64 - 27 * X ** 2))
    idP = sp.expand(200 * CURVP.subs(sg, (7 + 3 * X) / 10)
                    - (131 - 42 * X - 9 * X ** 2))
    idP2 = sp.expand(100 * (1 - ((7 + 3 * X) / 10) ** 2)
                     - 3 * (1 - X) * (17 + 3 * X))
    rep("P6c all substitution identities hold exactly",
        idM == 0 and idP == 0 and idP2 == 0)

    DS = corePoly.degree(sg)
    DTH = corePoly.degree(th)
    DT = corePoly.degree(t)
    DX = DS + 2 * DTH
    degs3 = [DX, DTH, DT]
    P2l = [2 ** cc for cc in range(DT + 1)]

    def build_TU(piece):
        """Integer tensor of core after the piece substitution.
        pieceM: sg = 3x/4,      th = u(64-27x^2)/64,      scale 4^DS 64^DTH
        pieceP: sg = (7+3x)/10, th = u(131-42x-9x^2)/200, scale 10^DS 200^DTH
        (t = 2 tau both).  Exact Fraction + Schwartz-Zippel guards below."""
        if piece == 'M':
            CP = {0: {0: 1}}
            for b in range(1, DTH + 1):
                prev = CP[b - 1]
                cur = {}
                for jj, cf in prev.items():
                    cur[jj] = cur.get(jj, 0) + 64 * cf
                    cur[jj + 2] = cur.get(jj + 2, 0) - 27 * cf
                CP[b] = cur
            SP = {a: {a: 3 ** a} for a in range(DS + 1)}      # (3x)^a
            PS = [4 ** a for a in range(DS + 1)]
            PC2 = [64 ** b for b in range(DTH + 1)]
        else:
            CP = {0: {0: 1}}
            for b in range(1, DTH + 1):
                prev = CP[b - 1]
                cur = {}
                for jj, cf in prev.items():
                    cur[jj] = cur.get(jj, 0) + 131 * cf
                    cur[jj + 1] = cur.get(jj + 1, 0) - 42 * cf
                    cur[jj + 2] = cur.get(jj + 2, 0) - 9 * cf
                CP[b] = cur
            SP = {0: {0: 1}}
            for a in range(1, DS + 1):
                prev = SP[a - 1]
                cur = {}
                for jj, cf in prev.items():
                    cur[jj] = cur.get(jj, 0) + 7 * cf
                    cur[jj + 1] = cur.get(jj + 1, 0) + 3 * cf
                SP[a] = cur
            PS = [10 ** a for a in range(DS + 1)]
            PC2 = [200 ** b for b in range(DTH + 1)]
        TU = {}
        for (a, b, cc), cf in Tcore.items():
            base = cf * PS[DS - a] * P2l[cc] * PC2[DTH - b]
            for ja, pa in SP[a].items():
                for jc, pc in CP[b].items():
                    key = (ja + jc, b, cc)
                    TU[key] = TU.get(key, 0) + base * pa * pc
        return {k: v for k, v in TU.items() if v != 0}

    PIECE = {}
    for piece, scale_pair, subs_fr in (
            ('M', (4, 64), lambda xv: (Fr(3, 4) * xv,
                                       1 - Fr(3, 4) * (Fr(3, 4) * xv) ** 2)),
            ('P', (10, 200), lambda xv: (Fr(7, 10) + Fr(3, 10) * xv,
                                         Fr(9, 10)
                                         - (Fr(7, 10) + Fr(3, 10) * xv) ** 2
                                         / 2))):
        print("-" * 78)
        print(f"P7  piece {piece}: substitution -> integer tensor + exact "
              f"cross-checks")
        TU = build_TU(piece)
        PIECE[piece] = TU
        print(f"     integer tensor: {len(TU)} nonzero, degs "
              f"({DX},{DTH},{DT})")
        sA, sB = scale_pair
        SCALE = Fr(sA) ** DS * Fr(sB) ** DTH
        random.seed(11)
        ok7 = True
        for _ in range(12):
            xv = Fr(random.randint(0, 16), 16)
            uv = Fr(random.randint(0, 16), 16)
            tv = Fr(random.randint(0, 16), 16)
            sgv, curv = subs_fr(xv)
            thv = uv * curv
            ttv = 2 * tv
            lhs = sum(Fr(cf) * xv ** a2 * uv ** b2 * tv ** c2
                      for (a2, b2, c2), cf in TU.items())
            rhs = SCALE * sum(Fr(cf) * sgv ** a2 * thv ** b2 * ttv ** c2
                              for (a2, b2, c2), cf in Tcore.items())
            ok7 &= (lhs == rhs)
        rep(f"P7 piece {piece}: exact Fraction cross-check at 12 points", ok7)
        # P7b: conclusive polynomial-identity check by Schwartz-Zippel over
        # F_q, q = 2^61 - 1 (prime).  The difference polynomial has total
        # degree <= DX+DTH+DT, so a nonzero difference survives 20 uniform
        # F_q points with probability <= ((DX+DTH+DT)/q)^20 ~ 1e-370.
        qP = (1 << 61) - 1
        random.seed(101)
        ok7b = True
        if piece == 'M':
            def sub_q(xq):
                sgq = 3 * xq % qP * pow(4, qP - 2, qP) % qP
                thq_num = (64 - 27 * xq * xq) % qP
                thq_den = pow(64, qP - 2, qP)
                return sgq, thq_num * thq_den % qP
        else:
            def sub_q(xq):
                sgq = (7 + 3 * xq) % qP * pow(10, qP - 2, qP) % qP
                thq_num = (131 - 42 * xq - 9 * xq * xq) % qP
                thq_den = pow(200, qP - 2, qP)
                return sgq, thq_num * thq_den % qP
        scale_q = pow(sA, DS, qP) * pow(sB, DTH, qP) % qP
        for _ in range(20):
            xq = random.randrange(qP)
            uq = random.randrange(qP)
            tq = random.randrange(qP)
            lhs = 0
            pxa = {a2: pow(xq, a2, qP) for a2 in set(k[0] for k in TU)}
            pub = {b2: pow(uq, b2, qP) for b2 in set(k[1] for k in TU)}
            ptc = {c2: pow(tq, c2, qP) for c2 in set(k[2] for k in TU)}
            for (a2, b2, c2), cf in TU.items():
                lhs = (lhs + cf * pxa[a2] % qP * pub[b2] % qP * ptc[c2]) % qP
            sgq, curq = sub_q(xq)
            thq = uq * curq % qP
            ttq = 2 * tq % qP
            psa = {a2: pow(sgq, a2, qP) for a2 in set(k[0] for k in Tcore)}
            psb = {b2: pow(thq, b2, qP) for b2 in set(k[1] for k in Tcore)}
            psc = {c2: pow(ttq, c2, qP) for c2 in set(k[2] for k in Tcore)}
            rhs = 0
            for (a2, b2, c2), cf in Tcore.items():
                rhs = (rhs + cf * psa[a2] % qP * psb[b2] % qP * psc[c2]) % qP
            ok7b &= (lhs % qP == rhs * scale_q % qP)
        rep(f"P7b piece {piece}: Schwartz-Zippel identity mod 2^61-1, "
            f"20 points", ok7b)

    # ---- P8/P9/P10 towers per piece ----------------------------------------
    for piece in ('M', 'P'):
        TU = PIECE[piece]
        closed_end = "x=1 (sg=3/4)" if piece == 'M' else "x=0 (sg=7/10)"
        open_end = "x=0 (p=4/5)" if piece == 'M' else "x=1 (p=0)"
        print("-" * 78)
        print(f"P8  piece {piece}: root box [0,1]^3, scaled Bernstein "
              f"coefficients")
        BB = dict(TU)
        for ax in (2, 1, 0):
            BB = bern_axis(BB, degs3, ax)
        total = (DX + 1) * (DTH + 1) * (DT + 1)
        negs = sum(1 for v in BB.values() if v < 0)
        nz = total - len(BB)
        print(f"     nonzero {len(BB)}, exact zeros {nz}, negatives {negs}")
        rep(f"P8 piece {piece}: NO negative Bernstein coeffs (core_u > 0 "
            f"on OPEN box)", negs == 0 and len(BB) > 0)

        print("-" * 78)
        print(f"P9  piece {piece}: faces {closed_end}, u=1, tau=1")
        print(f"    [{open_end} is the excluded degenerate limit of the "
              f"Gray family]")
        ok9 = True
        if piece == 'M':
            faces = ((0, "x=1", collapse_at1), (1, "u=1", collapse_at1),
                     (2, "tau=1", collapse_at1))
        else:
            faces = ((0, "x=0", restrict_at0), (1, "u=1", collapse_at1),
                     (2, "tau=1", collapse_at1))
        for ax, nm, fn in faces:
            Tf, degf = fn(TU, ax, degs3)
            Bf = dict(Tf)
            for ax2 in range(len(degf) - 1, -1, -1):
                Bf = bern_axis(Bf, degf, ax2)
            negf = sum(1 for v in Bf.values() if v < 0)
            okf = (negf == 0 and len(Bf) > 0)
            ok9 &= okf
            rep(f"P9 piece {piece} face {nm}: no negatives ({negf}), "
                f"nonzero ({len(Bf)}) => > 0 on open face", okf)

        print("-" * 78)
        print(f"P10 piece {piece}: edges (1-var no-negatives + nonzero) "
              f"and the corner")
        ok10 = True
        # the two x-end edges: (x-end, u=1) in tau and (x-end, tau=1) in u
        for keep, nm in ((2, "u=1) in tau"), (1, "tau=1) in u")):
            if piece == 'M':
                # x=1 end: collapse x (sum over all x-powers) -- value at x=1
                co = {}
                for mon, cf in TU.items():
                    co[mon[keep]] = co.get(mon[keep], 0) + cf
                enm = f"(x=1,{nm}"
            else:
                # x=0 end: restrict to x-power 0
                co = {}
                for mon, cf in TU.items():
                    if mon[0] != 0:
                        continue
                    co[mon[keep]] = co.get(mon[keep], 0) + cf
                enm = f"(x=0,{nm}"
            # careful: for keep=2 we must also fix the OTHER var at 1
            # (collapse) -- the dicts above already sum over both remaining
            # axes except `keep`, which is exactly evaluation at 1 there.
            d_ = degs3[keep]
            bb = [0] * (d_ + 1)
            for a2, cf in co.items():
                if cf == 0:
                    continue
                for b2 in range(a2, d_ + 1):
                    bb[b2] += comb(d_ - a2, b2 - a2) * cf
            negb = sum(1 for v in bb if v < 0)
            oke = (negb == 0 and any(v > 0 for v in bb))
            ok10 &= oke
            rep(f"P10 piece {piece} edge {enm}: no negatives => > 0 on "
                f"open edge", oke)
        # the (u=1,tau=1) edge in x
        co = {}
        for mon, cf in TU.items():
            co[mon[0]] = co.get(mon[0], 0) + cf
        bb = [0] * (DX + 1)
        for a2, cf in co.items():
            if cf == 0:
                continue
            for b2 in range(a2, DX + 1):
                bb[b2] += comb(DX - a2, b2 - a2) * cf
        negb = sum(1 for v in bb if v < 0)
        oke = (negb == 0 and any(v > 0 for v in bb))
        ok10 &= oke
        print(f"     (u=1,tau=1) end values: at x=0 -> "
              f"{float(bb[0]):.4g}, at x=1 -> {float(bb[-1]):.4g}")
        rep(f"P10 piece {piece} edge (u=1,tau=1): no negatives => > 0 on "
            f"open edge x in (0,1)", oke)
        if piece == 'M':
            corner = sum(cf for cf in TU.values())            # (1,1,1)
            cnm = "(1,1,1) [sg=3/4, th=curveM, t=2]"
        else:
            corner = sum(cf for (a2, b2, c2), cf in TU.items() if a2 == 0)
            cnm = "(0,1,1) [sg=7/10, th=curveP, t=2]"
        ok10 &= rep(f"P10d piece {piece} corner {cnm}: exact value > 0",
                    corner > 0)
        if piece == 'M':
            print("     [union: open box + faces x=1,u=1,tau=1 + edges")
            print("      (x=1,u=1),(x=1,tau=1),(u=1,tau=1) + corner (1,1,1)")
            print("      = (0,1]^3 = claim region of (A-M); x=0 (p=4/5),")
            print("      u=0 (D=0), tau=0 (s=0) are excluded open limits]")
        else:
            print("     [union: open box + faces x=0,u=1,tau=1 + edges")
            print("      (x=0,u=1),(x=0,tau=1),(u=1,tau=1) + corner (0,1,1)")
            print("      = [0,1) x (0,1] x (0,1] = claim region of (A-P);")
            print("      x=1 (p=0), u=0 (D=0), tau=0 (s=0) are excluded]")

    # ---- P11: certificate (B), exact Sturm, per piece ----------------------
    print("-" * 78)
    print("P11 certificate (B): D_c < curve * Dbar on each piece")
    N5, Dn5 = disc_pieces(AVAL)
    rep("P11a disc denominator == 2^18 (5D-4)^12  (> 0 for D < 4/5)",
        sp.expand(Dn5 - 2 ** 18 * (5 * D - 4) ** 12) == 0)
    dmin = min(m[0] for m in sp.Poly(N5, D).monoms())
    C5 = sp.expand(sp.cancel(N5 / D ** dmin))
    C50 = sp.expand(C5.subs(D, 0))
    id0 = sp.expand(C50 - 2 ** 18 * p ** 4 * (1 - p) ** 2 * (5 * p - 4) ** 4
                    * (p ** 2 - p + 1))
    rep("P11b N = D^4 C, C(0,p) == 2^18 p^4 (1-p)^2 (5p-4)^4 (p^2-p+1)",
        dmin == 4 and id0 == 0)
    rep("P11b' p^2-p+1 > 0 always (discriminant 1-4 = -3 < 0) "
        "=> disc > 0 as D->0+ on (0,4/5)",
        sp.discriminant(p ** 2 - p + 1, p) == -3)
    for piece, curve, (lo, hi), expstrip, expdeg, expconst in (
            ('M', CURVM, (sp.Integer(0), sp.Rational(3, 4)),
             {'sigma': 13, '1 - sigma': 12, 'sigma + 1': 2}, 29, 5 ** 12),
            ('P', CURVP, (sp.Rational(7, 10), sp.Integer(1)),
             {'sigma': 8, '1 - sigma': 13, 'sigma + 1': 2}, 33, 5 ** 22)):
        Xcurve = curve * DBAR
        expr = N5.subs({D: Xcurve, p: PSUB})
        num, den = sp.fraction(sp.cancel(sp.together(expr)))
        denf = sp.factor_list(den)
        okden = all(b.is_number or sp.simplify(b - (1 + sg ** 2)) == 0
                    for b, e in denf[1]) and denf[0] > 0
        rep(f"P11c piece {piece}: cleared denominator = {expconst} "
            f"(1+sg^2)^6 (positive)",
            okden and sp.simplify(den - expconst * (1 + sg ** 2) ** 6) == 0)
        Rp = sp.Poly(sp.expand(num), sg)
        strips = {}
        for f in (sg, 1 - sg, 1 + sg):
            fP = sp.Poly(f, sg)
            while True:
                q_, r_ = sp.div(Rp, fP)
                if r_.is_zero and not q_.is_zero:
                    Rp = q_
                    strips[f] = strips.get(f, 0) + 1
                else:
                    break
        okstrip = all(strips.get(sp.sympify(k) if isinstance(k, str) else k,
                                 0) == v
                      for k, v in {sg: expstrip['sigma'],
                                   1 - sg: expstrip['1 - sigma'],
                                   1 + sg: expstrip['sigma + 1']}.items())
        nroots = Rp.count_roots(lo, hi)
        mid = (lo + hi) / 2
        vl, vm, vh = Rp.eval(lo), Rp.eval(mid), Rp.eval(hi)
        print(f"     piece {piece}: numerator = sg^{expstrip['sigma']} "
              f"(1-sg)^{expstrip['1 - sigma']} (1+sg)^{expstrip['sigma + 1']}"
              f" * R{Rp.degree()}; Sturm roots in [{lo},{hi}]: {nroots}; "
              f"signs at (lo,mid,hi): {sp.sign(vl)},{sp.sign(vm)},"
              f"{sp.sign(vh)}")
        rep(f"P11d piece {piece}: strip multiplicities exact, residual "
            f"deg {expdeg}", okstrip and Rp.degree() == expdeg)
        rep(f"P11e piece {piece}: residual has NO roots in the piece "
            f"interval and is negative (Sturm)",
            nroots == 0 and vl < 0 and vm < 0 and vh < 0)
        rep(f"P11f piece {piece} => disc(curve*Dbar) < 0 on the piece; "
            f"with P11a/b/b' and first-root logic, D_c < curve*Dbar",
            okstrip and nroots == 0 and vm < 0)
    print("     [safe direction only: D_c = first sign change of disc in D")
    print("      (committed convention, lemma_7_34_dbar_upper_bound (i)-(iii));")
    print("      disc > 0 at D->0+ and disc < 0 on each curve place that")
    print("      first sign change strictly below each curve;  strip factors")
    print("      sg^a (1-sg)^b (1+sg)^c are strictly positive on the OPEN")
    print("      sg-interval of each piece -- sg=0 (p=pmax) and sg=1 (p=0)")
    print("      are excluded limits of the Gray family]")

    # ---- P12: exact anchors ON each curve ----------------------------------
    print("-" * 78)
    print("P12 exact anchors ON each curve at w = -1 (t = 2):")
    print("    pieceM: sg=1/2, th=curveM(1/2)=13/16, p=3/5,    D=13/100")
    print("    pieceP: sg=4/5, th=curveP(4/5)=29/50, p=36/125, D=58/5125")
    ok12 = True
    for piece, sga, xanch in (('M', sp.Rational(1, 2), Fr(2, 3)),
                              ('P', sp.Rational(4, 5), Fr(1, 3))):
        curve = CURVM if piece == 'M' else CURVP
        tha = curve.subs(sg, sga)
        pv = sp.Rational(4, 5) * (1 - sga ** 2)
        Dbv = sp.Rational(4, 5) * (1 - sga) ** 2 / (1 + sga ** 2)
        Dv = tha * Dbv
        ethv = Dv / ((A - 1) * (1 - Dv))

        def etaf_w(Wv):
            E_ = ethv * Wv
            return ((A - 1) * Dv * E_ + Dv - (A - 1) * E_) / \
                ((A - 1) * Dv * E_ + Dv - (A - 1))

        # NOTE (A=4 lesson): do NOT wrap these exact Rationals in
        # sp.nsimplify at float precision -- all quantities are exact
        # Rationals already; the P12a gate below asserts that.
        ew = etaf_w(-1)
        Tw = [[(1 - pv) if i == j else pv / (A - 1) for j in range(AVAL)]
              for i in range(AVAL)]
        Qw = sp.zeros(5, 5)
        for a_ in range(5):
            x_, xp, xq = STATES[REPS[a_]]
            for (y, yp, yq) in STATES:
                term = Tw[xp][yp] * Tw[xq][yq] / Tw[x_][y]
                if yp != y: term *= ew
                if yq != y: term *= ew
                Qw[a_, _pat((y, yp, yq))] += term

        def Cfw(Wv):
            E_ = ethv * Wv
            return ((A - 1) * Dv * E_ + Dv - (A - 1)) / (A * Dv - (A - 1))

        Crw = Cfw(-1) * Cfw(-1) / (Cfw(1) ** 2)
        Lw = 1 - sp.Rational(3, 2) * Dv * (1 - Dv) * 2
        chiw = ((Lw / Crw) * sp.eye(5) - Qw).det(method='berkowitz')
        ok12 &= rep(f"P12a piece {piece}: all anchor quantities exact "
                    f"Rationals", ew.is_Rational and Crw.is_Rational
                    and chiw.is_Rational)
        TU = PIECE[piece]
        mxa2 = max(m[0] for m in TU)
        pw = [Fr(1)] * (mxa2 + 1)
        for a2 in range(1, mxa2 + 1):
            pw[a2] = pw[a2 - 1] * xanch
        vc = sum(Fr(cf) * pw[a2] for (a2, b2, c2), cf in TU.items())
        print(f"     piece {piece}: chi(Lambda) = {float(chiw):.6e} (exact "
              f"rational), core_u({xanch},1,1) sign "
              f"{'+' if vc > 0 else '-'}")
        ok12 &= rep(f"P12 piece {piece}: exact chi > 0 at the anchor and "
                    f"sign matches core_u", chiw > 0 and vc > 0)

    # ---- P13: numeric separation scans (consistency) ------------------------
    print("-" * 78)
    print("P13 numeric separation per piece: D_c(bisect) < curve*Dbar and")
    print("    chi > 0 on the curve at t=2, 25 grid points each (mpmath")
    print("    consistency scan)")
    for piece, curvef, lo_, hi_ in (
            ('M', lambda s: 1.0 - 0.75 * s * s, 0.0, 0.75),
            ('P', lambda s: (9.0 - 5.0 * s * s) / 10.0, 0.70, 1.0)):
        ok13 = True
        for k in range(1, 26):
            sgv = lo_ + (hi_ - lo_) * k / 26.0
            pv2 = 0.8 * (1.0 - sgv * sgv)
            Db = 0.8 * (1.0 - sgv) ** 2 / (1.0 + sgv * sgv)
            dth = curvef(sgv)
            dc = float(Dc5_numeric(pv2))
            ok13 &= (dc / Db < dth)
            chc = chi_numeric(pv2, dth * Db, math.pi, 2.0)
            ok13 &= (chc > 0)
        rep(f"P13 piece {piece}: curve strictly above theta_c and chi > 0 "
            f"on the curve at t=2", ok13)

    # ---- P14: EXACT SPLICE GATE (the union logic) ---------------------------
    print("-" * 78)
    print("P14 EXACT SPLICE: overlap [7/10, 3/4] nonempty; on it")
    print("    min(curveM, curveP) = curveM >= the Gray edge (by (B-M))")
    ok14a = sp.Rational(7, 10) < sp.Rational(3, 4)
    dfe = sp.expand(CURVP - CURVM - (5 * sg ** 2 - 2) / 20)
    ok14b = (dfe == 0)
    vM7, vP7 = CURVM.subs(sg, sp.Rational(7, 10)), \
        CURVP.subs(sg, sp.Rational(7, 10))
    vM34, vP34 = CURVM.subs(sg, sp.Rational(3, 4)), \
        CURVP.subs(sg, sp.Rational(3, 4))
    ok14c = (vM7 == sp.Rational(253, 400) and vP7 == sp.Rational(131, 200)
             and vM7 < vP7 and vM34 == sp.Rational(37, 64)
             and vP34 == sp.Rational(99, 160) and vM34 < vP34)
    # curveP - curveM = (5sg^2-2)/20 is increasing in sg^2 and already > 0
    # at sg = 7/10 (5*(49/100) = 49/20 > 2), hence > 0 on ALL of the overlap.
    ok14d = (5 * sp.Rational(49, 100) - 2 > 0)
    rep("P14a 7/10 < 3/4: sg-coverage (0,3/4] u [7/10,1) = (0,1) seamless",
        ok14a)
    rep("P14b curveP - curveM == (5 sg^2 - 2)/20 exactly", ok14b)
    rep("P14c curveM(7/10)=253/400 < curveP(7/10)=262/400; "
        "curveM(3/4)=37/64 < curveP(3/4)=99/160", ok14c)
    rep("P14d 5 sg^2 - 2 >= 9/20 > 0 on [7/10,3/4] => min(curveM,curveP) "
        "= curveM on the whole overlap", ok14d)
    print("     [union logic, exactly as the A=3 smallp splice: on the")
    print("      overlap the pieceP slab NESTS OVER the pieceM slab")
    print("      (curveP > curveM), so the union region {th <= curveM,")
    print("      sg <= 3/4} u {th <= curveP, sg >= 7/10} contains Gray for")
    print("      ALL p in (0,4/5): each sg is covered by at least one piece")
    print("      whose own (B)-certificate places D_c strictly below its")
    print("      curve, and chi > 0 holds on each slab by its own (A)-")
    print("      certificate.  No discrete-p residual, no gap.]")

    # ---- P15: pinch consistency at sg -> 1 ---------------------------------
    print("-" * 78)
    print("P15 pinch consistency at sg = 0.995 (numeric): theta_c < curveP")
    print("    < theta_x  (om = 1-sg^2 units: th_c ~ 2/5 + 0.2403 om,")
    print("    curveP = 2/5 + om/2, th_x ~ 2/5 + 0.7204 om -- scout 8)")
    sgv = 0.995
    pv3 = 0.8 * (1.0 - sgv * sgv)
    Db3 = 0.8 * (1.0 - sgv) ** 2 / (1.0 + sgv * sgv)
    dth3 = (9.0 - 5.0 * sgv * sgv) / 10.0
    thc3 = float(Dc5_numeric(pv3)) / Db3
    lo, hi = dth3, 1.0
    okpocket = (chi_numeric(pv3, lo * Db3, math.pi, 2.0) > 0
                and chi_numeric(pv3, hi * Db3, math.pi, 2.0) < 0)
    if okpocket:
        for _ in range(40):
            midq = 0.5 * (lo + hi)
            if chi_numeric(pv3, midq * Db3, math.pi, 2.0) > 0:
                lo = midq
            else:
                hi = midq
    thx3 = 0.5 * (lo + hi)
    print(f"     theta_c = {thc3:.7f}  curveP = {dth3:.7f}  "
          f"theta_x = {thx3:.7f}")
    rep("P15 theta_c < curveP < theta_x at the pinch point",
        okpocket and thc3 < dth3 < thx3)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("CLAIM CERTIFIED (A-M): chi(Lambda) > 0 for sg in (0, 3/4],")
    print("th in (0, 1-(3/4)sg^2], t in (0,2]      [exact Bernstein tower]")
    print("CLAIM CERTIFIED (A-P): chi(Lambda) > 0 for sg in [7/10, 1),")
    print("th in (0, (9-5sg^2)/10], t in (0,2]     [exact Bernstein tower]")
    print("CLAIM CERTIFIED (B-M/B-P): D_c(p) < curve * Dbar(p) on each piece")
    print("[exact Sturm; D_c = first sign change of the alternating-product")
    print(" reduced-cubic discriminant, committed lemma_7_34_dbar_upper_bound")
    print(" convention; safe direction only]")
    print("SPLICE (P14, exact): 7/10 < 3/4 and curveP > curveM on the overlap")
    print("=> the A=5 Gray region for ALL p in (0, 4/5) is strictly inside")
    print("   the union of the two certified chi > 0 slabs -- TWO pieces,")
    print("   continuum grade end to end, no discrete-p residual (the A=3")
    print("   splice pattern; a single curve is provably too rigid at A=5).")
    return PASS


if __name__ == "__main__":
    main()

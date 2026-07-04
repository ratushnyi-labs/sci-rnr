#!/usr/bin/env python3
r"""
lemma_7_34_a4_cert5_on_gray.py
============================================================================
CERT5-ON-GRAY (BUG-009-D / Route B, A=4): exact certificate that the k=0
charpoly sign  chi(Lambda) > 0  holds on a GRAY-CONTAINING subdomain of the
rationalized superdomain -- the certificate that is FALSE on the full
superdomain (real Perron pocket, exact witness chi(Lambda) < 0 at
sg = 3/4, th = 1, t = 2, i.e. p = 21/64, D = Dbar = 3/100; check P5).

STATEMENT (all exact, A = 4).  Rationalized domain
    p = (3/4)(1 - sg^2),   Dbar = (3/4)(1-sg)^2/(1+sg^2),   D = th*Dbar,
    t = 1 - cos s in (0, 2],  sg in (0, 1),  th in (0, 1].
Separating curve (rational, degree 2 -- a SINGLE curve for the whole range):
    thetatilde4(sg) = 1 - sg^2/2.
 (A)  chi(Lambda) > 0  for all  sg in (0, 1),  th in (0, thetatilde4(sg)],
      t in (0, 2].
 (B)  D_c(p) < thetatilde4(sg) * Dbar(p)  for all  sg in (0, 1), i.e. for
      ALL p in (0, 3/4), where D_c is the A=4 Gray threshold = first sign
      change in D of the reduced-cubic discriminant disc(D,p) of the 3x3
      symmetry-reduced alternating product (the committed threshold object
      of lemma_7_34_dbar_upper_bound.py, steps (i)-(iii)).
 =>   the Gray region {0 < D <= D_c(p)} x {angles s in (0, pi]} for ALL
      p in (0, 3/4) lies strictly inside the region where (A) certifies
      chi(Lambda) > 0.

SCOPE (vs A = 3).  At A=3 the committed cert5-on-gray covers only
sg in (0, 25/32] and needs the separate Schur small/moderate-p certificate
(discrete-p grade) below p = 133/512.  At A=4 the SINGLE curve 1 - sg^2/2
threads the corridor (theta_c, theta_x) over the ENTIRE range sg in (0,1):
no splice, no discrete-p residual -- whole-Gray coverage at full continuum
grade in one certificate.  (Both sg-endpoints are degenerate corners --
sg=0 is p = pmax = 3/4 with D -> 3/4, sg=1 is p -> 0 -- and both are
excluded limits of the Gray family; the certificate is OPEN at both ends.)

SIGN LINKAGE (exact, no numerics needed).  With the t-domain assembly of
the committed A=3 pipeline (lemma_7_34_midband_chi_certificates.py) run at
A = 4:
    target := assemble(0) = (1+sg^2)^K * chat_t^5 * CrN_t^5 * chi(Lambda)
on the rationalized domain, where (check P3, exact expansions)
    chat_t = p(p-1) * B(D,t),      CrN_t = -B(D,t),
    B(D,t) = 2 D^2 t (1-D)(3-D) + (4D-3)^2  >  0   for 0 < D < 3/4
    (sum of nonnegative terms, second strictly positive),
hence  chat_t^5 CrN_t^5 = [p(1-p)]^5 B^10 > 0  on the claim region and
    sign(target) = sign(chi(Lambda))            (P4 = 140-pt numeric guard).
The identity is the charpoly expansion chi = sum (-1)^k e_k Lam^{5-k} with
e_k = Etld_k / chat^k, Lam = L*CrD_t/CrN_t, L = 1 - (3/2)D(1-D)t (target
rate R_4 >= 3/2, same L as A=3), valid once the grading gates P1 hold
(E_k symmetric at (-w)^k, chat at (-w)^1, mN = mD).

CERTIFICATE (A) STRUCTURE.  target = (1-sg)^7 * th * t * core  (exact
division, P6; all three stripped factors > 0 on the claim region).
Substitute
    sg = x,   th = u * (2 - x^2)/2,   t = 2 tau
and clear denominators (integer tensor, exact cross-check P7 + independent
Schwartz-Zippel identity check P7b).  Then on [0,1]^3 in (x, u, tau):
    P8   ALL scaled Bernstein coefficients >= 0, not all 0
         =>  core_u > 0 on the OPEN box (0,1)^3;
    P9   faces u=1, tau=1: 2-var Bernstein no-negatives + nonzero
         =>  core_u > 0 on the two open faces;
    P10  edge (u=1,tau=1): 1-var no-negatives + nonzero => positive on the
         open edge x in (0,1).
    Union = (0,1) x (0,1] x (0,1] exactly = the claim region of (A).
    (x=0 <=> p = 3/4 and x=1 <=> p = 0 are excluded limits of the Gray
    family; no x-faces are needed, unlike A=3 where the sg-range was
    CLOSED at 25/32.)  A zero-coefficient count is reported: zeros are
    allowed by the open-box+closure logic; the edge polynomial vanishes
    EXACTLY at both ends x=0 and x=1 (the two degenerate corners), which
    is why no closed-box certificate can exist and the faces/edge tower
    is used instead.

CERTIFICATE (B) (exact Sturm).  disc(D,p) = N(D,p) / (3^10 (3-4D)^12)
(identity P11a; the denominator is > 0 for D < 3/4, so sign(disc) =
sign(N) there).  N = D^4 * C(D,p) with
    C(0,p) = 3^9 p^4 (1-p)^2 (3-4p)^4 (3p^2-2p+3)     (identity P11b),
and 3p^2-2p+3 > 0 always (discriminant -32 < 0), so disc > 0 as D -> 0+
for every p in (0, 3/4) \ {roots at the endpoints only}.  On the curve
D = thetatilde4*Dbar, p = (3/4)(1-sg^2):  after clearing the positive
denominator 2^32 (1+sg^2)^6 and stripping sg^13 (1-sg)^13 (1+sg)^2 (all
positive on (0,1)), the residual degree-28 integer polynomial R28 has NO
roots in [0, 1] (Sturm) and is negative there (P11c).  A continuous
function of D positive as D -> 0+ and negative at D = thetatilde4*Dbar
has its FIRST sign change strictly below the curve; that first sign
change IS the A=4 threshold object D_c (committed convention:
lemma_7_34_dbar_upper_bound.py steps (i)-(iii), threshold mechanism
probe_7_34_beyond_gray_structure B1).  Only the SAFE direction is used:
D_c <= (any D where disc has already turned negative with disc > 0
before it), hence D_c < thetatilde4 * Dbar.  The identification of the
first disc sign change with the all-word Gray threshold inherits the
grade of the committed threshold mechanism, exactly as in the A=3
cert5-on-gray usage of the Lemma 7.34j quartic.

HONEST NOTES.
 *  The Perron pocket is REAL at A=4 (P5 exact witness at sg=3/4, th=1:
    numerically theta_x(sg=3/4, t=2) ~ 0.77 < 1) -- the restriction to
    th <= thetatilde4 is necessary, not an artifact.
 *  Near sg -> 1 the corridor pinches: at sg = 0.995 the admissible
    window for a in curves 1 - a sg^2 is (0.49747, 0.50251) (scout
    geography), and the curve a = 1/2 rides through the pinch; P14 checks
    theta_c < thetatilde4 < theta_x numerically at that point.  Near
    sg -> 0 the binding constraint is theta_c < thetatilde4, which holds
    with margin (theta_c = 1 - sg^2(1+o(1)) vs thetatilde4 = 1 - sg^2/2).
 *  chi(Lambda) > 0 is the k=0 member of the Budan-Fourier tower; the
    j >= 1 positive-side and all negative-side psi^(j) certificates hold
    on the FULL superdomain box at A=4 (scout run, root-box grade) and
    are packaged separately; interior-angle realness remains a separate
    open item (unchanged by this script).
 *  Numeric guards run at mpmath dps 40-60: double precision FLIPS signs
    of these polynomials near the degenerate corners (A=3 lesson,
    reproduced at A=4 in development; independent 400-pt hostile sweep
    found claim-region chi values down to ~5e-33).
 *  sp.nsimplify must NOT touch exact Rationals in the anchor P12: at
    ambient dps 40 it returns an only-approximately-equal radical-power
    form for Crw (caught in development; P12a gates against recurrence).

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.
Runtime: ~2-3 min (measured 102 s).  Machinery: committed A=3 pipeline
(scripts/verify/lemma_7_34_midband_chi_certificates.py build_Q/build_Cr/
Chebyshev t-domain assembly; lemma_7_34_midband_cert5_on_gray.py Bernstein
tower) run at AVAL = 4, and lemma_7_34_dbar_upper_bound.py disc_pieces.
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
# symbols and the rationalized domain (A = 4)
# ---------------------------------------------------------------------------
p, D, w, t = sp.symbols('p D w t')
sg, th, U, X = sp.symbols('sigma theta u x')
AVAL = 4
A = sp.Integer(AVAL)
PSUB = sp.Rational(3, 4) * (1 - sg ** 2)
DBAR = sp.Rational(3, 4) * (1 - sg) ** 2 / (1 + sg ** 2)
DSUB = th * DBAR
VS = (sg, th, t)
THT4 = 1 - sg ** 2 / 2                       # the single separating curve

# ---------------------------------------------------------------------------
# replica pattern-quotient machinery (committed pipeline, at A=4)
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
    """chi(Lambda) from the numeric eigenvalues of Q (A=4)."""
    Qn, Crn = Q_numeric(pv, Dv, sv)
    ev = mp.eig(Qn)[0]
    Lm = (1 - 1.5 * Dv * (1 - Dv) * tvf) / Crn
    out = mp.mpc(1)
    for lam in ev:
        out *= (Lm - lam)
    return float(mp.re(out))


# ---------------------------------------------------------------------------
# the A=4 threshold object: reduced-cubic discriminant of the alternating
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


def Dc4_numeric(pv, prec_bits=60):
    """A=4 Gray threshold via alternating-product collision bisection
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
# main
# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("CERT5-ON-GRAY (A=4): chi(Lambda) > 0 on th <= thetatilde4(sg) = 1-sg^2/2,")
    print("sg in (0,1) (ALL p in (0,3/4)), all angles; + D_c < thetatilde4*Dbar")
    print("(exact Sturm).  A = 4.  Exact rational/integer arithmetic end to end.")
    print("=" * 78)

    # ---- build the t-domain pieces (committed assembly, A=4) --------------
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
    Bpoly = 2 * D ** 2 * t * (1 - D) * (3 - D) + (4 * D - 3) ** 2
    ok3 = (sp.expand(chat - p * (p - 1) * Bpoly) == 0)
    ok3 &= (sp.expand(CrNt + Bpoly) == 0)
    ok3 &= (sp.expand(CrDt + (4 * D - 3) ** 2) == 0)
    rep("P3 chat_t = p(p-1)B, CrN_t = -B, B = 2D^2t(1-D)(3-D)+(4D-3)^2", ok3)
    rep("P3b => chat_t^5 CrN_t^5 = [p(1-p)]^5 B^10 > 0 for 0<p<1, D<3/4",
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
        pv = 0.75 * (1.0 - float(sv) ** 2)
        Dv = float(thv) * 0.75 * (1.0 - float(sv)) ** 2 / (1.0 + float(sv) ** 2)
        if Dv <= 1e-9 or Dv >= 0.7499:
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

    # ---- P5: pocket witness (why the restriction is necessary) ------------
    print("-" * 78)
    print("P5  pocket witness: target < 0 at (sg=3/4, th=1, t=2) [exact]")
    mxa = max(m[0] for m in Tg)
    pow34 = [Fr(1)] * (mxa + 1)
    for a2 in range(1, mxa + 1):
        pow34[a2] = pow34[a2 - 1] * Fr(3, 4)
    mxc = max(m[2] for m in Tg)
    pow2 = [Fr(2) ** c2 for c2 in range(mxc + 1)]
    wit = sum(cf * pow34[m[0]] * pow2[m[2]] for m, cf in Tg.items())
    print(f"     target = {float(wit):.6e}  (exact rational; "
          f"p = 21/64, D = Dbar = 3/100)")
    rep("P5 exact target < 0 at the pocket point (th=1 above the curve)",
        wit < 0)

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

    # ---- P7: u-substitution to (x,u,tau) on [0,1]^3, exact cross-check ----
    print("-" * 78)
    print("P7  substitution sg = x, th = u(2-x^2)/2, t = 2 tau")
    DS = corePoly.degree(sg)
    DTH = corePoly.degree(th)
    DT = corePoly.degree(t)
    DX = DS + 2 * DTH
    # (2 - x^2)^b expansion table (integer)
    POW = {0: {0: 1}}
    for b in range(1, DTH + 1):
        prev = POW[b - 1]
        cur = {}
        for jj, cf in prev.items():
            cur[jj] = cur.get(jj, 0) + 2 * cf
            cur[jj + 1] = cur.get(jj + 1, 0) - cf
        POW[b] = cur
    P2l = [2 ** cc for cc in range(DT + 1)]
    P2th = [2 ** b for b in range(DTH + 1)]
    # core(sg,th,t) * 2^DTH with the substitution: th^b contributes the
    # denominator 2^b ONLY (sg = x carries no scaling at A=4, unlike the
    # A=3 25/32 range compression), so the compensating power is 2^(DTH-b).
    TU = {}
    for (a, b, cc), cf in Tcore.items():
        base = cf * P2l[cc] * P2th[DTH - b]
        for jj, pc in POW[b].items():
            key = (a + 2 * jj, b, cc)
            TU[key] = TU.get(key, 0) + base * pc
    TU = {k: v for k, v in TU.items() if v != 0}
    print(f"     integer tensor: {len(TU)} nonzero, degs ({DX},{DTH},{DT})")
    SCALE = Fr(2) ** DTH
    random.seed(11)
    ok7 = True
    for _ in range(12):
        xv = Fr(random.randint(0, 16), 16)
        uv = Fr(random.randint(0, 16), 16)
        tv = Fr(random.randint(0, 16), 16)
        sgv = xv
        thv = uv * (2 - xv ** 2) / 2
        ttv = 2 * tv
        lhs = sum(Fr(cf) * xv ** a2 * uv ** b2 * tv ** c2
                  for (a2, b2, c2), cf in TU.items())
        rhs = SCALE * sum(Fr(cf) * sgv ** a2 * thv ** b2 * ttv ** c2
                          for (a2, b2, c2), cf in Tcore.items())
        ok7 &= (lhs == rhs)
    rep("P7 exact Fraction cross-check core_u == scaled core at 12 points", ok7)
    # P7b: conclusive polynomial-identity check by Schwartz-Zippel over F_q,
    # q = 2^61 - 1 (prime).  The difference polynomial has total degree
    # <= DX+DTH+DT, so a nonzero difference survives 20 uniform F_q points
    # with probability <= ((DX+DTH+DT)/q)^20 ~ 1e-370.
    qP = (1 << 61) - 1
    inv2 = pow(2, qP - 2, qP)
    scale_q = pow(2, DTH, qP)
    random.seed(101)
    ok7b = True
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
        sgq = xq
        thq = uq * ((2 - xq * xq) % qP) % qP * inv2 % qP
        ttq = 2 * tq % qP
        psa = {a2: pow(sgq, a2, qP) for a2 in set(k[0] for k in Tcore)}
        psb = {b2: pow(thq, b2, qP) for b2 in set(k[1] for k in Tcore)}
        psc = {c2: pow(ttq, c2, qP) for c2 in set(k[2] for k in Tcore)}
        rhs = 0
        for (a2, b2, c2), cf in Tcore.items():
            rhs = (rhs + cf * psa[a2] % qP * psb[b2] % qP * psc[c2]) % qP
        ok7b &= (lhs % qP == rhs * scale_q % qP)
    rep("P7b Schwartz-Zippel identity check mod 2^61-1, 20 uniform points",
        ok7b)

    # ---- P8: root-box scaled-Bernstein, no negatives ----------------------
    print("-" * 78)
    print("P8  root box [0,1]^3: scaled Bernstein coefficients")

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

    degs3 = [DX, DTH, DT]
    BB = dict(TU)
    for ax in (2, 1, 0):
        BB = bern_axis(BB, degs3, ax)
    total = (DX + 1) * (DTH + 1) * (DT + 1)
    negs = sum(1 for v in BB.values() if v < 0)
    nz = total - len(BB)
    print(f"     nonzero {len(BB)}, exact zeros {nz}, negatives {negs}")
    rep("P8 NO negative Bernstein coefficients (=> core_u > 0 on OPEN box)",
        negs == 0 and len(BB) > 0)

    # ---- P9: face certificates (direct, 2-var) ----------------------------
    print("-" * 78)
    print("P9  faces u=1, tau=1 (2-var Bernstein no-negatives + nonzero)")
    print("    [no x-faces: sg = 0 (p = 3/4) and sg = 1 (p = 0) are excluded")
    print("     limits of the Gray family -- the claim is OPEN in sg]")

    def collapse(T, ax, degs):
        out = {}
        for mon, cf in T.items():
            key = tuple(m for i, m in enumerate(mon) if i != ax)
            out[key] = out.get(key, 0) + cf     # value at var=1: sum coeffs
        return {k: v for k, v in out.items() if v != 0}, \
            [d_ for i, d_ in enumerate(degs) if i != ax]

    ok9 = True
    for ax, nm in ((1, "u=1"), (2, "tau=1")):
        Tf, degf = collapse(TU, ax, degs3)
        Bf = dict(Tf)
        for ax2 in range(len(degf) - 1, -1, -1):
            Bf = bern_axis(Bf, degf, ax2)
        negf = sum(1 for v in Bf.values() if v < 0)
        okf = (negf == 0 and len(Bf) > 0)
        ok9 &= okf
        rep(f"P9 face {nm}: no negatives ({negf}), nonzero ({len(Bf)}) "
            f"=> > 0 on open face", okf)

    # ---- P10: the (u=1, tau=1) edge ---------------------------------------
    print("-" * 78)
    print("P10 edge (u=1,tau=1): 1-var no-negatives + nonzero, open in x")
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
    ok10 = (negb == 0 and any(v > 0 for v in bb))
    print(f"     end values: at x=0 -> {bb[0]}, at x=1 -> {bb[-1]} "
          f"(both degenerate corners; exact zeros allowed, ends excluded)")
    rep("P10 edge (u=1,tau=1): no negatives => > 0 on open edge x in (0,1)",
        ok10)
    print("     [union: open box + 2 open faces + 1 open edge")
    print("      = (0,1) x (0,1] x (0,1] exactly = claim region of (A)]")

    # ---- P11: certificate (B), exact Sturm --------------------------------
    print("-" * 78)
    print("P11 certificate (B): D_c < thetatilde4 * Dbar on sg in (0, 1)")
    N4, Dn4 = disc_pieces(AVAL)
    rep("P11a disc denominator == 3^10 (3-4D)^12  (> 0 for D < 3/4)",
        sp.expand(Dn4 - 3 ** 10 * (3 - 4 * D) ** 12) == 0)
    dmin = min(m[0] for m in sp.Poly(N4, D).monoms())
    C4 = sp.expand(sp.cancel(N4 / D ** dmin))
    C40 = sp.expand(C4.subs(D, 0))
    id0 = sp.expand(C40 - 3 ** 9 * p ** 4 * (1 - p) ** 2 * (3 - 4 * p) ** 4
                    * (3 * p ** 2 - 2 * p + 3))
    rep("P11b N = D^4 C, C(0,p) == 3^9 p^4 (1-p)^2 (3-4p)^4 (3p^2-2p+3)",
        dmin == 4 and id0 == 0)
    rep("P11b' 3p^2-2p+3 > 0 always (discriminant 4-36 = -32 < 0) "
        "=> disc > 0 as D->0+ on (0,3/4)",
        sp.discriminant(3 * p ** 2 - 2 * p + 3, p) == -32)
    Xcurve = THT4 * DBAR
    expr = N4.subs({D: Xcurve, p: PSUB})
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    denf = sp.factor_list(den)
    okden = all(b.is_number or sp.simplify(b - (1 + sg ** 2)) == 0
                for b, e in denf[1]) and denf[0] > 0
    rep("P11c cleared denominator = 2^32 (1+sg^2)^6  (positive)",
        okden and sp.simplify(den - 2 ** 32 * (1 + sg ** 2) ** 6) == 0)
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
    okstrip = (strips.get(sg, 0) == 13 and strips.get(1 - sg, 0) == 13
               and strips.get(1 + sg, 0) == 2)
    nroots = Rp.count_roots(0, 1)
    v0, vh, v1 = Rp.eval(0), Rp.eval(sp.Rational(1, 2)), Rp.eval(1)
    print(f"     numerator = sg^13 (1-sg)^13 (1+sg)^2 * R28; deg R28 = "
          f"{Rp.degree()}; Sturm roots in [0,1]: {nroots}; signs at "
          f"(0,1/2,1): {sp.sign(v0)},{sp.sign(vh)},{sp.sign(v1)}")
    rep("P11d numerator = sg^13 (1-sg)^13 (1+sg)^2 R28 (positive strips), "
        "R28 deg 28", okstrip and Rp.degree() == 28)
    rep("P11e R28 has NO roots in [0,1] and is negative there (Sturm)",
        nroots == 0 and v0 < 0 and vh < 0 and v1 < 0)
    rep("P11f => disc(thetatilde4*Dbar) < 0 on (0,1); with P11a/b/b' and "
        "the first-root logic, D_c < thetatilde4*Dbar",
        okstrip and nroots == 0 and v0 < 0)
    print("     [safe direction only: D_c = first sign change of disc in D")
    print("      (committed convention, lemma_7_34_dbar_upper_bound (i)-(iii));")
    print("      disc > 0 at D->0+ and disc < 0 on the curve place that first")
    print("      sign change strictly below the curve]")

    # ---- P12: exact anchor -------------------------------------------------
    print("-" * 78)
    print("P12 exact anchor ON the curve: sg=1/2, th=thetatilde4(1/2)=7/8,")
    print("    p=9/16, D=21/160, w=-1 (t=2) vs core_u sign")
    pv = sp.Rational(9, 16)
    Dv = sp.Rational(7, 8) * sp.Rational(3, 20)
    ethv = Dv / ((A - 1) * (1 - Dv))

    def etaf_w(Wv):
        E_ = ethv * Wv
        return ((A - 1) * Dv * E_ + Dv - (A - 1) * E_) / \
            ((A - 1) * Dv * E_ + Dv - (A - 1))

    # NOTE (A=4 lesson, caught in development): do NOT wrap these exact
    # Rationals in sp.nsimplify -- at ambient mpmath dps 40 nsimplify
    # "identifies" Crw = 12737761/12390400 as an only-approximately-equal
    # radical-power form (1/7) 2^(162/193) 3^(190/579) 5^(319/579) 7^(43/579),
    # which is WRONG (not exactly equal) and makes the 5x5 det explode
    # symbolically.  All quantities here are exact Rationals already; the
    # P12 gate below asserts that.
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
    rep("P12a all anchor quantities exact Rationals (no float, no nsimplify)",
        ew.is_Rational and Crw.is_Rational and chiw.is_Rational)
    vc = sum(Fr(cf) * Fr(1, 2) ** a2 for (a2, b2, c2), cf in TU.items())
    print(f"     chi(Lambda) = {float(chiw):.6e} (exact rational), "
          f"core_u(1/2,1,1) sign {'+' if vc > 0 else '-'}")
    rep("P12 exact chi > 0 at the anchor and sign matches core_u",
        chiw > 0 and vc > 0)

    # ---- P13: numeric separation scan (consistency) ------------------------
    print("-" * 78)
    print("P13 numeric separation: D_c(bisect) < thetatilde4*Dbar and chi > 0")
    print("    on the curve at t=2, 25 grid points (mpmath consistency scan)")
    ok13 = True
    for k in range(1, 26):
        sgv = k / 26.0
        pv2 = 0.75 * (1.0 - sgv * sgv)
        Db = 0.75 * (1.0 - sgv) ** 2 / (1.0 + sgv * sgv)
        dth = 1.0 - sgv * sgv / 2.0
        dc = float(Dc4_numeric(pv2))
        ok13 &= (dc / Db < dth)
        chc = chi_numeric(pv2, dth * Db, math.pi, 2.0)
        ok13 &= (chc > 0)
    rep("P13 curve strictly above theta_c and chi > 0 on the curve at t=2",
        ok13)

    # ---- P14: pinch consistency at sg -> 1 ---------------------------------
    print("-" * 78)
    print("P14 pinch consistency at sg = 0.995 (numeric): theta_c <")
    print("    thetatilde4 < theta_x (the corridor pinches ~ (0.4975,0.5025)")
    print("    in 1 - a sg^2 units; the curve a = 1/2 rides through)")
    sgv = 0.995
    pv3 = 0.75 * (1.0 - sgv * sgv)
    Db3 = 0.75 * (1.0 - sgv) ** 2 / (1.0 + sgv * sgv)
    dth3 = 1.0 - sgv * sgv / 2.0
    thc3 = float(Dc4_numeric(pv3)) / Db3
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
    print(f"     theta_c = {thc3:.7f}  thetatilde4 = {dth3:.7f}  "
          f"theta_x = {thx3:.7f}")
    rep("P14 theta_c < thetatilde4 < theta_x at the pinch point",
        okpocket and thc3 < dth3 < thx3)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("CLAIM CERTIFIED (A): chi(Lambda) > 0 for sg in (0,1),")
    print("th in (0, 1 - sg^2/2], t in (0,2]      [exact Bernstein tower]")
    print("CLAIM CERTIFIED (B): D_c(p) < (1 - sg^2/2) * Dbar(p) there")
    print("[exact Sturm; D_c = first sign change of the alternating-product")
    print(" reduced-cubic discriminant, committed lemma_7_34_dbar_upper_bound")
    print(" convention; safe direction only]")
    print("=> the A=4 Gray region for ALL p in (0, 3/4) is strictly inside")
    print("   the certified chi > 0 region -- ONE curve, no small-p splice")
    print("   (A=3 needed the separate Schur certificate below p = 133/512).")
    return PASS


if __name__ == "__main__":
    main()

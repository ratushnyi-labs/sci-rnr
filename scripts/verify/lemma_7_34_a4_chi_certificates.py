#!/usr/bin/env python3
r"""
lemma_7_34_a4_chi_certificates.py
============================================================================
BUG-009-D / Route B mid-band program, A=4 REPLICATION: the four positive-side
charpoly-derivative sign certificates
  C1: chi''''(Lambda)/24 = 5 Lambda - e1                              > 0
  C2: chi'''(Lambda)/6   = 10 Lambda^2 - 4 e1 Lambda + e2             > 0
  C3: chi''(Lambda)/2    = 10 Lambda^3 - 6 e1 Lambda^2 + 3 e2 Lambda
                           - e3                                       > 0
  C4: chi'(Lambda)       = 5 Lambda^4 - 4 e1 Lambda^3 + 3 e2 Lambda^2
                           - 2 e3 Lambda + e4                         > 0
CERTIFIED for A=4 on the ENTIRE rationalized Gray superdomain -- exact
rational arithmetic end to end, root-box Bernstein + closed-face towers.
(e_k = elementary symmetric functions of the 5x5 pattern-quotient replica
matrix Q(eta, etabar); chi = charpoly of Q; Lambda = L/Cr with
L = 1 - (3/2) D (1-D) t, t = 1 - cos s in [0,2] -- L and the target
R_4(s) >= 3/2 are unchanged from A=3.  Realness of chi's coefficients on
|w| = 1 is mechanical here: the grading gate verifies E_k/(-w)^k is
w-symmetric, hence real on the circle.)

DOMAIN (A=4 rationalization; same key move as A=3 -- rationalize the Gray
region, don't box it; lemma_7_34_dbar_upper_bound conventions):
    q = pA/(A-1),  sigma = sqrt(1-q)   =>   p = (3/4)(1-sigma^2),
    Dbar = (3/4)(1-sigma)^2/(1+sigma^2)      -- EXACTLY RATIONAL,
    D = theta * Dbar,  theta in [0,1].
The (sigma,theta,t) box [0,1] x [0,1] x [0,2] covers ALL p in (0,3/4), ALL
D in (0, Dbar] (a superset of Gray: D_c < Dbar, slack 3%..81% over p --
scout S0c), ALL spectral angles s in (0,pi].  Cleared denominators are
const*(1+sigma^2)^k -- manifestly positive.

METHOD -- t-domain small-piece assembly (the C4/C5 path of the committed
A=3 script lemma_7_34_midband_chi_certificates.py, cross-validated there
against the direct rational path for certs 1-3):
  M = c Q is POLYNOMIAL in (p,D,w), c = p(1-p) G Gt,
  G = D^2 w + (D-3)(1-D), Gt = D^2 + (D-3)(1-D) w;  Berkowitz charpoly
  gives E_k = c^k e_k; the grading E_k = (-w)^k Etld_k(t) symmetrizes
  (Chebyshev in u = w + 1/w = 2 - 2t); assembling the domain-substituted
  pieces yields ONE polynomial per certificate,
      target_j = (1+sigma^2)^maxden (chat CrN)^{5-j} chi^(j)(Lambda).

EXACT MULTIPLIER GATE (P3a-P3g; an upgrade over the committed A=3 scripts,
per review recommendation -- there the global target sign was fixed only by
the numeric chain, with an auto-flip):
    chat = p(p-1) B,   CrN = -B,   CrD = -(4D-3)^2,   with
    B(D,t) = (4D-3)^2 + 2 t D^2 (1-D)(3-D)        [G Gt = w B on |w|=1]
  =>  chat CrN = p(1-p) B^2  EXACTLY (polynomial identities, gates
  P3a-P3d,P3g), and B > 0 on the claimed domain {sigma in (0,1)} because
    3 - 4D = 3[(1-theta)(1+sigma^2) + 2 theta sigma]/(1+sigma^2) > 0
  (gate P3e, exact identity) while B is box-nonnegative (gate P3f,
  root-box Bernstein).  Hence sign(target_j) == sign(chi^(j)(Lambda))
  EXACTLY -- the numeric sign chain (90 points, mpmath dps 80; float64
  flips signs on these polynomials) is demoted to a cross-check and NO
  sign flip is tolerated (the A=3 auto-flip is gone).  A value chain
  additionally checks target == (1+sg^2)^m (chat CrN)^{5-j} chi^(j)(Lambda)
  numerically to < 1e-20 relative at 8 points per target.

RESULT: all four certify AT THE ROOT BOX (no subdivision), all face towers
close (numbers from the recorded run of this script):
  C1 (j=4): target 367 terms,   strip (1-sg)^1, core 360,   min ~ 1.89
  C2 (j=3): target 2268 terms,  strip (1-sg)^2, core 2146,  min ~ 3.91e-3
  C3 (j=2): target 6784 terms,  strip (1-sg)^3, core 6496,  min ~ 4.49e-6
  C4 (j=1): target 14994 terms, strip (1-sg)^4, core 14206, min ~ 2.62e-9
Term counts for j <= 3 coincide EXACTLY with the A=3 assembly counts
(2268 / 6784 / 14994) -- the assembly's monomial support is A-independent
there; only the E1-driven j=4 target differs (367 vs 347 on the A=3
negative side: E1's support grows with A).  On the correct domain the
positivity is Bernstein-visible without any subdivision, exactly as at
A=3: the entire difficulty was the domain, not the polynomial.

WHAT IS NOT CLAIMED (honest scope, mirroring the corrected A=3 headers):
  * cert5 = chi(Lambda) > 0 on the FULL superdomain is FALSE at A=4 as at
    A=3: the root-box Bernstein of the j=0 target has min ~ -5.5e21 (scout
    run, 26210-term target) and the real-Perron pocket has an A=4
    analogue: theta_x(sigma) exists for sigma >~ 0.3 at t=2, pocket
    present for p/pmax in ~(0.01, 0.91) (scout 3 geography), mirroring
    A=3's (0.02, 0.92).  cert5-on-Gray needs a restricted separating
    curve (A=4 candidate: thetatilde4 = 1 - sigma^2/2, scout 6) -- a
    separate deliverable, exactly as the committed A=3 split
    (lemma_7_34_midband_cert5_on_gray.py).
  * realness of the Q spectrum on the superdomain is numeric-grade only
    (scout 4: 500-point scan, max |Im|/rho ~ 1.6e-28, min relative gap
    ~ 1.9e-9), and the A=4 collision margin is MUCH tighter than A=3's:
    min collision-surface/Dbar ratio ~ 1.005 (at sigma ~ 0.7, t = 2)
    vs 1.35 at A=3.  The Budan-Fourier consequence (no real eigenvalue
    >= Lambda on Gray, once cert5-on-Gray lands) is conditional on
    realness for interior angles, as at A=3.

CHECKS (default run, ~1.5 min; exact rational arithmetic unless marked):
  P*  pieces: M = cQ polynomial; charpoly grading (E_k symmetric at
      (-w)^k, chat at (-w)^1, mN == mD); P3a-g exact multiplier
      identities + B >= 0; all domain denominators const*(1+sigma^2)^k.
  C1..C4  per certificate: dps-80 numeric sign chain (90 pts, no flip
      allowed) + dps-80 value chain (8 pts, < 1e-20 relative); zero-face
      strip; EXACT root-box Bernstein via the Fraction tensor engine
      (cross-checked against the committed A=3 sympy engine on C1: same
      min/max/zero-count); face towers theta=1 / t=2 (2-var root-box
      Bernstein) + edge Sturm, closing the claimed domain
      {sigma in (0,1)} x {theta in (0,1]} x {t in (0,2]}.

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.
"""
import itertools
import random
import time
from fractions import Fraction
from math import comb

import sympy as sp
import mpmath as mp

mp.mp.dps = 30
PASS = True


def rep(name, ok):
    global PASS
    PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}")
    return ok


p, D, w, t = sp.symbols('p D w t')
sg, th = sp.symbols('sigma theta')
X = sp.symbols('x')
U = sp.symbols('u')
AVAL = 4
A = sp.Integer(AVAL)


def _pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4


eth = D / ((A - 1) * (1 - D))


def _etaf(W):
    E = eth * W
    return ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))


STATES = list(itertools.product(range(AVAL), repeat=3))
REPS = [None] * 5
for _k, _tr in enumerate(STATES):
    if REPS[_pat(_tr)] is None:
        REPS[_pat(_tr)] = _k


def build_Q():
    eta, etb = _etaf(w), _etaf(1 / w)
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(AVAL)]
         for i in range(AVAL)]
    Q = sp.zeros(5, 5)
    for a_ in range(5):
        x, xp, xq = STATES[REPS[a_]]
        for (y, yp, yq) in STATES:
            term = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x][y])
            if yp != y: term *= eta
            if yq != y: term *= etb
            Q[a_, _pat((y, yp, yq))] += term
    return Q


def build_Cr():
    def Cf(W):
        E = eth * W
        return ((A - 1) * D * E + D - (A - 1)) / (A * D - (A - 1))
    return sp.cancel(Cf(w) * Cf(1 / w) / (Cf(1) ** 2))


# A=4 structural pieces (eta = D(1-w)(1-D)/G exactly; G Gt = w B on |w|=1)
G = D ** 2 * w + (D - (AVAL - 1)) * (1 - D)
Gt = D ** 2 + (D - (AVAL - 1)) * (1 - D) * w
CFAC = sp.expand(p * (1 - p) * G * Gt)
BDT = (AVAL * D - (AVAL - 1)) ** 2 + 2 * t * D ** 2 * (1 - D) * (AVAL - 1 - D)


def laurent_sym_to_t(expr, K=60):
    """w-symmetric Laurent polynomial -> polynomial in t via Chebyshev in
    u = w + 1/w = 2 - 2t on |w| = 1.  Returns None if not symmetric."""
    ex = sp.expand(expr)
    npoly = sp.Poly(sp.expand(ex * w ** K), w)
    co = {}
    for (md,), c in npoly.terms():
        co[md - K] = co.get(md - K, 0) + c
    ks = sorted(co.keys())
    maxk = max(abs(k) for k in ks) if ks else 0
    C = {0: sp.Integer(2), 1: U}
    for k2 in range(2, maxk + 1):
        C[k2] = sp.expand(U * C[k2 - 1] - C[k2 - 2])
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
    return sp.expand(out.subs(U, 2 - 2 * t))


def sym_to_t_mid(expr, K=60):
    """Laurent in w -> (t-polynomial, central power m of (-w)); None if the
    central power is half-integer or the Laurent part is asymmetric."""
    ex = sp.expand(expr)
    poly = sp.Poly(sp.expand(ex * w ** K), w)
    degs = [md - K for (md,), cc in poly.terms()]
    mid = sp.Rational(min(degs) + max(degs), 2)
    if int(mid) != mid:
        return None, None
    r = laurent_sym_to_t(sp.cancel(ex / (-w) ** int(mid)), K)
    return r, (int(mid) if r is not None else None)


# rationalized Gray superdomain (A=4):  p = (3/4)(1-sg^2),  D = th*Dbar
PSUB = sp.Rational(3, 4) * (1 - sg ** 2)
DBAR = sp.Rational(3, 4) * (1 - sg) ** 2 / (1 + sg ** 2)
DSUB = th * DBAR
VS = (sg, th, t)
CANDS = [1 - sg, sg, th, t, 1 + sg, 1 + sg ** 2]


def strip_trivial(polyE):
    """Iterated exact division by box-nonnegative trivial factors."""
    stripped = {}
    P = sp.Poly(polyE, *VS)
    for f in CANDS:
        fP = sp.Poly(f, *VS)
        while True:
            q_, r_ = sp.div(P, fP, *VS)
            if r_ == 0 and not q_.is_zero:
                P = sp.Poly(q_, *VS)
                stripped[f] = stripped.get(f, 0) + 1
            else:
                break
    return P.as_expr(), stripped


def dom_piece(ex):
    """Substitute the rationalized domain; return (Poly, (1+sg^2)-exponent,
    numeric constant) with denominator exactly const*(1+sg^2)^k, else None."""
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


# ------------- sympy-Rational Bernstein engine (A=3 committed) -------------
def to_tensor(polyE, vs=VS):
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


def bernstein_root_box(polyE, box, vs=VS):
    """Exact Bernstein coefficients of polyE on box; returns (min, max,
    nzeros).  nzeros counts EXACT-ZERO coefficients (dropped from the sparse
    dict by axis_apply -- recovered from the dense count).  min>0 AND
    nzeros==0 certifies CLOSED-box strict positivity; min>0 with nzeros>0
    certifies f>0 on the OPEN box and f>=0 on its closure only (a zero
    vertex coefficient equals a zero of f at that box corner -- this
    distinction once masked a boundary zero in a subdivided run)."""
    T, degs = to_tensor(polyE, vs)
    total = 1
    for d in degs:
        total *= (d + 1)
    for ax, (lo, hi) in enumerate(box):
        d = degs[ax]
        h = hi - lo
        # v -> lo + h v:  v^a = sum_b C(a,b) lo^(a-b) h^b v^b.  The h power is
        # b, NOT a: the h**a version silently evaluates every lo != 0 sub-box
        # on the wrong region (this exact bug once made a false A=3 cert5
        # "certification" pass -- caught by two independent adversarial
        # re-derivations; lo = 0 root boxes are unaffected).
        M = {}
        for a in range(d + 1):
            for b in range(a + 1):
                M[(b, a)] = sp.binomial(a, b) * lo ** (a - b) * h ** b
        T = axis_apply(T, M, ax, d + 1)
    for ax, d in enumerate(degs):
        M = {}
        for b in range(d + 1):
            for a in range(b + 1):
                M[(b, a)] = sp.Rational(comb(b, a), comb(d, a))
        T = axis_apply(T, M, ax, d + 1)
    nzeros = total - len(T)
    vals = list(T.values())
    if not vals:
        return sp.Integer(0), sp.Integer(0), nzeros
    return min(vals), max(vals), nzeros


# ---------------- Fraction tensor engine (root-box Bernstein) --------------
# Exact rational arithmetic (fractions.Fraction), ~20x faster than the sympy
# path on these targets; cross-checked against bernstein_root_box on C1.
def poly_to_ftensor(polyE, vs=VS):
    P = sp.Poly(polyE, *vs)
    degs = [P.degree(v) for v in vs]
    T = {}
    for mon, c in zip(P.monoms(), P.coeffs()):
        T[mon] = Fraction(int(sp.numer(c)), int(sp.denom(c)))
    return T, degs


def bernstein_root_box_frac(T, degs, scales=None):
    """Exact Bernstein coefficients on the root box prod_i [0, scales[i]]
    (Fraction arithmetic).  Returns (min, max, nzeros); same semantics as
    bernstein_root_box (root boxes only: lo = 0 on every axis)."""
    if scales is not None:
        out = {}
        for m, c in T.items():
            f = c
            for i, s in enumerate(scales):
                if s != 1:
                    f *= Fraction(s) ** m[i]
            out[m] = f
        T = out
    for ax, d in enumerate(degs):
        out = {}
        for mon, c in T.items():
            a = mon[ax]
            base = Fraction(1, comb(d, a))
            for b in range(a, d + 1):
                f = comb(b, a) * base
                nm = list(mon); nm[ax] = b; nm = tuple(nm)
                out[nm] = out.get(nm, Fraction(0)) + f * c
        T = {k: v for k, v in out.items() if v != 0}
    total = 1
    for d in degs:
        total *= (d + 1)
    nzeros = total - len(T)
    vals = list(T.values())
    if not vals:
        return Fraction(0), Fraction(0), nzeros
    return min(vals), max(vals), nzeros


def ftensor_eval_mp(T, pt, dps=80):
    """mpmath evaluation of a Fraction tensor at exact rational points."""
    with mp.workdps(dps):
        vals = [mp.mpf(x.numerator) / x.denominator for x in pt]
        n = len(pt)
        mx = [max(m[i] for m in T) for i in range(n)]
        pows = [[vals[i] ** k for k in range(mx[i] + 1)] for i in range(n)]
        s = mp.mpf(0)
        for m, c in T.items():
            term = mp.mpf(c.numerator) / c.denominator
            for i in range(n):
                term *= pows[i][m[i]]
            s += term
        return s


def strip_1d(e, var):
    """Strip var^k and (1-var)^k factors from a univariate polynomial."""
    P = sp.Poly(sp.expand(e), var)
    for f in (var, 1 - var):
        fP = sp.Poly(f, var)
        while True:
            q_, r_ = sp.div(P, fP, var)
            if r_ == 0 and not q_.is_zero:
                P = sp.Poly(q_, var)
            else:
                break
    return P


def face_closure(core, name, repfn):
    """The 3-var no-negative-coefficients certificate gives f > 0 only on the
    OPEN box.  The claimed domain also includes the closed faces theta=1
    (D = Dbar) and t=2 (s = pi).  Certify those:
      * face restriction, 2-var Bernstein: no negative coefficients and not
        identically zero  =>  f > 0 on the OPEN face;
      * shared edge (theta=1, t=2), 1-var: exact Sturm after stripping the
        sigma / (1-sigma) endpoint factors  =>  f > 0 on sigma in (0,1).
    Union: open box + open faces + open edge = {sigma in (0,1)} x
    {theta in (0,1]} x {t in (0,2]} -- exactly the claimed domain."""
    ok = True
    f1 = sp.expand(core.subs(th, 1))
    T1, d1 = poly_to_ftensor(f1, vs=(sg, t))
    mn1, mx1, _ = bernstein_root_box_frac(T1, d1, scales=(1, 2))
    ok &= repfn(f"{name} face theta=1: no negative coeffs, face nonzero",
                mn1 >= 0 and mx1 > 0)
    f2 = sp.expand(core.subs(t, 2))
    T2, d2 = poly_to_ftensor(f2, vs=(sg, th))
    mn2, mx2, _ = bernstein_root_box_frac(T2, d2, scales=(1, 1))
    ok &= repfn(f"{name} face t=2: no negative coeffs, face nonzero",
                mn2 >= 0 and mx2 > 0)
    e = strip_1d(core.subs({th: 1, t: 2}), sg)
    ok &= repfn(f"{name} edge theta=1,t=2: Sturm 0 roots in [0,1], positive",
                e.count_roots(0, 1) == 0 and e.eval(sp.Rational(1, 2)) > 0)
    return ok


# ------------------------- numerics (A = 4, dps 80) ------------------------
def Q_numeric(pv, Dv, sv):
    """Numeric Q and Cr at the AMBIENT mpmath precision (call in workdps)."""
    thv = mp.log(Dv / ((AVAL - 1) * (1 - Dv)))
    E = mp.e ** (thv + 1j * sv)
    Cn = ((AVAL - 1) * Dv * E + Dv - (AVAL - 1)) / (AVAL * Dv - (AVAL - 1))
    et = ((AVAL - 1) * Dv * E + Dv - (AVAL - 1) * E) / ((AVAL - 1) * Dv * E + Dv - (AVAL - 1))
    C0 = ((AVAL - 1) * Dv * mp.e ** thv + Dv - (AVAL - 1)) / (AVAL * Dv - (AVAL - 1))
    Crn = abs(Cn / C0) ** 2
    Tn = [[(1 - pv) if i == j else pv / (AVAL - 1) for j in range(AVAL)]
          for i in range(AVAL)]
    etb = mp.conj(et)
    Qn = mp.zeros(5, 5)
    for a2 in range(5):
        x, xp, xq = STATES[REPS[a2]]
        for (y, yp, yq) in STATES:
            Qn[a2, _pat((y, yp, yq))] += (Tn[xp][yp] * Tn[xq][yq] / Tn[x][y]
                                          * (et if yp != y else 1)
                                          * (etb if yq != y else 1))
    return Qn, Crn


def chi_j_value(pfr, Dfr, tfr, j, neg=False):
    """chi^(j)(Lambda) (or psi^(j)(Lambda), psi(x) = -chi(-x), if neg) from
    the numeric eigenvalues of Q, at exact rational (p, D, t) inputs and the
    ambient mpmath precision (call inside mp.workdps)."""
    pv = mp.mpf(pfr.numerator) / pfr.denominator
    Dv = mp.mpf(Dfr.numerator) / Dfr.denominator
    tv = mp.mpf(tfr.numerator) / tfr.denominator
    sv = mp.acos(1 - tv)
    Qn, Crn = Q_numeric(pv, Dv, sv)
    ev = mp.eig(Qn)[0]
    if neg:
        ev = [-lam for lam in ev]
    Lm = (1 - mp.mpf(3) / 2 * Dv * (1 - Dv) * tv) / Crn
    cs = [mp.mpc(1)]
    for lam in ev:
        ns = [mp.mpc(0)] * (len(cs) + 1)
        for k, cc in enumerate(cs):
            ns[k] += cc
            ns[k + 1] -= cc * lam
        cs = ns
    for _ in range(j):
        n = len(cs) - 1
        cs = [cs[k] * (n - k) for k in range(n)]
    val = mp.mpc(0)
    for cc in cs:
        val = val * Lm + cc
    return mp.re(val)


def _frac(x):
    """sympy Rational -> Fraction."""
    return Fraction(int(x.p), int(x.q))


# ---------------------------------------------------------------------------
# shared pieces build + gates (used by this script and the negside script)
# ---------------------------------------------------------------------------
def build_pieces(repfn, prefix):
    """Build Q, Cr, the Berkowitz charpoly of M = cQ, the t-domain pieces,
    and run the shared gates: polynomiality, grading, the P3 exact
    multiplier identities, and the domain-denominator check.  Returns the
    piece dict or None on gate failure."""
    t0 = time.time()
    Q = build_Q()
    Cr = build_Cr()
    CrN, CrD = sp.fraction(Cr)
    L = 1 - sp.Rational(3, 2) * D * (1 - D) * t
    M = sp.zeros(5, 5)
    poly_ok = True
    for i in range(5):
        for j_ in range(5):
            n_, d_ = sp.fraction(sp.cancel(CFAC * Q[i, j_]))
            poly_ok = poly_ok and bool(d_.is_number)
            M[i, j_] = sp.expand(n_ / d_)
    repfn(f"{prefix} M = c*Q polynomial in (p,D,w)", poly_ok)
    if not poly_ok:
        return None
    E = [sp.expand(co * (-1) ** k)
         for k, co in enumerate(M.charpoly(X).all_coeffs())]
    Etld, mks = [sp.Integer(1)], [0]
    for k in range(1, 6):
        ek, mk = sym_to_t_mid(E[k])
        Etld.append(ek)
        mks.append(mk)
    chat, mc = sym_to_t_mid(CFAC)
    CrNt, mN = sym_to_t_mid(CrN)
    CrDt, mD = sym_to_t_mid(CrD)
    grading_ok = (None not in mks and mc == 1 and mN == mD
                  and all(mks[k] == k * mc for k in range(6)))
    repfn(f"{prefix} charpoly grading: E_k symmetric at (-w)^k, chat at (-w)^1",
          grading_ok)
    if not grading_ok:
        return None
    # ---- P3: exact multiplier identities (proof-grade sign bookkeeping) ----
    ok = repfn(f"{prefix} P3a EXACT: chat == p(p-1)*B",
               sp.expand(chat - p * (p - 1) * BDT) == 0)
    ok &= repfn(f"{prefix} P3b EXACT: CrN == -B",
                sp.expand(CrNt + BDT) == 0)
    ok &= repfn(f"{prefix} P3c EXACT: CrD == -(4D-3)^2",
                sp.expand(CrDt + (AVAL * D - (AVAL - 1)) ** 2) == 0)
    ok &= repfn(f"{prefix} P3d EXACT: chat*CrN == p(1-p)*B^2 (multiplier)",
                sp.expand(chat * CrNt - p * (1 - p) * BDT ** 2) == 0)
    ok &= repfn(f"{prefix} P3e EXACT: 3-4D == 3[(1-th)(1+sg^2)+2 th sg]/(1+sg^2)",
                sp.cancel(sp.together(
                    3 - 4 * DSUB
                    - 3 * ((1 - th) * (1 + sg ** 2) + 2 * th * sg) / (1 + sg ** 2))) == 0)
    PB = dom_piece(BDT)
    if PB is None:
        ok &= repfn(f"{prefix} P3f B on domain: positive denom, no negative "
                    f"Bernstein coeffs", False)
    else:
        TB, dB = poly_to_ftensor(PB[0].as_expr())
        mnB, mxB, _ = bernstein_root_box_frac(TB, dB, scales=(1, 1, 2))
        ok &= repfn(f"{prefix} P3f B on domain: positive denom, no negative "
                    f"Bernstein coeffs", PB[2] > 0 and mnB >= 0 and mxB > 0)
    ok &= repfn(f"{prefix} P3g EXACT: G*Gt == w*B under t -> 1-(w+1/w)/2",
                sp.cancel(G * Gt
                          - w * BDT.subs(t, 1 - (w + 1 / w) / 2)) == 0)
    if not ok:
        return None
    PC = {'L': dom_piece(L), 'CrN': dom_piece(CrNt), 'CrD': dom_piece(CrDt),
          'chat': dom_piece(chat)}
    for k in range(1, 6):
        PC[f'E{k}'] = dom_piece(Etld[k])
    ok_dom = repfn(f"{prefix} all pieces have const*(1+sigma^2)^k denominators",
                   all(v is not None for v in PC.values()))
    if not ok_dom:
        return None
    print(f"     [pieces built in {time.time() - t0:.1f}s]")
    return PC


def assemble(j, PC, neg=False):
    """(1+sg^2)^maxden (chat CrN)^{5-j} chi^(j)(Lambda) [psi^(j) if neg] as a
    (sigma,theta,t) polynomial; returns (expr, maxden).  By gates P3a-P3f the
    multiplier (chat CrN)^{5-j} = (p(1-p)B^2)^{5-j} is > 0 on the claimed
    domain, so sign(expr) == sign(chi^(j)(Lambda)) EXACTLY."""
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
        sgn = sp.Integer(1) if neg else sp.Integer((-1) ** k)
        num = num + term * sgn * fall / consts
    return num.as_expr(), maxden


def run_target(j, PC, neg, name, repfn, seed, engine_check=False):
    """Full per-certificate pipeline: assemble, dps-80 sign + value chains,
    zero-face strip, exact root-box Bernstein, face-closure towers.
    Returns True iff every gate passes."""
    t0 = time.time()
    tg, maxden = assemble(j, PC, neg=neg)
    nterms = len(sp.Poly(tg, *VS).terms())
    print(f"     assemble: {time.time() - t0:.1f}s; {nterms} terms; "
          f"maxden {maxden}")
    Ttens, _ = poly_to_ftensor(tg)
    random.seed(seed)
    t0 = time.time()
    agree = disagree = 0
    vc_ok, vc_n = True, 0
    for _ in range(90):
        svr = Fraction(random.randint(1, 31), 32)
        tvr = Fraction(random.randint(1, 32), 16)
        thr = Fraction(random.randint(1, 16), 16)
        pfr = _frac(PSUB.subs(sg, sp.Rational(svr.numerator, svr.denominator)))
        Dfr = _frac(DSUB.subs({sg: sp.Rational(svr.numerator, svr.denominator),
                               th: sp.Rational(thr.numerator, thr.denominator)}))
        if Dfr <= Fraction(1, 10 ** 9) or Dfr >= Fraction(3, 4):
            continue
        tgv = ftensor_eval_mp(Ttens, (svr, thr, tvr), dps=80)
        with mp.workdps(80):
            cn = chi_j_value(pfr, Dfr, tvr, j, neg=neg)
            if vc_n < 8:
                pv_ = mp.mpf(pfr.numerator) / pfr.denominator
                Dv_ = mp.mpf(Dfr.numerator) / Dfr.denominator
                tv_ = mp.mpf(tvr.numerator) / tvr.denominator
                sv_ = mp.mpf(svr.numerator) / svr.denominator
                Bv = ((4 * Dv_ - 3) ** 2
                      + 2 * tv_ * Dv_ ** 2 * (1 - Dv_) * (3 - Dv_))
                pred = ((1 + sv_ ** 2) ** maxden
                        * (pv_ * (1 - pv_) * Bv ** 2) ** (5 - j) * cn)
                dn = max(abs(tgv), abs(pred))
                if dn > 0:
                    vc_ok = vc_ok and (abs(tgv - pred) / dn < mp.mpf('1e-20'))
                    vc_n += 1
        if tgv == 0:
            continue
        if (cn > 0) == (tgv > 0):
            agree += 1
        else:
            disagree += 1
    repfn(f"{name} sign chain dps80 ({agree} agree, {disagree} flip; "
          f"P3: no flip allowed)", disagree == 0 and agree >= 50)
    repfn(f"{name} value chain == (1+sg^2)^m (chat CrN)^(5-j) chi^(j), "
          f"{vc_n} pts", vc_ok and vc_n >= 5)
    print(f"     chains: {time.time() - t0:.1f}s")
    if disagree:
        return False
    t0 = time.time()
    core, stripped = strip_trivial(tg)
    strip_ok = all(f in CANDS for f in stripped)
    print(f"     target {nterms} terms -> core "
          f"{len(sp.Poly(core, *VS).terms())} terms; stripped "
          f"{[(str(kk), v) for kk, v in stripped.items()]}; "
          f"{time.time() - t0:.1f}s")
    repfn(f"{name} zero-face strip uses box-nonnegative factors only",
          strip_ok)
    t0 = time.time()
    Ct, Cdegs = poly_to_ftensor(core)
    mn_b, mx_b, nz_b = bernstein_root_box_frac(Ct, Cdegs, scales=(1, 1, 2))
    print(f"     root-box Bernstein: min-nonzero ~ {float(mn_b):.6g}, "
          f"max ~ {float(mx_b):.6g}, exact zeros: {nz_b}; "
          f"{time.time() - t0:.1f}s; degs {Cdegs}")
    ok_open = repfn(f"{name} no negative Bernstein coeffs (OPEN-box f > 0)",
                    mn_b > 0)
    if engine_check:
        box = [(sp.Integer(0), sp.Integer(1)), (sp.Integer(0), sp.Integer(1)),
               (sp.Integer(0), sp.Integer(2))]
        mn_s, mx_s, nz_s = bernstein_root_box(core, box)
        okx = (_frac(mn_s) == mn_b and _frac(mx_s) == mx_b and nz_s == nz_b)
        repfn(f"{name} engine cross-check: sympy Bernstein == Fraction "
              f"Bernstein", okx)
        ok_open = ok_open and okx
    ok_face = face_closure(core, name, repfn)
    return ok_open and ok_face


def main():
    t00 = time.time()
    print("=" * 78)
    print("A=4 mid-band chi-derivative certificates C1..C4: chi'''', chi''', chi'',")
    print("chi' > 0 at Lambda = L/Cr on the full rationalized Gray superdomain:")
    print("sigma in [0,1] (all p < 3/4), theta in [0,1] (all D <= Dbar), t in [0,2]")
    print("-- exact root-box Bernstein positivity + closed-face towers")
    print("=" * 78)
    PC = build_pieces(rep, "P*")
    if PC is None:
        print("=" * 78)
        print("OVERALL -> FAIL")
        return
    for j, cname, desc in (
            (4, "C1", "chi''''(Lambda)/24 = 5 Lambda - e1 > 0"),
            (3, "C2", "chi'''(Lambda)/6 = 10 Lambda^2 - 4 e1 Lambda + e2 > 0"),
            (2, "C3", "chi''(Lambda)/2 = 10 Lambda^3 - 6 e1 Lambda^2 "
                      "+ 3 e2 Lambda - e3 > 0"),
            (1, "C4", "chi'(Lambda) = 5 Lambda^4 - 4 e1 Lambda^3 "
                      "+ 3 e2 Lambda^2 - 2 e3 Lambda + e4 > 0")):
        print("-" * 78)
        print(f"{cname}  {desc}")
        run_target(j, PC, neg=False, name=cname, repfn=rep, seed=23,
                   engine_check=(cname == "C1"))
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}   "
          f"({time.time() - t00:.0f}s)")
    print("CERTIFIED: chi', chi'', chi''', chi'''' > 0 at Lambda on the full")
    print("A=4 superdomain (root boxes + closed faces theta=1, t=2).  NOT")
    print("claimed: chi(Lambda) > 0 full-domain (FALSE at A=4 as at A=3 --")
    print("real-Perron pocket, scout j=0 Bernstein min ~ -5.5e21); realness of")
    print("the spectrum (numeric-grade only; A=4 collision margin 1.005 x Dbar")
    print("vs 1.35 at A=3).  Remaining for R_4(s) >= 3/2 on Gray: cert5-on-Gray")
    print("(restricted curve, candidate thetatilde4 = 1 - sigma^2/2),")
    print("negative-side certificates (companion script), interior-angle")
    print("realness.")


if __name__ == "__main__":
    main()

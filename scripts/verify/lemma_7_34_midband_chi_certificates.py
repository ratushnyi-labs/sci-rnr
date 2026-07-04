#!/usr/bin/env python3
r"""
lemma_7_34_midband_chi_certificates.py
============================================================================
BUG-009-D / Route B mid-band program: the first THREE of the five
charpoly-derivative sign certificates PROVEN for A=3 on the ENTIRE
rationalized Gray superdomain -- exact rational arithmetic end to end.

TARGETS (Lambda = L/Cr, L = 1 - (3/2) D(1-D) t, t = 1-cos s in [0,2]):
  cert1:  chi''''(Lambda)/24 = 5 Lambda - e1                      > 0
  cert2:  chi'''(Lambda)/6   = 10 Lambda^2 - 4 e1 Lambda + e2     > 0
  cert3:  chi''(Lambda)/2    = 10 L^3 - 6 e1 L^2 + 3 e2 L - e3    > 0
(e_k = elementary symmetric functions of the 5x5 pattern-quotient replica
matrix Q(eta, etabar); chi = charpoly of Q, real coefficients by the proven
swap-similarity P Q(eta,etabar) P = Q(etabar,eta).)

DOMAIN (the key move -- rationalize the Gray region, don't box it):
  Rectangular (p,D) boxes FAIL: any box containing the band overshoots the
  Gray region by up to ~7x in D at moderate p, entering territory where the
  certificates have no reason to hold (Bernstein subdivision stalls there).
  Instead parametrize by the PROVEN upper bound Dbar of the Gray threshold
  (lemma_7_34_dbar_upper_bound.py):
      q = pA/(A-1),  sigma = sqrt(1-q)   =>   p = (A-1)(1-sigma^2)/A,
      Dbar = ((A-1)/A) (1-sigma)^2/(1+sigma^2)     -- EXACTLY RATIONAL,
      D = theta * Dbar,  theta in [0,1].
  The (sigma,theta,t) box [0,1] x [0,1] x [0,2] then covers ALL p in
  (0, (A-1)/A), ALL D in (0, Dbar(p)] (a superset of Gray since D_c < Dbar),
  ALL spectral angles s in (0, pi].  Cleared denominators are powers of
  81(1+sigma^2) -- manifestly positive.

METHOD (per certificate):
  (i)   build the certificate as an exact rational function of (p,D,w) with
        w = e^{is} on the unit circle (t enters through L);
  (ii)  symmetrize numerator and denominator by their central (-w)^m power
        (the two m's MUST match -- realness bookkeeping), convert the
        w-symmetric Laurent parts to polynomials in t via the Chebyshev
        recurrence C_0=2, C_1=u, C_k = u C_{k-1} - C_{k-2}, u = 2-2t;
  (iii) substitute the rationalized domain, clear the positive denominators,
        multiply numerator-of-num by numerator-of-den (sign-tracked) to get
        ONE polynomial target(sigma,theta,t) with target>0 <=> cert>0;
  (iv)  STRIP ZERO-FACE FACTORS: the target vanishes identically on the
        sigma=1 face (p -> 0 degeneracy; cert k carries (1-sigma)^{2k}
        (1+sigma)^k) -- strict positivity certification is impossible for a
        box touching a zero face; divide these nonnegative factors out
        exactly and certify the strictly positive core;
  (v)   VALIDATE the sign chain numerically at ~90 random rational points
        against a direct mpmath spectral evaluation of the certificate
        (this step caught a real sign error in development: the cleared
        denominator carries p(p-1) < 0, flipping the target);
  (vi)  certify core > 0 by EXACT MULTIVARIATE BERNSTEIN COEFFICIENTS
        (b_beta = sum_{alpha<=beta} a_alpha * prod C(beta,alpha)/C(d,alpha);
        all coefficients > 0 => polynomial > 0 on the box; the test is
        complete under subdivision for strictly positive polynomials --
        unlike the shifted-monomial test, which certified NOTHING here).

RESULT: all three certificates certify AT THE ROOT BOX (no subdivision):
  cert1: 413-term target, strip (1-sigma)^2(1+sigma), 374-term core, True
  cert2: 2680-term target, strip (1-sigma)^4(1+sigma)^2, 2350-term core, True
  cert3: 8050-term target, strip (1-sigma)^6(1+sigma)^3, 7132-term core, True
On the correct domain the positivity is Bernstein-visible without any
subdivision -- the entire difficulty was the domain, not the polynomial.

STATUS -- HONEST CORRECTION (supersedes the first version of this header).
Certificates 1-4 (chi'''', chi''', chi'', chi') are CERTIFIED on the full
rationalized superdomain, all at lo=0 ROOT BOXES (sound).  The original
cert5 (chi(Lambda) > 0 on the FULL superdomain) claim was FALSE and is
RETRACTED -- caught by two independent adversarial re-derivations:
  * chi(Lambda) < 0 on a real interior pocket: at (p = 8/15, D = Dbar,
    w = -1) exact arithmetic gives chi(Lambda) = -0.000730647... < 0 with a
    FULLY REAL spectrum whose Perron root 0.542298 exceeds Lambda =
    0.534999.  The docstring's earlier "rho Cr/L = 1.013 = complex-pair
    modulus beyond D_c" was misattributed: it is a REAL Perron root
    crossing Lambda on a curve D_perron(p) strictly between D_c and Dbar
    (D_perron/D_c ~ 1.06-1.12); the pocket is {theta_x(p) < theta <= 1},
    present for p/pmax in ~(0.02, 0.92), t near 2.
  * the original "cert5 certified with 7 nodes" rested on an affine-
    transform bug in bernstein_root_box (h**a for h**b), which evaluated
    every lo != 0 sub-box on the wrong region.  Fixed below; root-box
    certifications (certs 1-4) never touched a lo != 0 box and stand.
WHAT REMAINS TRUE AND CERTIFIED: chi'''' , chi''', chi'', chi' > 0 at
Lambda on the whole superdomain, and chi(Lambda) > 0 HOLDS ON GRAY
(D <= D_c; normalized margin >= 1.8e-7 on a fine adversarial grid --
numeric grade, certification on a Gray-containing pocket-free subdomain
is the open task).  Realness of the spectrum, per the same adversarial
round: holds on the ENTIRE superdomain (collision surface >= 1.35 x Dbar,
min over p at t=2; at w = -1 realness is PROVEN exactly on the whole strip
via the charpoly split f1*f4 and full discriminant factorization).  D_c is
spectrally INVISIBLE in Q -- the earlier "disc -> 0 at D_c" scout note was
wrong.  Remaining for R_3(s) >= 3/2 on Gray:
  (a) cert5-on-Gray: certify chi(Lambda) > 0 on theta <= thetatilde(sigma),
      a rational curve separating D_c/Dbar from D_perron/Dbar (gap >=
      0.021 uniformly);
  (b) negative-side certificates (-1)^{k+1} chi^(k)(-Lambda) > 0 (Budan-
      Fourier at +Lambda bounds only the positive side; numerically
      max |lambda_min| / Lambda ~ 0.834, ~1.20x margin -- an earlier
      0.747/1.34x figure understated the worst case);
  (c) realness certificate for interior angles (disc >= 0; at w=-1 done).

CERTIFICATION RESULTS (exact arithmetic):
  cert1: 413-term target,   strip (1-sg)^2(1+sg),      root box, min ~ 1.09e3
  cert2: 2680-term target,  strip (1-sg)^4(1+sg)^2,    root box, min ~ 2.00e6
  cert3: 8050-term target,  strip (1-sg)^6(1+sg)^3,    root box
  cert4: 14994-term target, strip (1-sg)^4,            root box
  cert5: FALSE on the full superdomain (exact pocket witness, C5 below);
         true-on-Gray is numeric-grade pending the restricted certification

CHECKS (default run, ~8 min):
  C1  cert1 pipeline: den positivity, symmetrization parity, 90-pt numeric
      sign validation, root-box Bernstein certification.
  C2  same for cert2.
  C3  (env CERT3=1; ~1h)  same for cert3.  Verified in development twice:
      via the direct rational path AND the polynomial-charpoly fast path.
  C4  cert4 via the t-domain assembly: charpoly grading check (E_k
      w-symmetric at (-w)^k), piece positivity, 90-pt numeric validation
      against eigenvalue-built chi^(1), Bernstein certification.
  C5  cert5 CORRECTED STATUS: (C5a) exact-arithmetic pocket witness
      chi(Lambda) < 0 at (p=8/15, D=Dbar, w=-1) -- the retraction is
      reproducible; (C5b) chi(Lambda) > 0 on a Gray grid (D <= 0.99 D_c
      by bisected threshold, 240 points) -- numeric-grade Gray evidence.

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.
"""
import itertools
import os
import random
import math
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
AVAL = 3
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
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(AVAL)] for i in range(AVAL)]
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


U = sp.symbols('u')


def laurent_sym_to_t(expr, K=40):
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


# rationalized Gray superdomain (A=3):  p = (2/3)(1-sg^2),  D = th*Dbar
PSUB = sp.Rational(2, 3) * (1 - sg ** 2)
DBAR = sp.Rational(2, 3) * (1 - sg) ** 2 / (1 + sg ** 2)
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
        # on the wrong region (this exact bug once made a false cert5
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
    mn1, mx1, _ = bernstein_root_box(
        f1, [(sp.Integer(0), sp.Integer(1)), (sp.Integer(0), sp.Integer(2))],
        vs=(sg, t))
    ok &= repfn(f"{name} face theta=1: no negative coeffs, face nonzero",
                mn1 >= 0 and mx1 > 0)
    f2 = sp.expand(core.subs(t, 2))
    mn2, mx2, _ = bernstein_root_box(
        f2, [(sp.Integer(0), sp.Integer(1)), (sp.Integer(0), sp.Integer(1))],
        vs=(sg, th))
    ok &= repfn(f"{name} face t=2: no negative coeffs, face nonzero",
                mn2 >= 0 and mx2 > 0)
    e = strip_1d(core.subs({th: 1, t: 2}), sg)
    ok &= repfn(f"{name} edge theta=1,t=2: Sturm 0 roots in [0,1], positive",
                e.count_roots(0, 1) == 0 and e.eval(sp.Rational(1, 2)) > 0)
    return ok


def Q_numeric(pv, Dv, sv):
    thv = mp.log(Dv / ((AVAL - 1) * (1 - Dv)))
    E = mp.e ** (thv + 1j * sv)
    Cn = ((AVAL - 1) * Dv * E + Dv - (AVAL - 1)) / (AVAL * Dv - (AVAL - 1))
    et = ((AVAL - 1) * Dv * E + Dv - (AVAL - 1) * E) / ((AVAL - 1) * Dv * E + Dv - (AVAL - 1))
    C0 = ((AVAL - 1) * Dv * mp.e ** thv + Dv - (AVAL - 1)) / (AVAL * Dv - (AVAL - 1))
    Crn = abs(Cn / C0) ** 2
    Tn = [[(1 - pv) if i == j else pv / (AVAL - 1) for j in range(AVAL)] for i in range(AVAL)]
    etb = mp.conj(et)
    Qn = mp.zeros(5, 5)
    for a2 in range(5):
        x, xp, xq = STATES[REPS[a2]]
        for (y, yp, yq) in STATES:
            Qn[a2, _pat((y, yp, yq))] += (Tn[xp][yp] * Tn[xq][yq] / Tn[x][y]
                                          * (et if yp != y else 1) * (etb if yq != y else 1))
    return Qn, Crn


def run_certificate(cert_expr, numfun, name):
    """Full pipeline (i)-(vi); returns True iff every gate passes."""
    num, den = sp.fraction(sp.cancel(sp.together(cert_expr)))

    def sym_part(ex, K=40):
        poly = sp.Poly(sp.expand(sp.expand(ex) * w ** K), w)
        degs = [md - K for (md,), c in poly.terms()]
        mid = sp.Rational(min(degs) + max(degs), 2)
        if int(mid) != mid:
            return None, None
        return laurent_sym_to_t(sp.cancel(ex / (-w) ** int(mid)), K), int(mid)

    numT, mn = sym_part(num)
    denT, md_ = sym_part(den)
    ok_sym = numT is not None and denT is not None and mn == md_
    rep(f"{name} (ii) symmetrization, matching (-w)^m powers", ok_sym)
    if not ok_sym:
        return False

    def rz(exprT):
        return sp.fraction(sp.cancel(sp.together(exprT.subs({p: PSUB, D: DSUB}))))

    n2, d2 = rz(numT)
    dn2, dd2 = rz(denT)

    def pospow(ex):
        f = sp.factor_list(ex)
        ok = all(b.is_number or sp.simplify(b - (1 + sg ** 2)) == 0 for b, e in f[1])
        return ok, f[0]

    ok1, c1 = pospow(d2)
    ok2, c2 = pospow(dd2)
    rep(f"{name} (iii) cleared denominators = const*(1+sigma^2)^k", ok1 and ok2)
    if not (ok1 and ok2):
        return False
    sign = (1 if c1 > 0 else -1) * (1 if c2 > 0 else -1)
    tg = sp.expand(n2 * dn2 * sign)
    core, stripped = strip_trivial(tg)
    strip_ok = all(f in CANDS for f in stripped)
    print(f"     target {len(sp.Poly(tg, *VS).terms())} terms -> core "
          f"{len(sp.Poly(core, *VS).terms())} terms; stripped "
          f"{[(str(k), v) for k, v in stripped.items()]}")
    rep(f"{name} (iv) zero-face strip uses box-nonnegative factors only", strip_ok)

    random.seed(17)
    okall, npts = True, 0
    for _ in range(90):
        sv = sp.Rational(random.randint(1, 31), 32)
        tv = sp.Rational(random.randint(1, 32), 16)
        thv = sp.Rational(random.randint(1, 16), 16)
        pv = float(PSUB.subs(sg, sv))
        Dv = float(DSUB.subs({sg: sv, th: thv}))
        if Dv <= 1e-9 or Dv >= 0.666:
            continue
        tvf = float(tv)
        cn = numfun(pv, Dv, math.acos(1 - tvf), tvf)
        tgv = float(tg.subs({sg: sv, th: thv, t: tv}))
        okall &= (cn > 0) == (tgv > 0)
        npts += 1
    rep(f"{name} (v) numeric sign chain, {npts} random points", okall)
    if not okall:
        return False

    box = [(sp.Integer(0), sp.Integer(1)), (sp.Integer(0), sp.Integer(1)),
           (sp.Integer(0), sp.Integer(2))]
    mn_b, mx_b, nz_b = bernstein_root_box(core, box)
    print(f"     root-box Bernstein: min-nonzero ~ {float(mn_b):.6g}, "
          f"exact zeros: {nz_b}")
    ok_open = rep(f"{name} (vi) no negative Bernstein coeffs (OPEN-box f > 0)",
                  mn_b > 0)
    ok_face = face_closure(core, name, rep)
    return ok_open and ok_face


def main():
    print("=" * 78)
    print("Mid-band chi-derivative certificates 1-3, A=3, full rationalized Gray")
    print("superdomain: sigma in [0,1] (all p), theta in [0,1] (all D <= Dbar),")
    print("t in [0,2] (all angles) -- exact Bernstein positivity")
    print("=" * 78)
    Q = build_Q()
    Cr = build_Cr()
    L = 1 - sp.Rational(3, 2) * D * (1 - D) * t
    Lam = L / Cr
    s1 = sp.cancel(sp.together(sum(Q[i, i] for i in range(5))))

    def nf1(pv, Dv, sv, tvf):
        Qn, Crn = Q_numeric(pv, Dv, sv)
        return float(5 * (1 - 1.5 * Dv * (1 - Dv) * tvf) / Crn
                     - mp.re(sum(Qn[i, i] for i in range(5))))

    print("-" * 78)
    print("C1  cert1: chi''''(Lambda)/24 = 5 Lambda - e1 > 0")
    run_certificate(5 * Lam - s1, nf1, "C1")

    print("-" * 78)
    print("C2  cert2: chi'''(Lambda)/6 = 10 Lambda^2 - 4 e1 Lambda + e2 > 0")
    Q2 = Q * Q
    s2 = sp.cancel(sp.together(sum(Q2[i, i] for i in range(5))))

    def nf2(pv, Dv, sv, tvf):
        Qn, Crn = Q_numeric(pv, Dv, sv)
        QQ = Qn * Qn
        t1 = mp.re(sum(Qn[i, i] for i in range(5)))
        t2 = mp.re(sum(QQ[i, i] for i in range(5)))
        Lm = (1 - 1.5 * Dv * (1 - Dv) * tvf) / Crn
        return float(10 * Lm ** 2 - 4 * t1 * Lm + (t1 ** 2 - t2) / 2)

    run_certificate(10 * Lam ** 2 - 4 * s1 * Lam + (s1 ** 2 - s2) / 2, nf2, "C2")

    if os.environ.get("CERT3"):
        print("-" * 78)
        print("C3  cert3: chi''(Lambda)/2 = 10 L^3 - 6 e1 L^2 + 3 e2 L - e3 > 0 (~1h)")
        Q3 = Q2 * Q
        s3 = sp.cancel(sp.together(sum(Q3[i, i] for i in range(5))))
        e2 = (s1 ** 2 - s2) / 2
        e3 = (s1 ** 3 - 3 * s1 * s2 + 2 * s3) / 6

        def nf3(pv, Dv, sv, tvf):
            Qn, Crn = Q_numeric(pv, Dv, sv)
            QQ = Qn * Qn
            QQQ = QQ * Qn
            t1 = mp.re(sum(Qn[i, i] for i in range(5)))
            t2 = mp.re(sum(QQ[i, i] for i in range(5)))
            t3 = mp.re(sum(QQQ[i, i] for i in range(5)))
            ee2 = (t1 ** 2 - t2) / 2
            ee3 = (t1 ** 3 - 3 * t1 * t2 + 2 * t3) / 6
            Lm = (1 - 1.5 * Dv * (1 - Dv) * tvf) / Crn
            return float(10 * Lm ** 3 - 6 * t1 * Lm ** 2 + 3 * ee2 * Lm - ee3)

        run_certificate(10 * Lam ** 3 - 6 * s1 * Lam ** 2 + 3 * e2 * Lam - e3,
                        nf3, "C3")
    else:
        print("-" * 78)
        print("C3  skipped (set CERT3=1 to run, ~1h; certified in development via")
        print("    both the direct rational path and the polynomial-charpoly path)")

    run_certs_45(Q)

    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("CERTIFIED: chi', chi'', chi''', chi'''' > 0 at Lambda on the full")
    print("superdomain (root boxes).  RETRACTED: full-superdomain chi(Lambda) > 0")
    print("(real Perron pocket beyond D_perron(p) in (D_c, Dbar); exact witness")
    print("C5a).  ON GRAY chi(Lambda) > 0 holds (numeric C5b); remaining for")
    print("R_3 >= 3/2 on Gray: cert5-on-Gray certification, negative-side")
    print("certificates chi^(k)(-Lambda), realness for interior angles.")


# ---------------------------------------------------------------------------
# C4/C5: t-domain small-piece assembly (chi', chi)
# ---------------------------------------------------------------------------
X = sp.symbols('x')


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


def chi_j_numeric(pv, Dv, sv, tvf, j):
    """chi^(j)(Lambda) from the numeric eigenvalues of Q."""
    Qn, Crn = Q_numeric(pv, Dv, sv)
    ev = mp.eig(Qn)[0]
    Lm = (1 - 1.5 * Dv * (1 - Dv) * tvf) / Crn
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
    return float(mp.re(val))


def run_certs_45(Q):
    print("-" * 78)
    print("C4/C5  cert4 = chi'(Lambda) > 0 and cert5 = chi(Lambda) > 0 via the")
    print("       t-domain small-piece assembly (~6 min)")
    Cr = build_Cr()
    CrN, CrD = sp.fraction(Cr)
    L = 1 - sp.Rational(3, 2) * D * (1 - D) * t
    G = D ** 2 * w + (D - 2) * (1 - D)
    Gt = D ** 2 + (D - 2) * (1 - D) * w
    c = sp.expand(p * (1 - p) * G * Gt)
    M = sp.zeros(5, 5)
    for i in range(5):
        for j in range(5):
            n_, d_ = sp.fraction(sp.cancel(c * Q[i, j]))
            M[i, j] = sp.expand(n_ / d_)
    E = [sp.expand(co * (-1) ** k)
         for k, co in enumerate(M.charpoly(X).all_coeffs())]
    Etld, mks = [sp.Integer(1)], [0]
    for k in range(1, 6):
        ek, mk = sym_to_t_mid(E[k])
        Etld.append(ek)
        mks.append(mk)
    chat, mc = sym_to_t_mid(c)
    CrNt, mN = sym_to_t_mid(CrN)
    CrDt, mD = sym_to_t_mid(CrD)
    grading_ok = (None not in mks and mc == 1 and mN == mD
                  and all(mks[k] == k * mc for k in range(6)))
    rep("C4/C5 charpoly grading: E_k symmetric at (-w)^k, chat at (-w)^1",
        grading_ok)
    if not grading_ok:
        return
    PC = {'L': dom_piece(L), 'CrN': dom_piece(CrNt), 'CrD': dom_piece(CrDt),
          'chat': dom_piece(chat)}
    for k in range(1, 6):
        PC[f'E{k}'] = dom_piece(Etld[k])
    rep("C4/C5 all pieces have const*(1+sigma^2)^k denominators",
        all(v is not None for v in PC.values()))
    if not all(v is not None for v in PC.values()):
        return

    def assemble(j):
        """LD^{5-j} chat^{5-j} chi^(j)(Lambda) as a (sigma,theta,t) Poly."""
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
        return num.as_expr()

    for j, name in ((1, "C4"),):
        tg = assemble(j)
        random.seed(23)
        agree = disagree = 0
        for _ in range(90):
            sv = sp.Rational(random.randint(1, 31), 32)
            tv = sp.Rational(random.randint(1, 32), 16)
            thv = sp.Rational(random.randint(1, 16), 16)
            pv = float(PSUB.subs(sg, sv))
            Dv = float(DSUB.subs({sg: sv, th: thv}))
            if Dv <= 1e-9 or Dv >= 0.666:
                continue
            tvf = float(tv)
            cn = chi_j_numeric(pv, Dv, math.acos(1 - tvf), tvf, j)
            tgv = float(tg.subs({sg: sv, th: thv, t: tv}))
            if tgv == 0:
                continue
            if (cn > 0) == (tgv > 0):
                agree += 1
            else:
                disagree += 1
        rep(f"{name} numeric sign chain uniform ({agree} agree, {disagree} flip)",
            agree == 0 or disagree == 0)
        if disagree and not agree:
            tg = sp.expand(-tg)
        elif disagree:
            return
        core, stripped = strip_trivial(tg)
        strip_ok = all(f in CANDS for f in stripped)
        print(f"     target {len(sp.Poly(tg, *VS).terms())} terms -> core "
              f"{len(sp.Poly(core, *VS).terms())} terms; stripped "
              f"{[(str(kk), v) for kk, v in stripped.items()]}")
        rep(f"{name} zero-face strip uses box-nonnegative factors only", strip_ok)
        box = [(sp.Integer(0), sp.Integer(1)), (sp.Integer(0), sp.Integer(1)),
               (sp.Integer(0), sp.Integer(2))]
        mn_b, _, nz_b = bernstein_root_box(core, box)
        print(f"     min-nonzero ~ {float(mn_b):.6g}, exact zeros: {nz_b}")
        rep(f"{name} no negative Bernstein coeffs (OPEN-box f > 0)", mn_b > 0)
        face_closure(core, name, rep)

    # ---- C5: corrected status for cert5 = chi(Lambda) ----
    print("-" * 78)
    print("C5a  RETRACTION WITNESS: chi(Lambda) < 0 at (p=8/15, D=Dbar, w=-1),")
    print("     exact arithmetic (real Perron root above Lambda, beyond Gray)")
    pw = sp.Rational(8, 15)
    s2 = 1 - sp.Rational(3, 2) * pw
    Dw = sp.simplify(sp.Rational(2, 3) * (1 - sp.sqrt(s2)) ** 2 / (1 + s2))
    # build Q, Cr, Lambda exactly at w = -1 (eta is real there)
    ethw = Dw / ((A - 1) * (1 - Dw))

    def etaf_w(Wv):
        E = ethw * Wv
        return ((A - 1) * Dw * E + Dw - (A - 1) * E) / ((A - 1) * Dw * E + Dw - (A - 1))

    ew, ebw = sp.simplify(etaf_w(-1)), sp.simplify(etaf_w(-1))
    Tw = [[(1 - pw) if i == j else pw / (A - 1) for j in range(AVAL)] for i in range(AVAL)]
    Qw = sp.zeros(5, 5)
    for a_ in range(5):
        x, xp, xq = STATES[REPS[a_]]
        for (y, yp, yq) in STATES:
            term = Tw[xp][yp] * Tw[xq][yq] / Tw[x][y]
            if yp != y: term *= ew
            if yq != y: term *= ebw
            Qw[a_, _pat((y, yp, yq))] += term

    def Cfw(Wv):
        E = ethw * Wv
        return ((A - 1) * Dw * E + Dw - (A - 1)) / (A * Dw - (A - 1))

    Crw = sp.simplify(Cfw(-1) * Cfw(-1) / (Cfw(1) ** 2))
    Lw = 1 - sp.Rational(3, 2) * Dw * (1 - Dw) * 2
    Lamw = sp.simplify(Lw / Crw)
    chiw = sp.simplify((Lamw * sp.eye(5) - Qw).det())
    print(f"     chi(Lambda) = {float(chiw):.9g}  (exact sign: {sp.sign(chiw)})")
    rep("C5a chi(Lambda) < 0 at the pocket witness (retraction reproducible)",
        sp.sign(chiw) == -1)

    print("C5b  chi(Lambda) > 0 on Gray (D <= 0.99 D_c, bisected threshold),")
    print("     240-point adversarial grid -- numeric grade")
    ok5b, npts5b = True, 0
    for pf in (0.05, 0.2, 0.4, 0.6, 0.8, 0.95):
        pv = pf * 2.0 / 3.0
        dc = _Dc_numeric(pv)
        for fD in (0.3, 0.7, 0.9, 0.99):
            Dv = fD * dc
            for tvf in (0.25, 0.75, 1.25, 1.6, 1.85, 2.0):
                for extra in (0, 1):
                    tt = tvf if not extra else min(2.0, tvf + 0.05)
                    val = chi_j_numeric(pv, Dv, math.acos(1 - tt), tt, 0)
                    ok5b &= val > 0
                    npts5b += 1
    rep(f"C5b chi(Lambda) > 0 at all {npts5b} Gray grid points", ok5b)


def _Dc_numeric(pv):
    """Gray threshold D_c(p) for A=3: complexification of the 3x3
    symmetry-reduced alternating product (bisection, mpmath)."""
    Aq = 3

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

    lo, hi = mp.mpf('1e-8'), mp.mpf(Aq - 1) / Aq - mp.mpf('1e-9')
    for _ in range(60):
        mid = (lo + hi) / 2
        if imflag(mid) < mp.mpf('1e-30'):
            lo = mid
        else:
            hi = mid
    return float((lo + hi) / 2)




if __name__ == "__main__":
    main()

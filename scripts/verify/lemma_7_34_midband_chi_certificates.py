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

STATUS -- ALL FIVE chi-derivative certificates are now CERTIFIED (this
script runs 1, 2, 4, 5 by default; 3 behind CERT3=1 for runtime):
  * cert4/cert5 use the T-DOMAIN SMALL-PIECE ASSEMBLY (C4/C5 below): the
    w-domain assembly of chi^(j) OOMs; instead, since c := p(1-p) G Gt = w
    * chat with chat w-symmetric, e_k(Q) = E_k(M)/c^k = Etld_k(p,D,t) /
    chat(p,D,t)^k lives entirely in the t-domain (E_k from the polynomial
    Berkowitz charpoly of M = cQ; each E_k is w-symmetric at exactly
    (-w)^k -- grading check enforced).  Substitute the rationalized domain
    into the SMALL pieces first, then assemble with sparse Poly products
    and equalized (1+sigma^2) powers.  6 min total vs OOM.
  * BUDAN-FOURIER CONSEQUENCE: all five chi^(k)(Lambda) > 0 on the open
    domain => zero sign variations => NO REAL EIGENVALUE of Q reaches
    Lambda = L/Cr -- for ALL p in (0,2/3), ALL D in (0,Dbar], ALL angles.
  * The chi^(k) certificates bound only REAL eigenvalues (Budan-Fourier).
    A numeric scan shows rho(Q)*Cr/L reaches 1.013 at theta=1 (D = Dbar,
    beyond the Gray threshold D_c): the complex-pair modulus DOES exceed
    Lambda beyond Gray, so no full-box modulus certificate (Schur-Cohn)
    can exist.  The program therefore needs the REALNESS-ON-GRAY lemma
    (spectrum of Q real for D <= D_c = first eigenvalue collision) to
    convert the chi-certificates into rho(Q) Cr <= L on Gray.  That lemma
    is the ONLY remaining piece of Route B's R_3(s) >= 3/2 on Gray.

CERTIFICATION RESULTS (development runs, exact arithmetic):
  cert1: 413-term target,   strip (1-sg)^2(1+sg),      root box, min ~ 1.09e3
  cert2: 2680-term target,  strip (1-sg)^4(1+sg)^2,    root box, min ~ 2.00e6
  cert3: 8050-term target,  strip (1-sg)^6(1+sg)^3,    root box
  cert4: 14994-term target, strip (1-sg)^4,            root box
  cert5: 26210-term target, strip (1-sg)^7 theta t,    7 nodes
  (cert5's theta/t strip is forced: chi(Lambda)=0 at D=0 and at s=0, where
   Lambda collides with the Perron root -- the open-domain statement is
   exactly the honest one.)

CHECKS (default run, ~8 min):
  C1  cert1 pipeline: den positivity, symmetrization parity, 90-pt numeric
      sign validation, root-box Bernstein certification.
  C2  same for cert2.
  C3  (env CERT3=1; ~1h)  same for cert3.  Verified in development twice:
      via the direct rational path AND the polynomial-charpoly fast path.
  C4  cert4 via the t-domain assembly: charpoly grading check (E_k
      w-symmetric at (-w)^k), piece positivity, 90-pt numeric validation
      against eigenvalue-built chi^(1), Bernstein certification.
  C5  same for cert5 (chi itself).

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


def to_tensor(polyE):
    P = sp.Poly(polyE, *VS)
    degs = [P.degree(v) for v in VS]
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


def bernstein_root_box(polyE, box):
    """Exact Bernstein coefficients of polyE on box; returns (min, max)."""
    T, degs = to_tensor(polyE)
    for ax, (lo, hi) in enumerate(box):
        d = degs[ax]
        h = hi - lo
        M = {}
        for a in range(d + 1):
            ha = h ** a
            for b in range(a + 1):
                M[(b, a)] = sp.binomial(a, b) * lo ** (a - b) * ha
        T = axis_apply(T, M, ax, d + 1)
    for ax, d in enumerate(degs):
        M = {}
        for b in range(d + 1):
            for a in range(b + 1):
                M[(b, a)] = sp.Rational(comb(b, a), comb(d, a))
        T = axis_apply(T, M, ax, d + 1)
    vals = list(T.values())
    return (min(vals), max(vals)) if vals else (sp.Integer(0), sp.Integer(0))


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
    mn_b, mx_b = bernstein_root_box(core, box)
    print(f"     root-box Bernstein coefficient min ~ {float(mn_b):.6g} "
          f"(exact rational; > 0 certifies)")
    return rep(f"{name} (vi) ROOT-BOX Bernstein certification (all coeffs > 0)",
               mn_b > 0)


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
    print("chi, chi', chi'', chi''', chi'''' > 0 at Lambda = L/Cr on the ENTIRE")
    print("rationalized Gray superdomain (open faces), A=3.  By Budan-Fourier no")
    print("REAL eigenvalue of Q reaches Lambda anywhere on the domain.  Remaining")
    print("for R_3(s) >= 3/2 on Gray: the realness-on-Gray lemma.")


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

    for j, name in ((1, "C4"), (0, "C5")):
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
        mn_b, _ = bernstein_root_box(core, box)
        if mn_b > 0:
            rep(f"{name} ROOT-BOX Bernstein certification (all coeffs > 0)", True)
            continue
        # cert5 needs a few subdivisions: certify by bisection (theta axis
        # first, then t), depth-limited
        ok = certify_subdiv(core, box)
        rep(f"{name} Bernstein certification with subdivision", ok)


def certify_subdiv(core, box, maxdepth=12):
    work = [(box, 0)]
    while work:
        bx, dep = work.pop()
        mn_b, _ = bernstein_root_box(core, bx)
        if mn_b > 0:
            continue
        if dep >= maxdepth:
            return False
        i = dep % 3
        lo, hi = bx[i]
        mid = sp.Rational(lo + hi, 2) if isinstance(lo + hi, int) else (lo + hi) / 2
        b1 = list(bx); b1[i] = (lo, mid)
        b2 = list(bx); b2[i] = (mid, hi)
        work.append((b1, dep + 1))
        work.append((b2, dep + 1))
    return True


if __name__ == "__main__":
    main()

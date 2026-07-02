#!/usr/bin/env python3
r"""
lemma_7_34_e3_beta_sign.py
============================================================================
BUG-009-D / Route B: monotonicity-in-s at its FIRST NONTRIVIAL ORDER O(D^2) --
PROVEN for A = 2, 3, 4, 5 via the one-variable sign certificate beta_A(p) <= 0.

CONTEXT. The flatness lemma (probe_7_34_monotone_s_leading_D.py) shows R_A(s;D) is
exactly flat in s through O(D): e_1 = 2t, e_2 = t q_A(p). The genuine s-dependence
starts at O(D^2) via e_3. THIS lemma closes that order:

STRUCTURE (exact, A=2..5): the t-polynomial ladder continues --
    e_3(t,p) = alpha_A(p) t + beta_A(p) t^2      (degree 2, zero constant term),
so, with R = (1-g)/(D(1-D)t) = [2 + q_A D]/(1-D) + [e_3/t] D^2/(1-D) + O(D^3),
    dR_A/dt = beta_A(p) * D^2/(1-D) + O(D^3):
monotone non-increasing in s at O(D^2)  <=>  beta_A(p) <= 0.   (ONE VARIABLE!)

CLOSED FORMS (derived by third-order D-direct rank-1 PT; verified here):
    beta_2 = -4 (2p-1)^2 / (p^2 (1-p)^2)                        -- MANIFESTLY <= 0
    beta_3 = -3 (3p-2)^2 (9p^3-21p^2+12p-4) / (4 p^2 (p-1))
    beta_4 = -4 (4p-3)^2 (256p^4-896p^3+1080p^2-576p+135) / (243 p^2 (1-p)^2)
    beta_5 = -(5p-4)^2 (375p^4-1350p^3+1710p^2-960p+224) / (128 p^2 (1-p)^2)

SIGN PROOFS (per A, on p in (0, (A-1)/A)):
    A=2: minus a square over a square -- immediate.
    A=3: denominator 4p^2(p-1) < 0, so beta_3 <= 0 <=> the cubic
         c(p) = 9p^3-21p^2+12p-4 <= 0 on (0, 2/3): c has NO real root in [0, 2/3]
         (Sturm count = 0) and c(0) = -4 < 0 -- hence c < 0 throughout, and
         beta_3 = -3(3p-2)^2 c(p)/(4p^2(p-1)) = -(3(3p-2)^2/(4p^2(1-p))) * (-c) <= 0.
    A=4: quartic q4 = 256p^4-896p^3+1080p^2-576p+135 has NO root in [0, 3/4]
         (Sturm) and q4(0) = 135 > 0 => q4 > 0 => beta_4 <= 0.
    A=5: quartic q5 = 375p^4-1350p^3+1710p^2-960p+224: same argument on [0, 4/5].

CONSEQUENCE. Combined with the endpoint theorems, the small-D Route B picture for
A>=3 is now: FLAT in s at O(D), MONOTONE NON-INCREASING in s at O(D^2) (this lemma),
endpoint >= 3/2 with the universal two-term margin -- i.e. the monotone-reduction's
both sub-lemmas hold at their leading nontrivial orders for every tested alphabet.
Residual: all-orders-in-D (large D up to D_c) monotonicity; universal-A beta formula.

UNIVERSAL FORMULA (derived by the same third-order PT with A SYMBOLIC, series-ring
in D + Laurent-safe Chebyshev w->t conversion; verified against all per-A anchors):
    beta(A,p) = -4 (1 - A(1-p))^2 P(A,p) / (p^2 (A-1)^5 (1-p)^2),
    P(A,p) = 2A^4(1-p)^4 - A^3(4p^4-20p^3+39p^2-32p+9)
             - A^2(8p^3-33p^2+40p-15) - A(6p^2-16p+11) + 3,
and P is EXACTLY the lambda_3 numerator polynomial (beta proportional to lam_3 --
internal consistency across independent derivations). Prefactor and denominator are
sign-definite, so UNIVERSAL beta<=0  <=>  P(A,p) >= 0 on {A>=2, 0<p<(A-1)/A} --
one two-variable polynomial positivity (proven per-A below; universal proof is the
remaining step).

CHECKS:
  X1  e_3 is degree 2 in t with zero constant term (A=2..5, exact).
  X2  beta_A(p) matches the closed forms above (exact).
  X3  sign proofs: A=2 manifest; A=3,4,5 via Sturm root-count 0 on the closed
      interval + endpoint sign (exact rational arithmetic).
  X4  assembled: dR/dt at O(D^2) = beta_A(p)/(1-D) * D^2 <= 0 (numeric spot-check
      against the raw replica derivative at A=3, p=0.2, D=0.4 D_c).
  X5  universal formula: subs A=2..5 reproduces the per-A closed forms (exact).
  X6  out-of-sample: A=6 via a FRESH per-A derivation matches the universal
      formula (A=6 was not used in constructing it).
  X7  sign grid: P(A,p) > 0 at exact rational points for A=6..12 (universal
      positivity proof = remaining step).

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~3-5 min.
"""
import itertools
import sympy as sp
import mpmath as mp
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

# --- import the flatness-probe machinery (derive returns E1t,E2t,E3t exact) ---
import importlib.util as ilu
_spec = ilu.spec_from_file_location(
    "flat", "/Users/para/work/rnr/scripts/verify/probe_7_34_monotone_s_leading_D.py")
_flat = ilu.module_from_spec(_spec)
import sys as _sys; _sys.modules["flat"] = _flat
_spec.loader.exec_module(_flat)
p = _flat.p; tsym = _flat.tsym
RES = _flat.RES

BETA_CLOSED = {
    2: -4 * (2 * p - 1)**2 / (p**2 * (p - 1)**2),
    3: -3 * (3 * p - 2)**2 * (9 * p**3 - 21 * p**2 + 12 * p - 4) / (4 * p**2 * (p - 1)),
    4: -4 * (4 * p - 3)**2 * (256 * p**4 - 896 * p**3 + 1080 * p**2 - 576 * p + 135) / (243 * p**2 * (p - 1)**2),
    5: -(5 * p - 4)**2 * (375 * p**4 - 1350 * p**3 + 1710 * p**2 - 960 * p + 224) / (128 * p**2 * (p - 1)**2),
}

def X1X2():
    print("-" * 78); print("X1/X2  e_3 = alpha t + beta t^2 (deg 2, no t^0); beta matches closed forms")
    ok1 = True; ok2 = True; betas = {}
    for A, (_, _, E3t) in RES.items():
        e3 = sp.expand(sp.cancel(E3t))
        deg = sp.Poly(e3, tsym).degree()
        c0 = sp.simplify(e3.coeff(tsym, 0))
        beta = sp.cancel(e3.coeff(tsym, 2)); betas[A] = beta
        ok1 = ok1 and (deg == 2) and (c0 == 0)
        match = sp.simplify(beta - BETA_CLOSED[A]) == 0
        ok2 = ok2 and match
        print(f"     A={A}: deg={deg} t^0={c0}  beta matches closed form: {match}")
    rep("X1 e_3 degree 2 in t, zero constant (ladder continues)", ok1)
    rep("X2 beta_A(p) closed forms verified", ok2)
    return betas

def X3():
    print("-" * 78); print("X3  sign proofs beta_A(p) <= 0 on (0,(A-1)/A)")
    ok = True
    # A=2: manifest
    print("     A=2: -4(2p-1)^2/(p^2(1-p)^2) -- minus square/square, manifest")
    # A=3: cubic < 0 on [0, 2/3]
    c3 = sp.Poly(9 * p**3 - 21 * p**2 + 12 * p - 4, p)
    n3 = c3.count_roots(0, sp.Rational(2, 3))
    s3 = c3.eval(0)
    ok3 = (n3 == 0) and (s3 < 0)
    print(f"     A=3: cubic roots in [0,2/3]: {n3} (want 0); c(0)={s3}<0 => cubic<0 => beta_3<=0: {ok3}")
    # A=4: quartic > 0 on [0, 3/4]
    q4 = sp.Poly(256 * p**4 - 896 * p**3 + 1080 * p**2 - 576 * p + 135, p)
    n4 = q4.count_roots(0, sp.Rational(3, 4))
    s4 = q4.eval(0)
    ok4 = (n4 == 0) and (s4 > 0)
    print(f"     A=4: quartic roots in [0,3/4]: {n4} (want 0); q(0)={s4}>0 => q>0 => beta_4<=0: {ok4}")
    # A=5: quartic > 0 on [0, 4/5]
    q5 = sp.Poly(375 * p**4 - 1350 * p**3 + 1710 * p**2 - 960 * p + 224, p)
    n5 = q5.count_roots(0, sp.Rational(4, 5))
    s5 = q5.eval(0)
    ok5 = (n5 == 0) and (s5 > 0)
    print(f"     A=5: quartic roots in [0,4/5]: {n5} (want 0); q(0)={s5}>0 => q>0 => beta_5<=0: {ok5}")
    ok = ok3 and ok4 and ok5
    return rep("X3 beta_A(p) <= 0 PROVEN (A=2 manifest; A=3,4,5 Sturm + endpoint)", ok)

def X4():
    print("-" * 78); print("X4  spot-check vs raw replica: dR/dt ~ beta_A D^2/(1-D) at A=3")
    mp.mp.dps = 40
    def _pat(tr):
        x, u, up = tr
        if x == u == up: return 0
        if x == u and u != up: return 1
        if x == up and u != up: return 2
        if u == up and x != u: return 3
        return 4
    def gA(A, pv, Dv, s):
        A = int(A); pv = mp.mpf(pv); Dv = mp.mpf(Dv); s = mp.mpf(s)
        th = mp.log(Dv / ((A - 1) * (1 - Dv))); E = mp.e**(th + 1j * s)
        C = ((A - 1) * Dv * E + Dv - (A - 1)) / (A * Dv - (A - 1))
        eta = ((A - 1) * Dv * E + Dv - (A - 1) * E) / ((A - 1) * Dv * E + Dv - (A - 1))
        C0 = ((A - 1) * Dv * mp.e**th + Dv - (A - 1)) / (A * Dv - (A - 1)); Cr = abs(C / C0)**2
        T = [[(1 - pv) if i == j else pv / (A - 1) for j in range(A)] for i in range(A)]
        etb = mp.conj(eta)
        st = list(itertools.product(range(A), repeat=3)); reps = [None] * 5
        for k, tr in enumerate(st):
            if reps[_pat(tr)] is None: reps[_pat(tr)] = k
        Q = mp.zeros(5, 5)
        for a_ in range(5):
            i = reps[a_]; x, xp, xq = st[i]
            for j, (y, yp, yq) in enumerate(st):
                Q[a_, _pat((y, yp, yq))] += T[xp][yp] * T[xq][yq] / T[x][y] * (eta if yp != y else 1) * (etb if yq != y else 1)
        ev, _ = mp.eig(Q); return Cr * max(abs(e) for e in ev)
    A = 3; pv = 0.2
    # D_c(3, 0.2) approx via known scale p^2/8
    Dv = mp.mpf('0.4') * mp.mpf(pv)**2 / 8
    def R(s):
        t = 1 - mp.cos(s)
        return (1 - gA(A, pv, Dv, s)) / (Dv * (1 - Dv) * t)
    s0 = mp.mpf('1.5'); h = mp.mpf('1e-3')
    dRdt_num = (R(s0 + h) - R(s0 - h)) / (mp.cos(s0 - h) - mp.cos(s0 + h))  # dR/dt = dR/ds / (ds->dt)
    beta_val = mp.mpf(str(sp.N(BETA_CLOSED[3].subs(p, sp.Rational(1, 5)), 30)))
    pred = beta_val * Dv**2 / (1 - Dv)
    rel = abs(dRdt_num - pred) / abs(pred)
    ok = rel < mp.mpf('0.15')  # O(D^3) corrections allowed
    print(f"     numeric dR/dt = {mp.nstr(dRdt_num,6)}  predicted beta D^2/(1-D) = {mp.nstr(pred,6)}  rel={mp.nstr(rel,3)}")
    return rep("X4 raw-replica derivative matches beta_A D^2/(1-D) (15% at O(D^3))", ok)

BETA_UNIV = None
def _beta_univ():
    global BETA_UNIV
    if BETA_UNIV is None:
        A = sp.symbols('A')
        Ppoly = (2*A**4*p**4 - 8*A**4*p**3 + 12*A**4*p**2 - 8*A**4*p + 2*A**4
                 - 4*A**3*p**4 + 20*A**3*p**3 - 39*A**3*p**2 + 32*A**3*p - 9*A**3
                 - 8*A**2*p**3 + 33*A**2*p**2 - 40*A**2*p + 15*A**2
                 - 6*A*p**2 + 16*A*p - 11*A + 3)
        BETA_UNIV = (-4*(A*p - A + 1)**2 * Ppoly / (p**2*(A-1)**5*(p-1)**2), A, Ppoly)
    return BETA_UNIV

def X5():
    print("-" * 78); print("X5  universal formula reproduces per-A closed forms")
    bu, A, _ = _beta_univ()
    ok = True
    for Av, cf in BETA_CLOSED.items():
        here = sp.simplify(bu.subs(A, Av) - cf) == 0
        ok = ok and here
        print(f"     A={Av}: universal == per-A closed form: {here}")
    return rep("X5 universal beta(A,p) matches A=2..5 anchors (exact)", ok)

def X6():
    print("-" * 78); print("X6  out-of-sample A=6: fresh per-A derivation vs universal formula")
    E1t6, E2t6, E3t6 = _flat.derive(6)
    beta6 = sp.cancel(sp.expand(sp.cancel(E3t6)).coeff(tsym, 2))
    bu, A, _ = _beta_univ()
    ok = sp.simplify(beta6 - bu.subs(A, 6)) == 0
    print(f"     A=6 fresh derivation == universal formula: {ok}")
    return rep("X6 out-of-sample A=6 match (universal formula validated)", ok)

def X7():
    print("-" * 78); print("X7  P(A,p) > 0 grid A=6..12 (universal positivity = remaining step)")
    _, A, Ppoly = _beta_univ()
    ok = True
    for Av in range(6, 13):
        for k in range(1, 20):
            pv = sp.Rational(k, 20)
            if pv >= sp.Rational(Av - 1, Av): continue
            val = Ppoly.subs({A: Av, p: pv})
            if val <= 0: ok = False; print(f"     VIOLATION A={Av} p={pv}: {val}")
    print(f"     all grid points positive: {ok}")
    return rep("X7 P(A,p)>0 on grid A=6..12 (sign extends; universal proof pending)", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("beta_A(p) <= 0: monotonicity-in-s at O(D^2) PROVEN (A=2..5)")
    print("=" * 78)
    X1X2(); X3(); X4(); X5(); X6(); X7()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("Both monotone-reduction sub-lemmas now hold at their leading nontrivial")
    print("orders for A=2..5: flat at O(D), monotone at O(D^2) (this lemma), endpoint")
    print(">= 3/2 (universal two-term corner + all-D). Residuals: all-orders-in-D;")
    print("universal-A beta formula.")

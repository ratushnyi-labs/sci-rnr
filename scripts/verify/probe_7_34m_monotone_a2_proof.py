#!/usr/bin/env python3
r"""
probe_7_34m_monotone_a2_proof.py
============================================================================
BUG-009-D / Route B, A=2: MONOTONICITY of R_2(s) in s -- PROVEN in closed form.

This closes the second sub-lemma of the monotone-reduction (the first, the endpoint
R_2(pi)>=3/2, is proven in probe_7_34m_endpoint_sharp_a2.py). Together they give the
SHARP A=2 replica inequality R_2(s) >= 3/2 for all s (vs the committed softer 11/16),
via the non-direct route "R_2 monotone in s + endpoint at s=pi" -- no convexity of the
Perron branch, no joint (s,D) bound.

CLOSED FORM (A=2, t = 1-cos s in (0,2], a = D(1-D), kappa^2=(1-2p)^2/(p^2(1-p)^2)).
From the exact factorization the replica floor is g_2 = (1/2)[F + sqrt(F^2+G)] with
    F = 1 - 2 a t,     G = 8 kappa^2 a^2 t (1 - 4a + 2 a^2 t) / (1-4a)^2,
and R_2(s) = (1-g_2)/(a t) = [(1+2at) - sqrt(F^2+G)] / (2 a t). (Verified vs the
replica g_A to 1e-36.) Writing H := F^2+G = 1 + B t + C t^2, elementary calculus gives
    sign(dR_2/dt) = -sign( 2 sqrt(H) - (2 + B t) )
(dR_2/dt = [(2+Bt) - 2 sqrt(H)]/(4 a t^2 sqrt(H)); monotone-decreasing iff
2 sqrt(H) >= 2 + B t). The squaring step below is safe ONE-SIDEDLY:
4C - B^2 >= 0 gives 2 sqrt(H) >= |2 + B t| >= 2 + B t regardless of the sign
of 2+Bt (and in fact 2+Bt > 0 on the Gray domain),
and since 2H - t H' = 2 + B t exactly, the monotone condition 2 sqrt(H) >= 2 + B t
squares (both sides >=0 in the binding regime) to t^2 (4C - B^2) >= 0, i.e.

    dR_2/dt <= 0  for all t   <=>   4C - B^2 >= 0.

THE DISCRIMINANT (exact):  4C - B^2 = 8 a m (1 - (3 + kappa^2) a),  m = 8 kappa^2 a^2/(1-4a)^2,
so   4C - B^2 >= 0   <=>   a <= 1/(3 + kappa^2).

ON THE GRAY REGION  a <= a_c = D_c(1-D_c) = p^2/(4(1-p)^2)  (D_c = 1/2 - sqrt(1-2p)/(2(1-p))),
the binding is a = a_c, and 1 - (3+kappa^2) a_c >= 0 is equivalent to the polynomial
inequality
    4(1-p)^4 - 3 p^2 (1-p)^2 - (1-2p)^2  =  p^4 - 10 p^3 + 17 p^2 - 12 p + 3  >= 0,
which, in u := 1/2 - p in (0,1/2), is
    (16 u^4 + 128 u^3 + 56 u^2 + 32 u + 1) / 16  >  0     (ALL coefficients positive).
Hence a_c < 1/(3+kappa^2) STRICTLY on p in (0,1/2), so 4C-B^2 > 0, so dR_2/dt < 0 on
the whole Gray region -- MONOTONICITY PROVEN.

CONCLUSION (A=2, closed form, both sub-lemmas proven):
    R_2(s) >= R_2(pi) >= 3/2   for all s in (0,pi], on the entire A=2 Gray region.
The sharp constant 3/2 (tight at s=pi, p->0, D->D_c). A>=3 remains open (its endpoint
does not reduce to a 2x2 block; same wall as Route B's general-A charpoly effort).

CHECKS:
  T1  closed-form R_2(t) equals the replica g_2 (via gA), t=1-cos s, to 1e-30.
  T2  (sympy) H=F^2+G has H(0)=1 and 4C-B^2 = 8 a m (1-(3+kappa^2)a) exactly.
  T3  (sympy) boundary poly = p^4-10p^3+17p^2-12p+3, and in u=1/2-p equals
      (16u^4+128u^3+56u^2+32u+1)/16 -- all positive coeffs => >0 on (0,1/2) (PROVEN).
  T4  a_c < 1/(3+kappa^2) strictly on the Gray region (=> 4C-B^2>0 => monotone).
  T5  assembled: dR_2/dt < 0 on the Gray grid (numeric confirm of the proof).

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~1 min.
"""
import itertools
import mpmath as mp
mp.mp.dps = 40
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

def _pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def gA(A, p, D, s):
    A = int(A); p = mp.mpf(p); D = mp.mpf(D); s = mp.mpf(s)
    th = mp.log(D / ((A - 1) * (1 - D))); E = mp.e**(th + 1j * s)
    C = ((A - 1) * D * E + D - (A - 1)) / (A * D - (A - 1))
    eta = ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))
    C0 = ((A - 1) * D * mp.e**th + D - (A - 1)) / (A * D - (A - 1)); Cr = abs(C / C0)**2
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]; etb = mp.conj(eta)
    st = list(itertools.product(range(A), repeat=3)); npat = 5 if A >= 3 else 4; reps = [None] * npat
    for k, tr in enumerate(st):
        pp = _pat(tr)
        if reps[pp] is None: reps[pp] = k
    Q = mp.zeros(npat, npat)
    for a_ in range(npat):
        i = reps[a_]; x, xp, xq = st[i]
        for j, (y, yp, yq) in enumerate(st):
            Q[a_, _pat((y, yp, yq))] += T[xp][yp] * T[xq][yq] / T[x][y] * (eta if yp != y else 1) * (etb if yq != y else 1)
    ev, _ = mp.eig(Q); return Cr * max(abs(e) for e in ev)

def kap2(p): p = mp.mpf(p); return (1 - 2 * p)**2 / (p**2 * (1 - p)**2)
def Dc2(p): p = mp.mpf(p); return mp.mpf('0.5') - mp.sqrt(1 - 2 * p) / (2 * (1 - p))

def R2_closed(p, D, t):
    a = D * (1 - D); k2 = kap2(p); F = 1 - 2 * a * t
    G = 8 * k2 * a**2 * t * (1 - 4 * a + 2 * a**2 * t) / (1 - 4 * a)**2
    return ((1 + 2 * a * t) - mp.sqrt(F**2 + G)) / (2 * a * t)

PS = ['0.02', '0.05', '0.1', '0.2', '0.35', '0.45']

def T1():
    print("-" * 78); print("T1  closed-form R_2(t) == replica g_2 (t=1-cos s), all s")
    worst = mp.mpf(0)
    for p in PS:
        dc = Dc2(p)
        for fr in ('0.3', '0.9'):
            D = mp.mpf(fr) * dc
            for sd in (15, 55, 100, 150, 179):
                s = mp.mpf(sd) / 180 * mp.pi; t = 1 - mp.cos(s)
                worst = max(worst, abs(R2_closed(p, D, t) - (1 - gA(2, p, D, s)) / (D * (1 - D) * t)))
    print(f"     worst |closed-form - replica| = {mp.nstr(worst, 3)}")
    return rep("T1 closed-form R_2(t) matches replica g_2", worst < mp.mpf('1e-30'))

def T2():
    print("-" * 78); print("T2  (sympy) H(0)=1 and 4C-B^2 = 8 a m (1-(3+kappa^2)a)")
    import sympy as sp
    a, t, k2 = sp.symbols('a t k2', positive=True)
    m = 8 * k2 * a**2 / (1 - 4 * a)**2
    H = sp.expand((1 - 2 * a * t)**2 + m * ((1 - 4 * a) * t + 2 * a**2 * t**2))
    Hp = sp.Poly(H, t); C = Hp.coeff_monomial(t**2); B = Hp.coeff_monomial(t); c0 = Hp.coeff_monomial(1)
    ok0 = sp.simplify(c0 - 1) == 0
    disc = sp.simplify(4 * C - B**2); target = 8 * a * m * (1 - (3 + k2) * a)
    okd = sp.simplify(disc - target) == 0
    print(f"     H(0)={sp.simplify(c0)}   4C-B^2 == 8 a m (1-(3+k2)a): {okd}")
    return rep("T2 monotone discriminant 4C-B^2 = 8 a m (1-(3+kappa^2)a)", bool(ok0 and okd))

def T3():
    print("-" * 78); print("T3  (sympy) boundary poly p^4-10p^3+17p^2-12p+3 = (16u^4+..+1)/16 (u=1/2-p), all + coeffs")
    import sympy as sp
    p = sp.symbols('p', positive=True); u = sp.symbols('u', positive=True)
    poly = sp.expand(4 * (1 - p)**4 - 3 * p**2 * (1 - p)**2 - (1 - 2 * p)**2)
    okform = sp.simplify(poly - (p**4 - 10 * p**3 + 17 * p**2 - 12 * p + 3)) == 0
    polyu = sp.expand(poly.subs(p, sp.Rational(1, 2) - u))
    coeffs = sp.Poly(polyu * 16, u).all_coeffs()  # [16,128,56,32,1]
    allpos = all(cc > 0 for cc in coeffs)
    print(f"     poly == p^4-10p^3+17p^2-12p+3: {okform}   16*poly(u)= {coeffs} (all>0: {allpos})")
    return rep("T3 boundary poly manifestly >0 on (0,1/2) (positive-coeff SOS-like)", bool(okform and allpos))

def T4():
    print("-" * 78); print("T4  a_c < 1/(3+kappa^2) strictly on the Gray region (=> 4C-B^2>0 => monotone)")
    ok = True
    for p in PS:
        a_c = mp.mpf(p)**2 / (4 * (1 - mp.mpf(p))**2); thr = 1 / (3 + kap2(p))
        ok = ok and a_c < thr
        print(f"     p={p}: a_c={mp.nstr(a_c,5)}  1/(3+kappa^2)={mp.nstr(thr,5)}  a_c<thr:{a_c<thr}")
    return rep("T4 a_c < 1/(3+kappa^2) on Gray (monotonicity certificate holds)", ok)

def T5():
    print("-" * 78); print("T5  assembled: dR_2/dt < 0 on the Gray grid (numeric confirm of the proof)")
    worst = -mp.inf
    for p in PS:
        dc = Dc2(p)
        for fr in [mp.mpf(f) for f in ('0.1', '0.5', '0.9', '0.999')]:
            D = mp.mpf(fr) * dc
            for j in range(1, 100):
                t = mp.mpf(2) * j / 100
                worst = max(worst, mp.diff(lambda tt: R2_closed(p, D, tt), t))
    print(f"     max dR_2/dt over Gray grid = {mp.nstr(worst,3)} (<0 => monotone decreasing)")
    return rep("T5 dR_2/dt < 0 on the Gray region (R_2 monotone; min at s=pi)", worst < 0)

if __name__ == "__main__":
    print("=" * 78)
    print("BUG-009-D / Route B (A=2): monotonicity of R_2(s) -- PROVEN in closed form")
    print("=" * 78)
    T1(); T2(); T3(); T4(); T5()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("A=2 SHARP inequality R_2(s)>=3/2 for all s is PROVEN (monotone + endpoint),")
    print("stronger than the committed 11/16. A>=3 endpoint remains open (no 2x2 reduction).")

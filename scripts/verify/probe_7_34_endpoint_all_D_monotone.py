#!/usr/bin/env python3
r"""
probe_7_34_endpoint_all_D_monotone.py
============================================================================
BUG-009-D / Route B: the endpoint infimum over the WHOLE Gray region sits at the
D_c corner -- R_A(pi; D) is strictly DECREASING in D on (0, D_c], so
    inf_{0<D<=D_c} R_A(pi; D) = R_A(pi; D_c),
and the (proven, universal) two-term corner theorem
    R_A(pi)|_{D_c} - 3/2 = [(A-2)/(A-1)] p - [(4A^2+24A-65)/(8(A-1)^2)] p^2 + O(p^3)
then bounds the endpoint on ALL of the Gray region (small p): the all-D endpoint
sub-lemma reduces to the already-derived corner.

MECHANISM. R_A(pi; D) = 2 - E_A(p) D + F_A(p) D^2 + ... with E_A = [2(A-1)/p^2](1+...)
POSITIVE and large, while the curvature correction obeys 2 F_A D <= 2 F_A D_c with
D_c ~ p^2/(4(A-1)): the ratio 2 F D_c / E -> 0 as p -> 0 (F's leading order is O(1/p^2),
same as E, so the ratio is O(D_c) = O(p^2)). Hence R'(D) = -E (1 + O(p^2)) < 0
throughout the Gray interval for small p -- monotone decrease, minimum at D_c.
(At A=2 this is also visible in the exact closed form.)

STATUS. The D-monotonicity is verified numerically across the full grid (V1, p<=0.3)
and justified at leading order via the symbolic ratio (V3); the corner value is a
THEOREM (universal two-term margin). SCOPE CORRECTION (hostile sweep, see
probe_7_34_corner_pocket_refutation.py): D-monotonicity FAILS in a bounded pocket
near the uniform corner p >= ~0.98 (A-1)/A for every A>=3 -- this probe's mechanism
is REGIONAL (clean for p <= 0.97 (A-1)/A, which contains this probe's entire grid);
in the pocket the target R>=3/2 holds by direct margin (>=1.74) instead. A fully
rigorous all-orders-in-D proof on the clean region needs remainder control beyond
4 PT orders -- same class as the O(p^3)/uniformity gap, flagged honestly.

CHECKS:
  V1  dR_A(pi;D)/dD < 0 numerically across A=2..6, p in {0.02,0.1,0.3}, D over
      20 points spanning (0.02, 0.999) D_c  (central differences, mpmath 40 dps).
  V2  consequence: min_D R_A(pi;D) = R_A(pi;D_c) >= 3/2 on the grid; the all-D
      endpoint inherits the corner margin.
  V3  symbolic (per-A, exact): the decrease-dominance ratio 2 F_A(p) D_c / E_A(p)
      is O(p^2) with a finite leading coefficient (computed exactly for A=3,4,5),
      so R'(D) < 0 on (0,D_c] at leading order in p.

Deps: mpmath, sympy.  Python: /Users/para/.venvs/rnr/bin/python.  ~3-6 min.
"""
import itertools
import mpmath as mp
import sympy as sp
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

def Dc_num(A, pv):
    import numpy as np
    if A == 2:
        return float(mp.mpf('0.5') - mp.sqrt(1 - 2 * mp.mpf(pv)) / (2 * (1 - mp.mpf(pv))))
    def im(Dv):
        b0 = Dv / (A - 1); K = np.full((A, A), b0); np.fill_diagonal(K, 1 - Dv); Ki = np.linalg.inv(K)
        Tn = np.full((A, A), pv / (A - 1)); np.fill_diagonal(Tn, 1 - pv)
        By = lambda y: np.array([[Ki[y, a] * Tn[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-14, (A - 1) / A - 1e-7
    for _ in range(90):
        m = (lo + hi) / 2; lo, hi = (m, hi) if im(m) < 1e-13 else (lo, m)
    return lo

def R(A, p, D):
    return (1 - gA(A, p, D, mp.pi)) / (D * (1 - D) * 2)

def V1V2():
    print("-" * 78)
    print("V1/V2  dR/dD < 0 across the Gray interval; min at D_c; corner >= 3/2")
    ok1 = True; ok2 = True
    for A in (2, 3, 4, 5, 6):
        for pv in (0.02, 0.1, 0.3):
            if pv >= (A - 1) / A - 1e-6: continue
            dc = Dc_num(A, pv)
            fracs = [0.02 + 0.979 * k / 19 for k in range(20)]
            vals = [R(A, pv, f * dc) for f in fracs]
            mono = all(vals[k + 1] < vals[k] + mp.mpf('1e-12') for k in range(len(vals) - 1))
            corner = vals[-1]
            ok1 = ok1 and mono
            ok2 = ok2 and (min(vals) >= corner - mp.mpf('1e-10')) and (corner >= mp.mpf('1.5') - mp.mpf('1e-9'))
            print(f"     A={A} p={pv}: monotone-decreasing:{mono}  min=corner(0.999Dc)={mp.nstr(corner,7)} >=3/2:{corner>=1.5}")
    rep("V1 R_A(pi;D) strictly decreasing in D on (0,D_c] (grid)", ok1)
    rep("V2 min_D R = corner value >= 3/2 (all-D endpoint inherits corner margin)", ok2)

def V3():
    print("-" * 78)
    print("V3  symbolic dominance: 2 F_A(p) D_c / E_A(p) = O(p^2) exactly (A=3,4,5)")
    ok = True
    for Aval in (3, 4, 5):
        A = sp.Integer(Aval)
        p, et, D = sp.symbols('p eta_s D')
        T = [[(1 - p) if i == j else p / (A - 1) for j in range(Aval)] for i in range(Aval)]
        st = list(itertools.product(range(Aval), repeat=3)); reps = [None] * 5
        for k, tr in enumerate(st):
            if reps[_pat(tr)] is None: reps[_pat(tr)] = k
        Q0 = sp.zeros(5, 5); Q1 = sp.zeros(5, 5); Q2 = sp.zeros(5, 5)
        for a_ in range(5):
            x, xp, xq = st[reps[a_]]
            for (y, yp, yq) in st:
                t = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x][y])
                ne = (1 if yp != y else 0) + (1 if yq != y else 0)
                [Q0, Q1, Q2][ne][a_, _pat((y, yp, yq))] += t
        Q0 = Q0.applyfunc(sp.cancel); Q1 = Q1.applyfunc(sp.cancel); Q2 = Q2.applyfunc(sp.cancel)
        v = Q0[:, 0]
        xs = {0: v}; lams = {}
        for k in range(1, 5):
            r = Q1 * xs[k - 1] + (Q2 * xs[k - 2] if k >= 2 else sp.zeros(5, 1))
            lam_k = sp.cancel(r[0]); xk = r - lam_k * v
            for j in range(1, k): xk = xk - lams[j] * xs[k - j]
            xs[k] = xk.applyfunc(sp.cancel); lams[k] = lam_k
        etaD = 2 * D * (1 - D) / (A * D - 2 * D**2 - (A - 1))
        CrD = ((A * D - 2 * D**2 - (A - 1)) / (A * D - (A - 1)))**2
        etas = sp.series(etaD, D, 0, 4).removeO()
        Crs = sp.series(CrD, D, 0, 4).removeO()
        g = sp.expand(Crs * (1 + sum(sp.expand(lams[k] * etas**k) for k in range(1, 5))))
        omg = sp.expand(1 - g)
        e2 = sp.cancel(sp.expand(omg).coeff(D, 2)); e3 = sp.cancel(sp.expand(omg).coeff(D, 3))
        E = sp.cancel(-(e2 + 4) / 2); F = sp.cancel((4 + e2 + e3) / 2)
        Dc = sp.Rational(1, 4 * (Aval - 1)) * p**2
        ratio = sp.cancel(2 * F * Dc / E)
        lead = sp.limit(ratio / p**2, p, 0)
        ok_here = (lead not in (sp.oo, -sp.oo, sp.zoo)) and sp.simplify(lead) != 0
        ok = ok and ok_here
        print(f"     A={Aval}: 2FD_c/E = [{sp.nsimplify(lead)}] p^2 + O(p^3)  (finite nonzero: {ok_here})")
    return rep("V3 decrease-dominance ratio = O(p^2) exactly => R'(D)<0 at small p", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("BUG-009-D endpoint: all-D reduction -- R_A(pi;D) decreasing => infimum at D_c corner")
    print("=" * 78)
    V1V2(); V3()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("All-D endpoint inherits the (proven) two-term corner margin at small p.")
    print("Honest residual: all-orders-in-D monotonicity for large p (same class as the")
    print("O(p^3)/uniformity remainder).")

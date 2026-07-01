#!/usr/bin/env python3
"""
probe_7_34m_monotone_route.py
============================================================================
BUG-009-D non-direct route: replace the (hard) CONVEXITY of G(u)=g_A(arccos u)
with the (simpler, first-order) MONOTONICITY of the replica floor ratio
    R_A(s) = [1 - g_A(s)] / [D(1-D)(1-cos s)]   on s in (0, pi].

WHY THIS IS A GAP-SPLIT (project rule: non-direct paths, split gaps smaller).
The achievability of the operational RD-dispersion (Thm 7.34m / Rem 7.34e')
needs the s=pi BINDING: min_{s in (0,pi]} R_A(s) attained at s=pi. Remark 7.34m'
attacks this via the SECOND-order condition "G(u):=g_A(arccos u) convex on
[-1,1]" (equivalently g4 >= v/12 at the endpoint) -- the sole remaining
conjecture, only proven at leading order in D (c_A>0), all A, on the Gray bulk.

A monotone non-increasing R_A(s) is a strictly WEAKER, FIRST-order sufficient
condition for the same s=pi binding: if R_A' <= 0 on (0,pi] then the minimum is
at the right endpoint s=pi automatically, with NO appeal to convexity of the
Perron branch. First-order sign conditions are typically far easier to certify
(one transfer-matrix derivative) than second-order ones, so IF R_A is monotone
this splits BUG-009-D into a smaller, cleaner sub-lemma.

THIS PROBE tests, numerically and honestly, whether that non-direct route is
even AVAILABLE: is R_A(s) monotone non-increasing on (0,pi], for A=2..4, several
p, and D across the whole Gray region (0, D_c)? It reports:
  * whether the grid-minimum of R_A sits at s=pi (the binding, model-free);
  * whether R_A is monotone non-increasing (the non-direct sufficient condition);
  * R_A(pi) and R_A(0+) endpoints (R_A(0+)=2v/(D(1-D)) analytically).
This is EXPLORATION that either green-lights the monotone sub-lemma or refutes
it (in which case convexity stays the only known route). It does NOT by itself
prove achievability -- it maps whether the simpler route exists.

RELATION TO ROUTE B (the active BUG-009-D lead).  Route B's certificate is the
s-UNIFORM spectral inequality  g_A(s) <= 1 - (3/2) D(1-D)(1-cos s)  on the Gray
region, i.e. exactly  R_A(s) >= 3/2  for all s (proven there directly via a
characteristic-polynomial localization for every s, with the binding corner
s=pi, D=D_c, p->0 checked as C3).  Monotonicity DECOMPOSES that single uniform
inequality into two smaller sub-lemmas -- a "split the gap smaller" reduction:
     (i)  R_A(s) monotone non-increasing in s   (first-order, s-direction);
     (ii) R_A(pi) >= 3/2                          (endpoint, D-direction only).
(i)+(ii) => R_A(s) >= 3/2 for all s.  (i) is exactly the statement that s=pi is
the binding corner -- so this supplies the conceptual JUSTIFICATION for why
Route B need only check the s=pi corner, and an alternative proof path that
avoids a joint (s,D) bound.  (ii) is tight: R_A(pi) -> 3/2 as A=2, p->0, D->D_c.

CHECKS (PASS = the non-direct monotone-reduction route is available on the grid):
  M1  grid-min of R_A(s) is at s=pi (binding holds pointwise), all params.
  M2  R_A(s) monotone non-increasing on (0,pi] (first-order route available).
  M3  consistency: R_A(0+) ~ 2v/(D(1-D)) (v the per-site posterior variance).
  M4  endpoint sub-lemma R_A(pi) >= 3/2 up to D_c (tight at A=2,p->0,D->D_c);
      together with M2 this reconstructs Route B's R_A(s)>=3/2 uniform bound.

Deps: mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~1-3 min.
"""
import itertools
import mpmath as mp
mp.mp.dps = 30
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<64} {'PASS' if ok else 'FAIL'}"); return ok

def _pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def gA(A, p, D, s):
    """g_A(s) = |C_s/C_0|^2 rho(Q_A(s)) via the (<=5x5) S_A pattern quotient."""
    A = int(A); p = mp.mpf(p); D = mp.mpf(D); s = mp.mpf(s)
    th = mp.log(D / ((A - 1) * (1 - D))); E = mp.e**(th + 1j * s); den = A * D - (A - 1)
    C = ((A - 1) * D * E + D - (A - 1)) / den
    eta = ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))
    C0 = ((A - 1) * D * mp.e**th + D - (A - 1)) / den; Cr = abs(C / C0)**2
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]
    etb = mp.conj(eta)
    st = list(itertools.product(range(A), repeat=3))
    npat = 5 if A >= 3 else 4; reps = [None] * npat
    for k, tr in enumerate(st):
        pp = _pat(tr)
        if reps[pp] is None: reps[pp] = k
    Q = mp.zeros(npat, npat)
    for a in range(npat):
        i = reps[a]; x, xp, xq = st[i]
        for j, (y, yp, yq) in enumerate(st):
            Q[a, _pat((y, yp, yq))] += T[xp][yp] * T[xq][yq] / T[x][y] * (eta if yp != y else 1) * (etb if yq != y else 1)
    ev, _ = mp.eig(Q); return Cr * max(abs(e) for e in ev)

def Dc(A, pv):
    import numpy as np
    def im(Dv):
        b0 = Dv / (A - 1); K = np.full((A, A), b0); np.fill_diagonal(K, 1 - Dv); Ki = np.linalg.inv(K)
        Tn = np.full((A, A), pv / (A - 1)); np.fill_diagonal(Tn, 1 - pv)
        By = lambda y: np.array([[Ki[y, a] * Tn[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-12, (A - 1) / A - 1e-7
    for _ in range(60):
        m = (lo + hi) / 2; lo, hi = (m, hi) if im(m) < 1e-11 else (lo, m)
    return lo

def R_of_s(A, p, D, s):
    return (1 - gA(A, p, D, s)) / (D * (1 - D) * (1 - mp.cos(s)))

PARAMS = []
for A in (2, 3, 4):
    for p in ('0.1', '0.25'):
        if mp.mpf(p) >= mp.mpf(A - 1) / A: continue
        dc = Dc(A, float(p))
        for frac in ('0.15', '0.4', '0.7', '0.92'):
            PARAMS.append((A, p, mp.mpf(frac) * dc, frac))

def M1_M2():
    print("-" * 78)
    print("M1/M2  grid-min at s=pi (binding)  &  R_A(s) monotone non-increasing")
    Mgrid = 60
    ok_bind = True; ok_mono = True
    worst_mono = mp.mpf(0); worst_bindgap = mp.mpf(0)
    for (A, p, D, frac) in PARAMS:
        sj = [mp.pi * j / Mgrid for j in range(1, Mgrid + 1)]  # (0, pi], includes pi
        Rj = [R_of_s(A, p, D, s) for s in sj]
        # M1: min at s=pi (last node)
        argmin = max(range(len(Rj)), key=lambda k: -Rj[k])
        at_pi = (argmin == len(Rj) - 1)
        bindgap = Rj[argmin] - Rj[-1]  # <=0 if pi is min; magnitude of violation if not
        # M2: monotone non-increasing: R_{k+1} - R_k <= tol
        incr = max((Rj[k + 1] - Rj[k]) for k in range(len(Rj) - 1))
        mono = incr <= mp.mpf('1e-9')
        ok_bind = ok_bind and at_pi
        ok_mono = ok_mono and mono
        worst_mono = max(worst_mono, incr)
        worst_bindgap = min(worst_bindgap, bindgap)
        tag = 'OK' if (at_pi and mono) else ('MIN@interior' if not at_pi else 'NON-MONO')
        print(f"     A={A} p={p} D={frac}Dc: R(pi)={mp.nstr(Rj[-1],5)} R(0+)~={mp.nstr(Rj[0],5)} "
              f"min@pi={at_pi} mono={mono} maxDelta={mp.nstr(incr,2)} [{tag}]")
    rep("M1 grid-min of R_A(s) at s=pi (binding holds pointwise)", ok_bind)
    rep(f"M2 R_A(s) monotone non-increasing (worst up-step={mp.nstr(worst_mono,2)})", ok_mono)
    return ok_bind, ok_mono

def M3():
    print("-" * 78)
    print("M3  endpoint consistency  R_A(0+) ~ 2 v / (D(1-D))   (v = per-site variance)")
    ok = True
    for (A, p, D, frac) in PARAMS[:4]:
        v = -mp.diff(lambda s: gA(A, p, D, s), 0, 2) / 2
        pred = 2 * v / (D * (1 - D))
        s0 = mp.mpf('1e-3'); got = R_of_s(A, p, D, s0)
        relerr = abs(got - pred) / abs(pred)
        ok = ok and relerr < mp.mpf('1e-3')
        print(f"     A={A} p={p} D={frac}Dc: R(0+)={mp.nstr(got,6)} 2v/(D(1-D))={mp.nstr(pred,6)} rel={mp.nstr(relerr,2)}")
    return rep("M3 R_A(0+) matches 2v/(D(1-D)) analytic endpoint", ok)

def M4():
    print("-" * 78)
    print("M4  endpoint sub-lemma  R_A(pi) >= 3/2  up to D_c (tight at A=2,p->0,D->D_c)")
    ok = True
    for A in (2, 3, 4, 5):
        for p in ('0.05', '0.15'):
            if mp.mpf(p) >= mp.mpf(A - 1) / A: continue
            dc = Dc(A, float(p))
            vals = [ (1 - gA(A, p, mp.mpf(f) * dc, mp.pi)) / (mp.mpf(f) * dc * (1 - mp.mpf(f) * dc) * 2)
                     for f in ('0.5', '0.9', '0.999') ]
            here = all(v >= mp.mpf('1.5') - mp.mpf('1e-9') for v in vals)
            ok = ok and here
            print(f"     A={A} p={p}: R(pi)@0.5,0.9,0.999 Dc = " +
                  ", ".join(mp.nstr(v, 7) for v in vals) + f"  >=3/2:{here}")
    return rep("M4 endpoint R_A(pi) >= 3/2 up to D_c (with M2 => R_A(s)>=3/2 all s)", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("BUG-009-D non-direct route probe: monotonicity of R_A(s) vs convexity of G(u)")
    print("=" * 78)
    ob, om = M1_M2()
    M3()
    M4()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    if ob and not om:
        print("NOTE: binding holds (min@pi) but R_A NOT globally monotone -> monotone route")
        print("      unavailable as-is; a weaker single-crossing/quasi-monotone lemma may still work.")

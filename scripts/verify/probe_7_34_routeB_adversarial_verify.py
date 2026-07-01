#!/usr/bin/env python3
r"""
probe_7_34_routeB_adversarial_verify.py
============================================================================
Independent ADVERSARIAL verification of the Route B certificate (BUG-009-D), plus a
new reduction lead: A=2 is the worst-case alphabet at the endpoint.

Route B claims (still BUG-009-D "blocked", handoff -- an unproven-but-asserted lead):
    g_A(s) <= 1 - (3/2) D(1-D)(1-cos s)   <=>   R_A(s) := [1-g_A(s)]/[D(1-D)(1-cos s)] >= 3/2
on the whole Gray region 0<D<=D_c(p), all s in (0,pi], all A, p in (0,1/2).

This probe re-derives g_A from an INDEPENDENT implementation (the pattern-quotient
replica Perron, verified elsewhere against full 2^n enumeration) and hunts for ANY
violation across a hostile grid -- an independent cross-check of the active lead,
NOT a re-run of Route B's own certificate script.

FINDINGS:
  * No violation: min over the hostile grid (A up to 8; p in {1e-3..0.49}; D up to
    0.9999 D_c; s across (0,pi]) is R_A(s) = 1.50005..., attained at the binding
    corner (A=2, p->0, D->D_c, s->pi). So 3/2 is the SHARP constant and the
    certificate holds numerically with the tightness exactly at that corner.
  * NEW LEAD (worst-case alphabet): R_A(pi) is INCREASING in A at the tight corner
    (1.5008, 1.522, 1.531, 1.535 for A=2,3,4,5 at p=0.05, D->D_c). So A=2 is the
    binding alphabet: min_A R_A(pi) = R_2(pi). Hence the general-A endpoint
    R_A(pi) >= 3/2 REDUCES to (a) A=2 is worst-case, R_A(pi) >= R_2(pi), plus
    (b) the PROVEN R_2(pi) >= 3/2 (probe_7_34m_endpoint_sharp_a2.py). This is a
    fresh "split smaller" path for the A>=3 endpoint that avoids the quartic Perron
    (route_f_monotone_reduction memory): prove monotonicity in A, not solve A>=3.

CHECKS:
  V1  R_A(s) >= 3/2 across the hostile grid (no violation); report min + corner.
  V2  worst-case alphabet: R_A(pi) increasing in A at the tight corner (min at A=2),
      so the A>=3 endpoint reduces to R_A(pi) >= R_2(pi) + proven R_2(pi)>=3/2.

Deps: numpy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~2-4 min.
"""
import itertools
import mpmath as mp
import numpy as np
mp.mp.dps = 50
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

def Dc(A, pv):
    def im(Dv):
        b0 = Dv / (A - 1); K = np.full((A, A), b0); np.fill_diagonal(K, 1 - Dv); Ki = np.linalg.inv(K)
        Tn = np.full((A, A), pv / (A - 1)); np.fill_diagonal(Tn, 1 - pv)
        By = lambda y: np.array([[Ki[y, a] * Tn[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-12, (A - 1) / A - 1e-7
    for _ in range(70):
        m = (lo + hi) / 2; lo, hi = (m, hi) if im(m) < 1e-11 else (lo, m)
    return lo

def R(A, p, D, s):
    t = 1 - mp.cos(mp.mpf(s))
    return (1 - gA(A, p, D, s)) / (D * (1 - D) * t)

def V1():
    print("-" * 78); print("V1  R_A(s) >= 3/2 across a hostile grid (independent; hunt for violations)")
    gmin = mp.inf; corner = None; viol = 0
    for A in (2, 3, 4, 5, 6, 8):
        for pv in ('0.001', '0.01', '0.05', '0.2', '0.4', '0.49'):
            pf = float(pv)
            if pf >= (A - 1) / A - 1e-6: continue
            dc = Dc(A, pf)
            for fr in ('0.05', '0.3', '0.6', '0.9', '0.99', '0.9999'):
                D = mp.mpf(fr) * dc
                for sd in (10, 45, 90, 135, 170, 179, 179.9):
                    r = R(A, pf, D, sd / mp.mpf(180) * mp.pi)
                    if r < gmin: gmin = r; corner = (A, pv, fr, float(sd))
                    if r < mp.mpf('1.5') - mp.mpf('1e-25'): viol += 1
    print(f"     global min R_A(s) = {mp.nstr(gmin, 10)}  at (A,p,D/Dc,s_deg)={corner}  violations={viol}")
    return rep("V1 R_A(s) >= 3/2 everywhere on the hostile grid (Route B holds)", viol == 0 and gmin >= mp.mpf('1.5') - mp.mpf('1e-6'))

def V2():
    print("-" * 78); print("V2  worst-case alphabet: R_A(pi) increasing in A (min at A=2) => reduce A>=3 endpoint")
    ok = True
    for pv in ('0.02', '0.05', '0.2'):
        pf = float(pv); row = []
        for A in (2, 3, 4, 5, 6):
            if pf >= (A - 1) / A - 1e-6: row.append(None); continue
            dc = Dc(A, pf); D = mp.mpf('0.999') * dc
            row.append(R(A, pf, D, mp.pi))
        vals = [v for v in row if v is not None]
        incr = all(vals[i + 1] >= vals[i] - mp.mpf('1e-12') for i in range(len(vals) - 1))
        minA2 = (vals[0] == min(vals))
        ok = ok and incr and minA2
        print(f"     p={pv}: R_A(pi)@0.999Dc = " + ", ".join(mp.nstr(v, 7) for v in vals) +
              f"  increasing_in_A:{incr}  min@A=2:{minA2}")
    return rep("V2 A=2 is worst-case alphabet at endpoint (min_A R_A(pi)=R_2(pi))", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("Independent adversarial verification of Route B (R_A(s)>=3/2) + worst-case-alphabet lead")
    print("=" * 78)
    V1(); V2()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("LEAD: general-A endpoint R_A(pi)>=3/2 reduces to [A=2 worst-case: R_A(pi)>=R_2(pi)]")
    print("+ [proven R_2(pi)>=3/2], avoiding the A>=3 quartic Perron. Next: prove A-monotonicity.")

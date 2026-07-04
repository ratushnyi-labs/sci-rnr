#!/usr/bin/env python3
r"""
probe_7_34_route_b_corner.py
============================================================================
ROUTE B step 7 -- the BINDING CORNER. The worst c in g_A<=1-c D(1-D)(1-cos s)
is at s=pi, D->0, p->0, value ~1.50-1.51 (slightly increasing in A). Pin down
the EXACT limit there: compute the small-D expansion of 1-g_A(pi) and of
(1-g_A(s)) at general s, leading order in D, as a function of (A,p,u=cos s).

If 1-g_A(s) = c_A(p,u) D(1-D)(1-cos s) + O(D^2(...)) with c_A(p,u) bounded
below by a clean constant (e.g. 3/2 at the corner), we get the floor.

Use mpmath high-precision Richardson in D at fixed small D to extract the
leading coefficient kappa(A,p,s) = (1-g_A(s))/(D(1-cos s)) as D->0, and study
its inf over s in (0,pi], p in (0,(A-1)/A), and over A.

Also: closed-form at s=pi. At s=pi, eta is REAL (arg=pi from cw3): eta(pi)<0.
Compute g_A(pi) symbolically in D,p,A via the pattern quotient (sympy) at
leading order in D.
"""
import itertools
import math
import numpy as np
import mpmath as mp
import sympy as sp
mp.mp.dps = 50

def pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def build_Q_num(A, pv, eta, etb):
    T = [[(1 - pv) if i == j else pv / (A - 1) for j in range(A)] for i in range(A)]
    st = list(itertools.product(range(A), repeat=3)); npat = 5 if A >= 3 else 4
    reps = [None] * npat
    for k, tr in enumerate(st):
        if reps[pat(tr)] is None: reps[pat(tr)] = k
    Q = mp.zeros(npat, npat)
    for a in range(npat):
        x, xp, xq = st[reps[a]]
        for (y, yp, yq) in st:
            t = T[xp][yp] * T[xq][yq] / T[x][y]
            if yp != y: t *= eta
            if yq != y: t *= etb
            Q[a, pat((y, yp, yq))] += t
    return Q

def eta_of(A, pv, Dv, sv):
    Dv = mp.mpf(Dv); sv = mp.mpf(sv); th = mp.log(Dv / ((A - 1) * (1 - Dv))); E = mp.e**(th + 1j * sv)
    den = (A - 1) * Dv * E + Dv - (A - 1)
    return ((A - 1) * Dv * E + Dv - (A - 1) * E) / den

def perron_mp(Q):
    ev, _ = mp.eig(Q); return max(ev, key=lambda e: abs(e))

def main():
    print("Leading coefficient kappa(A,p,s) = lim_{D->0} (1-g_A(s))/(D(1-cos s)).")
    print("(D(1-D)~D as D->0.) Scan s and p; report inf over s for each (A,p).")
    for A in (2, 3, 4, 5, 6):
        for p in (0.01, 0.05, 0.1, 0.2):
            if p >= (A - 1) / A - 1e-6: continue
            # extract leading coeff at very small D via two D values + Richardson
            kap_min = mp.inf; smin = None
            for s in [mp.pi * k / 16 for k in range(1, 17)]:
                vals = []
                for D in (mp.mpf('1e-7'), mp.mpf('5e-8')):
                    eta = eta_of(A, p, float(D), float(s)); etb = mp.conj(eta)
                    g = mp.re(perron_mp(build_Q_num(A, p, eta, etb)))
                    vals.append((1 - g) / (D * (1 - mp.cos(s))))
                # Richardson (linear in D): kappa = 2*v(D/2) - v(D) approx; here ratios already ~const
                kap = 2 * vals[1] - vals[0]
                if kap < kap_min:
                    kap_min = kap; smin = s
            print(f"  A={A} p={p}: inf_s kappa = {mp.nstr(kap_min, 8)} at s={mp.nstr(smin,5)} "
                  f"(={mp.nstr(smin/mp.pi,4)} pi)")
    print()
    print("=" * 72)
    print("At s=pi exactly, p->0, D->0: closed leading value of kappa(A).")
    for A in (2, 3, 4, 5, 6, 10, 20):
        rows = []
        for p in (mp.mpf('1e-3'), mp.mpf('1e-4')):
            D = mp.mpf('1e-8')
            eta = eta_of(A, p, float(D), float(mp.pi)); etb = mp.conj(eta)
            g = mp.re(perron_mp(build_Q_num(A, p, eta, etb)))
            kap = (1 - g) / (D * 2)  # 1-cos(pi)=2
            rows.append(kap)
        print(f"  A={A}: kappa(pi) at p->0 = {mp.nstr(rows[0],8)}, {mp.nstr(rows[1],8)} "
              f"(p=1e-3,1e-4)")

if __name__ == "__main__":
    main()

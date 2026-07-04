#!/usr/bin/env python3
r"""
probe_7_34_route_b_corner_limit.py
============================================================================
The BINDING CORNER mechanism (why the floor 3/2 is a finite limit). The 2nd-
order drop coefficient c2/(A-1)^2 blows up as p->0 (c2 ~ kappa^2 ~ 1/p^2), but
D_c -> 0 like p^2/4, so the PRODUCT c2 * D_c stays FINITE -- exactly the A=2
Psi_2(D_c)<5/16 mechanism generalized. We:
  (1) confirm c2(A,p) ~ C(A)/p^2 as p->0 (so c2 p^2 -> C(A) finite);
  (2) compute the exact corner product (c2/(A-1)^2)*D_c at p->0 -> finite L_A;
  (3) the leading-order floor prediction: 1-rho/(2Dt) = 1 - (c2/(A-1)^2)D + ...,
      and the WORST (smallest) value on Gray is at D=D_c, giving
      floor_lead(A) = 2*(1 - (c2/(A-1)^2) D_c) at p->0. Confirm ->3/2 region.
  (4) identify the closed form of c2: test c2 =? kappa^2 * f(A,p) /something, with
      kappa^2=(1-2p)^2/(p^2(1-p)^2) (the A=2 constant). Print c2 / kappa^2.
"""
import itertools
import math
import numpy as np
import mpmath as mp
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

def perron(Q):
    ev, _ = mp.eig(Q); return max(ev, key=lambda e: abs(e))

def c2_extract(A, p):
    h = mp.mpf('1e-4')
    def r(e1, e2): return mp.re(perron(build_Q_num(A, p, e1, e2)))
    return (r(h, h) - r(h, -h) - r(-h, h) + r(-h, -h)) / (4 * h * h)

def Dc_num(A, pv):
    def im(Dv):
        b0 = Dv / (A - 1); K = np.full((A, A), b0); np.fill_diagonal(K, 1 - Dv); Ki = np.linalg.inv(K)
        Tm = np.full((A, A), pv / (A - 1)); np.fill_diagonal(Tm, 1 - pv)
        By = lambda y: np.array([[Ki[y, a] * Tm[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-12, (A - 1) / A - 1e-7
    for _ in range(60):
        m = (lo + hi) / 2; lo, hi = (m, hi) if im(m) < 1e-10 else (lo, m)
    return lo

def kappa2(p):
    return (1 - 2 * p)**2 / (p**2 * (1 - p)**2)

def main():
    print("(1)&(4) c2 ~ C(A)/p^2; c2 / kappa^2 as p->0 (kappa^2=(1-2p)^2/(p^2(1-p)^2)):")
    for A in (2, 3, 4, 5):
        for p in (mp.mpf('0.1'), mp.mpf('0.02'), mp.mpf('0.005')):
            c2 = c2_extract(A, p)
            k2 = kappa2(p)
            print(f"  A={A} p={float(p)}: c2={mp.nstr(c2,8)}  c2*p^2={mp.nstr(c2*p*p,7)}  "
                  f"c2/kappa^2={mp.nstr(c2/k2,7)}")
        print()
    print("=" * 72)
    print("(2)&(3) corner: drop = (c2/(A-1)^2)*D_c at small p, and floor_lead=2(1-drop).")
    print("This is the leading-order prediction of the worst ratio (->3/2 region).")
    for A in (2, 3, 4, 5, 6):
        for p in (mp.mpf('0.05'), mp.mpf('0.01'), mp.mpf('0.003')):
            dc = mp.mpf(Dc_num(A, float(p)))
            c2 = c2_extract(A, p)
            drop = (c2 / (A - 1)**2) * dc
            print(f"  A={A} p={float(p)}: D_c={mp.nstr(dc,5)}  (c2/(A-1)^2)D_c={mp.nstr(drop,6)}  "
                  f"floor_lead=2(1-drop)={mp.nstr(2*(1-drop),7)}")
        print()

if __name__ == "__main__":
    main()

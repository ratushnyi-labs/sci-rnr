#!/usr/bin/env python3
r"""
probe_7_34_route_b_hostile_corner.py
============================================================================
ADVERSARIAL corner drill: the 2nd-order-TRUNCATED prediction floor_lead dipped
slightly below 3/2 at p=0.003, A>=5 (probe corner_limit). That is a truncation
artifact (higher orders help, since the binding is at the LARGEST D where the
2nd-order Taylor is least accurate). Verify the FULL (exact) g_A floor
(1-g_A)/(D(1-D)(1-cos s)) stays >= 3/2 at the hostile corner s=pi, D=D_c, tiny
p, large A. If it ever dips below 3/2, the constant must be lowered (still >0).
"""
import itertools
import math
import numpy as np
import mpmath as mp
mp.mp.dps = 60

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

def Cratio(A, pv, Dv, sv):
    Dv = mp.mpf(Dv); sv = mp.mpf(sv); th = mp.log(Dv / ((A - 1) * (1 - Dv))); E = mp.e**(th + 1j * sv)
    den = A * Dv - (A - 1)
    C = ((A - 1) * Dv * E + Dv - (A - 1)) / den
    C0 = ((A - 1) * Dv * mp.e**th + Dv - (A - 1)) / den
    return C / C0

def Dc_num(A, pv):
    def im(Dv):
        b0 = Dv / (A - 1); K = np.full((A, A), b0); np.fill_diagonal(K, 1 - Dv); Ki = np.linalg.inv(K)
        Tm = np.full((A, A), pv / (A - 1)); np.fill_diagonal(Tm, 1 - pv)
        By = lambda y: np.array([[Ki[y, a] * Tm[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-12, (A - 1) / A - 1e-7
    for _ in range(80):
        m = (lo + hi) / 2; lo, hi = (m, hi) if im(m) < 1e-11 else (lo, m)
    return lo

def perron(Q):
    ev, _ = mp.eig(Q); return max(ev, key=lambda e: abs(e))

def gA(A, p, D, s):
    P = abs(Cratio(A, p, D, s))**2
    rho = mp.re(perron(build_Q_num(A, p, eta_of(A, p, float(D), float(s)),
                                    mp.conj(eta_of(A, p, float(D), float(s))))))
    return P * rho

def main():
    print("FULL exact (1-g_A)/(D(1-D)t) at the hostile corner s=pi, D=D_c, tiny p, large A.")
    print("(2nd-order truncation dipped <3/2 at p=0.003,A>=5; check the EXACT value.)")
    worst = mp.inf; worst_at = None
    for A in (2, 3, 4, 5, 6, 8, 10, 12):
        for p in (mp.mpf('0.05'), mp.mpf('0.01'), mp.mpf('0.003'), mp.mpf('0.001'),
                  mp.mpf('0.0003')):
            if p >= mp.mpf(A - 1) / A: continue
            dc = mp.mpf(Dc_num(A, float(p)))
            if dc < 1e-13: continue
            # near-corner: D=Dc and also a few D below, s near pi and a couple below
            for fracD in (mp.mpf('1.0'), mp.mpf('0.95'), mp.mpf('0.7')):
                D = fracD * dc
                for sv in (mp.pi, mp.pi * mp.mpf('0.9'), mp.pi * mp.mpf('0.5')):
                    t = 1 - mp.cos(sv)
                    r = (1 - gA(A, p, D, sv)) / (D * (1 - D) * t)
                    if r < worst:
                        worst = r; worst_at = (A, float(p), float(fracD), float(sv))
            print(f"  A={A} p={float(p)}: D_c={mp.nstr(dc,4)} ratio(s=pi,D=Dc)="
                  f"{mp.nstr((1-gA(A,p,dc,mp.pi))/(dc*(1-dc)*2),9)}")
        print()
    print("=" * 72)
    print(f"GLOBAL worst exact ratio over hostile corner = {mp.nstr(worst,9)} at "
          f"A={worst_at[0]} p={worst_at[1]} fracD={worst_at[2]} s={worst_at[3]:.3f}")
    print(f"  >= 3/2 ?  {worst >= mp.mpf('1.5')}")
    print(f"  >= 1   ?  {worst >= mp.mpf('1.0')}  (any positive c works for GAP-1)")

if __name__ == "__main__":
    main()

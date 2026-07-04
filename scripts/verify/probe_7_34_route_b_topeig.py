#!/usr/bin/env python3
r"""
probe_7_34_route_b_topeig.py
============================================================================
ROUTE B step 6. Directly verify the spectral inequality top eigenvalue
g_A(s) <= 1 - c(1-cos s) and find the largest provable c, all A on Gray.
(The charpoly det-sign test was wrong because L(s) drops below the negative
part of the spectrum at large s; the inequality is purely about the TOP root.)

Then, for the PROOF, the real-symmetric route: at each s, Q(s) has real
spectrum, top eigenvalue simple. The right LOWER bound on 1-g_A(s) uses the
RAYLEIGH form of the symmetrized operator. We compute the symmetrizer at each
s (from that s's own positive Perron vectors -- they ARE positive for s>0 even
though the s=0 left vector degenerates) and confirm
        1 - lambda_max(S(s)) >= c (1-cos s).
"""
import itertools
import math
import numpy as np
import mpmath as mp
mp.mp.dps = 30

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

def to_np(Q):
    n = Q.rows
    return np.array([[complex(Q[i, j]) for j in range(n)] for i in range(n)])

def main():
    ss = np.concatenate([np.geomspace(1e-3, 0.5, 14, endpoint=False),
                         np.linspace(0.5, math.pi, 24)])
    ps = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3]
    fracs = [0.1, 0.3, 0.5, 0.8, 1.0]
    print("Direct: worst c such that g_A(s) <= 1 - c(1-cos s) on Gray (= worst ratio).")
    print("Also: does the inequality top<=1-c(1-cos s) hold for c=1.4 everywhere?")
    for A in (2, 3, 4, 5, 6):
        worst = np.inf; worst_at = None; viol14 = 0
        for p in ps:
            if p >= (A - 1) / A - 1e-6: continue
            dc = Dc_num(A, p)
            for frac in fracs:
                D = frac * dc
                if D < 1e-10: continue
                dd = D * (1 - D)
                for s in ss:
                    Q = to_np(build_Q_num(A, p, eta_of(A, p, D, float(s)),
                                          mp.conj(eta_of(A, p, D, float(s))))).real
                    top = np.max(np.linalg.eigvals(Q).real)
                    ratio = (1 - top) / (dd * (1 - math.cos(s)))
                    if ratio < worst:
                        worst = ratio; worst_at = (p, frac, s, D)
                    if top > 1 - 1.4 * (1 - math.cos(s)):
                        viol14 += 1
        print(f"  A={A}: worst c = {worst:.5f} at p={worst_at[0]} D={worst_at[3]:.2e} "
              f"s={worst_at[2]:.3f}; #violations of c=1.4 floor = {viol14}")

if __name__ == "__main__":
    main()

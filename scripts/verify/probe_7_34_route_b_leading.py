#!/usr/bin/env python3
r"""
probe_7_34_route_b_leading.py
============================================================================
ROUTE B step 8. The leading coefficient is EXACTLY 2:
   1 - g_A(s) = 2 D (1-cos s) + O(D^2) * (something(A,p,s)),
uniform in A, p, s. Confirm "exactly 2" to high precision via Richardson, and
isolate the O(D^2) coefficient's SIGN (it must be NEGATIVE so the finite-D
ratio drops from 2 toward 1.5, never above 2 -- giving the floor c=3/2 with
margin, OR if the 2nd-order term can push the ratio BELOW some c we need to
bound it).

CRITICAL for the proof: if 1-g_A(s) = 2 D(1-cos s) - (2nd order) and the
2nd-order piece is itself <= (1/2) * 2 D(1-cos s) on Gray (D<=D_c), then
1-g_A >= (3/2) D(1-cos s) and we are done with c=3/2.

So we need:  the O(D^2) deviation R(A,p,s,D) satisfies
   2 D(1-cos s) - (1-g_A(s)) <= (1/2) D(1-cos s)  on Gray,
i.e.  (1-g_A(s)) >= (3/2) D(1-cos s).  This is EXACTLY the worst-ratio>=1.5
finding. We now check it directly as the target inequality and measure margin,
AND check the leading-2 claim to pin the proof structure.
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

def perron_mp(Q):
    ev, _ = mp.eig(Q); return max(ev, key=lambda e: abs(e))

def gAval(A, p, D, s):
    eta = eta_of(A, p, float(D), float(s)); etb = mp.conj(eta)
    return mp.re(perron_mp(build_Q_num(A, p, eta, etb)))

def main():
    print("(L1) Leading coeff = lim_{D->0} (1-g_A)/(D(1-cos s)) -> 2 exactly?")
    print("     Richardson with D=1e-6,5e-7,2.5e-7 (account D(1-D), so use D not D(1-D)).")
    for A in (2, 3, 5):
        for p in (0.05, 0.2):
            for s in (mp.mpf('0.5'), mp.mpf('1.5'), mp.pi):
                Ds = [mp.mpf('1e-6'), mp.mpf('5e-7'), mp.mpf('2.5e-7')]
                rs = [(1 - gAval(A, p, D, s)) / (D * (1 - D) * (1 - mp.cos(s))) for D in Ds]
                # extrapolate linearly: r(D) ~ k0 + k1 D
                k0 = rs[2] + (rs[2] - rs[1])  # Richardson for halving
                print(f"  A={A} p={p} s={float(s):.3f}: ratios={[mp.nstr(r,9) for r in rs]} "
                      f"-> leading {mp.nstr(k0,9)}")
    print()
    print("=" * 72)
    print("(L2) TARGET inequality 1-g_A(s) >= (3/2) D(1-D)(1-cos s) on FULL Gray.")
    print("     min over fine grid of (1-g_A)/(D(1-D)(1-cos s)); must be >= 1.5.")
    ss = np.concatenate([np.geomspace(1e-3, 0.5, 18, endpoint=False),
                         np.linspace(0.5, math.pi, 30)])
    ps = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4]
    fracs = np.linspace(0.05, 1.0, 12)
    for A in (2, 3, 4, 5, 6, 8):
        worst = mp.inf; at = None
        for p in ps:
            if p >= (A - 1) / A - 1e-6: continue
            dc = Dc_num(A, p)
            for frac in fracs:
                D = mp.mpf(float(frac)) * dc
                if D < 1e-11: continue
                dd = D * (1 - D)
                for s in ss:
                    g = gAval(A, p, D, mp.mpf(float(s)))
                    r = (1 - g) / (dd * (1 - mp.cos(mp.mpf(float(s)))))
                    if r < worst:
                        worst = r; at = (p, float(frac), float(s))
        flag = "OK c>=3/2" if worst >= mp.mpf('1.5') - mp.mpf('1e-6') else "BELOW 3/2!"
        print(f"  A={A}: min ratio over Gray = {mp.nstr(worst,8)}  [{flag}] at "
              f"p={at[0]} frac={at[1]:.2f} s={at[2]:.3f}")

if __name__ == "__main__":
    main()

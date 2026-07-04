#!/usr/bin/env python3
r"""
probe_7_34_route_b_quadratic.py
============================================================================
ROUTE B step 9 -- generalize the A=2 quadratic-Perron + Jensen-floor proof.
A=2 (tex 12928): g(s)=1/2[F+sqrt(F^2+G)], F=|C_s/C0|^2|1+eta|^2 = 1-2D(1-D)t
(JENSEN FLOOR identity, t=1-cos s), G=4kappa^2|eta|^2|C_s/C0|^4>=0; then
1-g = (4(1-F)-G)/(2[(2-F)+sqrt(F^2+G)]) >= (11/16)D(1-D)t.

GENERALIZE. The affineness probe shows Q_A(eta,0) has RANK 2 -- so the Perron
likely sits in a 2x2 effective block, root of a QUADRATIC
   lambda^2 - tr_A lambda + det_A = 0,  g_A = 1/2[tr_A + sqrt(tr_A^2 - 4 det_A)].
Multiply by the prefactor |C_s/C0|^2. We test:
  (Q1) Does g_A(s)=rho(Q)/... equal the larger root of a quadratic whose
       coefficients we can read off? Check: is g_A a root of lambda^2 - b lambda
       + c with b,c the trace & "2x2 det" of the rank-2 part? Fit b,c from the
       two largest eigenvalues lam1=g_A, lam2 and verify lam1+lam2, lam1 lam2
       match a clean (A,p,D,s) form.
  (Q2) The JENSEN FLOOR for general A: is there an F_A with
       F_A = |prefactor|^2 * (trace-like) = 1 - 2 D(1-D) t  (same identity)?
       Equivalently does lam1+lam2 (the two Perron-block eigenvalues), times the
       prefactor, equal the Jensen floor 1-2D(1-D)t at leading order?
This pins whether the A=2 proof TRANSFERS structurally.
"""
import itertools
import math
import numpy as np
import mpmath as mp
mp.mp.dps = 40

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
    for _ in range(60):
        m = (lo + hi) / 2; lo, hi = (m, hi) if im(m) < 1e-10 else (lo, m)
    return lo

def main():
    print("(Q1) The eta/etb operator Q is NOT prefactor-normalized: g_A = |C_s/C0|^2 rho(Q).")
    print("     The 2 largest eigenvalues of Q (call mu1>=|mu2|). Test g_A=|C|^2 mu1 and")
    print("     whether mu1 is the larger root of a quadratic lam^2 - b lam + c with")
    print("     b = mu1+mu2, c = mu1 mu2, AND whether the OTHER eigenvalues are tiny (O(D^2)).")
    for A in (3, 4, 5):
        p = 0.2; D = 0.6 * Dc_num(A, p)
        for s in (0.5, 1.5, math.pi):
            Q = build_Q_num(A, p, eta_of(A, p, D, s), mp.conj(eta_of(A, p, D, s)))
            ev, _ = mp.eig(Q)
            evs = sorted(ev, key=lambda e: -abs(e))
            C2 = abs(Cratio(A, p, D, s))**2
            g = C2 * mp.re(evs[0])
            rest = [abs(e) for e in evs[2:]]
            print(f"  A={A} s={s:.2f}: mu=[{', '.join(mp.nstr(e,6) for e in evs[:3])}...]  "
                  f"|C|^2 mu1 = g_A = {mp.nstr(g,8)}  max|rest|={mp.nstr(max(rest) if rest else mp.mpf(0),4)}")
        print()

    print("=" * 72)
    print("(Q2) JENSEN FLOOR: is |C_s/C0|^2 (mu1+mu2) = 1 - 2D(1-D)t at LEADING order?")
    print("     (the A=2 identity F=1-2D(1-D)t is the source of the leading-2 coefficient.)")
    for A in (2, 3, 4, 5):
        p = 0.2
        for frac in (0.3, 1.0):
            D = frac * Dc_num(A, p)
            for s in (0.5, 1.5, math.pi):
                t = 1 - math.cos(s)
                Q = build_Q_num(A, p, eta_of(A, p, D, s), mp.conj(eta_of(A, p, D, s)))
                ev, _ = mp.eig(Q); evs = sorted(ev, key=lambda e: -abs(e))
                C2 = abs(Cratio(A, p, D, s))**2
                trace2 = C2 * (evs[0] + evs[1])  # prefactor times 2x2 trace
                jensen = 1 - 2 * D * (1 - D) * t
                print(f"  A={A} D={D:.4f} s={s:.2f}: |C|^2(mu1+mu2)={mp.nstr(mp.re(trace2),9)}  "
                      f"Jensen 1-2D(1-D)t={mp.nstr(jensen,9)}  diff={mp.nstr(mp.re(trace2)-jensen,3)}")
        print()

if __name__ == "__main__":
    main()

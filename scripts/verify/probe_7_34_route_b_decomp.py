#!/usr/bin/env python3
r"""
probe_7_34_route_b_decomp.py
============================================================================
ROUTE B step 10 -- decompose g_A = |C_s/C0|^2 * rho(Q) and locate the source
of the leading-2 coefficient, to build the general-A proof.

Plan. Write g_A(s) = P(s) * rho(s), P=|C_s/C0|^2, rho=rho(Q(eta,etb)).
1 - g_A = (1-P) + P(1-rho) = (1-P rho).
We expand BOTH factors in D (eta=O(D)). The known A=2 identity says the LEADING
piece is the prefactor-trace combination = Jensen floor. For general A, test the
SPLIT:
  (a) prefactor alone: 1 - P(s) = ?  (closed form from Cratio)
  (b) Perron defect:  1 - rho(s) = ?  (first/2nd order perturbation)
and see which carries the 2D(1-D)t.

Actually cleaner: P(s) = |C_s/C0|^2 has a CLOSED FORM. Compute it and its
small-D expansion. And rho(s) = 1 + (perturbation). The product's 1-g_A leading
term must be 2D(1-D)t. Decompose to see the additive split.

ALSO: test the GENERALIZED quadratic. Q's char poly may factor as
(stuff)*(lambda^2 - b lambda + c) with the Perron the root of the quadratic.
Read b=trace of 2x2 Perron block, c. Then g_A=|C|^2 * (1/2)(b+sqrt(b^2-4c)).
Check if |C|^2 * b and |C|^2^2 * (b^2-4c) have clean forms (the A=2 F and G).
"""
import itertools
import math
import numpy as np
import mpmath as mp
mp.mp.dps = 45

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

def perron_mp(Q):
    ev, _ = mp.eig(Q); return max(ev, key=lambda e: abs(e))

def main():
    print("Split 1-g_A = (1-P) + P(1-rho).  P=|C_s/C0|^2, rho=rho(Q).")
    print("Both as multiples of D(1-D)t -- which carries the leading 2?")
    for A in (2, 3, 5):
        p = 0.2
        for frac in (0.2, 1.0):
            D = frac * Dc_num(A, p)
            for s in (0.8, math.pi):
                t = 1 - math.cos(s); base = D * (1 - D) * t
                P = abs(Cratio(A, p, D, s))**2
                rho = mp.re(perron_mp(build_Q_num(A, p, eta_of(A, p, D, s), mp.conj(eta_of(A, p, D, s)))))
                g = P * rho
                print(f"  A={A} D={D:.4f} s={s:.2f}: (1-P)/base={float((1-P)/base):.4f}  "
                      f"P(1-rho)/base={float(P*(1-rho)/base):.4f}  (1-g)/base={float((1-g)/base):.4f}")
        print()

    print("=" * 72)
    print("eta small-D structure: eta_s ~ ? Confirm eta = -alpha(1-w) form (w=e^is).")
    print("A=2 had eta_s=D(1-D)(w-1)/((1-D)^2-D^2 w). General A leading: eta~ -(w-1)*?")
    A = 3; p = 0.2
    for D in (mp.mpf('1e-4'), mp.mpf('1e-5')):
        for s in (mp.mpf('0.5'), mp.mpf('1.5'), mp.pi):
            w = mp.e**(1j * s)
            eta = eta_of(A, p, float(D), float(s))
            # candidate leading: eta ~ D/(A-1)*(w-1)*const ... measure eta/(D(w-1))
            ratio = eta / (D * (w - 1))
            print(f"  A={A} D={mp.nstr(D,2)} s={float(s):.2f}: eta/(D(w-1)) = {mp.nstr(ratio,8)}")
    print()

    print("=" * 72)
    print("Char-poly factor: does Q factor with Perron in a 2x2 block? Print the two")
    print("largest eigs mu1,mu2 (NOT a +-pair for A>=3) and the next pair (mu3~-mu2-ish?).")
    A = 3; p = 0.2; D = Dc_num(A, p)
    for s in (0.5, 1.5, math.pi):
        Q = build_Q_num(A, p, eta_of(A, p, D, s), mp.conj(eta_of(A, p, D, s)))
        ev, _ = mp.eig(Q); evs = sorted(ev, key=lambda e: -abs(e))
        print(f"  A=3 s={s:.2f}: eigs = {[mp.nstr(e,5) for e in evs]}")
        # the 'paired' antisym structure: mu2,mu3 near +-x; mu1 perron; mu4,mu5 tiny
        print(f"          mu2+mu3={mp.nstr(evs[1]+evs[2],4)}  mu2*mu3={mp.nstr(evs[1]*evs[2],4)}")

if __name__ == "__main__":
    main()

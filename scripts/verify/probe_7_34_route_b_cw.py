#!/usr/bin/env python3
r"""
probe_7_34_route_b_cw.py
============================================================================
ROUTE B scratch: direct replica SPECTRAL INEQUALITY for general A via
Collatz-Wielandt. We want a POSITIVE test vector w (entrywise > 0) such that
the NORMALIZED replica operator M(s) = W_A(s)/gamma^2 satisfies, on the
pattern-quotient,
        rho(M(s)) <= max_i (M(s) w)_i / w_i  <= 1 - c (1-cos s)
with an explicit c>0 for all A on the Gray region. If found, GAP-1 closes
general-A achievability WITHOUT the curvature positivity c_A>0.

Subtlety: W_A(s) is COMPLEX (G_s has e^{is}), so Collatz-Wielandt does not
apply to W_A directly. But rho(W_A(s)) is REAL (the per-site replica rate;
g_A(s)=rho/gamma^2 is a real eigenvalue). The pattern-quotient operator Q in
the affineness probe (eta,etb parametrization, with etb=conj(eta)) is the right
object: it is NOT real either. We need a real nonnegative matrix with the same
Perron value to apply CW. Strategy: test CW on |Q(s)| (entrywise modulus) which
dominates rho(Q(s)) by Perron-Frobenius (rho(Q) <= rho(|Q|)). If rho(|Q(s)|) <=
1 - c(1-cos s) we are done a fortiori. Check whether that is even true (the
modulus bound may be too lossy). Otherwise work with the real symmetrized form.

This script DIAGNOSES which object to bound and tries test vectors.
"""
import itertools
import numpy as np
import mpmath as mp

# ---- operator builders (from probe_7_34l / probe_7_34m) --------------------
def Tmat(A, p):
    T = np.full((A, A), p / (A - 1)); np.fill_diagonal(T, 1 - p); return T

def Kinv(A, D):
    b0 = D / (A - 1); a0 = 1 - D - b0
    return (np.eye(A) - b0 * np.ones((A, A))) / a0

def lam_star(A, D):
    return math.log((1 - D) * (A - 1) / D)

import math

def Dc_num(A, pv):
    def im(Dv):
        Ki = Kinv(A, Dv); Tn = Tmat(A, pv)
        By = lambda y: np.array([[Ki[y, a] * Tn[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0))
        t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-12, (A - 1) / A - 1e-7
    for _ in range(70):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if im(mid) < 1e-11 else (lo, mid)
    return lo

def gamma_dual(A, p, D, n=4):
    lam = lam_star(A, D)
    ws = list(itertools.product(range(A), repeat=n))
    T = Tmat(A, p)
    PX = {}
    for w in ws:
        v = 1.0 / A
        for k in range(1, n):
            v *= T[w[k - 1], w[k]]
        PX[w] = v
    Ki = Kinv(A, D)
    PY = {}
    for y in ws:
        acc = 0.0
        for x in ws:
            f = 1.0
            for k in range(n):
                f *= Ki[y[k], x[k]]
            acc += f * PX[x]
        PY[y] = acc
    def dH(x, y): return sum(1 for a, b in zip(x, y) if a != b)
    ratios = []
    for x in ws:
        M = sum(PY[y] * math.exp(-lam * dH(x, y)) for y in ws)
        ratios.append(M / PX[x])
    return float(np.mean(ratios)) ** (1.0 / n)

def W_operator(A, p, D, s):
    T = Tmat(A, p); Ki = Kinv(A, D); nu = math.exp(-lam_star(A, D))
    w = np.exp(1j * s)
    def G(x, xp):
        return sum(Ki[y, xp] * (1.0 if y == x else nu * w) for y in range(A))
    Gm = np.array([[G(x, xp) for xp in range(A)] for x in range(A)])
    states = list(itertools.product(range(A), repeat=3))
    SI = {st: i for i, st in enumerate(states)}; nS = len(states)
    W = np.zeros((nS, nS), complex)
    for (xa, xb, xc) in states:
        i = SI[(xa, xb, xc)]
        for (x, xp, xq) in states:
            j = SI[(x, xp, xq)]
            W[i, j] = (1 / T[xa, x]) * T[xb, xp] * T[xc, xq] * Gm[x, xp] * np.conj(Gm[x, xq])
    return W

# pattern-quotient (5x5 A>=3, 4x4 A=2), A-independent dim, carries Perron
def pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def quotient_from_W(A, p, D, s):
    """build the pattern-quotient of the NORMALIZED W/gamma^2 directly."""
    T = Tmat(A, p); Ki = Kinv(A, D); nu = math.exp(-lam_star(A, D))
    w = np.exp(1j * s)
    def G(x, xp):
        return sum(Ki[y, xp] * (1.0 if y == x else nu * w) for y in range(A))
    Gm = np.array([[G(x, xp) for xp in range(A)] for x in range(A)])
    states = list(itertools.product(range(A), repeat=3))
    npat = 5 if A >= 3 else 4
    reps = [None] * npat
    for k, tr in enumerate(states):
        if reps[pat(tr)] is None:
            reps[pat(tr)] = k
    Q = np.zeros((npat, npat), complex)
    for a in range(npat):
        xa, xb, xc = states[reps[a]]
        for (x, xp, xq) in states:
            val = (1 / T[xa, x]) * T[xb, xp] * T[xc, xq] * Gm[x, xp] * np.conj(Gm[x, xq])
            Q[a, pat((x, xp, xq))] += val
    return Q

def perron(M):
    ev = np.linalg.eigvals(M)
    return ev[np.argmax(np.abs(ev))]

# ---------------------------------------------------------------------------
def main():
    print("Verify quotient Perron == W Perron, and g_A=rho/gamma^2 real:")
    for A in (2, 3, 4, 5):
        for p in (0.1, 0.2):
            D = 0.7 * Dc_num(A, p)
            for s in (0.5, 1.5, math.pi):
                W = W_operator(A, p, D, s)
                Q = quotient_from_W(A, p, D, s)
                rW = perron(W); rQ = perron(Q)
                print(f"  A={A} p={p} s={s:.2f}: rho(W)={rW.real:.6f}{rW.imag:+.2e}i "
                      f"rho(Q)={rQ.real:.6f}{rQ.imag:+.2e}i  |diff|={abs(rW-rQ):.1e}")
        print()

    print("=" * 70)
    print("Now test: is rho(Q(s))/gamma^2 = g_A(s) and the floor ratio:")
    ss = np.concatenate([np.geomspace(1e-2, 0.5, 15, endpoint=False),
                         np.linspace(0.5, math.pi, 30)])
    for A in (2, 3, 4, 5, 6):
        worst = np.inf
        for p in (0.02, 0.05, 0.1, 0.2, 0.3):
            if p >= (A - 1) / A - 1e-6:
                continue
            for frac in (0.3, 0.6, 1.0):
                D = frac * Dc_num(A, p)
                if D < 1e-9: continue
                gam = gamma_dual(A, p, D)
                dd = D * (1 - D)
                for s in ss:
                    Q = quotient_from_W(A, p, D, float(s))
                    g = perron(Q).real / gam ** 2
                    ratio = (1 - g) / (dd * (1 - math.cos(s)))
                    worst = min(worst, ratio)
        print(f"  A={A}: worst (1-g)/(D(1-D)(1-cos s)) over Gray = {worst:.4f}")

if __name__ == "__main__":
    main()

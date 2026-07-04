#!/usr/bin/env python3
r"""
probe_7_34_route_b_cw2.py
============================================================================
ROUTE B step 2: find the RIGHT real nonnegative object for Collatz-Wielandt,
and a test vector w giving g_A(s) <= 1 - c(1-cos s).

KEY STRUCTURAL FACT to exploit. g_A(s) = rho(Q(s))/gamma^2. There is a
well-defined gamma^2-normalization built into the eta/etb parametrization of
the affineness probe (Q_A(eta,etb) there ALREADY equals W/gamma^2 in the
pattern quotient: g_A(0)=1, i.e. Q_A(eta_0,etb_0) has Perron exactly 1 at s=0).
So work with the eta/etb operator from probe_7_34m: at s=0 it has Perron 1.

We seek a real nonneg matrix R(s) with rho(Q(s)) <= rho(R(s)) and a positive
w with (R w)_i/w_i <= 1 - c(1-cos s). Candidates:
  (1) Entrywise modulus |Q(s)|: rho(Q) <= rho(|Q|) always (PF). Test if the
      modulus inequality is tight enough (it may overshoot 1).
  (2) The eta/etb operator at s=0 (REAL, Perron 1) plus a real curvature bound:
      since Q(s) is analytic in u=cos s on the symmetric (etb=conj eta) locus,
      maybe a real similarity makes it real.

This script: (a) build Q via eta/etb (matches probe_7_34m exactly, Perron-1 at
s=0); (b) compute rho(|Q(s)|) and the floor ratio; (c) try CW test vectors on
Q(s) treated via the real 2nd-moment: actually rho(Q) is real so we can use the
variational form for the dominant eigenvalue of the symmetric part? No -- non-normal.
Instead: directly compute the CW bound max_i Re[(Q w)_i/w_i] is NOT a bound.

Correct tool: rho(Q) is a SIMPLE real eigenvalue (Perron of the s=0 limit,
continued). For a matrix with a simple real dominant eigenvalue lambda and
positive right/left Perron vectors at the base point, lambda admits the
Rayleigh-type bound via the SIMILARITY to a real matrix when Q is similar to a
real nonneg matrix. CHECK: is Q(s) similar (by a real positive diagonal /
the s=0 Perron vectors) to a real matrix? Test the "DB-symmetrized" form
S = diag(l)^{1/2-like}... Let's just empirically test several real surrogates.
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

def perron_mp(Q):
    ev, _ = mp.eig(Q)
    return max(ev, key=lambda e: abs(e))

def to_np(Q):
    n = Q.rows
    return np.array([[complex(Q[i, j]) for j in range(n)] for i in range(n)])

def perron_np(M):
    ev = np.linalg.eigvals(M)
    return ev[np.argmax(np.abs(ev))]

def main():
    print("Check Q(eta,etb) has Perron 1 at s=0 (so g_A = rho(Q) directly, no gamma):")
    for A in (2, 3, 4, 5):
        for p in (0.1, 0.2):
            D = 0.7 * Dc_num(A, p)
            eta0 = eta_of(A, p, D, 0.0); etb0 = mp.conj(eta0)
            r0 = perron_mp(build_Q_num(A, p, eta0, etb0))
            print(f"  A={A} p={p}: rho(Q) at s=0 = {mp.nstr(r0, 12)}")
    print()

    print("=" * 72)
    print("Floor ratio with g_A = rho(Q(s)) directly (NO gamma normalization):")
    ss = np.concatenate([np.geomspace(1e-2, 0.5, 12, endpoint=False),
                         np.linspace(0.5, math.pi, 25)])
    for A in (2, 3, 4, 5, 6):
        worst = np.inf; worst_at = None
        for p in (0.02, 0.05, 0.1, 0.2, 0.3):
            if p >= (A - 1) / A - 1e-6: continue
            for frac in (0.3, 0.6, 1.0):
                D = frac * Dc_num(A, p)
                if D < 1e-9: continue
                dd = D * (1 - D)
                for s in ss:
                    eta = eta_of(A, p, D, float(s)); etb = mp.conj(eta)
                    g = float(mp.re(perron_mp(build_Q_num(A, p, eta, etb))))
                    ratio = (1 - g) / (dd * (1 - math.cos(s)))
                    if ratio < worst:
                        worst = ratio; worst_at = (p, frac, s)
        print(f"  A={A}: worst (1-g)/(D(1-D)(1-cos s)) = {worst:.4f} at p={worst_at[0]} "
              f"frac={worst_at[1]} s={worst_at[2]:.3f}")
    print()

    print("=" * 72)
    print("Collatz-Wielandt surrogate test. rho(Q) is REAL & simple (Perron-1 at s=0).")
    print("Try real-nonneg surrogate R = |Q| (entrywise modulus); CW with several w:")
    for A in (2, 3, 4, 5):
        p = 0.2; D = Dc_num(A, p)
        # s=0 Perron right/left vectors (real)
        eta0 = eta_of(A, p, D, 0.0); etb0 = mp.conj(eta0)
        Q0 = to_np(build_Q_num(A, p, eta0, etb0))
        ev0, vr0 = np.linalg.eig(Q0)
        k0 = np.argmax(np.abs(ev0))
        r0_vec = np.abs(vr0[:, k0].real); r0_vec /= r0_vec.sum()
        evl0, vl0 = np.linalg.eig(Q0.T)
        l0_vec = np.abs(vl0[:, np.argmax(np.abs(evl0))].real); l0_vec /= l0_vec.sum()
        for s in (0.3, 1.0, 2.0, math.pi):
            eta = eta_of(A, p, D, s); etb = mp.conj(eta)
            Q = to_np(build_Q_num(A, p, eta, etb))
            absQ = np.abs(Q)
            rho_absQ = perron_np(absQ).real
            rho_Q = perron_np(Q).real
            # CW upper bounds for absQ with several w
            def cw(M, w):
                return float(np.max((M @ w) / w))
            ones = np.ones(Q.shape[0])
            cw_ones = cw(absQ, ones)
            cw_r0 = cw(absQ, r0_vec)
            cw_l0 = cw(absQ, l0_vec)
            dd = D * (1 - D); base = dd * (1 - math.cos(s))
            print(f"  A={A} s={s:.2f}: rho(Q)={rho_Q:.5f} rho(|Q|)={rho_absQ:.5f}  "
                  f"CW(|Q|,1)={cw_ones:.5f} CW(|Q|,r0)={cw_r0:.5f} CW(|Q|,l0)={cw_l0:.5f}")
            print(f"          (1-rho(|Q|))/base = {(1-rho_absQ)/base:.3f}  "
                  f"(1-CW_r0)/base={(1-cw_r0)/base:.3f}")
        print()

if __name__ == "__main__":
    main()

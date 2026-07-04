#!/usr/bin/env python3
r"""
probe_7_34_route_b_sym.py
============================================================================
ROUTE B step 4 -- the SYMMETRIZATION route (replaces CW, which the phase-
cancellation kills). DISCOVERY (cw3): the spectrum of Q(s) is ENTIRELY REAL
for all A,p,D,s on the Gray symmetric locus. So Q(s) is similar to a real
matrix; if it is similar to a real SYMMETRIC S(s) by a fixed positive diagonal
(reversibility), then
        g_A(s) = rho(Q(s)) = lambda_max(S(s)),
and the spectral inequality g_A(s) <= 1 - c(1-cos s) becomes the PSD statement
        (1 - c(1-cos s)) I  -  S(s)  >=  0     (Loewner),
which is a clean, A-uniform, checkable certificate (no CW test vector needed,
no curvature c_A>0 expansion needed).

This script:
  (S1) Confirm the spectrum of Q(s) is real (all A) -- the enabling fact.
  (S2) Find a positive diagonal D_diag with D_diag^{-1} Q D_diag SYMMETRIC
       (detailed balance). Test via the s=0 right Perron r0 and left l0:
       the canonical symmetrizer for a matrix with positive Perron data is
       D_sym = diag(sqrt(r0_i / l0_i)) (when Q is reversible w.r.t. l0).
       Check symmetry defect of D_sym^{-1} Q D_sym.
  (S3) If symmetrizable, compute lambda_max(S(s)) == rho(Q) and the floor ratio
       (1 - lambda_max)/(D(1-D)(1-cos s)) for all A -- same ~1.5 number.
  (S4) The PSD certificate: smallest eigenvalue of (1 - c*(1-cos s)) I - S(s)
       for c = 1 (a safe floor below the worst ratio 1.5); check >= 0 for all
       A,p,D,s on Gray. If yes -> g_A <= 1 - (1-cos s) -> GAP-1 closes general A.
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
    print("(S2) Symmetrizability test. Build D_sym from s=0 Perron r0,l0; check")
    print("     symmetry defect of D_sym^{-1} Q(s) D_sym at general s.")
    for A in (2, 3, 4, 5):
        p = 0.2; D = Dc_num(A, p)
        Q0 = to_np(build_Q_num(A, p, eta_of(A, p, D, 0.0), mp.conj(eta_of(A, p, D, 0.0))))
        ev0, vr = np.linalg.eig(Q0); k0 = np.argmax(ev0.real)
        r0 = vr[:, k0].real; r0 = r0 / r0[0] * np.sign(r0[0])
        evl, vl = np.linalg.eig(Q0.T); kl = np.argmax(evl.real)
        l0 = vl[:, kl].real; l0 = l0 / l0[0] * np.sign(l0[0])
        # reversibility symmetrizer: S = diag(d) Q diag(d)^{-1} symmetric where d_i^2 = l0_i/r0_i
        # (standard: if l0_i Q_ij = ... ) -- test the canonical d_i = sqrt(l0_i/r0_i)
        ok_rev = np.all(l0 > 0) and np.all(r0 > 0)
        if ok_rev:
            d = np.sqrt(np.abs(l0 / r0))
            for s in (0.5, 1.0, 2.0, math.pi):
                Q = to_np(build_Q_num(A, p, eta_of(A, p, D, s), mp.conj(eta_of(A, p, D, s)))).real
                S = np.diag(d) @ Q @ np.diag(1.0 / d)
                defect = np.max(np.abs(S - S.T))
                lam_max = np.max(np.linalg.eigvalsh((S + S.T) / 2))
                print(f"  A={A} s={s:.2f}: sym-defect(D_sym^-1 Q D_sym)={defect:.2e}  "
                      f"lam_max(sym part)={lam_max:.6f}")
        else:
            print(f"  A={A}: r0 or l0 not all positive (r0={r0}, l0={l0})")
        print()

if __name__ == "__main__":
    main()

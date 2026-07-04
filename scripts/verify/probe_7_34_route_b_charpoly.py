#!/usr/bin/env python3
r"""
probe_7_34_route_b_charpoly.py
============================================================================
ROUTE B step 5 -- the CHARACTERISTIC-POLYNOMIAL route. Since the spectrum of
Q(s) is REAL (cw3) and the Perron value g_A(s)=rho(Q(s)) is the LARGEST real
eigenvalue, the spectral inequality
        g_A(s) <= L(s) := 1 - c (1 - cos s)
is EQUIVALENT to: the characteristic polynomial chi(lambda)=det(lambda I-Q(s))
has NO root in (L(s), +infty), i.e. chi(L(s)) and the leading sign agree AND
chi is monotone increasing beyond L(s). Concretely, since chi(lambda) ->
+infty and the largest root is g_A, we have g_A <= L(s) IFF chi(lambda) > 0 for
all lambda >= L(s). A SUFFICIENT, checkable condition: chi(L(s)) > 0 AND
chi'(L(s)) > 0 AND chi is convex-increasing on [L(s),infty) (Newton: if
chi(L)>0 and chi'(L)>0 and all higher derivs of one sign... or simply: L(s) >
g_A iff chi(L)>0 given L is to the right of all roots).

PLAN. Work with the pattern quotient (<=5x5). Use the s=0 deflation: at s=0,
chi_0(lambda) = lambda^{n-1}(lambda-1) (Perron 1, rest 0 -- the gap-1 fact). At
s>0, eta=O(D) perturbs. The cleanest CERTIFICATE for g_A <= L(s):
        det(L(s) I - Q(s)) > 0  AND  L(s) lies to the right of ALL eigenvalues.
The second clause: since the NEXT eigenvalues are O(D) (bounded away below L for
small D) we check rho_2 := 2nd largest |eig| < L(s). Then det(L I - Q)>0 with L
above the 2nd eigenvalue forces the top eigenvalue < L (sign-count argument).

We test, for c=1 and also the empirical c, across the Gray region all A:
  (P1) g_A(s) = top real eigenvalue; floor ratio (recap).
  (P2) det(L(s) I - Q(s)) with L=1-(1-cos s): sign (>0 needed), all A.
  (P3) 2nd-largest eigenvalue < L(s) (so det>0 => top<L): all A.
  (P4) MONOTONE certificate: define h(lambda)=det(lambda I - Q(s)). Verify
       h(L)>0 and the largest root g_A < L by direct comparison (sanity).
The pair (P2 & P3) is a rigorous sufficient condition; we measure its margin.
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
    ss = np.concatenate([np.geomspace(1e-2, 0.5, 10, endpoint=False),
                         np.linspace(0.5, math.pi, 18)])
    ps = [0.02, 0.05, 0.1, 0.2, 0.3]
    fracs = [0.2, 0.5, 0.8, 1.0]
    for c in (1.0, 1.4):
        print("=" * 72)
        print(f"Certificate with c={c}:  L(s)=1-c(1-cos s)")
        for A in (2, 3, 4, 5, 6):
            min_det = np.inf; min_gap = np.inf; cert_ok = True; worst_ratio = np.inf
            for p in ps:
                if p >= (A - 1) / A - 1e-6: continue
                dc = Dc_num(A, p)
                for frac in fracs:
                    D = frac * dc
                    if D < 1e-9: continue
                    dd = D * (1 - D)
                    for s in ss:
                        L = 1 - c * (1 - math.cos(s))
                        Q = to_np(build_Q_num(A, p, eta_of(A, p, D, float(s)),
                                              mp.conj(eta_of(A, p, D, float(s))))).real
                        ev = np.sort(np.linalg.eigvals(Q).real)[::-1]
                        top = ev[0]; second = ev[1]
                        det_LQ = np.prod(L - ev)  # = det(L I - Q)
                        ratio = (1 - top) / (dd * (1 - math.cos(s)))
                        worst_ratio = min(worst_ratio, ratio)
                        # sufficient cert: second < L  AND det(L I - Q) > 0  => top < L
                        gap = L - second
                        ok = (gap > 0) and (det_LQ > 0)
                        cert_ok = cert_ok and ok and (top <= L)
                        min_det = min(min_det, det_LQ)
                        min_gap = min(min_gap, gap)
            print(f"  A={A}: worst (1-g)/base={worst_ratio:.4f}  min det(L-Q)={min_det:.3e}  "
                  f"min(L-2nd eig)={min_gap:.4f}  cert(top<=L all)={cert_ok}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
r"""
probe_7_34_route_b_cw3.py
============================================================================
ROUTE B step 3. The modulus surrogate |Q| FAILS (rho(|Q|)>1): the decay
g_A<1 is a PHASE-CANCELLATION effect, invisible to entrywise-positive CW.

So we need either (a) a REAL form of Q(s) whose Perron equals rho(Q) and is a
genuine sub-stochastic-type contraction, or (b) abandon pure CW and prove the
direct inequality 1 - g_A(s) >= c(1-cos s) via the analytic structure.

This script explores:
  (1) The symmetric locus etb = conj(eta). Is Q(s) similar to a REAL matrix
      via a fixed (s-independent or simple) similarity? Check if D^{-1} Q D is
      real for D = diag built from s=0 Perron data, or if Q is "essentially
      real" up to a unitary phase.
  (2) The s=0 LEFT Perron vector l0 and RIGHT r0: form the rank-1 deflation and
      look at the SECOND-order (in s) behavior of rho directly via perturbation:
        rho(s) = 1 + l0^T Q'(0) r0 * (deriv) ... but rho(0)=1, rho'(0)=0 (real
      even function of s), and rho''(0) = the curvature = -2 vbar. We want a
      LOWER bound on -rho''-type but globally.
  (3) The HERMITIAN part route: for the dominant eigenvalue of a matrix with
      simple real Perron, is rho(Q) <= rho_max of the field-of-values? No.
      Instead use: g_A(s) = (1/A) * something? Check the 2nd-moment quadratic
      form: g_A relates to E|Phi|^2. Actually rho(Q(s)) = lim (E2)^{1/n}; and
      E2 = sum over configs of a GIBBS weight. The Gibbs/transfer structure
      might give 1-g via a VARIANCE that is manifestly >= c(1-cos s).
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

def perron_np(M):
    ev = np.linalg.eigvals(M)
    return ev[np.argmax(np.abs(ev))]

def main():
    print("Examine the eta(s) structure: eta and etb=conj(eta) on the symmetric locus.")
    print("Print eta for a few s; note |eta| and arg, and the structure of Q entries.")
    A = 3; p = 0.2; D = Dc_num(A, p)
    for s in (0.0, 0.5, 1.0, 2.0, math.pi):
        eta = eta_of(A, p, D, s)
        print(f"  s={s:.3f}: eta={complex(eta):.5f}  |eta|={abs(complex(eta)):.5f}  "
              f"arg={math.atan2(eta.imag, eta.real):.4f}")
    print()

    print("=" * 72)
    print("Try: is Q(s) similar to a REAL matrix via diag phase? Check if there is")
    print("a diagonal U=diag(e^{i theta_k}) with U^{-1} Q U real. Compute the")
    print("'realness defect' after optimal diagonal phase alignment per row/col.")
    # A simpler invariant: rho(Q) is real. Check eigenvalues of Q -- are they all
    # real, or complex? If the SPECTRUM is real, Q is similar to a real matrix.
    for A in (2, 3, 4, 5):
        p = 0.2; D = Dc_num(A, p)
        for s in (0.5, 1.0, 2.0, math.pi):
            eta = eta_of(A, p, D, s); etb = mp.conj(eta)
            Q = to_np(build_Q_num(A, p, eta, etb))
            ev = np.linalg.eigvals(Q)
            max_im = np.max(np.abs(ev.imag))
            ev_sorted = ev[np.argsort(-np.abs(ev))]
            print(f"  A={A} s={s:.2f}: spectrum max|Im(eig)|={max_im:.2e}  "
                  f"eigs={[f'{e.real:.4f}{e.imag:+.1e}i' for e in ev_sorted]}")
        print()

if __name__ == "__main__":
    main()

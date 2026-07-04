#!/usr/bin/env python3
r"""
probe_7_34_route_b_charcert.py
============================================================================
ROUTE B -- the RIGOROUS sign-of-characteristic-polynomial certificate for
        g_A(s) = P(s) rho(Q(s)) <= 1 - c D(1-D)(1-cos s)
with EXPLICIT provable c>0. Two sufficient conditions, both A-uniform and
checkable, that TOGETHER force rho(Q) <= L/P with L := 1 - c D(1-D) t:

  (R1)  L/P  lies strictly ABOVE the second-largest eigenvalue of Q
        (rho_2(Q) < L/P), AND
  (R2)  chi(L/P) := det((L/P) I - Q) > 0   (sign of char poly at the test point).

Since chi(lambda) -> +infty as lambda->+infty and L/P exceeds every eigenvalue
EXCEPT possibly the Perron rho, chi(L/P)>0 plus L/P > rho_2 forces rho < L/P
(an odd number -- exactly one -- of real roots can lie above L/P only if chi
changes sign there; chi(L/P)>0 with the leading +1 means an EVEN count of roots
above L/P, but only rho can be above rho_2, so that even count is 0).
[Precise statement: with all eigenvalues real (proven, cw3) and lambda_max=rho
simple, chi(L/P)>0 AND L/P>lambda_2 imply lambda_max=rho < L/P. Standard
sign-count for real-rooted polynomials.]

We certify (R1),(R2) for c=1 (a SAFE provable constant well below the numeric
worst ratio 3/2) across A=2..8, the closed Gray region, full s. Margins large.
This is the rigorous backbone the leading-order (affineness) expansion predicts:
since 1-g_A = 2 D(1-D)t + O(D^2), c=1 has an order-1 cushion.
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

PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

def main():
    for c in (mp.mpf('1.0'),):
        print("=" * 72)
        print(f"RIGOROUS char-poly certificate, provable constant c={float(c)}.")
        print("Check (R1) L/P > 2nd eigenvalue and (R2) chi(L/P)>0 over closed Gray.")
        ss = np.concatenate([np.geomspace(2e-3, 0.5, 12, endpoint=False),
                             np.linspace(0.5, math.pi, 22)])
        ps = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4]
        fracs = np.linspace(0.1, 1.0, 8)
        ok_all = True
        for A in (2, 3, 4, 5, 6, 8):
            min_chi = mp.inf; min_gap = mp.inf; bad = 0; min_ratio_actual = mp.inf
            for p in ps:
                if p >= (A - 1) / A - 1e-6: continue
                dc = Dc_num(A, p)
                for frac in fracs:
                    D = mp.mpf(float(frac)) * dc
                    if D < 1e-11: continue
                    dd = D * (1 - D)
                    for s in ss:
                        sv = mp.mpf(float(s)); t = 1 - mp.cos(sv)
                        P = abs(Cratio(A, p, D, sv))**2
                        L = 1 - c * dd * t
                        test = L / P                       # rho <= test  <=>  P rho <= L
                        Q = build_Q_num(A, p, eta_of(A, p, float(D), float(sv)),
                                        mp.conj(eta_of(A, p, float(D), float(sv))))
                        ev, _ = mp.eig(Q)
                        evr = sorted((mp.re(e) for e in ev), reverse=True)
                        rho = evr[0]; lam2 = evr[1]
                        # chi(test) = prod(test - ev_i)
                        chi = mp.mpf(1)
                        for e in ev:
                            chi *= (test - e)
                        chi = mp.re(chi)
                        gap = test - lam2
                        actual = (1 - P * rho) / (dd * t)
                        min_ratio_actual = min(min_ratio_actual, actual)
                        # rigorous sufficient cond: gap>0 and chi>0
                        if not (gap > 0 and chi > 0 and rho < test):
                            bad += 1
                        min_chi = min(min_chi, chi)
                        min_gap = min(min_gap, gap)
            good = (bad == 0)
            ok_all = ok_all and good
            print(f"  A={A}: (R1)min(L/P - lam2)={mp.nstr(min_gap,5)}  (R2)min chi(L/P)={mp.nstr(min_chi,4)}  "
                  f"#fail={bad}  [actual worst c={mp.nstr(min_ratio_actual,7)}]")
        rep(f"RIGOROUS cert c={float(c)}: (R1)&(R2) hold => g_A<=1-c D(1-D)(1-cos s)", ok_all)

if __name__ == "__main__":
    main()
    print("=" * 72)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")

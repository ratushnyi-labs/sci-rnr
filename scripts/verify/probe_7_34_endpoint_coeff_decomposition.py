#!/usr/bin/env python3
r"""
probe_7_34_endpoint_coeff_decomposition.py
============================================================================
BUG-009-D / Route B: the small-p endpoint coefficient (A-2)/(A-1) DECOMPOSED into
two independently-derivable first-order perturbation coefficients.

CONTEXT. probe_7_34_endpoint_smallp_asymptotic.py established the leading margin
    R_A(pi)|_{D=D_c} - 3/2  ~  (A-2)/(A-1) * p     (A>=3; A=2 drops to p^2/8).
This probe SPLITS that coefficient into smaller, separately-provable pieces (the
"split the gap smaller" step toward a closed-form proof).

THE DECOMPOSITION. Write the two expansions at s=pi, p->0:
    slope:      R_A(pi; D) = 2 - E_A(p) D + O(D^2),
                E_A(p) = [2(A-1)/p^2] * (1 + eps_1(A) p + O(p^2)),
    threshold:  D_c(A,p) = [p^2/(4(A-1))] * (1 + delta_1(A) p + O(p^2)).
Then E*D_c -> 1/2 forces the leading margin structure, and
    R_A(pi)|_{D_c} - 3/2 = -[(eps_1 + delta_1)/2] p + O(p^2)
(the O(D^2) term contributes only at O(p^2): D_c^2 ~ p^4 against a 1/p^2-scale
coefficient -- consistent with A=2, where eps_1+delta_1=0 and the margin is p^2/8).

MEASURED COEFFICIENTS (all clean rationals; Richardson in D then in p, 60 dps):
    delta_1(A) = 2                      UNIVERSAL in A  (D_c = p^2/(4(A-1)) (1+2p+...))
    eps_1(A)   = -2 - 2 (A-2)/(A-1)     ( -2, -3, -10/3, -7/2, -18/5 for A=2..6 )
    identity:  -(eps_1 + delta_1)/2 = (A-2)/(A-1)   EXACT for every A
               (vanishes at A=2 => the O(p^2) tight case).

WHY THIS IS PROGRESS. Each piece is a FIRST-ORDER perturbation coefficient:
  * delta_1 = 2 comes from the Gray-threshold discriminant expansion (at A=2 it is
    read off the closed form D_c = 1/2 - sqrt(1-2p)/(2(1-p)) = (p^2/4)(1+2p+...);
    its A-universality says the threshold's relative p-correction is alphabet-free);
  * eps_1 is the O(p) correction to the O(D) slope of the replica Perron at s=pi --
    a first-order eigenvalue-perturbation object at D=0 (NO root-solving; the
    quartic wall of the exact A>=3 Perron is irrelevant to it).
A closed-form proof of the two coefficients would upgrade the A>=3 leading-order
endpoint closure from numerics to proof; both are strictly smaller sub-problems than
the original margin.

CHECKS:
  X1  delta_1(A) = 2 for A=2..6 (universal; and matches the A=2 closed form).
  X2  eps_1(A) = -2 - 2(A-2)/(A-1) for A=2..6.
  X3  identity -(eps_1+delta_1)/2 = (A-2)/(A-1) reproduces the measured margin
      coefficient (cross-check vs probe_7_34_endpoint_smallp_asymptotic.py).

Deps: numpy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~2-4 min.
"""
import itertools
import mpmath as mp
import numpy as np
mp.mp.dps = 60
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

def _pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def gA(A, p, D, s):
    A = int(A); p = mp.mpf(p); D = mp.mpf(D); s = mp.mpf(s)
    th = mp.log(D / ((A - 1) * (1 - D))); E = mp.e**(th + 1j * s)
    C = ((A - 1) * D * E + D - (A - 1)) / (A * D - (A - 1))
    eta = ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))
    C0 = ((A - 1) * D * mp.e**th + D - (A - 1)) / (A * D - (A - 1)); Cr = abs(C / C0)**2
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]; etb = mp.conj(eta)
    st = list(itertools.product(range(A), repeat=3)); npat = 5 if A >= 3 else 4; reps = [None] * npat
    for k, tr in enumerate(st):
        pp = _pat(tr)
        if reps[pp] is None: reps[pp] = k
    Q = mp.zeros(npat, npat)
    for a_ in range(npat):
        i = reps[a_]; x, xp, xq = st[i]
        for j, (y, yp, yq) in enumerate(st):
            Q[a_, _pat((y, yp, yq))] += T[xp][yp] * T[xq][yq] / T[x][y] * (eta if yp != y else 1) * (etb if yq != y else 1)
    ev, _ = mp.eig(Q); return Cr * max(abs(e) for e in ev)

def Dc(A, pv):
    def im(Dv):
        b0 = Dv / (A - 1); K = np.full((A, A), b0); np.fill_diagonal(K, 1 - Dv); Ki = np.linalg.inv(K)
        Tn = np.full((A, A), pv / (A - 1)); np.fill_diagonal(Tn, 1 - pv)
        By = lambda y: np.array([[Ki[y, a] * Tn[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-14, (A - 1) / A - 1e-7
    for _ in range(120):
        m = (lo + hi) / 2; lo, hi = (m, hi) if im(m) < 1e-15 else (lo, m)
    return mp.mpf(lo)

AS = range(2, 7)
PGRID = ('2e-3', '1e-3', '5e-4')

def delta1(A):
    vals = []
    for p in PGRID:
        pf = mp.mpf(p); dc = Dc(A, float(pf))
        vals.append((dc / (pf**2 / (4 * (A - 1))) - 1) / pf)
    return 2 * vals[-1] - vals[-2]  # Richardson in p

def eps1(A):
    vals = []
    for p in PGRID:
        pf = mp.mpf(p)
        Ds = [pf**2 / (4 * (A - 1)) / 50, pf**2 / (4 * (A - 1)) / 100]
        Es = []
        for D in Ds:
            R = (1 - gA(A, float(pf), D, mp.pi)) / (D * (1 - D) * 2)
            Es.append((2 - R) / D)
        E0 = 2 * Es[1] - Es[0]  # Richardson in D
        vals.append((E0 / (2 * (A - 1) / pf**2) - 1) / pf)
    return 2 * vals[-1] - vals[-2]  # Richardson in p

def main():
    print("-" * 78); print("X1  delta_1(A) = 2 universal (threshold p-correction, alphabet-free)")
    okd = True; d1 = {}
    for A in AS:
        d = delta1(A); d1[A] = d
        here = abs(d - 2) < mp.mpf('0.02'); okd = okd and here
        print(f"     A={A}: delta_1 = {mp.nstr(d,6)}  (=2: {here})")
    rep("X1 delta_1 = 2 for all A (universal)", okd)

    print("-" * 78); print("X2  eps_1(A) = -2 - 2(A-2)/(A-1)")
    oke = True; e1 = {}
    for A in AS:
        e = eps1(A); e1[A] = e
        tgt = -2 - 2 * mp.mpf(A - 2) / (A - 1)
        here = abs(e - tgt) < mp.mpf('0.02'); oke = oke and here
        print(f"     A={A}: eps_1 = {mp.nstr(e,6)}  target {mp.nstr(tgt,6)}  match:{here}")
    rep("X2 eps_1 = -2 - 2(A-2)/(A-1) for all A", oke)

    print("-" * 78); print("X3  identity: -(eps_1+delta_1)/2 = (A-2)/(A-1) (margin coefficient recovered)")
    oki = True
    for A in AS:
        lhs = -(e1[A] + d1[A]) / 2; tgt = mp.mpf(A - 2) / (A - 1)
        here = abs(lhs - tgt) < mp.mpf('0.02'); oki = oki and here
        print(f"     A={A}: -(eps_1+delta_1)/2 = {mp.nstr(lhs,6)}  (A-2)/(A-1) = {mp.nstr(tgt,6)}  match:{here}")
    rep("X3 decomposition reproduces the (A-2)/(A-1) margin coefficient", oki)

if __name__ == "__main__":
    print("=" * 78)
    print("BUG-009-D endpoint coefficient decomposition: (A-2)/(A-1) = -(eps_1+delta_1)/2")
    print("=" * 78)
    main()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("Both pieces are first-order perturbation coefficients (no root-solving):")
    print("delta_1=2 from the threshold discriminant; eps_1 from D=0 eigenvalue perturbation")
    print("of W_A(pi). Proving them in closed form upgrades the A>=3 leading-order endpoint")
    print("closure from numerics to proof.")

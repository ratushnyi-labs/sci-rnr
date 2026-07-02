#!/usr/bin/env python3
r"""
probe_7_34_endpoint_c2_numeric.py
============================================================================
BUG-009-D / Route B: the O(p^2) endpoint-margin coefficient c_2(A) -- numerically
pinned to clean rationals; A=2 VALIDATES the machinery against the proven 1/8.

CONTEXT. The leading-order endpoint theorem (all A, derived) is
    R_A(pi)|_{D_c} - 3/2 = [(A-2)/(A-1)] p + c_2(A) p^2 + O(p^3).
This probe pins c_2(A). Precision bottleneck is D_c: float64 bisection cannot
resolve O(p^2) corrections, so D_c is bisected on the 3x3 symmetry-reduced cubic
DISCRIMINANT (lemma_7_34_dc_smallp_delta1 structure) in mpmath at 60 dps; A=2 uses
the exact closed form D_c = 1/2 - sqrt(1-2p)/(2(1-p)).

RESULTS (Richardson over p = 4e-3, 2e-3, 1e-3):
    A=2:  c_2 -> 0.1250013 = 1/8   -- recovers the PROVEN A=2 value (validation:
          the entire expansion machinery closes the loop end-to-end);
    A=3..6:  c_2 -> -43/32, -95/72, -155/128, -223/200  (6-7 digits), i.e.
          denominators 8(A-1)^2, numerators quadratic with second difference 8:
    CONJECTURE (A>=3):   c_2(A) = -(4A^2 + 24A - 65) / (8 (A-1)^2).
(A=2 does NOT lie on the A>=3 family -- expected: the linear coefficient vanishes
there and the pattern space differs, 4 vs 5 classes.)

CONSEQUENCE. c_2(A) < 0 for A>=3 (range -1.1..-1.35): the second-order term
erodes the linear margin, with positivity of the two-term expansion up to
p* = [(A-2)/(A-1)]/|c_2| (e.g. A=3: p* ~ 0.37) -- comfortably covering the
small-p regime; consistent with the adversarial all-D verification (no violations).

STATUS: c_2(A) values are NUMERIC (rational fit conjectured); the exact derivation
needs (i) delta_2 (next order of the discriminant ansatz-series, mechanical) and
(ii) third-order rank-1 perturbation theory for the slope (the Q0 P_perp = 0
collapse extends; m3 has no new resolvent structure). A=2's 1/8 is already proven.

CHECKS:
  N1  A=2: c_2 -> 1/8 (machinery validation against the proven closed form).
  N2  A=3..6: c_2 matches -(4A^2+24A-65)/(8(A-1)^2) within Richardson tolerance.
  N3  positivity window: (A-2)/(A-1) p + c_2 p^2 > 0 for p <= 0.05, A=3..6
      (the sharpened two-term margin stays positive in the tight regime).

Deps: mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~3-6 min.
"""
import itertools
import mpmath as mp
mp.mp.dps = 60
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

# ---- high-precision D_c via the 3x3 symmetry-reduced cubic discriminant ----
def disc3(A, p, D):
    lam2 = 1 - D * A / (A - 1)
    al = mp.mpf(1) / A + (1 - mp.mpf(1) / A) / lam2
    be = mp.mpf(1) / A - (mp.mpf(1) / A) / lam2
    td = 1 - p; to = p / (A - 1)
    def Bred(y):
        Ki = {0: (al if y == 0 else be), 1: (be if y == 0 else al), 2: be}
        Trow = {0: {0: td, 1: to, 2: (A - 2) * to},
                1: {0: to, 1: td, 2: (A - 2) * to},
                2: {0: to, 1: to, 2: td + (A - 3) * to}}
        M = mp.zeros(3, 3)
        for i in range(3):
            for j in range(3):
                M[i, j] = Ki[i] * Trow[i][j]
        return M
    M = Bred(1) * Bred(0)
    c2 = M[0, 0] + M[1, 1] + M[2, 2]
    c1 = (M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]) + (M[0, 0] * M[2, 2] - M[0, 2] * M[2, 0]) + (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1])
    c0 = mp.det(M)
    a = -c2; b = c1; c = -c0
    return 18 * a * b * c - 4 * a**3 * c + a**2 * b**2 - 4 * b**3 - 27 * c**2

def Dc_hp(A, p):
    lo = p**2 / (8 * (A - 1)); hi = p**2 / (2 * (A - 1))
    flo = disc3(A, p, lo)
    for _ in range(200):
        m = (lo + hi) / 2
        if disc3(A, p, m) * flo > 0: lo = m
        else: hi = m
    return (lo + hi) / 2

def Dc2(p): return mp.mpf('0.5') - mp.sqrt(1 - 2 * p) / (2 * (1 - p))

# ---- replica margin ----
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

def margin(A, p):
    dc = Dc2(p) if A == 2 else Dc_hp(A, mp.mpf(p))
    return (1 - gA(A, float(p), dc, mp.pi)) / (dc * (1 - dc) * 2) - mp.mpf('1.5')

def c2_measure(A):
    c1 = mp.mpf(A - 2) / (A - 1)
    vals = []
    for p in ('4e-3', '2e-3', '1e-3'):
        pf = mp.mpf(p)
        vals.append((margin(A, pf) - c1 * pf) / pf**2)
    r1 = 2 * vals[1] - vals[0]; r2 = 2 * vals[2] - vals[1]
    return 2 * r2 - r1

def c2_conj(A): return -mp.mpf(4 * A**2 + 24 * A - 65) / (8 * (A - 1)**2)

def N1():
    print("-" * 78); print("N1  A=2: c_2 -> 1/8 (validation vs the proven closed form)")
    c = c2_measure(2)
    ok = abs(c - mp.mpf(1) / 8) < mp.mpf('5e-4')
    print(f"     c_2(2) = {mp.nstr(c,7)}  vs 1/8 = 0.125  match:{ok}")
    return rep("N1 A=2 recovers the proven 1/8 (machinery validated)", ok)

def N2():
    print("-" * 78); print("N2  A=3..6: c_2 = -(4A^2+24A-65)/(8(A-1)^2) (conjectured rational family)")
    ok = True
    for A in (3, 4, 5, 6):
        c = c2_measure(A); tgt = c2_conj(A)
        here = abs(c - tgt) < mp.mpf('2e-3'); ok = ok and here
        print(f"     A={A}: c_2 = {mp.nstr(c,7)}  conj {mp.nstr(tgt,7)} "
              f"(= -{4*A**2+24*A-65}/{8*(A-1)**2})  match:{here}")
    return rep("N2 c_2 matches -(4A^2+24A-65)/(8(A-1)^2), A=3..6 (CONJECTURE)", ok)

def N3():
    print("-" * 78); print("N3  positivity window: (A-2)/(A-1) p + c_2 p^2 > 0 for p <= 0.05")
    ok = True
    for A in (3, 4, 5, 6):
        c1 = mp.mpf(A - 2) / (A - 1); c2 = c2_conj(A)
        pstar = c1 / abs(c2)
        worst = min(c1 * mp.mpf(p) + c2 * mp.mpf(p)**2 for p in ('0.001', '0.01', '0.05'))
        here = worst > 0 and pstar > mp.mpf('0.05'); ok = ok and here
        print(f"     A={A}: p* = {mp.nstr(pstar,4)}  two-term margin at p<=0.05 min = {mp.nstr(worst,4)} (>0: {here})")
    return rep("N3 sharpened two-term margin positive through the tight regime", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("BUG-009-D endpoint O(p^2) coefficient c_2(A): numeric pin + rational conjecture")
    print("=" * 78)
    N1(); N2(); N3()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("Next: derive c_2(A) exactly -- delta_2 (discriminant order p^14) + third-order")
    print("rank-1 PT (m3; the Q0 P_perp=0 collapse extends). A=2's 1/8 already proven.")

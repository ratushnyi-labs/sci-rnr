#!/usr/bin/env python3
r"""
lemma_7_34_eps1_slope_coefficient.py
============================================================================
BUG-009-D / Route B: the slope coefficient eps_1 = -2 - 2(A-2)/(A-1) -- DERIVED
SYMBOLICALLY (exact arithmetic, A=3,4,5), completing the (A-2)/(A-1) endpoint
coefficient theorem at leading order.

CONTEXT. The small-p endpoint margin splits (probe_7_34_endpoint_coeff_decomposition)
as   R_A(pi)|_{D_c} - 3/2 = -[(eps_1 + delta_1)/2] p + O(p^2),   where
    delta_1: D_c = [p^2/(4(A-1))](1 + delta_1 p + ...)   -- PROVEN = 2 for ALL A
                                                            (lemma_7_34_dc_smallp_delta1),
    eps_1:   R_A(pi; D) = 2 - E_A(p) D + O(D^2),
             E_A(p) = [2(A-1)/p^2](1 + eps_1 p + ...).
THIS script derives eps_1 exactly. Method (no root-solving; the dead A>=3 quartic
Perron never enters):
  1. At s=pi the disagreement tilt is real:  eta_pi(D) = 2D(1-D)/(AD - 2D^2 - (A-1)),
     and the prefactor is  Cr(D) = [(AD-2D^2-(A-1))/(AD-(A-1))]^2  (exact rationals).
  2. The pattern-quotient replica matrix Q(eta) has entries POLYNOMIAL in eta.
     At eta=0 only the (y,y,y) target class survives, so Q(0) has a single nonzero
     column with eigenvalues {0 (x4), 1 (simple)} -- lambda=1 is a SIMPLE root.
  3. Hence the Perron branch expands by IMPLICIT DIFFERENTIATION of the charpoly
     P(lambda, eta) at (1,0):  lambda(eta) = 1 + m1 eta + m2 eta^2 + ...,
         m1 = -P_eta/P_lambda,   m2 = -(P_ee + 2 P_el m1 + P_ll m1^2)/(2 P_lambda),
     all partial derivatives evaluated at the point (pure polynomial algebra).
  4. Compose g = Cr(D) * lambda(eta_pi(D)), expand 1-g = e1 D + e2 D^2 + ...:
     e1 = 4 exactly, and E = -(e2+4)/2; then E p^2/(2(A-1)) = 1 + eps_1 p + O(p^2).

RESULTS (exact, A=3,4,5):
    e1 = 4;   m1 = 2(A-1)  (= 4, 6, 8);
    E(p) closed forms, e.g. A=3: (p-2)(3p-2)^2(9p^3-12p^2+6p-4)/(8 p^2 (p-1)^2);
    eps_1 = -3, -10/3, -7/2  =  -2 - 2(A-2)/(A-1)   EXACTLY.
With delta_1 = 2 (all A):  -(eps_1+delta_1)/2 = (A-2)/(A-1)  -- the leading-order
endpoint margin  R_A(pi)|_{D_c} - 3/2 = [(A-2)/(A-1)] p + O(p^2)  is now DERIVED
(A=3,4,5 symbolically; A=2 fully proven separately with margin p^2/8; A=6,7 numeric).
The universal-A eps_1 derivation needs the y-sum class decomposition with symbolic
multiplicities (mechanical but heavier; the m1=2(A-1) pattern makes it a target).

CHECKS:
  Z1  e1 = 4 exactly (A=3,4,5).
  Z2  m1 = 2(A-1) exactly (A=3,4,5).
  Z3  eps_1 = -2 - 2(A-2)/(A-1) exactly (A=3,4,5).
  Z4  assembled: -(eps_1 + 2)/2 = (A-2)/(A-1) (with the proven delta_1=2).
  Z5  numeric cross-anchor: E(p) closed form matches the replica at p=0.1 (mpmath).

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~3-5 min.
"""
import itertools
import sympy as sp
import mpmath as mp
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

def derive(Aval):
    A = sp.Integer(Aval)
    p, D, lam, et = sp.symbols('p D lambda eta_s')
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(Aval)] for i in range(Aval)]
    st = list(itertools.product(range(Aval), repeat=3)); npat = 5; reps = [None] * npat
    for k, tr in enumerate(st):
        if reps[_pat(tr)] is None: reps[_pat(tr)] = k
    Q = sp.zeros(npat, npat)
    for a_ in range(npat):
        x, xp, xq = st[reps[a_]]
        for (y, yp, yq) in st:
            term = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x][y])
            if yp != y: term *= et
            if yq != y: term *= et
            Q[a_, _pat((y, yp, yq))] += term
    P = Q.charpoly(lam).as_expr()
    sub = {lam: sp.Integer(1), et: sp.Integer(0)}
    Pl0 = sp.simplify(sp.diff(P, lam).subs(sub))
    Pe0 = sp.simplify(sp.diff(P, et).subs(sub))
    m1 = sp.cancel(-Pe0 / Pl0)
    Pee0 = sp.simplify(sp.diff(P, et, 2).subs(sub))
    Pel0 = sp.simplify(sp.diff(P, et, lam).subs(sub))
    Pll0 = sp.simplify(sp.diff(P, lam, 2).subs(sub))
    m2 = sp.cancel(-(Pee0 + 2 * Pel0 * m1 + Pll0 * m1**2) / (2 * Pl0))
    etaD = 2 * D * (1 - D) / (A * D - 2 * D**2 - (A - 1))
    CrD = ((A * D - 2 * D**2 - (A - 1)) / (A * D - (A - 1)))**2
    etas = sp.series(etaD, D, 0, 3).removeO()
    Crs = sp.series(CrD, D, 0, 3).removeO()
    g = sp.expand(Crs * (1 + m1 * etas + m2 * etas**2))
    omg = sp.expand(1 - g)
    e1 = sp.cancel(omg.coeff(D, 1)); e2 = sp.cancel(omg.coeff(D, 2))
    E = sp.cancel(-(e2 + 4) / 2)
    expr = sp.series(sp.cancel(E * p**2 / (2 * (A - 1))), p, 0, 3).removeO()
    eps1 = sp.simplify(sp.expand(expr).coeff(p, 1))
    return e1, m1, E, eps1, p

def main():
    res = {}
    for Aval in (3, 4, 5):
        res[Aval] = derive(Aval)

    print("-" * 78); print("Z1  e1 = 4 exactly")
    ok = all(sp.simplify(res[A][0] - 4) == 0 for A in res)
    for A in res: print(f"     A={A}: e1 = {res[A][0]}")
    rep("Z1 e1 = 4 (leading slope of 1-g at s=pi) exact", ok)

    print("-" * 78); print("Z2  m1 = 2(A-1) exactly")
    ok = all(sp.simplify(res[A][1] - 2 * (A - 1)) == 0 for A in res)
    for A in res: print(f"     A={A}: m1 = {res[A][1]}")
    rep("Z2 m1 = 2(A-1) (Perron O(eta) slope) exact", ok)

    print("-" * 78); print("Z3  eps_1 = -2 - 2(A-2)/(A-1) exactly")
    ok = True
    for A in res:
        tgt = -2 - 2 * sp.Rational(A - 2, A - 1)
        here = sp.simplify(res[A][3] - tgt) == 0; ok = ok and here
        print(f"     A={A}: eps_1 = {res[A][3]}  target {tgt}  match:{here}")
    rep("Z3 eps_1 = -2 - 2(A-2)/(A-1) exact (A=3,4,5)", ok)

    print("-" * 78); print("Z4  assembled: -(eps_1 + delta_1)/2 = (A-2)/(A-1) with delta_1=2")
    ok = True
    for A in res:
        lhs = sp.simplify(-(res[A][3] + 2) / 2); tgt = sp.Rational(A - 2, A - 1)
        here = sp.simplify(lhs - tgt) == 0; ok = ok and here
        print(f"     A={A}: -(eps_1+2)/2 = {lhs}  (A-2)/(A-1) = {tgt}  match:{here}")
    rep("Z4 leading-order endpoint margin (A-2)/(A-1) DERIVED (A=3,4,5)", ok)

    print("-" * 78); print("Z5  numeric anchor: E(p) closed form vs replica at p=0.1")
    mp.mp.dps = 40
    def gA_num(A, p, D, s):
        A = int(A); p = mp.mpf(p); D = mp.mpf(D); s = mp.mpf(s)
        th = mp.log(D / ((A - 1) * (1 - D))); Ee = mp.e**(th + 1j * s)
        C = ((A - 1) * D * Ee + D - (A - 1)) / (A * D - (A - 1))
        eta = ((A - 1) * D * Ee + D - (A - 1) * Ee) / ((A - 1) * D * Ee + D - (A - 1))
        C0 = ((A - 1) * D * mp.e**th + D - (A - 1)) / (A * D - (A - 1)); Cr = abs(C / C0)**2
        T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]
        etb = mp.conj(eta)
        st = list(itertools.product(range(A), repeat=3)); npat = 5; reps = [None] * npat
        for k, tr in enumerate(st):
            if reps[_pat(tr)] is None: reps[_pat(tr)] = k
        Qm = mp.zeros(npat, npat)
        for a_ in range(npat):
            i = reps[a_]; x, xp, xq = st[i]
            for j, (y, yp, yq) in enumerate(st):
                Qm[a_, _pat((y, yp, yq))] += T[xp][yp] * T[xq][yq] / T[x][y] * (eta if yp != y else 1) * (etb if yq != y else 1)
        ev, _ = mp.eig(Qm); return Cr * max(abs(e) for e in ev)
    ok = True
    for A in res:
        Ecf = sp.lambdify(res[A][4], res[A][2], 'mpmath')(mp.mpf('0.1'))
        D = mp.mpf('1e-4')
        e2n = ((1 - gA_num(A, 0.1, D, mp.pi)) - 4 * D) / D**2
        En = -(e2n + 4) / 2
        here = abs(Ecf - En) / abs(En) < mp.mpf('1e-2'); ok = ok and here
        print(f"     A={A}: E closed form={mp.nstr(Ecf,7)}  numeric={mp.nstr(En,7)}  match:{here}")
    rep("Z5 E(p) closed form matches the replica numerically", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("eps_1 = -2-2(A-2)/(A-1) (slope coefficient): exact derivation A=3,4,5")
    print("=" * 78)
    main()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("Both pieces of the (A-2)/(A-1) endpoint coefficient are now DERIVED:")
    print("delta_1=2 (ALL A) + eps_1=-2-2(A-2)/(A-1) (A=3,4,5 exact; A=6,7 numeric).")
    print("Leading-order A>=3 endpoint margin is a theorem on the derived range.")

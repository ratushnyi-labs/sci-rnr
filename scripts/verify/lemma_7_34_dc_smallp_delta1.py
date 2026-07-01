#!/usr/bin/env python3
r"""
lemma_7_34_dc_smallp_delta1.py
============================================================================
BUG-009-D / Route B: the threshold coefficient delta_1 = 2 -- DERIVED SYMBOLICALLY
(exact rational arithmetic) for A = 2, 3, 4, 5; numerically confirmed A = 6, 7.

CONTEXT. probe_7_34_endpoint_coeff_decomposition.py split the small-p endpoint margin
coefficient as (A-2)/(A-1) = -(eps_1+delta_1)/2, with delta_1 the relative p-correction
of the Gray threshold,
    D_c(A,p) = [p^2/(4(A-1))] * (1 + delta_1 p + O(p^2)),
measured delta_1 = 2 for all A. THIS script derives delta_1 = 2 exactly:

  * A=2: read off the closed form D_c = 1/2 - sqrt(1-2p)/(2(1-p)) = (p^2/4)(1+2p+O(p^2)).
  * A>=3 (A=3,4,5 here): the alternating-word product B_1 B_0 reduces, by the symmetry
    of the A-ary symmetric channel/source, to a 3x3 block on the classes {0},{1},{rest}
    (K^{-1} has two distinct entries al/be; T likewise). The Gray threshold is where the
    class-0/1 eigenvalue pair of this cubic turns complex, i.e. the CUBIC DISCRIMINANT
    vanishes. Substituting the ansatz D = W p^2 (1 + d1 p) and series-expanding the
    discriminant in p (exact sympy arithmetic):
       - the leading order is p^12 (eigenvalue scale ~ p^2, disc ~ scale^6), and its
         coefficient vanishes iff W = 1/(4(A-1))  [or the degenerate W=0], recovering
         the known leading threshold;
       - the p^13 coefficient, evaluated at W* = 1/(4(A-1)), vanishes iff d1 = 2.
    Both steps are exact polynomial algebra -- no numerics, no root-solving of the
    (dead) A>=3 quartic Perron.

CONCLUSION: delta_1 = 2 is PROVEN for A = 2..5 and A-universality is confirmed
numerically to A=7 (the fully-symbolic-A derivation is the same computation with A a
parameter; it is heavy but contains no new mathematical content). One of the two
pieces of the (A-2)/(A-1) coefficient is now closed-form; the remaining piece is
eps_1 = -2 - 2(A-2)/(A-1) (D=0 eigenvalue perturbation of W_A(pi)).

CHECKS:
  Y1  A=2: series of the closed-form D_c gives (p^2/4)(1 + 2p + O(p^2)) exactly.
  Y2  A=3,4,5: cubic-discriminant leading order p^12 vanishes iff W=1/(4(A-1)) (exact).
  Y3  A=3,4,5: p^13 coefficient at W* vanishes iff d1=2 (exact).
  Y4  A=6,7: numeric confirmation delta_1 -> 2 (bisected D_c, Richardson in p).

Deps: sympy, numpy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~2-3 min.
"""
import sympy as sp
import numpy as np
import mpmath as mp
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

def Y1():
    print("-" * 78); print("Y1  A=2 closed form: D_c = (p^2/4)(1+2p+O(p^2))")
    p = sp.symbols('p', positive=True)
    Dc = sp.Rational(1, 2) - sp.sqrt(1 - 2 * p) / (2 * (1 - p))
    ser = sp.series(Dc, p, 0, 5).removeO()
    lead = sp.simplify(ser - (p**2 / 4 + p**3 / 2))
    ok = sp.simplify(sp.series(lead, p, 0, 4).removeO()) == 0
    print(f"     series(D_c) = {sp.expand(ser)}  -> p^2/4 + p^3/2 + O(p^4): {ok}")
    return rep("Y1 A=2: delta_1 = 2 from the closed form (exact)", ok)

def _derive(Aval):
    A = sp.Integer(Aval)
    p = sp.symbols('p', positive=True)
    W, d1 = sp.symbols('W d1')
    D = W * p**2 * (1 + d1 * p)
    lam2 = 1 - D * A / (A - 1)
    al = sp.Rational(1, Aval) + (1 - sp.Rational(1, Aval)) / lam2
    be = sp.Rational(1, Aval) - sp.Rational(1, Aval) / lam2
    td = 1 - p; to = p / (A - 1)
    def Bred(y):
        Ki = {0: (al if y == 0 else be), 1: (be if y == 0 else al), 2: be}
        Trow = {0: {0: td, 1: to, 2: (A - 2) * to},
                1: {0: to, 1: td, 2: (A - 2) * to},
                2: {0: to, 1: to, 2: td + (A - 3) * to}}
        M = sp.zeros(3, 3)
        for i in range(3):
            for j in range(3):
                M[i, j] = Ki[i] * Trow[i][j]
        return M
    M = Bred(1) * Bred(0)
    M = M.applyfunc(lambda e: sp.cancel(sp.together(e)))
    c2 = sp.cancel(M.trace())
    c1 = sp.cancel(sum(M[[i for i in r], [j for j in r]].det() for r in ([0, 1], [0, 2], [1, 2])))
    c0 = sp.cancel(M.det())
    a = -c2; b = c1; c = -c0
    disc = 18 * a * b * c - 4 * a**3 * c + a**2 * b**2 - 4 * b**3 - 27 * c**2
    ser = sp.expand(sp.series(disc, p, 0, 17).removeO())
    poly = sp.Poly(ser, p)
    orders = sorted(m[0] for m in poly.monoms())
    lead = sp.factor(poly.coeff_monomial(p**orders[0]))
    solW = [sp.simplify(s) for s in sp.solve(sp.Eq(lead, 0), W)]
    Wstar = sp.Rational(1, 4 * (Aval - 1))
    nxt = poly.coeff_monomial(p**orders[1])
    nxtW = sp.factor(sp.simplify(nxt.subs(W, Wstar)))
    sold1 = [sp.simplify(s) for s in sp.solve(sp.Eq(nxtW, 0), d1)]
    return orders[0], solW, orders[1], sold1, Wstar

def Y2Y3():
    print("-" * 78); print("Y2/Y3  A=3,4,5: exact discriminant derivation of W*=1/(4(A-1)) and d1=2")
    ok2 = True; ok3 = True
    for Aval in (3, 4, 5):
        o0, solW, o1, sold1, Wstar = _derive(Aval)
        w_ok = (o0 == 12) and (Wstar in solW)
        d_ok = (sold1 == [2])
        ok2 = ok2 and w_ok; ok3 = ok3 and d_ok
        print(f"     A={Aval}: lead p^{o0} W-roots={solW} (W*={Wstar}: {w_ok});  p^{o1} => d1={sold1} ({d_ok})")
    rep("Y2 leading order vanishes iff W=1/(4(A-1)) (A=3,4,5, exact)", ok2)
    rep("Y3 next order at W* forces d1=2 (A=3,4,5, exact)", ok3)
    return ok2 and ok3

def Y4():
    print("-" * 78); print("Y4  A=6,7 numeric confirmation delta_1 -> 2")
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
    ok = True
    for A in (6, 7):
        vals = []
        for p in ('2e-3', '1e-3', '5e-4'):
            pf = mp.mpf(p); dc = Dc(A, float(pf))
            vals.append((dc / (pf**2 / (4 * (A - 1))) - 1) / pf)
        d = 2 * vals[-1] - vals[-2]
        here = abs(d - 2) < mp.mpf('0.02'); ok = ok and here
        print(f"     A={A}: delta_1 (Richardson) = {mp.nstr(d,5)}  (=2: {here})")
    return rep("Y4 A=6,7 numeric delta_1=2 (universality supported)", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("delta_1 = 2 (Gray-threshold small-p correction): exact derivation A=2..5")
    print("=" * 78)
    Y1(); Y2Y3(); Y4()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("First of the two (A-2)/(A-1) pieces is closed-form. Remaining: eps_1 =")
    print("-2-2(A-2)/(A-1) via D=0 eigenvalue perturbation of W_A(pi).")

#!/usr/bin/env python3
r"""
probe_7_34_route_b_pert.py
============================================================================
ROUTE B step 11 -- the SECOND-ORDER PERRON EXPANSION (the heart of the proof).
DISCOVERY (decomp): 1-g_A is carried entirely by the Perron defect 1-rho(Q),
with eta_s = (D/(A-1))(w-1)(1+O(D)), w=e^is, etb=conj(eta).

By boundary affineness rho(eta,0)=1+(A-1)eta, rho(0,etb)=1+(A-1)etb EXACTLY, so
the Taylor expansion of rho(eta,etb) about (0,0) has NO pure powers eta^k,etb^k
(k>=2). Hence to 2nd order:
        rho(eta,etb) = 1 + (A-1)(eta+etb) + c2 * eta*etb + O(3).
We compute c2 = d^2 rho/(d eta d etb) at 0 SYMBOLICALLY (sympy) in A,p. Then on
the diagonal etb=conj(eta) with eta=(D/(A-1))(w-1):
  eta+etb = 2 Re(eta) = (D/(A-1))(2 cos s - 2) = -(2D/(A-1)) t     (t=1-cos s)
  eta*etb = |eta|^2 = (D/(A-1))^2 |w-1|^2 = (D/(A-1))^2 * 2t
So
  rho = 1 - 2 D t + c2 (D/(A-1))^2 2t + O(D^3) = 1 - 2 D t [ 1 - c2 D/(A-1)^2 ] + ...
=> 1 - rho ~ 2 D t (leading), matching the measured leading-2 EXACTLY (note D
vs D(1-D): 1-D correction is part of O(D^2)). The SIGN/size of c2 controls the
2nd-order drop from 2 toward 1.5.

CHECKS:
  (1) symbolic c2(A,p) = mixed 2nd partial of rho at 0; print closed form.
  (2) verify leading 1-rho = 2 D t (the (A-1)(eta+etb) term) -- the JENSEN-floor
      analog comes from rho's LINEAR term + the exact eta leading form. CLEAN.
  (3) the full 1-g_A >= (3/2) D(1-D) t reduces to bounding c2*(D/(A-1)^2) and the
      prefactor + higher orders. Compute the exact 2nd-order coefficient of
      (1-rho)/(D t) and confirm it stays in [1.5, 2].
"""
import itertools
import sympy as sp

def pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def build_Q_sym(A, p, eta, etb):
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]
    st = list(itertools.product(range(A), repeat=3)); npat = 5 if A >= 3 else 4
    reps = [None] * npat
    for k, tr in enumerate(st):
        if reps[pat(tr)] is None: reps[pat(tr)] = k
    Q = sp.zeros(npat, npat)
    for a in range(npat):
        x, xp, xq = st[reps[a]]
        for (y, yp, yq) in st:
            t = T[xp][yp] * T[xq][yq] / T[x][y]
            if yp != y: t *= eta
            if yq != y: t *= etb
            Q[a, pat((y, yp, yq))] += t
    return Q

def perron_series(A, pval):
    """rho(eta,etb) to 2nd order about (0,0) via implicit char poly. Return
    (linear coeff of eta, mixed coeff c2)."""
    eta, etb, lam = sp.symbols('eta etb lambda')
    p = sp.Rational(pval) if not isinstance(pval, sp.Expr) else pval
    Q = build_Q_sym(A, p, eta, etb)
    n = Q.rows
    chi = (lam * sp.eye(n) - Q).det()
    chi = sp.expand(chi)
    # Perron branch lam = 1 + a1*(eta+etb) + c2 eta etb + (pure terms which vanish).
    # We KNOW pure terms vanish (affineness); fit a1, c2 by implicit differentiation.
    # Substitute the ansatz and match orders.
    a1, c2, a20, a02 = sp.symbols('a1 c2 a20 a02')
    ans = 1 + a1 * (eta + etb) + a20 * eta**2 + a02 * etb**2 + c2 * eta * etb
    sub = chi.subs(lam, ans)
    sub = sp.expand(sub)
    # collect to 2nd order in (eta,etb)
    eqs = []
    # order (1,0): coeff of eta
    eqs.append(sub.coeff(eta, 1).coeff(etb, 0).subs({eta: 0, etb: 0}))
    eqs.append(sub.coeff(etb, 1).coeff(eta, 0).subs({eta: 0, etb: 0}))
    eqs.append(sub.coeff(eta, 2).coeff(etb, 0).subs({eta: 0, etb: 0}))
    eqs.append(sub.coeff(etb, 2).coeff(eta, 0).subs({eta: 0, etb: 0}))
    eqs.append(sub.coeff(eta, 1).coeff(etb, 1).subs({eta: 0, etb: 0}))
    solu = sp.solve(eqs, [a1, c2, a20, a02], dict=True)
    return solu

def main():
    print("Symbolic 2nd-order Perron expansion rho=1+a1(eta+etb)+c2 eta etb+pure(=0):")
    for A in (2, 3, 4):
        for pv in ('1/5', '1/4'):
            sol = perron_series(A, pv)
            if sol:
                s = sol[0]
                a1 = sp.simplify(s[sp.Symbol('a1')])
                c2 = sp.simplify(s[sp.Symbol('c2')])
                a20 = sp.simplify(s[sp.Symbol('a20')])
                a02 = sp.simplify(s[sp.Symbol('a02')])
                print(f"  A={A} p={pv}: a1={a1} (should be A-1={A-1}); c2={c2}; "
                      f"pure a20={a20}, a02={a02} (should be 0)")
            else:
                print(f"  A={A} p={pv}: no solution")
    print()
    print("Now symbolic in p (A fixed): a1 and c2 as functions of p.")
    for A in (2, 3, 4):
        p = sp.Symbol('p', positive=True)
        sol = perron_series(A, p)
        if sol:
            s = sol[0]
            a1 = sp.simplify(s[sp.Symbol('a1')])
            c2 = sp.simplify(s[sp.Symbol('c2')])
            print(f"  A={A}: a1={a1}")
            print(f"        c2={c2}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
probe_7_34k_aary_gray.py
============================================================================
Lemma 7.34k (A-ary symmetric Gray threshold -- the general alphabet law).
Generalises binary (Lemma 7.34f) and ternary (Lemma 7.34j) to every A>=2.

Chain: A-ary symmetric Markov, T(i,i)=1-p, T(i,j)=p/(A-1) (j!=i), pi uniform.
Hamming distortion; backward test channel = A-ary symmetric channel
K(y|x)=1-D (y=x), D/(A-1) (y!=x), so K=(1-D-b0)I+b0 J, b0=D/(A-1),
K^{-1}=(1/(1-D-b0))(I - b0 J)  [since the J-eigenvalue of K is 1].
Deconvolution P_{Y*}=(K^{-1})^{(x)n}P_X; D_c^(A)(p)=sup{D: P_{Y*}>=0 for all n}.

Checks:
  K1  K^{-1} exact (K K^{-1}=I) for A=2..6.                       [sympy]
  K2  SPECTATOR DECOUPLING: the period-2 product M=B_1 B_0 (B_y[a,b]=
      K^{-1}[y,a]T(a,b)) has char poly = (spectator linear)^{A-2} x (cubic)
      for A>=3 (and a quadratic for A=2): the A-2 non-participating symbols
      form a symmetric block, leaving a 3-dim relevant block.        [sympy, A=2..5]
  K3  threshold degree: D_c^(A) is the complex-onset of the relevant
      cubic's dominant eigenpair = smallest root of its discriminant, an
      A-dependent QUARTIC in D for all A>=3 (clean sqrt for A=2).     [sympy, A=2..4]
  K4  small-p law: D_c^(A)(p) = p^2/(4(A-1)) + O(p^3) (exact leading coeff
      for A=2..4 by the D=p^2 u substitution; numeric all A).          [sympy + numeric]
  K5  binding word = the 2-symbol ALTERNATING word: min P_{Y*} over ALL
      A^n words just below D_c^(A) is >=0 (no other word binds first).  [exact-ish, A=3,4,5, n<=7]
  K6  D_c^(A)(p) strictly DECREASING in A (Gray region shrinks with the
      alphabet).                                                       [numeric]
Deps: sympy, numpy.
"""
import itertools
import numpy as np
import sympy as sp

PASS = True
def rep(name, ok):
    global PASS
    PASS = PASS and ok
    print(f"  {name:<60} {'PASS' if ok else 'FAIL'}")
    return ok


def _sym_setup(A):
    p, D, lam = sp.symbols('p D lam')
    b0 = sp.Rational(1, A-1)*D
    a0 = 1 - D - b0
    I = sp.eye(A); J = sp.ones(A, A)
    Kinv = (1/a0)*(I - b0*J)          # a0 + A*b0 = 1
    T = sp.Matrix(A, A, lambda i, j: (1-p) if i == j else p/(A-1))
    B = lambda y: sp.Matrix(A, A, lambda a, b: Kinv[y, a]*T[a, b])
    return p, D, lam, Kinv, T, B


def K1():
    print("-"*78); print("K1  A-ary K^{-1} exact, A=2..6")
    ok = True
    for A in range(2, 7):
        D = sp.symbols('D'); b0 = sp.Rational(1, A-1)*D; a0 = 1-D-b0
        I = sp.eye(A); J = sp.ones(A, A)
        K = a0*I + b0*J; Kinv = (1/a0)*(I - b0*J)
        ok = ok and sp.simplify(K*Kinv - I) == sp.zeros(A, A)
    return rep("K1 K K^{-1}=I (A=2..6)", ok)


def K2():
    print("-"*78); print("K2  spectator decoupling: char poly M = (cubic) x (linear)^{A-3}, A>=3")
    ok = True
    for A in range(2, 6):
        p, D, lam, Kinv, T, B = _sym_setup(A)
        M = B(1)*B(0)
        facs = sp.factor_list(sp.factor(M.charpoly(lam).as_expr()), lam)[1]
        degs = sorted(sp.Poly(f, lam).degree() for f, m in facs for _ in range(m))
        if A == 2:
            good = degs == [2]                 # binary: 2-dim relevant block
        else:
            # one cubic relevant block + (A-3) linear spectators (the A-2 spectator
            # symbols split into 1 coupled average mode + A-3 decoupled modes)
            good = degs.count(3) == 1 and degs.count(1) == (A-3) and sum(degs) == A
        print(f"     A={A}: eigen-factor degrees {degs}  ({'2-dim relevant (binary)' if A==2 else '3-dim cubic relevant + %d spectators'%(A-3)})")
        ok = ok and good
    return rep("K2 spectator decoupling (3-dim relevant block for A>=3)", ok)


def K3():
    print("-"*78); print("K3  D_c^(A) = smallest root of the relevant-block discriminant: QUARTIC for A>=3")
    ok = True
    for A in range(2, 5):
        p, D, lam, Kinv, T, B = _sym_setup(A)
        M = B(1)*B(0)
        facs = sp.factor_list(sp.factor(M.charpoly(lam).as_expr()), lam)[1]
        # relevant block = the highest-degree factor (cubic for A>=3, quadratic for A=2)
        rel = max((f for f, m in facs), key=lambda f: sp.Poly(f, lam).degree())
        disc = sp.factor(sp.together(sp.discriminant(sp.Poly(rel, lam), lam)))
        num, den = sp.fraction(disc)
        binding = max((g for g, _ in sp.factor_list(num, D)[1]),
                      key=lambda g: sp.Poly(g, D).degree())
        deg = sp.Poly(binding, D).degree()
        expect = 2 if A == 2 else 4
        print(f"     A={A}: binding factor degree in D = {deg} (expect {expect})")
        ok = ok and deg == expect
    return rep("K3 binding factor: quadratic (A=2) / quartic (A=3,4)", ok)


def K4():
    print("-"*78); print("K4  small-p law D_c^(A) = p^2/(4(A-1)) + O(p^3)")
    ok = True
    for A in range(2, 5):
        p, D, lam, Kinv, T, B = _sym_setup(A)
        M = B(1)*B(0)
        facs = sp.factor_list(sp.factor(M.charpoly(lam).as_expr()), lam)[1]
        rel = max((f for f, m in facs), key=lambda f: sp.Poly(f, lam).degree())
        num = sp.fraction(sp.together(sp.discriminant(sp.Poly(rel, lam), lam)))[0]
        u = sp.symbols('u', positive=True)
        sub = sp.expand(num.subs(D, p**2*u))
        po = sp.Poly(sub, p); low = min(m[0] for m in po.monoms())
        roots = sp.solve(po.coeff_monomial(p**low), u)
        ok = ok and sp.Rational(1, 4*(A-1)) in roots
        print(f"     A={A}: small-p u-root(s)={roots}, 1/(4(A-1))={sp.Rational(1,4*(A-1))}")
    # numeric all A
    def Dc_num(A, pv):
        def im(Dv):
            b0 = Dv/(A-1); K = np.full((A, A), b0); np.fill_diagonal(K, 1-Dv); Ki = np.linalg.inv(K)
            Tn = np.full((A, A), pv/(A-1)); np.fill_diagonal(Tn, 1-pv)
            By = lambda y: np.array([[Ki[y, a]*Tn[a, b] for b in range(A)] for a in range(A)])
            ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
            return max(abs(t[0].imag), abs(t[1].imag))
        lo, hi = 1e-12, (A-1)/A-1e-7
        for _ in range(70):
            mid = (lo+hi)/2; lo, hi = (mid, hi) if im(mid) < 1e-11 else (lo, mid)
        return lo
    num_ok = all(abs(Dc_num(A, 1e-3)/1e-6/(1/(4*(A-1))) - 1) < 5e-3 for A in (2, 3, 4, 5, 6, 8))
    rep("K4a exact small-p coeff 1/(4(A-1)) (A=2,3,4)", ok)
    return rep("K4b numeric D_c/p^2 -> 1/(4(A-1)) (A=2..8)", num_ok) and ok


def _Dc_num(A, pv):
    def im(Dv):
        b0 = Dv/(A-1); K = np.full((A, A), b0); np.fill_diagonal(K, 1-Dv); Ki = np.linalg.inv(K)
        Tn = np.full((A, A), pv/(A-1)); np.fill_diagonal(Tn, 1-pv)
        By = lambda y: np.array([[Ki[y, a]*Tn[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-12, (A-1)/A-1e-7
    for _ in range(60):
        mid = (lo+hi)/2; lo, hi = (mid, hi) if im(mid) < 1e-11 else (lo, mid)
    return lo


def K5():
    print("-"*78); print("K5  binding word = alternating: min P_Y* over ALL A^n words just below D_c is >=0")
    ok = True
    for A in (3, 4, 5):
        dc = _Dc_num(A, 0.2); worst = 0.0
        for fac in (0.99, 0.999):
            b0 = dc*fac/(A-1); K = np.full((A, A), b0); np.fill_diagonal(K, 1-dc*fac); Ki = np.linalg.inv(K)
            Tn = np.full((A, A), 0.2/(A-1)); np.fill_diagonal(Tn, 0.8)
            By = lambda y: np.array([[Ki[y, a]*Tn[a, b] for b in range(A)] for a in range(A)])
            for w in itertools.product(range(A), repeat=7):
                v = np.ones(A)/A
                for y in w:
                    v = v @ By(y)
                worst = min(worst, v.sum())
        print(f"     A={A}: D_c={dc:.6f}; min P_Y* over {A}^7 words at 0.99/0.999 D_c = {worst:+.2e}")
        ok = ok and worst >= -1e-12
    return rep("K5 alternating-onset = the all-n threshold (no other word binds first)", ok)


def K6():
    print("-"*78); print("K6  D_c^(A)(p) strictly decreasing in A")
    vals = [_Dc_num(A, 0.2) for A in range(2, 9)]
    ok = all(vals[i] > vals[i+1] for i in range(len(vals)-1))
    print(f"     D_c^(A)(0.2), A=2..8: {[round(v,6) for v in vals]}")
    return rep("K6 monotone decreasing in A (Gray region shrinks with alphabet)", ok)


if __name__ == "__main__":
    print("="*78)
    print("Lemma 7.34k (A-ary symmetric Gray threshold) -- verifier")
    print("="*78)
    K1(); K2(); K3(); K4(); K5(); K6()
    print("="*78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")

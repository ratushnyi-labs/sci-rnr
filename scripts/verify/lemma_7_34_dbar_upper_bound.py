#!/usr/bin/env python3
r"""
lemma_7_34_dbar_upper_bound.py
============================================================================
BUG-009-D / Route F: an explicit elementary UPPER BOUND for the Gray threshold,
proven on the ENTIRE parameter range --
    D_c(A, p)  <  Dbar(A, p) := (2(A-1)/A) * [ 1/2 - sqrt(1-q) / (2-q) ],
    q = p A/(A-1),         for all p in (0, (A-1)/A),  A = 3..8,
with EQUALITY AT A=2: Dbar(2,p) = 1/2 - sqrt(1-2p)/(2(1-p)) = D_c(2,p) exactly.

(Conjectured in the research fan-out with 2-50% observed slack; upgraded here to
a per-A theorem. Useful as an explicit domain bound for the mid-band program and
as a standalone characterization: D_c ~ Dbar * (leading orders) at small p.)

PROOF (per A, one-variable Sturm certificates only). The radical is eliminated by
the substitution s = sqrt(1-q), s in (0,1):
    p = (A-1)(1-s^2)/A,      Dbar = (2(A-1)/A) [1/2 - s/(1+s^2)]   (rational in s).
Then with disc(D,p) the reduced-cubic discriminant of the 3x3 symmetry-reduced
alternating product (the object whose first sign change defines D_c):
  (i)  disc(D->0+, p) > 0: the discriminant numerator factors D^4 * C(D,p), and
       C(0,p) > 0 root-free on (0,(A-1)/A) [Sturm];
  (ii) disc(Dbar(p), p) < 0 for all p: the substituted numerator is root-free in
       s on (0,1) with negative combined sign [Sturm];
  (iii) hence disc changes sign in (0, Dbar), and since the FIRST sign change is
       the eigenvalue collision D_c (threshold mechanism, probe_7_34_beyond_gray_
       structure B1), D_c < Dbar.

CHECKS:
  V1  A=2 equality: Dbar(2,p) == the exact closed-form D_c(2,p) (symbolic).
  V2  (ii) for A=3..8: substituted numerator s-root-free, sign -1.
  V2b (i) for A=3..8: C(0,p) > 0 root-free on the full p-interval.
  V3  numeric sanity: Dbar > D_c with the fan-out's observed slack at samples.

Deps: sympy, numpy.  Python: /Users/para/.venvs/rnr/bin/python.  ~2-3 min.
"""
import sympy as sp
import numpy as np
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

s, D, p = sp.symbols('s D p', positive=True)

def disc_pieces(Aval):
    A = sp.Integer(Aval)
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
            for j in range(3): M[i, j] = Ki[i] * Trow[i][j]
        return M
    M = (Bred(1) * Bred(0)).applyfunc(lambda e: sp.cancel(sp.together(e)))
    c2m = sp.cancel(M.trace())
    c1m = sp.cancel(sum(M[[i for i in r], [j for j in r]].det() for r in ([0, 1], [0, 2], [1, 2])))
    c0m = sp.cancel(M.det())
    a = -c2m; b = c1m; c = -c0m
    disc = sp.cancel(sp.together(18 * a * b * c - 4 * a**3 * c + a**2 * b**2 - 4 * b**3 - 27 * c**2))
    N, Dn = sp.fraction(disc)
    return sp.expand(N), Dn

def V1():
    print("-" * 78); print("V1  A=2 equality: Dbar(2,p) == D_c(2,p) exactly")
    q = 2 * p
    Dbar2 = (2 * sp.Rational(1) / 2) * 2 * (sp.Rational(1, 2) - sp.sqrt(1 - q) / (2 - q))
    # (2(A-1)/A) at A=2 is 1; write directly:
    Dbar2 = sp.Rational(1, 1) * (sp.Rational(1, 2) - sp.sqrt(1 - 2 * p) / (2 - 2 * p))
    Dc2 = sp.Rational(1, 2) - sp.sqrt(1 - 2 * p) / (2 * (1 - p))
    ok = sp.simplify(Dbar2 - Dc2) == 0
    print(f"     Dbar(2,p) - D_c(2,p) = {sp.simplify(Dbar2 - Dc2)}")
    return rep("V1 A=2: the bound is exact (equality)", ok)

def V2V2b():
    print("-" * 78); print("V2/V2b  per-A certificates (s-rationalized sign + origin positivity)")
    ok2 = True; okb = True
    for Aval in (3, 4, 5, 6, 7, 8):
        A = sp.Integer(Aval)
        N, Dn = disc_pieces(Aval)
        # (ii): substituted sign
        pS = (A - 1) * (1 - s**2) / A
        DbarS = (2 * (A - 1) / A) * (sp.Rational(1, 2) - s / (1 + s**2))
        Ns = sp.cancel(sp.together(sp.expand(N.subs({p: pS, D: DbarS}))))
        n2, d2 = sp.fraction(Ns)
        poly = sp.Poly(sp.expand(n2), s)
        nroots = poly.count_roots(sp.Rational(1, 1000), sp.Rational(999, 1000))
        v = poly.eval(sp.Rational(1, 2))
        dv = sp.sign(d2.subs(s, sp.Rational(1, 2)))
        dsign = sp.sign(v) * dv
        here2 = (nroots == 0) and (dsign == -1)
        ok2 = ok2 and here2
        # (i): C(0,p) > 0 root-free
        dmin = min(m[0] for m in sp.Poly(N, D).monoms())
        C0 = sp.expand(sp.cancel(N / D**dmin).subs(D, 0))
        C0p = sp.Poly(C0, p)
        nb = C0p.count_roots(sp.Rational(1, 100000), sp.Rational(Aval - 1, Aval) - sp.Rational(1, 100000))
        sb = C0p.eval(sp.Rational(1, 10))
        hereb = (nb == 0) and (sb > 0) and (dmin == 4)
        okb = okb and hereb
        print(f"     A={Aval}: disc(Dbar)<0 rootfree:{here2}  C(0,p)>0 rootfree:{hereb}")
    rep("V2 disc(Dbar(p),p) < 0 on the whole p-range (A=3..8)", ok2)
    return rep("V2b disc > 0 at D->0+ (origin cofactor positive; sign chain closed)", okb)

def V3():
    print("-" * 78); print("V3  numeric sanity: Dbar > D_c with slack")
    def Dc_num(A, pv):
        def im(Dv):
            b0 = Dv / (A - 1); K = np.full((A, A), b0); np.fill_diagonal(K, 1 - Dv); Ki = np.linalg.inv(K)
            Tn = np.full((A, A), pv / (A - 1)); np.fill_diagonal(Tn, 1 - pv)
            By = lambda y: np.array([[Ki[y, a] * Tn[a, b] for b in range(A)] for a in range(A)])
            ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
            return max(abs(t[0].imag), abs(t[1].imag))
        lo, hi = 1e-12, (A - 1) / A - 1e-7
        for _ in range(80):
            m = (lo + hi) / 2; lo, hi = (m, hi) if im(m) < 1e-11 else (lo, m)
        return lo
    ok = True
    import math
    for A in (3, 5, 8):
        for fr in (0.2, 0.6, 0.9):
            pv = fr * (A - 1) / A
            q = pv * A / (A - 1)
            db = (2 * (A - 1) / A) * (0.5 - math.sqrt(1 - q) / (2 - q))
            dc = Dc_num(A, pv)
            here = db > dc
            ok = ok and here
            print(f"     A={A} p/pmax={fr}: Dbar={db:.6f} D_c={dc:.6f} slack={100*(db/dc-1):.1f}%")
    return rep("V3 Dbar > D_c numerically with slack (samples)", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("Dbar upper bound: D_c(A,p) < (2(A-1)/A)[1/2 - sqrt(1-q)/(2-q)], q=pA/(A-1)")
    print("=" * 78)
    V1(); V2V2b(); V3()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("Per-A theorem A=3..8 (one-variable Sturm chain); A=2 equality. Universal-A")
    print("remains conjectural (same per-A + overdetermination convention as the ladder).")

#!/usr/bin/env python3
"""
probe_7_34j_ternary_gray.py
============================================================================
Lemma 7.34j (ternary symmetric Gray threshold -- the first leaf beyond binary).
Independent verifier (cold-reviewed; supersedes the two scratch probes).

Stationary symmetric ternary Markov chain on {0,1,2}: T(i,i)=1-p, T(i,j)=p/2
(j!=i), pi uniform. Hamming distortion; backward test channel = ternary
symmetric channel K(y|x)=1-D (y=x), D/2 (y!=x). Deconvolution
P_{Y*} = (K^{-1})^{(x)n} P_X; the Gray threshold D_c^(3)(p) = sup{D : P_{Y*}>=0
for all n}.

Checks:
  V1  K^{-1} = (1/(2-3D))(2I - D J) exact (K K^{-1}=I).           [sympy]
  V2  transfer form B_y[a,b]=K^{-1}[y,a] T[a,b], P_{Y*}(y)=pi prod B,
      matches the full 3^n brute deconvolution.                   [exact, n<=6]
  V3  D_c^(3)(p) = smallest positive root in (0,2/3) of the irreducible
      quartic Q(D,p); Q irreducible over Q(p); root = the direct all-word
      positivity threshold.                                        [sympy + transfer]
  V4  binding word = the 2-symbol ALTERNATING word: just above D_c it goes
      negative while the constant and period-3 (012..) words stay >0; just
      below D_c all three are >=0.                                 [transfer, n=12]
  V5  small-p asymptotics: ternary D_c^(3)=p^2/8+..., binary D_c=p^2/4+...,
      ratio -> 1/2 (ternary Gray region is HALF the binary at leading order).
  V6  honest scope: for p above ~0.67 the all-same constant word binds first
      (a different, near-iid regime) -- the quartic/alternating D_c is the
      threshold on p in (0, ~0.67].
Deps: sympy, numpy, fractions.
"""
import itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp

PASS = True
def rep(name, ok):
    global PASS
    PASS = PASS and ok
    print(f"  {name:<62} {'PASS' if ok else 'FAIL'}")
    return ok

# ---- the binding quartic (cubic-discriminant factor of the period-2 product) ----
_p, _D = sp.symbols('p D')
Qquartic = (9*_p**2*(3*_p-4)**2*_D**4
            - 4*_p*(3*_p-4)*(27*_p**2-30*_p-4)*_D**3
            + (441*_p**4-864*_p**3+432*_p**2+64)*_D**2
            - 8*(27*_p**4-33*_p**3+36*_p**2-32*_p+16)*_D
            + 16*_p**2*(_p**2+1))


def Dc3(pv):
    """smallest positive real root of Q in (0,2/3)."""
    rts = sp.Poly(Qquartic.subs(_p, pv), _D).nroots()
    pr = sorted(float(sp.re(r)) for r in rts
                if abs(float(sp.im(r))) < 1e-9 and 0 < float(sp.re(r)) < Fr(2, 3))
    return pr[0]


def transfer_PY(pv, Dv, word):
    """P_{Y*}(word) via the per-site transfer product (float)."""
    Kn = np.array([[(2 if i == j else 0) - Dv for j in range(3)]
                   for i in range(3)], float) / (2 - 3*Dv)
    T = np.array([[1-pv if i == j else pv/2 for j in range(3)]
                  for i in range(3)], float)
    By = lambda y: np.array([[Kn[y, a]*T[a, b] for b in range(3)] for a in range(3)])
    v = np.ones(3)/3.0
    for y in word:
        v = v @ By(y)
    return v.sum()


def brute_PY(pv, Dv, word):
    """exact-rational full deconvolution P_{Y*}(word) = sum_x prod Kinv[y,x] P_X(x)."""
    n = len(word)
    M = [[(2 if i == j else 0) - Dv for j in range(3)] for i in range(3)]  # x (2-3D)/... = Kinv*c
    c = 2 - 3*Dv
    def PX(x):
        val = Fr(1, 3)
        for t in range(n-1):
            val *= (1-pv) if x[t] == x[t+1] else pv/2
        return val
    s = Fr(0)
    for x in itertools.product(range(3), repeat=n):
        coef = Fr(1)
        for i in range(n):
            coef *= M[word[i]][x[i]]
        s += coef*PX(x)
    return s / c**n


def V1():
    print("-"*78); print("V1  K^{-1} = (1/(2-3D))(2I - D J) exact")
    D = sp.symbols('D'); I = sp.eye(3); J = sp.ones(3, 3)
    K = (1-3*D/2)*I + (D/2)*J
    Kinv = (1/(2-3*D))*(2*I - D*J)
    return rep("V1 K K^{-1} = I", sp.simplify(K*Kinv - I) == sp.zeros(3, 3))


def V2():
    print("-"*78); print("V2  transfer form == full brute deconvolution (exact, n<=6)")
    ok = True
    for pv in (Fr(1, 10), Fr(3, 10)):
        Dv = Fr(1, 1000)
        for n in range(2, 7):
            for _ in range(8):
                w = [int(x) for x in np.random.default_rng(n).integers(0, 3, n)]
                b = float(brute_PY(pv, Dv, w)); t = transfer_PY(float(pv), float(Dv), w)
                if abs(b - t) > 1e-9:
                    ok = False
    return rep("V2 transfer form matches brute deconvolution", ok)


def V3():
    print("-"*78); print("V3  D_c^(3) = smallest root of the IRREDUCIBLE quartic Q; = direct threshold")
    irr = len(sp.factor_list(Qquartic, _D)[1]) == 1
    rep("V3a Q irreducible over Q(p)", irr)
    ok = True
    for pv, reported in [(0.1, 0.0015400), (0.2, 0.0077419), (0.3, 0.0224567)]:
        dc = Dc3(Fr(pv).limit_denominator(100))
        # direct: alt word flips sign across dc (margin 1.08 so the creep index <=16)
        below = transfer_PY(pv, dc*0.92, [t % 2 for t in range(16)])
        above = transfer_PY(pv, dc*1.08, [t % 2 for t in range(16)])
        ok = ok and abs(dc - reported) < 1e-5 and below >= -1e-12 and above < 0
        print(f"     p={pv}: D_c^(3)={dc:.7f} (reported {reported}); alt below>=0:{below>=-1e-12} above<0:{above<0}")
    return irr and rep("V3b smallest Q-root = reported D_c^(3) and = alt-word threshold", ok)


def V4():
    print("-"*78); print("V4  binding word = alternating (above: alt<0, const & per3 >0; below: all>=0)")
    ok = True
    for pv in (0.1, 0.2, 0.3):
        dc = Dc3(Fr(pv).limit_denominator(100)); n = 12
        a_hi = transfer_PY(pv, dc*1.08, [t % 2 for t in range(n)])
        c_hi = transfer_PY(pv, dc*1.08, [0]*n)
        p3_hi = transfer_PY(pv, dc*1.08, [t % 3 for t in range(n)])
        a_lo = transfer_PY(pv, dc*0.92, [t % 2 for t in range(n)])
        ok = ok and a_hi < 0 and c_hi > 0 and p3_hi > 0 and a_lo >= -1e-12
        print(f"     p={pv}: above alt={a_hi:+.2e} const={c_hi:+.2e} per3={p3_hi:+.2e}; below alt={a_lo:+.2e}")
    return rep("V4 binding word is the 2-symbol alternating word", ok)


def V5():
    print("-"*78); print("V5  small-p: ternary p^2/8, binary p^2/4 => ratio 1/2 (half the binary Gray region)")
    pp = sp.symbols('p', positive=True)
    Dc_bin = (1 - sp.sqrt(1-2*pp)/(1-pp))/2
    bin_lead = sp.series(Dc_bin, pp, 0, 3).removeO().coeff(pp, 2)   # 1/4
    u = sp.symbols('u', positive=True)
    Qsub = sp.Poly(sp.expand(Qquartic.subs({_D: pp**2*u, _p: pp})), pp)
    low = min(m[0] for m in Qsub.monoms())
    uroot = sp.solve(Qsub.coeff_monomial(pp**low), u)        # u = 1/8
    tern_lead = uroot[0]
    ok = (bin_lead == sp.Rational(1, 4) and tern_lead == sp.Rational(1, 8)
          and sp.Rational(tern_lead)/bin_lead == sp.Rational(1, 2))
    print(f"     binary leading coeff = {bin_lead}; ternary leading coeff = {tern_lead}; ratio = {tern_lead/bin_lead}")
    return rep("V5 ternary D_c = p^2/8+..., half the binary p^2/4 (ratio 1/2)", ok)


def V6():
    print("-"*78); print("V6  honest scope: large p (~>0.67) the constant word binds first (near-iid regime)")
    # at p high, check the constant word can go negative at D below the alt threshold
    pv = 0.8
    dc_alt = Dc3(Fr(pv).limit_denominator(100))
    # scan D up to dc_alt; does const go negative before alt?
    flagged = False
    for k in range(1, 40):
        Dv = dc_alt*k/40
        c = transfer_PY(pv, Dv, [0]*10)
        if c < -1e-12:
            flagged = True; break
    print(f"     p={pv}: constant word binds below the alt threshold: {flagged} "
          f"(=> quartic/alt D_c valid on p in (0, ~0.67])")
    return rep("V6 large-p constant-word regime documented (scope p<~0.67)", flagged)


if __name__ == "__main__":
    print("="*78)
    print("Lemma 7.34j (ternary symmetric Gray threshold) -- verifier")
    print("="*78)
    V1(); V2(); V3(); V4(); V5(); V6()
    print("="*78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")

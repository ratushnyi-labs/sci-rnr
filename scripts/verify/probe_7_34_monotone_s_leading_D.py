#!/usr/bin/env python3
r"""
probe_7_34_monotone_s_leading_D.py
============================================================================
BUG-009-D / Route B: FLATNESS-IN-s AT LEADING ORDER IN D -- and the honest
localization of the open monotonicity-in-s to O(D^2).

CORRECTED SCOPE (post adversarial review). The original version of this probe
claimed a "leading-D monotonicity certificate d/dt[e_2/t] <= 0". That certificate
is VACUOUS: e_2(t,p) is EXACTLY LINEAR IN t (verified symbolically here, W2), so
e_2/t is t-free and its t-derivative is identically zero. The truthful content at
this order is a FLATNESS LEMMA:
    1 - g_A(s; D) = e_1(t,p) D + e_2(t,p) D^2 + O(D^3),  t = 1-cos s,
    e_1(t,p) = 2t  EXACTLY (all A)   and   e_2(t,p) = t * q_A(p)  EXACTLY,
so
    R_A(s; D) = 2 + [2 + q_A(p)] D + O(D^2)   -- EXACTLY FLAT in s through O(D).
CONSEQUENCE (honest): at leading order in D the whole-s-interval inequality
R_A(s) >= 3/2 carries NO information beyond the endpoint theorem (any single s
determines the O(D) level). Genuine monotonicity-in-s for A >= 3 is decided at
O(D^2) and beyond -- i.e. by e_3(t,p)/t -- which is checked NUMERICALLY here
(favorable sign: e_3/t decreasing in t) but remains UNPROVEN. The A=2 case is
separately fully proven (exact, all D).

CHECKS:
  W1  e_1(t,p) = 2t exactly for A=2..5 (flat D->0 limit).
  W2  e_2(t,p)/t is t-FREE (exact symbolic): the flatness lemma at O(D); the
      old "certificate" d/dt[e_2/t] == 0 identically (vacuous as a monotonicity
      statement -- this check documents WHY).
  W3  A=2 anchor: q_2(p) = e_2/t matches the exact closed form's O(D) slope.
  W4  small-p poles: q_A(p) = -2(A-1)/p^2 + (8A-12)/p + O(1) (matches E_A at t=2).
  W5  the OPEN object: e_3(t,p)/t strictly decreasing in t on a numeric grid
      (A=3,4; the favorable direction) -- labeled NUMERIC, not proven.

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~2-3 min.
"""
import itertools
import sympy as sp
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

p, D, w, tsym = sp.symbols('p D w t')

def _pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def serD(e, n=4):
    return sp.expand(sp.cancel(sp.series(sp.together(e), D, 0, n).removeO()))

def derive(Aval):
    """returns (e1(t), e2(t)) exact for alphabet Aval"""
    A = sp.Integer(Aval)
    eth = D / ((A - 1) * (1 - D))
    def etaf(W):
        E = eth * W
        return ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))
    def Cf(W):
        E = eth * W
        return ((A - 1) * D * E + D - (A - 1)) / (A * D - (A - 1))
    eta = etaf(w); etb = etaf(1 / w)
    Cr = sp.cancel(Cf(w) * Cf(1 / w) / (Cf(1) * Cf(1)))
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(Aval)] for i in range(Aval)]
    st = list(itertools.product(range(Aval), repeat=3))
    npat = 5 if Aval >= 3 else 4; reps = [None] * npat
    for k, tr in enumerate(st):
        if reps[_pat(tr)] is None: reps[_pat(tr)] = k
    Q = sp.zeros(npat, npat)
    for a_ in range(npat):
        x, xp, xq = st[reps[a_]]
        for (y, yp, yq) in st:
            term = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x][y])
            if yp != y: term *= eta
            if yq != y: term *= etb
            Q[a_, _pat((y, yp, yq))] += term
    Q = Q.applyfunc(serD)
    Q0 = Q.applyfunc(lambda e: sp.expand(e).coeff(D, 0))
    QI = Q.applyfunc(lambda e: sp.expand(e).coeff(D, 1))
    QII = Q.applyfunc(lambda e: sp.expand(e).coeff(D, 2))
    QIII = Q.applyfunc(lambda e: sp.expand(e).coeff(D, 3))
    v = Q0[:, 0]
    l1 = sp.cancel((QI * v)[0])
    x1 = (QI * v - l1 * v)
    l2 = sp.cancel((QII * v)[0] + (QI * x1)[0])
    x2 = (QII * v + QI * x1 - l1 * x1 - l2 * v)
    l3 = sp.cancel((QIII * v)[0] + (QII * x1)[0] + (QI * x2)[0])
    g = serD(sp.expand(serD(Cr) * (1 + l1 * D + l2 * D**2 + l3 * D**3)))
    omg = sp.expand(1 - g)
    e1 = sp.cancel(sp.expand(omg).coeff(D, 1))
    e2 = sp.cancel(sp.expand(omg).coeff(D, 2))
    e3 = sp.cancel(sp.expand(omg).coeff(D, 3))
    wt = (1 - tsym) + sp.sqrt((1 - tsym)**2 - 1)   # w with w+1/w = 2(1-t)
    E1t = sp.simplify(sp.radsimp(sp.cancel(sp.together(e1).subs(w, wt))))
    E2t = sp.simplify(sp.radsimp(sp.cancel(sp.together(e2).subs(w, wt))))
    E3t = sp.simplify(sp.radsimp(sp.cancel(sp.together(e3).subs(w, wt))))
    return E1t, E2t, E3t

RES = {A: derive(A) for A in (2, 3, 4, 5)}

def W1():
    print("-" * 78); print("W1  e_1(t,p) = 2t exactly (flat D->0 limit), A=2..5")
    ok = True
    for A, (E1t, _, _) in RES.items():
        here = sp.simplify(E1t - 2 * tsym) == 0
        ok = ok and here
        print(f"     A={A}: e_1 = {sp.simplify(E1t)}  (=2t: {here})")
    return rep("W1 e_1 = 2t for every A (R -> 2 flat in s as D->0)", ok)

def W2W3():
    print("-" * 78); print("W2/W3  FLATNESS: e_2/t is t-FREE (symbolic); A=2 anchor")
    ok2 = True; ok3 = True
    for A, (_, E2t, _) in RES.items():
        q = sp.cancel(E2t / tsym)
        tfree = sp.simplify(sp.diff(q, tsym)) == 0
        ok2 = ok2 and tfree
        print(f"     A={A}: e_2/t t-free: {tfree}  (q_{A}(p) = {sp.sstr(sp.simplify(q))[:60]})")
    # A=2 anchor: q_2(p) from the exact closed form R_2 = 2 - E D + ...:
    # E = -(e2+4)/2 evaluated at t=2 must match the endpoint slope; equivalently
    # q_2 = e_2/t at A=2 equals the exact O(D)-slope coefficient of (1-g_2)/t.
    a, t = sp.symbols('a t', positive=True)
    kap2 = (1 - 2 * p)**2 / (p**2 * (1 - p)**2)
    F = 1 - 2 * a * t
    G = 8 * kap2 * a**2 * t * (1 - 4 * a + 2 * a**2 * t) / (1 - 4 * a)**2
    omg_exact = sp.expand(((1 + 2 * a * t) - sp.sqrt(F**2 + G)) / 2)   # 1-g = (1/2)[(2-F)-sqrt(F^2+G)]
    # e2^exact = coefficient of a^2 in series (a ~ D at leading order differs by (1-D): compare at O(D^2) carefully)
    ser_a = sp.series(omg_exact, a, 0, 3).removeO()
    e2_exact_in_a = sp.cancel(sp.expand(ser_a).coeff(a, 2) / t)
    q2 = sp.cancel(RES[2][1] / tsym)
    # a = D(1-D): omg(D) = e1a*a + e2a*a^2 = e1a*D + (e2a - e1a)*D^2 => q_2 = e2a - e1a where e1a=2 (t-normalized)
    diffq = sp.simplify(q2 - (e2_exact_in_a - 2))
    ok3 = diffq == 0
    print(f"     A=2 anchor: q_2 == closed-form (e2_a - 2): {ok3}")
    rep("W2 e_2 = t * q_A(p) exactly (flatness at O(D); old certificate vacuous)", ok2)
    return rep("W3 A=2 anchor matches the exact closed form", ok3)

def W4():
    print("-" * 78); print("W4  small-p structure: e_2/t = -2(A-1)/p^2 + (8A-12)/p + O(1)")
    ok = True
    for A, (_, E2t, _) in RES.items():
        f = sp.cancel(E2t / tsym + 2 * (A - 1) / p**2 - (8 * A - 12) / p)
        lim = sp.limit(f.subs(tsym, sp.Rational(1, 2)), p, 0)
        here = lim.is_finite is True or (lim not in (sp.oo, -sp.oo, sp.zoo))
        ok = ok and here
        print(f"     A={A}: residual p->0 limit (t=1/2) = {sp.nsimplify(lim)}  finite:{here}")
    return rep("W4 leading small-p poles match -2(A-1)/p^2 + (8A-12)/p", ok)

def W5():
    print("-" * 78); print("W5  OPEN object: d/dt[e_3/t] < 0 on a numeric grid (NUMERIC, not proven)")
    ok = True
    for A in (3, 4):
        E3t = RES[A][2]
        obj = sp.cancel(sp.diff(sp.cancel(E3t / tsym), tsym))
        worst = None
        for pv in (sp.Rational(1, 10), sp.Rational(3, 10), sp.Rational(1, 2)):
            if pv >= sp.Rational(A - 1, A): continue
            for k in range(1, 11):
                tv = sp.Rational(k, 5)
                val = sp.N(obj.subs({p: pv, tsym: tv}), 30)
                if val > 0: worst = (pv, tv, val)
        here = worst is None
        ok = ok and here
        print(f"     A={A}: d/dt[e_3/t] < 0 on grid: {here}" + (f"  VIOLATION {worst}" if worst else ""))
    return rep("W5 e_3/t decreasing in t (favorable direction; NUMERIC ONLY)", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("Flatness-in-s at leading D (e_2 = t q_A(p)); the open piece is e_3 (A>=3)")
    print("=" * 78)
    W1(); W2W3(); W4(); W5()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("HONEST SUMMARY: R_A(s;D) is exactly flat in s through O(D), so the leading-D")
    print("s-interval statement reduces to the endpoint theorem (which IS proven). Genuine")
    print("A>=3 monotonicity-in-s is decided at O(D^2) (e_3/t decreasing: numeric only,")
    print("unproven). A=2 is separately fully proven (exact, all D).")

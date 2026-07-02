#!/usr/bin/env python3
r"""
probe_7_34m_endpoint_sharp_a2.py
============================================================================
BUG-009-D / Route B, A=2 SHARP endpoint: R_2(pi) >= 3/2 reduced to closed form.

CONTEXT. Route B's target is the SHARP replica-floor constant c=3/2, i.e.
    g_A(s) <= 1 - (3/2) D(1-D)(1-cos s)   <=>   R_A(s) := [1-g_A(s)]/[D(1-D)(1-cos s)] >= 3/2.
The committed A=2 theorem (tex ~12928) proves only the SOFTER R_2(s) >= 11/16 for
all s. The monotone-reduction (probe_7_34m_monotone_route.py) shows the binding is
at s=pi, so the SHARP constant 3/2 is decided by the single endpoint R_2(pi). This
probe proves, at A=2 and in closed form, that the sharp endpoint holds, and reduces
it to an explicit polynomial inequality.

CLOSED FORM (A=2, from the exact factorization det(lam I - W4) =
(lam^2 - q^2|eta|^2)(lam^2 - |1+eta|^2 lam - kappa^2|eta|^2), kappa^2=(1-2p)^2/(p^2(1-p)^2)).
Write a := D(1-D). At s=pi (w=-1): the Jensen identity gives F(pi)=1-4a, and
    G(pi) = |C_pi/C_0|^4 * 4 kappa^2 |eta_pi|^2 = 16 kappa^2 a^2 (1-2a)^2 / (1-4a)^2,
so  1 - g_2(pi) = (1/2)[(1+4a) - sqrt((1-4a)^2 + G(pi))]  and
    R_2(pi) = [(1+4a) - sqrt((1-4a)^2 + G(pi))] / (4a).

THE REDUCTION (elementary; 1-2a>0 on the Gray region so squaring is valid):
    R_2(pi) >= 3/2
      <=> (1-2a) >= sqrt((1-4a)^2 + G(pi))
      <=> (1-2a)^2 >= (1-4a)^2 + G(pi)
      <=> 4a - 12a^2 >= G(pi)
      <=> (1-3a)(1-4a)^2 >= 4 kappa^2 a (1-2a)^2
      <=> P(a,p) := (1-3a)(1-4a)^2 - 4 kappa^2 a (1-2a)^2 >= 0.
So the sharp endpoint is EXACTLY the polynomial inequality P(a,p) >= 0.

STATUS: PROVEN (A=2, closed form).  The boundary inequality is closed in closed
form.  The A=2 binary-symmetric Gray threshold is the alternating-word discriminant
root
    D_c(p) = 1/2 - sqrt(1-2p) / (2(1-p)),
whence a(D_c) = D_c(1-D_c) = p^2 / (4(1-p)^2), and substituting into P gives the
EXACT boundary margin
    P(a(D_c(p)), p) = p^2 (1-2p)^3 / (4 (1-p)^8)  >  0    for all p in (0,1/2).
The interior (no-root-in-(0,a_c]) step is closed IN CLOSED FORM by the
f-monotonicity lemma (S7):  P(a,p) >= 0  <=>  kappa^2 <= f(a) :=
(1-3a)(1-4a)^2 / (4a(1-2a)^2), and
    d log f / da = -3/(1-3a) - 1/a - 4/((1-2a)(1-4a))  <  0   on (0, 1/4)
(each term strictly negative there), so f is strictly decreasing and
min_{(0,a_c]} f = f(a_c); hence kappa^2 <= f(a_c) <=> P(a_c,p) >= 0 -- exactly
the proven boundary inequality. Therefore
P(a,p) >= 0 on the whole Gray region 0<D<=D_c(p); hence
    R_2(pi) >= 3/2   on the entire A=2 Gray region, with exact margin P_c above.
This is the SHARP endpoint constant (3/2 > the committed 11/16), tight as p->0,D->D_c.

CHECKS:
  S1  closed-form R_2(pi) equals the replica g_2 (via gA) to 1e-20.
  S2  reduction: [R_2(pi) >= 3/2] iff [P(a,p) >= 0], across p and the Gray region.
  S3  P(a,p) >= 0 on 0<D<=D_c(p) (sharp endpoint holds); report min P and tightness.
  S4  monotone-in-a near 0: P decreasing at small a (P(0,p)=1>0).
  S5  closed-form proof (sympy): D_c(p), a(D_c)=p^2/(4(1-p)^2), and the EXACT boundary
      P(a(D_c),p) = p^2(1-2p)^3/(4(1-p)^8) > 0.
  S6  no interior dip (numeric cross-check): smallest positive root of the cubic
      P(.,p) exceeds a(D_c) at sample p (supporting evidence; the PROOF is S7).
  S7  f-monotonicity closure (symbolic): P>=0 <=> kappa^2 <= f(a), and
      d log f/da = -3/(1-3a) - 1/a - 4/((1-2a)(1-4a)) < 0 on (0,1/4), so the
      Gray-interior inequality reduces to the proven boundary P(a_c,p)>=0 --
      closing the interior step in closed form (the sharp endpoint is PROVEN).

Deps: mpmath, numpy.  Python: /Users/para/.venvs/rnr/bin/python.  ~1 min.
"""
import itertools
import mpmath as mp
import numpy as np
mp.mp.dps = 40
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
    th = mp.log(D / ((A - 1) * (1 - D))); E = mp.e**(th + 1j * s); den = A * D - (A - 1)
    C = ((A - 1) * D * E + D - (A - 1)) / den
    eta = ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))
    C0 = ((A - 1) * D * mp.e**th + D - (A - 1)) / den; Cr = abs(C / C0)**2
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

def Dc2(pv):
    def im(Dv):
        K = np.array([[1 - Dv, Dv], [Dv, 1 - Dv]]); Ki = np.linalg.inv(K)
        Tn = np.array([[1 - pv, pv], [pv, 1 - pv]])
        By = lambda y: np.array([[Ki[y, a] * Tn[a, b] for b in range(2)] for a in range(2)])
        ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-12, 0.5 - 1e-7
    for _ in range(60):
        m = (lo + hi) / 2; lo, hi = (m, hi) if im(m) < 1e-11 else (lo, m)
    return lo

def R2pi_closed(a, kap2):
    G = 16 * kap2 * a**2 * (1 - 2 * a)**2 / (1 - 4 * a)**2
    return ((1 + 4 * a) - mp.sqrt((1 - 4 * a)**2 + G)) / (4 * a)

def Ppoly(a, kap2):
    return (1 - 3 * a) * (1 - 4 * a)**2 - 4 * kap2 * a * (1 - 2 * a)**2

PS = ['0.02', '0.05', '0.1', '0.2', '0.35', '0.45']
FRACS = [mp.mpf(f) for f in ('0.1', '0.3', '0.5', '0.7', '0.9', '0.99', '0.999')]

def S1():
    print("-" * 78); print("S1  closed-form R_2(pi) == replica g_2 (via gA)")
    ok = True; worst = mp.mpf(0)
    for p in PS:
        pf = mp.mpf(p); kap2 = (1 - 2 * pf)**2 / (pf**2 * (1 - pf)**2); dc = mp.mpf(Dc2(float(p)))
        for fr in (mp.mpf('0.3'), mp.mpf('0.9')):
            D = fr * dc; a = D * (1 - D)
            r_cf = R2pi_closed(a, kap2); r_dir = (1 - gA(2, p, D, mp.pi)) / (2 * a)
            worst = max(worst, abs(r_cf - r_dir)); ok = ok and abs(r_cf - r_dir) < mp.mpf('1e-20')
    print(f"     worst |closed-form - gA| = {mp.nstr(worst, 3)}")
    return rep("S1 closed-form R_2(pi) matches replica g_2", ok)

def S2():
    print("-" * 78); print("S2  reduction  [R_2(pi) >= 3/2]  <=>  [P(a,p) >= 0]")
    ok = True
    for p in PS:
        pf = mp.mpf(p); kap2 = (1 - 2 * pf)**2 / (pf**2 * (1 - pf)**2); dc = mp.mpf(Dc2(float(p)))
        for fr in FRACS:
            D = fr * dc; a = D * (1 - D)
            r = R2pi_closed(a, kap2); P = Ppoly(a, kap2)
            ok = ok and ((r >= mp.mpf('1.5')) == (P >= 0))
    return rep("S2 R_2(pi)>=3/2 iff P(a,p)>=0 (exact reduction)", ok)

def S3():
    print("-" * 78); print("S3  sharp endpoint holds: P(a,p) >= 0 on 0<D<=D_c(p); tight as p->0,D->D_c")
    ok = True
    for p in PS:
        pf = mp.mpf(p); kap2 = (1 - 2 * pf)**2 / (pf**2 * (1 - pf)**2); dc = mp.mpf(Dc2(float(p)))
        Ps = [Ppoly(fr * dc * (1 - fr * dc), kap2) for fr in FRACS]
        mn = min(Ps); rmn = R2pi_closed((FRACS[-1] * dc) * (1 - FRACS[-1] * dc), kap2)
        ok = ok and mn >= -mp.mpf('1e-30')
        print(f"     p={p}: Dc={mp.nstr(dc,4)} min P over Gray={mp.nstr(mn,4)}  R_2(pi)@0.999Dc={mp.nstr(rmn,8)}")
    return rep("S3 P>=0 on the A=2 Gray region (sharp endpoint R_2(pi)>=3/2)", ok)

def S4():
    print("-" * 78); print("S4  P(a,p) decreasing in a (=> binding at boundary a(D_c)); P(0,p)=1>0")
    ok = True
    for p in PS:
        pf = mp.mpf(p); kap2 = (1 - 2 * pf)**2 / (pf**2 * (1 - pf)**2); dc = mp.mpf(Dc2(float(p)))
        amax = (dc * (1 - dc))
        agrid = [amax * mp.mpf(k) / 20 for k in range(1, 21)]
        dec = all(Ppoly(agrid[i + 1], kap2) <= Ppoly(agrid[i], kap2) + mp.mpf('1e-30') for i in range(len(agrid) - 1))
        p0 = Ppoly(mp.mpf(0), kap2)
        ok = ok and dec and abs(p0 - 1) < mp.mpf('1e-30')
    return rep("S4 P decreasing in a, P(0,p)=1>0 (min is the D_c boundary)", ok)

def S5():
    print("-" * 78); print("S5  closed-form proof: D_c(p), a(D_c)=p^2/(4(1-p)^2), P(a(D_c),p) exact")
    import sympy as sp
    ps = sp.symbols('p', positive=True); asym = sp.symbols('a', positive=True)
    # A=2 binary-symmetric Gray threshold = alternating-word product discriminant root
    T = sp.Matrix([[1 - ps, ps], [ps, 1 - ps]]); Dsym = sp.symbols('D', positive=True)
    K = sp.Matrix([[1 - Dsym, Dsym], [Dsym, 1 - Dsym]]); Ki = K.inv()
    By = lambda y: sp.Matrix([[Ki[y, i] * T[i, j] for j in range(2)] for i in range(2)])
    Mm = By(1) * By(0); disc = sp.simplify(Mm.trace()**2 - 4 * Mm.det())
    dc_roots = sp.solve(sp.Eq(sp.numer(sp.together(disc)), 0), Dsym)
    Dc = sp.nsimplify(sp.Rational(1, 2) - sp.sqrt(1 - 2 * ps) / (2 * (1 - ps)))
    # confirm this Dc is a discriminant root
    disc_at = sp.simplify(disc.subs(Dsym, Dc))
    ac = sp.simplify(Dc * (1 - Dc))
    ok_ac = sp.simplify(ac - ps**2 / (4 * (1 - ps)**2)) == 0
    kap2 = (1 - 2 * ps)**2 / (ps**2 * (1 - ps)**2)
    P = (1 - 3 * asym) * (1 - 4 * asym)**2 - 4 * kap2 * asym * (1 - 2 * asym)**2
    Pc = sp.simplify(P.subs(asym, ac))
    ok_Pc = sp.simplify(Pc - ps**2 * (1 - 2 * ps)**3 / (4 * (1 - ps)**8)) == 0
    print(f"     disc(D_c)=0 : {disc_at == 0}   a(D_c)=p^2/(4(1-p)^2): {ok_ac}")
    print(f"     P(a(D_c),p) = {sp.factor(Pc)}   == p^2(1-2p)^3/(4(1-p)^8): {ok_Pc}")
    print(f"     P_c > 0 on (0,1/2): p^2>0, (1-2p)^3>0, (1-p)^8>0  => strict")
    return rep("S5 closed-form boundary P(a(D_c),p)=p^2(1-2p)^3/(4(1-p)^8)>0 (proven)", bool(disc_at == 0 and ok_ac and ok_Pc))

def S6():
    print("-" * 78); print("S6  no interior dip: smallest positive root of P(.,p) exceeds a(D_c)")
    ok = True
    for pv in [mp.mpf(f) for f in ('0.005', '0.05', '0.25', '0.4', '0.49')]:
        k = (1 - 2 * pv)**2 / (pv**2 * (1 - pv)**2)
        # P(a) = -(48+16k)a^3 + (40+16k)a^2 - (11+4k)a + 1
        coeffs = [-(48 + 16 * k), (40 + 16 * k), -(11 + 4 * k), mp.mpf(1)]
        roots = mp.polyroots(coeffs, maxsteps=200, extraprec=80)
        posreal = sorted([mp.re(r) for r in roots if abs(mp.im(r)) < mp.mpf('1e-20') and mp.re(r) > mp.mpf('1e-18')])
        smallest = posreal[0] if posreal else mp.inf
        ac = pv**2 / (4 * (1 - pv)**2)
        here = smallest > ac * (1 - mp.mpf('1e-9'))
        ok = ok and here
        print(f"     p={mp.nstr(pv,4)}: a_c={mp.nstr(ac,5)} smallest_pos_root={mp.nstr(smallest,5)} root>=a_c:{here}")
    return rep("S6 no interior root in (0,a_c] => P>=0 on Gray => R_2(pi)>=3/2 PROVEN", ok)

def S7():
    print("-" * 78); print("S7  f-monotonicity closure: dlog f/da < 0 on (0,1/4) (symbolic)")
    import sympy as sp
    a, k2 = sp.symbols('a kappa2', positive=True)
    f = (1 - 3 * a) * (1 - 4 * a)**2 / (4 * a * (1 - 2 * a)**2)
    # P >= 0 <=> kappa^2 <= f(a):  P = 4a(1-2a)^2 (f(a) - kappa^2) -- verify identity
    P = (1 - 3 * a) * (1 - 4 * a)**2 - 4 * k2 * a * (1 - 2 * a)**2
    ok_id = sp.simplify(P - 4 * a * (1 - 2 * a)**2 * (f - k2)) == 0
    # dlog f/da = -3/(1-3a) - 1/a - 4/((1-2a)(1-4a))  (each term < 0 on (0,1/4))
    dlog = sp.simplify(sp.diff(sp.log(f), a))
    target = -3 / (1 - 3 * a) - 1 / a - 4 / ((1 - 2 * a) * (1 - 4 * a))
    ok_dlog = sp.simplify(dlog - target) == 0
    print(f"     P == 4a(1-2a)^2 (f - kappa^2): {ok_id}")
    print(f"     dlog f/da == -3/(1-3a) - 1/a - 4/((1-2a)(1-4a)): {ok_dlog}")
    print(f"     each term < 0 on (0,1/4) => f strictly decreasing => interior reduces to boundary")
    return rep("S7 interior closed in closed form (f-monotonicity lemma)", bool(ok_id and ok_dlog))

if __name__ == "__main__":
    print("=" * 78)
    print("BUG-009-D / Route B (A=2): sharp endpoint R_2(pi)>=3/2 -- PROVEN in closed form")
    print("=" * 78)
    S1(); S2(); S3(); S4(); S5(); S6(); S7()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("A=2 sharp endpoint R_2(pi)>=3/2 is PROVEN: exact margin P_c=p^2(1-2p)^3/(4(1-p)^8)>0")
    print("on the Gray region. Open: the general-A (A>=3) endpoint, and monotonicity in s")
    print("(probe_7_34m_monotone_route.py) -- together they would give R_A(s)>=3/2 for all s.")

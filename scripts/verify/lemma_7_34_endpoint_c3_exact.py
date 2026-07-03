#!/usr/bin/env python3
r"""
lemma_7_34_endpoint_c3_exact.py
============================================================================
BUG-009-D / Route B: the endpoint ladder extended to THIRD order, universal in A --
    delta_3(A) = (8A^3 - 32A^2 + 53A - 32) / (2 (A-1)^3),
    c_3(A)     = -(18A^2 - 153A + 232) / (8 (A-1)^3),
completing the three-term corner expansion
    R_A(pi)|_{D_c} - 3/2 = [(A-2)/(A-1)] p - [(4A^2+24A-65)/(8(A-1)^2)] p^2
                           - [(18A^2-153A+232)/(8(A-1)^3)] p^3 + O(p^4),
with D_c = [p^2/(4(A-1))] (1 + 2p + delta_2 p^2 + delta_3 p^3 + O(p^4)).
(Note c_3 > 0 for A <= 6 (c_3(7) = -43/1728 < 0) -- the p^3 term HELPS the
margin at small alphabets -- and c_3(2) = 1/4.)

PROVENANCE. Derived in the user-authorized research fan-out (agent artifacts in
scratchpad/research_p3_remainder/), then INDEPENDENTLY verified through three
routes before packaging: (a) the agent's symbolic-A derivations; (b) the agent's
numeric checks; (c) a from-scratch high-precision residual pipeline (this file's
Y4: mpmath 70 dps, discriminant-bisected D_c, margins Richardson-extrapolated) --
all agreeing to 7 significant figures, plus the exact A=2 anchor delta_3(2) = 5
computed by hand from the closed form D_c = 1/2 - sqrt(1-2p)/(2(1-p)).

GRADING REFINEMENT (design input for the ladder): lam_k ~ p^{-2 floor(k/2)}
(even-graded poles), so ladder term k enters the margin at order p^{2 ceil(k/2)-2}.
Hence lam_5 first contributes at O(p^4) and c_3 is carried ENTIRELY by lam_1..lam_4
-- the same KMAX=4 assembly as c_2, pushed one p-order deeper. (The earlier
"margin order p^k needs lam_{k+2}" rule is the crude even-k version of this.)

CHECKS:
  Y1  A=2 anchor: delta_3(2) = 5 from the exact closed form (symbolic series).
  Y2  delta_3 exact per-A (A=3..10) via the p^15 discriminant coefficient with
      d3-linear solve -- 8 exact points overdetermine the degree-3 numerator of
      the universal form (symbolic-A formulations exceed practical runtime; the
      research artifact holds a direct symbolic-A derivation).
  Y3  c_3 UNIVERSAL: the series-ring KMAX=4 assembly (pole-shifted coefficient
      dicts, one order deeper than the c_2 computation) yields margin p^3
      coefficient = c_3(A) for every A (exact, A symbolic).
  Y4  independent numerics: delta_3 residual (A=3,4,5) and c_3 margin residual
      (A=3,4) from the raw replica + discriminant-bisected D_c match the
      universal formulas (Richardson, 70 dps).

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~5-8 min.
"""
import itertools
import sympy as sp
import mpmath as mp
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

A, p, D, w, W = sp.symbols('A p D w W')
DELTA2 = (12 * A**2 - 28 * A + 21) / (4 * (A - 1)**2)
DELTA3 = (8 * A**3 - 32 * A**2 + 53 * A - 32) / (2 * (A - 1)**3)
C3 = -(18 * A**2 - 153 * A + 232) / (8 * (A - 1)**3)

def _pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def Y1():
    print("-" * 78); print("Y1  A=2 anchor: delta_3(2) = 5 from the exact closed form")
    p2 = sp.symbols('p2', positive=True)
    Dc2 = sp.Rational(1, 2) - sp.sqrt(1 - 2 * p2) / (2 * (1 - p2))
    ser = sp.series(Dc2 / (p2**2 / 4), p2, 0, 4).removeO()
    d3 = sp.simplify(sp.expand(ser).coeff(p2, 3))
    ok = (sp.simplify(d3 - 5) == 0) and (sp.simplify(DELTA3.subs(A, 2) - 5) == 0)
    print(f"     closed-form coefficient = {d3}; universal formula at A=2 = {sp.simplify(DELTA3.subs(A,2))}")
    return rep("Y1 delta_3(2) = 5 (exact anchor)", ok)

def Y2():
    print("-" * 78); print("Y2  delta_3 per-A exact (A=3..10) == universal rational form")
    # NOTE: three symbolic-A formulations of this discriminant step (raw rational,
    # series-ring, series-ring + W*-substituted + d3-linear) all exceeded practical
    # runtime (>30 min); the per-A route below is exact rational arithmetic in
    # seconds per alphabet. Universality evidence: (a) delta_3(A)*2(A-1)^3 is by
    # construction a polynomial in A of degree 3 whose four coefficients are
    # overdetermined by the 8 exact per-A values below (+ the A=2 anchor of Y1);
    # (b) the research-artifact symbolic-A derivation (research_p3_remainder/)
    # obtained the same closed form directly.
    NS = 17
    def serp(e):
        return sp.expand(sp.cancel(sp.series(sp.together(e), p, 0, NS).removeO()))
    def d3_of_A(Aval):
        Av = sp.Integer(Aval)
        Wstar = sp.Rational(1, 4 * (Aval - 1))
        d2v = sp.Rational(12 * Aval**2 - 28 * Aval + 21, 4 * (Aval - 1)**2)
        d3 = sp.symbols('d3v')
        def disc_p15(d3val):
            Dans = Wstar * p**2 * (1 + 2 * p + d2v * p**2 + d3val * p**3)
            lam2 = 1 - Dans * Av / (Av - 1)
            al = sp.Rational(1, Aval) + (1 - sp.Rational(1, Aval)) / lam2
            be = sp.Rational(1, Aval) - sp.Rational(1, Aval) / lam2
            td = 1 - p; to = p / (Av - 1)
            def Bred(y):
                Ki = {0: (al if y == 0 else be), 1: (be if y == 0 else al), 2: be}
                Trow = {0: {0: td, 1: to, 2: (Av - 2) * to},
                        1: {0: to, 1: td, 2: (Av - 2) * to},
                        2: {0: to, 1: to, 2: td + (Av - 3) * to}}
                M = sp.zeros(3, 3)
                for i in range(3):
                    for j in range(3):
                        M[i, j] = Ki[i] * Trow[i][j]
                return M
            M = (Bred(1) * Bred(0)).applyfunc(serp)
            def smul(x, y): return serp(sp.expand(x * y))
            c2m = serp(M.trace())
            def det2(i, j): return serp(M[i, i] * M[j, j] - M[i, j] * M[j, i])
            c1m = serp(det2(0, 1) + det2(0, 2) + det2(1, 2))
            c0m = serp(sp.expand(M[0, 0] * (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1])
                                 - M[0, 1] * (M[1, 0] * M[2, 2] - M[1, 2] * M[2, 0])
                                 + M[0, 2] * (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0])))
            a = -c2m; b = c1m; c = -c0m
            a2 = smul(a, a); a3 = smul(a2, a); b2 = smul(b, b); b3 = smul(b2, b); cc = smul(c, c)
            disc = serp(18 * smul(smul(a, b), c) - 4 * smul(a3, c) + smul(a2, b2) - 4 * b3 - 27 * cc)
            return sp.cancel(sp.expand(disc).coeff(p, 15))
        f0 = disc_p15(sp.Integer(0)); f1 = disc_p15(sp.Integer(1))
        return sp.cancel(-f0 / (f1 - f0))
    ok = True
    for Aval in range(3, 11):
        got = d3_of_A(Aval)
        tgt = sp.Rational(8 * Aval**3 - 32 * Aval**2 + 53 * Aval - 32, 2 * (Aval - 1)**3)
        here = sp.simplify(got - tgt) == 0
        ok = ok and here
        print(f"     A={Aval}: delta_3 = {got}  == universal form: {here}")
    return rep("Y2 delta_3 exact A=3..10 (8-point overdetermination of the cubic)", ok)

def Y3():
    print("-" * 78); print("Y3  c_3 universal via series-ring KMAX=4 assembly, one order deeper")
    et = sp.symbols('eta_s')
    NORD = 8
    def ser(e):
        return sp.expand(sp.cancel(sp.series(e, p, 0, NORD).removeO()))
    REPS = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (0, 1, 2)]
    FRESH = [3, 4, 5]
    def Tsym(a_, b_): return (1 - p) if a_ == b_ else p / (A - 1)
    Q0 = sp.zeros(5, 5); Q1 = sp.zeros(5, 5); Q2 = sp.zeros(5, 5)
    for a_, (x, xp, xq) in enumerate(REPS):
        for (y, yp, yq) in itertools.product(range(6), repeat=3):
            fr = [l for l in dict.fromkeys([y, yp, yq]) if l in FRESH]
            if fr != FRESH[:len(fr)]: continue
            mult = sp.Integer(1)
            for i in range(len(fr)): mult *= (A - 3 - i)
            t = mult * Tsym(xp, yp) * Tsym(xq, yq) / Tsym(x, y)
            ne = (1 if yp != y else 0) + (1 if yq != y else 0)
            [Q0, Q1, Q2][ne][a_, _pat((y, yp, yq))] += t
    Q0 = Q0.applyfunc(lambda e: sp.cancel(sp.together(e)))
    Q1 = Q1.applyfunc(lambda e: sp.expand(sp.cancel(sp.together(e))))
    Q2 = Q2.applyfunc(lambda e: sp.expand(sp.cancel(sp.together(e))))
    v = Q0[:, 0]
    xs = {0: v}; lams = {}
    for k in range(1, 5):
        r = (Q1 * xs[k - 1] + (Q2 * xs[k - 2] if k >= 2 else sp.zeros(5, 1))).applyfunc(ser)
        lam_k = ser(r[0]); xk = (r - lam_k * v)
        for j in range(1, k): xk = xk - lams[j] * xs[k - j]
        xs[k] = xk.applyfunc(ser); lams[k] = lam_k
    def laurent_c(e, pole, hi):
        ee = sp.expand(sp.cancel(sp.expand(e) * p**pole))
        return {k - pole: sp.cancel(ee.coeff(p, k)) for k in range(0, hi + pole + 1)}
    Wc = sp.Rational(1, 4) / (A - 1)
    Dc = Wc * p**2 * (1 + 2 * p + DELTA2 * p**2 + DELTA3 * p**3)
    etaD = 2 * D * (1 - D) / (A * D - 2 * D**2 - (A - 1))
    CrD = ((A * D - 2 * D**2 - (A - 1)) / (A * D - (A - 1)))**2
    etac_c = laurent_c(sp.series(etaD.subs(D, Dc), p, 0, 9).removeO(), 0, 8)
    Crc_c = laurent_c(sp.series(CrD.subs(D, Dc), p, 0, 7).removeO(), 0, 6)
    def conv(a_, b_, hi):
        out = {}
        for i, ai in a_.items():
            if ai == 0: continue
            for j, bj in b_.items():
                if i + j > hi or bj == 0: continue
                out[i + j] = sp.cancel(out.get(i + j, 0) + ai * bj)
        return out
    eta_pow = {1: etac_c}
    for k in (2, 3, 4): eta_pow[k] = conv(eta_pow[k - 1], etac_c, 12)
    lam_c = {k: laurent_c(lams[k], 4, 3) for k in range(1, 5)}
    lamc_c = {j: (sp.Integer(1) if j == 0 else sp.Integer(0)) for j in range(0, 6)}
    for k in range(1, 5):
        t = conv(lam_c[k], eta_pow[k], 5)
        for j in range(0, 6): lamc_c[j] = sp.cancel(lamc_c[j] + t.get(j, 0))
    Cl = conv(Crc_c, lamc_c, 5)
    Nc = {j: sp.cancel((1 if j == 0 else 0) - Cl.get(j, 0)) for j in range(0, 6)}
    inv_c = laurent_c(sp.series(1 / (2 * Dc * (1 - Dc)), p, 0, 4).removeO(), 2, 3)
    marg = conv(Nc, inv_c, 3)
    m0 = sp.cancel(marg.get(0, 0) - sp.Rational(3, 2))
    m1 = sp.simplify(marg.get(1, 0)); m2 = sp.simplify(marg.get(2, 0)); m3 = sp.simplify(marg.get(3, 0))
    ok0 = sp.simplify(m0) == 0
    ok1 = sp.simplify(m1 - (A - 2) / (A - 1)) == 0
    ok2 = sp.simplify(m2 + (4 * A**2 + 24 * A - 65) / (8 * (A - 1)**2)) == 0
    ok3 = sp.simplify(m3 - C3) == 0
    print(f"     p^0={sp.simplify(m0)}  p^1 ok:{ok1}  p^2 ok:{ok2}")
    print(f"     c_3(A) = {sp.factor(m3)}")
    print(f"     == -(18A^2-153A+232)/(8(A-1)^3): {ok3}")
    return rep("Y3 c_3 universal (series-ring assembly, exact)", bool(ok0 and ok1 and ok2 and ok3))

def Y4():
    print("-" * 78); print("Y4  independent numerics (raw replica + bisected D_c, 70 dps)")
    mp.mp.dps = 70
    def gA(Av, pv, Dv, s):
        Av = int(Av); pv = mp.mpf(pv); Dv = mp.mpf(Dv); s = mp.mpf(s)
        th = mp.log(Dv / ((Av - 1) * (1 - Dv))); E = mp.e**(th + 1j * s)
        C = ((Av - 1) * Dv * E + Dv - (Av - 1)) / (Av * Dv - (Av - 1))
        eta = ((Av - 1) * Dv * E + Dv - (Av - 1) * E) / ((Av - 1) * Dv * E + Dv - (Av - 1))
        C0 = ((Av - 1) * Dv * mp.e**th + Dv - (Av - 1)) / (Av * Dv - (Av - 1)); Cr = abs(C / C0)**2
        T = [[(1 - pv) if i == j else pv / (Av - 1) for j in range(Av)] for i in range(Av)]
        etb = mp.conj(eta)
        st = list(itertools.product(range(Av), repeat=3)); reps = [None] * 5
        for k, tr in enumerate(st):
            if reps[_pat(tr)] is None: reps[_pat(tr)] = k
        Q = mp.zeros(5, 5)
        for a_ in range(5):
            i = reps[a_]; x, xp, xq = st[i]
            for j, (y, yp, yq) in enumerate(st):
                Q[a_, _pat((y, yp, yq))] += T[xp][yp] * T[xq][yq] / T[x][y] * (eta if yp != y else 1) * (etb if yq != y else 1)
        ev, _ = mp.eig(Q); return Cr * max(abs(e) for e in ev)
    def disc3n(Av, pv, Dv):
        lam2 = 1 - Dv * Av / (Av - 1)
        al = mp.mpf(1) / Av + (1 - mp.mpf(1) / Av) / lam2
        be = mp.mpf(1) / Av - (mp.mpf(1) / Av) / lam2
        td = 1 - pv; to = pv / (Av - 1)
        def Bred(y):
            Ki = {0: (al if y == 0 else be), 1: (be if y == 0 else al), 2: be}
            Trow = {0: {0: td, 1: to, 2: (Av - 2) * to},
                    1: {0: to, 1: td, 2: (Av - 2) * to},
                    2: {0: to, 1: to, 2: td + (Av - 3) * to}}
            M = mp.zeros(3, 3)
            for i in range(3):
                for j in range(3): M[i, j] = Ki[i] * Trow[i][j]
            return M
        M = Bred(1) * Bred(0)
        c2 = M[0, 0] + M[1, 1] + M[2, 2]
        c1 = (M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]) + (M[0, 0] * M[2, 2] - M[0, 2] * M[2, 0]) + (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1])
        c0 = mp.det(M)
        a_ = -c2; b_ = c1; c_ = -c0
        return 18 * a_ * b_ * c_ - 4 * a_**3 * c_ + a_**2 * b_**2 - 4 * b_**3 - 27 * c_**2
    def Dc_hp(Av, pv):
        lo = pv**2 / (8 * (Av - 1)); hi = pv**2 / (2 * (Av - 1))
        flo = disc3n(Av, pv, lo)
        for _ in range(240):
            m = (lo + hi) / 2
            if disc3n(Av, pv, m) * flo > 0: lo = m
            else: hi = m
        return (lo + hi) / 2
    ok = True
    for Av in (3, 4, 5):
        d2 = mp.mpf(12 * Av**2 - 28 * Av + 21) / (4 * (Av - 1)**2)
        d3t = mp.mpf(8 * Av**3 - 32 * Av**2 + 53 * Av - 32) / (2 * (Av - 1)**3)
        vals = []
        for pv in ('2e-3', '1e-3', '5e-4'):
            pf = mp.mpf(pv); dc = Dc_hp(Av, pf)
            vals.append((dc / (pf**2 / (4 * (Av - 1))) - 1 - 2 * pf - d2 * pf**2) / pf**3)
        r = 2 * vals[-1] - vals[-2]
        here = abs(r - d3t) < mp.mpf('1e-2'); ok = ok and here
        print(f"     A={Av}: delta_3 -> {mp.nstr(r,8)} (target {mp.nstr(d3t,8)}): {here}")
    for Av in (3, 4):
        c1n = mp.mpf(Av - 2) / (Av - 1)
        c2n = -mp.mpf(4 * Av**2 + 24 * Av - 65) / (8 * (Av - 1)**2)
        c3t = -mp.mpf(18 * Av**2 - 153 * Av + 232) / (8 * (Av - 1)**3)
        vals = []
        for pv in ('4e-3', '2e-3', '1e-3'):
            pf = mp.mpf(pv); dc = Dc_hp(Av, pf)
            marg = (1 - gA(Av, pf, dc, mp.pi)) / (dc * (1 - dc) * 2) - mp.mpf('1.5')
            vals.append((marg - c1n * pf - c2n * pf**2) / pf**3)
        r1 = 2 * vals[1] - vals[0]; r2 = 2 * vals[2] - vals[1]; rr = 2 * r2 - r1
        here = abs(rr - c3t) < mp.mpf('2e-2'); ok = ok and here
        print(f"     A={Av}: c_3 -> {mp.nstr(rr,8)} (target {mp.nstr(c3t,8)}): {here}")
    return rep("Y4 independent numeric residuals match (70 dps, Richardson)", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("Endpoint ladder third order: delta_3(A), c_3(A) universal (exact)")
    print("=" * 78)
    Y1(); Y2(); Y3(); Y4()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("Three-term corner margin, every A. Next: the per-A EFFECTIVE remainder bound")
    print("(research artifact; packaging separately) and conjectured delta_4/c_4.")

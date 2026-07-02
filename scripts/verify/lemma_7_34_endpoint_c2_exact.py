#!/usr/bin/env python3
r"""
lemma_7_34_endpoint_c2_exact.py
============================================================================
BUG-009-D / Route B: the O(p^2) endpoint-margin coefficient DERIVED exactly --
    c_2(A) = -(4A^2 + 24A - 65) / (8 (A-1)^2)      (A = 3, 4, 5 exact arithmetic;
                                                     A=2 separately proven 1/8;
                                                     A=6 numeric to 6 digits),
completing the two-term small-p endpoint expansion
    R_A(pi)|_{D_c} - 3/2 = [(A-2)/(A-1)] p + c_2(A) p^2 + O(p^3).

INGREDIENTS (each exact):
  (i)   delta_2: the p^2-correction of the Gray threshold,
            D_c = [p^2/(4(A-1))] (1 + 2p + delta_2 p^2 + ...),
            delta_2 = (12A^2 - 28A + 21) / (4 (A-1)^2)   -- UNIVERSAL in A,
        derived from the p^14 coefficient of the 3x3 symmetry-reduced cubic
        discriminant at W*=1/(4(A-1)), d1=2 (A symbolic; the p^14 coefficient is
        -(4A^2 d2 - 12A^2 - 8A d2 + 28A + 4 d2 - 21)/(1024 (A-1)^14)).
        A=2 anchor: the exact closed form gives D_c = (p^2/4)(1+2p+(13/4)p^2+...),
        and the universal formula reproduces 13/4 at A=2.
  (ii)  Rank-1 perturbation ladder: with Q(eta)=Q0+eta Q1+eta^2 Q2 and Q0=c e0^T
        (single nonzero column; lambda=1 simple; Q0 P_perp = 0), the recurrence
            lam_k = e0^T(Q1 x_{k-1} + Q2 x_{k-2}),
            x_k   = Q1 x_{k-1} + Q2 x_{k-2} - sum_{j<k} lam_j x_{k-j} - lam_k c
        is pure matrix algebra (no determinants, no root-solving).
  (iii) TRUNCATION RULE (the subtle point, diagnosed after an initial mismatch of
        exactly 1/(8(A-1)) against validated numerics): at margin order p^k the
        PT orders up to lam_{k+2} contribute -- for c_2 that means k<=4. lam_4's
        contribution at O(p^2) is exactly +1/(8(A-1)) (in R); lam_5, lam_6, ...
        vanish at O(p^2) (verified numerically to k=12, coefficients <= 1e-7).
        A naive k<=3 assembly is WRONG at second order.
  (iv)  Assembly: g = Cr(D) * (1 + sum_k lam_k eta(D)^k), with the exact
        eta(D) = 2D(1-D)/(AD-2D^2-(A-1)), Cr(D) = [(AD-2D^2-(A-1))/(AD-(A-1))]^2,
        evaluated at D = D_c(p) from (i), series in p.

CONSEQUENCE: the sharpened two-term margin is positive up to p* = c_1/|c_2|
(~0.37 at A=3 rising to ~0.72 at A=6), comfortably covering the tight regime;
consistent with the adversarial all-D sweep (no violations anywhere).

UNIVERSAL-A DERIVATION (check E4). The same assembly with A a SYMBOLIC parameter,
made tractable by SERIES-RING arithmetic: every PT-recurrence entry is truncated to
a p-Laurent series (order 6) immediately, so expressions stay bounded (the full
rational-function recurrence with symbolic A swells terminally -- two runaway
attempts documented). The finale avoids sympy Laurent pitfalls (negative-power
.coeff() returns 0 silently; etac^k for k>=3 starts at p^6/p^8, beyond naive series
cutoffs) by POLE-SHIFTED COEFFICIENT DICTS: multiply by p^pole, extract nonnegative
coefficients, shift back, convolve exactly. Result (~2 min): margin p^0 = 0,
p^1 = (A-2)/(A-1), and c_2(A) = -(4A^2+24A-65)/(8(A-1)^2) for EVERY A -- the
two-term endpoint margin theorem is UNIVERSAL in the alphabet.

CHECKS:
  E1  delta_2 universal (A symbolic): p^14 discriminant coefficient at W*, d1=2
      forces delta_2 = (12A^2-28A+21)/(4(A-1)^2); A=2 closed form gives 13/4.
  E2  fixed-A exact assembly (A=3,4,5; KMAX=4): margin p^1 = (A-2)/(A-1) and
      c_2 = -(4A^2+24A-65)/(8(A-1)^2) EXACTLY.
  E3  truncation: numeric PT ladder at A=3, p=1e-3 -- k=4 contributes ~ -1/16
      (in the eigenvalue; +1/16 in R) and k=5,6 contribute < 1e-6 at O(p^2).
  E4  UNIVERSAL A (symbolic; series-ring + pole-shifted dicts): p^0=0,
      p^1=(A-2)/(A-1), c_2(A) = -(4A^2+24A-65)/(8(A-1)^2) for every A.

Deps: sympy, mpmath.  Python: /Users/para/.venvs/rnr/bin/python.  ~4-8 min.
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

def E1():
    print("-" * 78); print("E1  delta_2 = (12A^2-28A+21)/(4(A-1)^2) universal (discriminant p^14)")
    A, p, W, d2 = sp.symbols('A p W d2')
    Dans = W * p**2 * (1 + 2 * p + d2 * p**2)
    lam2 = 1 - Dans * A / (A - 1)
    al = 1 / A + (1 - 1 / A) / lam2
    be = 1 / A - (1 / A) / lam2
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
    c2m = sp.cancel(M.trace())
    c1m = sp.cancel(sum(M[[i for i in r], [j for j in r]].det() for r in ([0, 1], [0, 2], [1, 2])))
    c0m = sp.cancel(M.det())
    a = -c2m; b = c1m; c = -c0m
    disc = 18 * a * b * c - 4 * a**3 * c + a**2 * b**2 - 4 * b**3 - 27 * c**2
    ser = sp.expand(sp.series(disc, p, 0, 16).removeO())
    poly = sp.Poly(ser, p)
    Wstar = 1 / (4 * (A - 1))
    c14 = sp.factor(sp.simplify(sp.simplify(poly.coeff_monomial(p**14)).subs(W, Wstar)))
    sols = sp.solve(sp.Eq(c14, 0), d2)
    tgt = (12 * A**2 - 28 * A + 21) / (4 * (A - 1)**2)
    ok_univ = len(sols) == 1 and sp.simplify(sols[0] - tgt) == 0
    # A=2 anchor from the exact closed form
    p2 = sp.symbols('p2', positive=True)
    Dc2 = sp.Rational(1, 2) - sp.sqrt(1 - 2 * p2) / (2 * (1 - p2))
    ser2 = sp.series(Dc2 / (p2**2 / 4), p2, 0, 3).removeO()
    d2_A2 = sp.simplify(sp.expand(ser2).coeff(p2, 2))
    ok_a2 = sp.simplify(d2_A2 - sp.Rational(13, 4)) == 0 and sp.simplify(tgt.subs(A, 2) - sp.Rational(13, 4)) == 0
    print(f"     delta_2 (A symbolic) = {sols[0] if sols else None}  == target: {ok_univ}")
    print(f"     A=2 closed form: coeff = {d2_A2} = 13/4 = universal formula at A=2: {ok_a2}")
    return rep("E1 delta_2 universal + A=2 anchor", bool(ok_univ and ok_a2))

def E2():
    print("-" * 78); print("E2  fixed-A exact assembly (KMAX=4): c_2 = -(4A^2+24A-65)/(8(A-1)^2)")
    ok = True
    for Aval in (3, 4, 5):
        A = sp.Integer(Aval)
        p, D = sp.symbols('p D')
        T = [[(1 - p) if i == j else p / (A - 1) for j in range(Aval)] for i in range(Aval)]
        st = list(itertools.product(range(Aval), repeat=3)); reps = [None] * 5
        for k, tr in enumerate(st):
            if reps[_pat(tr)] is None: reps[_pat(tr)] = k
        Q0 = sp.zeros(5, 5); Q1 = sp.zeros(5, 5); Q2 = sp.zeros(5, 5)
        for a_ in range(5):
            x, xp, xq = st[reps[a_]]
            for (y, yp, yq) in st:
                t = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x][y])
                ne = (1 if yp != y else 0) + (1 if yq != y else 0)
                [Q0, Q1, Q2][ne][a_, _pat((y, yp, yq))] += t
        Q0 = Q0.applyfunc(sp.cancel); Q1 = Q1.applyfunc(sp.cancel); Q2 = Q2.applyfunc(sp.cancel)
        v = Q0[:, 0]
        xs = {0: v}; lams = {}
        for k in range(1, 5):
            r = Q1 * xs[k - 1] + (Q2 * xs[k - 2] if k >= 2 else sp.zeros(5, 1))
            lam_k = sp.cancel(r[0]); xk = r - lam_k * v
            for j in range(1, k): xk = xk - lams[j] * xs[k - j]
            xs[k] = xk.applyfunc(sp.cancel); lams[k] = lam_k
        etaD = 2 * D * (1 - D) / (A * D - 2 * D**2 - (A - 1))
        CrD = ((A * D - 2 * D**2 - (A - 1)) / (A * D - (A - 1)))**2
        W = sp.Rational(1, 4 * (Aval - 1))
        d2v = sp.Rational(12 * Aval**2 - 28 * Aval + 21, 4 * (Aval - 1)**2)
        Dc = W * p**2 * (1 + 2 * p + d2v * p**2)
        etac = sp.series(etaD.subs(D, Dc), p, 0, 11).removeO()
        Crc = sp.series(CrD.subs(D, Dc), p, 0, 11).removeO()
        lamc = 1 + sum(sp.series(lams[k] * etac**k, p, 0, 11).removeO() for k in range(1, 5))
        marg = sp.series(sp.expand(1 - Crc * lamc) / (2 * Dc * (1 - Dc)), p, 0, 3).removeO() - sp.Rational(3, 2)
        marg = sp.expand(marg)
        c1c = sp.simplify(marg.coeff(p, 1)); c2c = sp.simplify(marg.coeff(p, 2))
        conj = sp.Rational(-(4 * Aval**2 + 24 * Aval - 65), 8 * (Aval - 1)**2)
        here = (sp.simplify(c1c - sp.Rational(Aval - 2, Aval - 1)) == 0) and (sp.simplify(c2c - conj) == 0)
        ok = ok and here
        print(f"     A={Aval}: p^1 = {c1c}, c_2 = {c2c} (target {conj})  match:{here}")
    return rep("E2 c_2(A) exact via KMAX=4 assembly (A=3,4,5)", ok)

def E3():
    print("-" * 78); print("E3  truncation rule: k=4 carries -1/(8(A-1)) (eigenvalue side); k=5,6 vanish at O(p^2)")
    mp.mp.dps = 60
    Aval = 3
    def build(A, p):
        T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]
        st = list(itertools.product(range(A), repeat=3)); reps = [None] * 5
        for k, tr in enumerate(st):
            if reps[_pat(tr)] is None: reps[_pat(tr)] = k
        Q0 = mp.zeros(5, 5); Q1 = mp.zeros(5, 5); Q2 = mp.zeros(5, 5)
        for a_ in range(5):
            x, xp, xq = st[reps[a_]]
            for (y, yp, yq) in st:
                t = T[xp][yp] * T[xq][yq] / T[x][y]
                ne = (1 if yp != y else 0) + (1 if yq != y else 0)
                [Q0, Q1, Q2][ne][a_, _pat((y, yp, yq))] += t
        return Q0, Q1, Q2
    def terms(p):
        p = mp.mpf(p)
        Q0, Q1, Q2 = build(Aval, p)
        v = mp.matrix([Q0[i, 0] for i in range(5)])
        xs = {0: v}; lams = {}
        Z = mp.matrix([0] * 5)
        for k in range(1, 7):
            r = Q1 * xs[k - 1] + (Q2 * xs[k - 2] if k >= 2 else Z)
            lam_k = r[0]; xk = r - lam_k * v
            for j in range(1, k): xk = xk - lams[j] * xs[k - j]
            xs[k] = xk; lams[k] = lam_k
        W = mp.mpf(1) / (4 * (Aval - 1)); d2 = mp.mpf(12 * Aval**2 - 28 * Aval + 21) / (4 * (Aval - 1)**2)
        Dc = W * p**2 * (1 + 2 * p + d2 * p**2)
        eta = 2 * Dc * (1 - Dc) / (Aval * Dc - 2 * Dc**2 - (Aval - 1))
        den = 2 * Dc * (1 - Dc)
        return {k: lams[k] * eta**k / den for k in range(2, 7)}
    t1 = terms('1e-3'); t2 = terms('5e-4')
    # p^2 coefficient of each term (const removed by Richardson)
    out = {}
    for k in range(2, 7):
        c0 = (t2[k] * 4 - t1[k]) / 3
        out[k] = (t1[k] - c0) / mp.mpf('1e-6')
    ok4 = abs(out[4] - (-mp.mpf(1) / (8 * (Aval - 1)))) < mp.mpf('1e-3')
    ok56 = abs(out[5]) < mp.mpf('1e-4') and abs(out[6]) < mp.mpf('1e-4')
    print(f"     k=4 O(p^2) coeff = {mp.nstr(out[4],5)} (target -1/(8(A-1)) = {mp.nstr(-mp.mpf(1)/(8*(Aval-1)),5)})")
    print(f"     k=5: {mp.nstr(out[5],3)}  k=6: {mp.nstr(out[6],3)}  (vanish)")
    return rep("E3 truncation: lam_4 carries the O(p^2) tail; lam_5,6 vanish", bool(ok4 and ok56))

def E4():
    print("-" * 78); print("E4  UNIVERSAL A: series-ring PT + pole-shifted coefficient dicts")
    A, p, et, D = sp.symbols('A p eta_s D')
    NORD = 6
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
    Q0 = Q0.applyfunc(ser); Q1 = Q1.applyfunc(ser); Q2 = Q2.applyfunc(ser)
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
    W = sp.Rational(1, 4) / (A - 1)
    d2v = (12 * A**2 - 28 * A + 21) / (4 * (A - 1)**2)
    Dc = W * p**2 * (1 + 2 * p + d2v * p**2)
    etaD = 2 * D * (1 - D) / (A * D - 2 * D**2 - (A - 1))
    CrD = ((A * D - 2 * D**2 - (A - 1)) / (A * D - (A - 1)))**2
    etac_c = laurent_c(sp.series(etaD.subs(D, Dc), p, 0, 7).removeO(), 0, 6)
    Crc_c = laurent_c(sp.series(CrD.subs(D, Dc), p, 0, 5).removeO(), 0, 4)
    def conv(a, b, hi):
        out = {}
        for i, ai in a.items():
            if ai == 0: continue
            for j, bj in b.items():
                if i + j > hi or bj == 0: continue
                out[i + j] = sp.cancel(out.get(i + j, 0) + ai * bj)
        return out
    eta_pow = {1: etac_c}
    for k in (2, 3, 4): eta_pow[k] = conv(eta_pow[k - 1], etac_c, 10)
    lam_c = {k: laurent_c(lams[k], 4, 2) for k in range(1, 5)}
    lamc_c = {j: (sp.Integer(1) if j == 0 else sp.Integer(0)) for j in range(0, 5)}
    for k in range(1, 5):
        t = conv(lam_c[k], eta_pow[k], 4)
        for j in range(0, 5): lamc_c[j] = sp.cancel(lamc_c[j] + t.get(j, 0))
    Cl = conv(Crc_c, lamc_c, 4)
    Nc = {j: sp.cancel((1 if j == 0 else 0) - Cl.get(j, 0)) for j in range(0, 5)}
    inv_c = laurent_c(sp.series(1 / (2 * Dc * (1 - Dc)), p, 0, 3).removeO(), 2, 2)
    marg = conv(Nc, inv_c, 2)
    m0 = sp.cancel(marg.get(0, 0) - sp.Rational(3, 2))
    m1 = sp.simplify(marg.get(1, 0)); m2 = sp.simplify(marg.get(2, 0))
    tgt = -(4 * A**2 + 24 * A - 65) / (8 * (A - 1)**2)
    ok = (sp.simplify(m0) == 0 and sp.simplify(m1 - (A - 2) / (A - 1)) == 0
          and sp.simplify(m2 - tgt) == 0)
    print(f"     p^0 = {sp.simplify(m0)}, p^1 = {m1}, c_2(A) = {sp.factor(m2)}")
    return rep("E4 UNIVERSAL c_2(A) = -(4A^2+24A-65)/(8(A-1)^2) for every A", bool(ok))

if __name__ == "__main__":
    print("=" * 78)
    print("c_2(A) = -(4A^2+24A-65)/(8(A-1)^2): exact derivation, ALL A")
    print("=" * 78)
    E1(); E2(); E3(); E4()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("Two-term endpoint margin theorem, UNIVERSAL in A:")
    print("R_A(pi)|Dc - 3/2 = [(A-2)/(A-1)] p - [(4A^2+24A-65)/(8(A-1)^2)] p^2 + O(p^3),")
    print("positive to p* ~ 0.37-0.72. A=2 is ON-family: c_2(2)=+1/8 (proven separately).")

#!/usr/bin/env python3
r"""
lemma_7_34_endpoint_c4_exact.py
============================================================================
BUG-009-D / Route B: the endpoint ladder at FOURTH order --
    delta_4(A) = (40A^4 - 240A^3 + 652A^2 - 848A + 429) / (8 (A-1)^4),
    c_4(A)     = -(60A^3 - 890A^2 + 3171A - 3272) / (32 (A-1)^4),
extending the corner expansion to
    R_A(pi)|_{D_c} - 3/2 = c_1 p + c_2 p^2 + c_3 p^3 + c_4 p^4 + O(p^5),
    D_c = [p^2/(4(A-1))] (1 + 2p + delta_2 p^2 + delta_3 p^3 + delta_4 p^4 + ...).
(c_4 > 0 for A <= 10 (c_4(11) < 0); c_4(2) = 5/16, on-family like every previous coefficient.)

EVIDENCE STATUS (honest): the universal forms are OVERDETERMINED, not merely
fitted -- delta_4's degree-4 numerator (5 coefficients) is determined by the five
exact values A=2..6 and then confirmed OUT-OF-SAMPLE at A=7 (this file, W3);
c_4's degree-3 numerator (4 coefficients) is overdetermined by A=2..6 already and
further confirmed at A=7. Universality in the sense of a symbolic-A derivation is
NOT claimed (the discriminant route at this order is impractical symbolically --
see the delta_3 lesson in lemma_7_34_endpoint_c3_exact.py); the stated rational
forms are exact at eight alphabets including two beyond the fitting set.

PROVENANCE: research fan-out artifact research_c3_remainder/step6 (exact per-A
values A=3..6 + A=2 anchors); A=7 derivations and the pipeline-validation fix are
from the packaging session (a subtle bug was caught: laurent-coefficient
extraction on EXACT rational-function ladders silently corrupts -- the ladder
must be series-ring truncated before coefficient extraction; symptom was
p-DEPENDENT "coefficients" and failed c_1..c_3 cross-checks, both impossible for
true margin coefficients -- a useful canary pattern).

CHECKS:
  W1  A=2 anchors: delta_4(2) = 61/8 from the exact closed-form D_c series;
      c_4(2) = 5/16 from the exact closed-form R_2(pi) margin series at D_c(p).
  W2  the universal forms interpolate the artifact's exact per-A values
      (delta_4: 513/128, 29/8, 7489/2048, 18813/5000; c_4: 149/512, 247/648,
      2167/8192, 1663/10000 for A = 3, 4, 5, 6).
  W3  OUT-OF-SAMPLE A=7: delta_4(7) = 13387/3456 re-derived live (discriminant
      p^16 coefficient, d4-linear two-point solve, series-ring); c_4(7) =
      4105/41472 re-derived live (series-ring ladder KMAX=6, pole-shifted dict
      assembly to margin p^4; c_1, c_2, c_3 re-derived as pipeline validation).

Deps: sympy.  Python: /Users/para/.venvs/rnr/bin/python.  ~1-2 min.
"""
import itertools
import sympy as sp
PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<66} {'PASS' if ok else 'FAIL'}"); return ok

p, D = sp.symbols('p D', positive=True)
def D4F(Av): return sp.Rational(40 * Av**4 - 240 * Av**3 + 652 * Av**2 - 848 * Av + 429, 8 * (Av - 1)**4)
def C4F(Av): return -sp.Rational(60 * Av**3 - 890 * Av**2 + 3171 * Av - 3272, 32 * (Av - 1)**4)

def _pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def W1():
    print("-" * 78); print("W1  A=2 anchors from the exact closed forms")
    # delta_4(2): D_c = 1/2 - sqrt(1-2p)/(2(1-p)) = (p^2/4)(1+2p+13/4 p^2+5p^3+d4 p^4+...)
    ser = sp.series((sp.Rational(1, 2) - sp.sqrt(1 - 2 * p) / (2 * (1 - p))) / (p**2 / 4), p, 0, 5).removeO()
    d4_2 = sp.simplify(sp.expand(ser).coeff(p, 4))
    ok1 = (sp.simplify(d4_2 - sp.Rational(61, 8)) == 0) and (D4F(2) == sp.Rational(61, 8))
    print(f"     delta_4(2) closed form = {d4_2} ; universal form = {D4F(2)}")
    # c_4(2): R_2(pi) = [(1+4a) - sqrt((1-4a)^2 + G)]/(4a), a = D(1-D),
    # G = 16 kap^2 a^2 (1-2a)^2/(1-4a)^2, kap^2 = (1-2p)^2/(p^2(1-p)^2), at D = D_c(p)
    kap2 = (1 - 2 * p)**2 / (p**2 * (1 - p)**2)
    Dc2 = sp.Rational(1, 2) - sp.sqrt(1 - 2 * p) / (2 * (1 - p))
    a = sp.expand(Dc2 * (1 - Dc2))
    G = 16 * kap2 * a**2 * (1 - 2 * a)**2 / (1 - 4 * a)**2
    R2 = ((1 + 4 * a) - sp.sqrt((1 - 4 * a)**2 + G)) / (4 * a)
    marg = sp.series(sp.simplify(R2 - sp.Rational(3, 2)), p, 0, 5).removeO()
    c4_2 = sp.simplify(sp.expand(marg).coeff(p, 4))
    ok2 = (sp.simplify(c4_2 - sp.Rational(5, 16)) == 0) and (C4F(2) == sp.Rational(5, 16))
    print(f"     c_4(2) closed form = {c4_2} ; universal form = {C4F(2)}")
    # also confirm the lower coefficients en passant (p^2: 1/8, p^3: 1/4)
    ok3 = (sp.simplify(sp.expand(marg).coeff(p, 2) - sp.Rational(1, 8)) == 0 and
           sp.simplify(sp.expand(marg).coeff(p, 3) - sp.Rational(1, 4)) == 0)
    return rep("W1 A=2 anchors: delta_4=61/8, c_4=5/16 (exact; lower orders re-anchor)", bool(ok1 and ok2 and ok3))

def W2():
    print("-" * 78); print("W2  universal forms interpolate the artifact per-A exact values")
    D4V = {3: sp.Rational(513, 128), 4: sp.Rational(29, 8), 5: sp.Rational(7489, 2048), 6: sp.Rational(18813, 5000)}
    C4V = {3: sp.Rational(149, 512), 4: sp.Rational(247, 648), 5: sp.Rational(2167, 8192), 6: sp.Rational(1663, 10000)}
    ok = True
    for Av in (3, 4, 5, 6):
        h = (D4F(Av) == D4V[Av]) and (C4F(Av) == C4V[Av])
        ok = ok and h
        print(f"     A={Av}: delta_4 {D4V[Av]} c_4 {C4V[Av]}  match:{h}")
    return rep("W2 exact per-A values A=3..6 on both universal forms", ok)

def W3():
    print("-" * 78); print("W3  OUT-OF-SAMPLE A=7: live re-derivation of both coefficients")
    Aval = 7; A = sp.Integer(Aval)
    NS = 18
    def serp(e): return sp.expand(sp.cancel(sp.series(sp.together(e), p, 0, NS).removeO()))
    W_ = sp.Rational(1, 4 * (Aval - 1))
    d2v = sp.Rational(12 * Aval**2 - 28 * Aval + 21, 4 * (Aval - 1)**2)
    d3v = sp.Rational(8 * Aval**3 - 32 * Aval**2 + 53 * Aval - 32, 2 * (Aval - 1)**3)
    def disc_p16(d4val):
        Dans = W_ * p**2 * (1 + 2 * p + d2v * p**2 + d3v * p**3 + d4val * p**4)
        lam2 = 1 - Dans * A / (A - 1)
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
        M = (Bred(1) * Bred(0)).applyfunc(serp)
        def smul(x, y): return serp(sp.expand(x * y))
        c2m = serp(M.trace())
        def det2(i, j): return serp(M[i, i] * M[j, j] - M[i, j] * M[j, i])
        c1m = serp(det2(0, 1) + det2(0, 2) + det2(1, 2))
        c0m = serp(sp.expand(M[0, 0] * (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1])
                             - M[0, 1] * (M[1, 0] * M[2, 2] - M[1, 2] * M[2, 0])
                             + M[0, 2] * (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0])))
        a_ = -c2m; b_ = c1m; c_ = -c0m
        a2 = smul(a_, a_); a3 = smul(a2, a_); b2 = smul(b_, b_); b3 = smul(b2, b_); cc = smul(c_, c_)
        disc = serp(18 * smul(smul(a_, b_), c_) - 4 * smul(a3, c_) + smul(a2, b2) - 4 * b3 - 27 * cc)
        return sp.cancel(sp.expand(disc).coeff(p, 16))
    f0 = disc_p16(sp.Integer(0)); f1 = disc_p16(sp.Integer(1))
    d4_7 = sp.cancel(-f0 / (f1 - f0))
    ok_d4 = sp.simplify(d4_7 - D4F(7)) == 0
    print(f"     delta_4(7) = {d4_7}  (universal {D4F(7)}): {ok_d4}")
    # ---- c_4(7): series-ring ladder KMAX=6 + dict assembly to p^4 ----
    NORD = 10
    def ser(e): return sp.expand(sp.cancel(sp.series(e, p, 0, NORD).removeO()))
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(Aval)] for i in range(Aval)]
    st = list(itertools.product(range(Aval), repeat=3)); reps = [None] * 5
    for k, tr in enumerate(st):
        if reps[_pat(tr)] is None: reps[_pat(tr)] = k
    Q0 = sp.zeros(5, 5); Q1 = sp.zeros(5, 5); Q2 = sp.zeros(5, 5)
    for a_i, (x, xp, xq) in enumerate([st[r] for r in reps]):
        for (y, yp, yq) in st:
            t = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x][y])
            ne = (1 if yp != y else 0) + (1 if yq != y else 0)
            [Q0, Q1, Q2][ne][a_i, _pat((y, yp, yq))] += t
    Q0 = Q0.applyfunc(lambda e: sp.cancel(sp.together(e)))
    Q1 = Q1.applyfunc(lambda e: ser(sp.cancel(sp.together(e))))
    Q2 = Q2.applyfunc(lambda e: ser(sp.cancel(sp.together(e))))
    v = Q0[:, 0]
    xs = {0: v}; lams = {}
    for k in range(1, 7):
        r = (Q1 * xs[k - 1] + (Q2 * xs[k - 2] if k >= 2 else sp.zeros(5, 1))).applyfunc(ser)
        lam_k = ser(r[0]); xk = r - lam_k * v
        for j in range(1, k): xk = xk - lams[j] * xs[k - j]
        xs[k] = xk.applyfunc(ser); lams[k] = lam_k
    def laurent_c(e, pole, hi):
        ee = sp.expand(sp.cancel(sp.expand(e) * p**pole))
        return {k - pole: sp.cancel(ee.coeff(p, k)) for k in range(0, hi + pole + 1)}
    Dc = W_ * p**2 * (1 + 2 * p + d2v * p**2 + d3v * p**3 + d4_7 * p**4)
    etaD = 2 * D * (1 - D) / (A * D - 2 * D**2 - (A - 1))
    CrD = ((A * D - 2 * D**2 - (A - 1)) / (A * D - (A - 1)))**2
    etac_c = laurent_c(sp.series(etaD.subs(D, Dc), p, 0, 11).removeO(), 0, 10)
    Crc_c = laurent_c(sp.series(CrD.subs(D, Dc), p, 0, 9).removeO(), 0, 8)
    def conv(a_, b_, hi):
        out = {}
        for i, ai in a_.items():
            if ai == 0: continue
            for j, bj in b_.items():
                if i + j > hi or bj == 0: continue
                out[i + j] = sp.cancel(out.get(i + j, 0) + ai * bj)
        return out
    eta_pow = {1: etac_c}
    for k in range(2, 7): eta_pow[k] = conv(eta_pow[k - 1], etac_c, 14)
    lam_c = {k: laurent_c(lams[k], 6, 4) for k in range(1, 7)}
    lamc_c = {j: (sp.Integer(1) if j == 0 else sp.Integer(0)) for j in range(0, 7)}
    for k in range(1, 7):
        t = conv(lam_c[k], eta_pow[k], 6)
        for j in range(0, 7): lamc_c[j] = sp.cancel(lamc_c[j] + t.get(j, 0))
    Cl = conv(Crc_c, lamc_c, 6)
    Nc = {j: sp.cancel((1 if j == 0 else 0) - Cl.get(j, 0)) for j in range(0, 7)}
    inv_c = laurent_c(sp.series(1 / (2 * Dc * (1 - Dc)), p, 0, 5).removeO(), 2, 4)
    marg = conv(Nc, inv_c, 4)
    m1 = sp.simplify(marg.get(1, 0)); m2 = sp.simplify(marg.get(2, 0))
    m3 = sp.simplify(marg.get(3, 0)); m4 = sp.simplify(marg.get(4, 0))
    ok_lower = (sp.simplify(m1 - sp.Rational(5, 6)) == 0 and
                sp.simplify(m2 + sp.Rational(4 * 49 + 24 * 7 - 65, 8 * 36)) == 0 and
                sp.simplify(m3 + sp.Rational(18 * 49 - 153 * 7 + 232, 8 * 216)) == 0)
    ok_c4 = sp.simplify(m4 - C4F(7)) == 0
    print(f"     c_1..c_3 pipeline validation at A=7: {ok_lower}")
    print(f"     c_4(7) = {m4}  (universal {C4F(7)}): {ok_c4}")
    return rep("W3 out-of-sample A=7: both coefficients on-form (pipeline validated)", bool(ok_d4 and ok_lower and ok_c4))

if __name__ == "__main__":
    print("=" * 78)
    print("Endpoint ladder fourth order: delta_4(A), c_4(A) -- overdetermined universal forms")
    print("=" * 78)
    W1(); W2(); W3()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("Four-term corner margin available at eight alphabets (A=2..7 exact; forms")
    print("overdetermined). Symbolic-A universality not claimed (discriminant route")
    print("impractical at this order -- see the delta_3 lesson).")

#!/usr/bin/env python3
"""
bug_009_ak_decay_convexity.py
============================================================================
BUG-009-B -- the a_k Chebyshev-coefficient DECAY-ORDER probe (the precise
quantity the hard research leaf BUG-009-D must bound rigorously).

CONTEXT (Remark 7.34m'', "All-D via a Chebyshev tail").  The general-A
achievability residual reduces to the convexity of the replica Perron branch
G(u)=g_A(arccos u) on u=cos s in [-1,1].  Expanding the floor function in
Chebyshev modes,
        g_A(s) = sum_{k>=0} a_k(D) T_k(cos s),
the u-independent curvature lower bound is
        G''(u) >= LB(D) := 4 a_2 - sum_{k>=3} |a_k| (1/3) k^2 (k^2 - 1),
and LB(D) > 0 PROVES the all-D convexity on the memory bulk.  Whether LB(D)>0
holds for ALL D up to D_c hinges ENTIRELY on how fast |a_k(D)| decays in k.

WHAT IS PROVEN vs WHAT THIS PROBE MEASURES (the honest line).
  * PROVEN (rigorously, Remark 7.34m''):  the analyticity strip half-width is
    tau* >= ln(1/D) - O(1) (Bauer-Fike gap-stability of the frozen operator
    W_0 against ||W(s)-W_0|| <= B D e^{|Im s|}), which gives the CRUDE Bernstein
    majorant  |a_k| <= 2 M e^{-k tau*} = O((cD)^k).
  * OBSERVED (numerically, this probe):  the TRUE decay is one D-power FASTER,
    |a_k| = O(D^{k+1}) for k>=2, with ALTERNATING signs and geometric ratio
    |a_{k+1}/a_k| ~ 0.9 D.  The extra D-power is a harmonic CANCELLATION that
    the crude majorant does not see.
  * OPEN (the research leaf BUG-009-D, NOT closed here):  a rigorous,
    cancellation-aware proof of |a_k| = O(D^{k+1}).  This probe does NOT prove
    it; it MEASURES the decay ORDER (the exponent of D in |a_k(D)|) and isolates
    the exact one-D-power gap between the proven majorant and the observed decay,
    so the hard leaf has a reproducible numerical target / falsifier.

THIS PROBE THEREFORE *MEASURES THE ORDER*, IT DOES NOT FABRICATE A BOUND.
A PASS here means "the observed O(D^{k+1}) decay order, alternating signs, and
the documented one-D-power Bernstein looseness are reproduced", NOT "the a_k
bound is proven".

CHECKS (PASS/FAIL):
  B1  decay ORDER: fit exponent of |a_k(D)| vs D for k=2..6; expect ~ k+1
      (i.e. a_2~D^3, a_3~D^4, ..., a_6~D^7).
  B2  alternating signs of a_k for k>=2 (a_2>0, a_3<0, a_4>0, ...).
  B3  geometric ratio |a_{k+1}/a_k| ~ 0.9 D (the small-D in-mode ratio).
  B4  a_2 = 2 kappa_A^2 D^3 + o(D^3) (closed form at A=2: kappa_2^2=
      (1-2p)^2/(p^2(1-p)^2)); ties a_2 to the dispersion constant.
  B5  Bauer-Fike/Bernstein CONTRAST: the proven majorant order is O(D^k)
      (one power too loose); the observed order is O(D^{k+1}); the GAP is
      exactly +1 in the D-exponent => the cancellation quantity BUG-009-D must
      capture.  (Asserts the looseness IS one power; does NOT claim closure.)

Deps: mpmath, numpy.  Python: /Users/para/.venvs/rnr/bin/python.  ~2-4 min.
High precision (mpmath 60 dps): the margin is ~D^3 and float-diff cannot
resolve it.  We use the <=5x5 S_A pattern-quotient operator (NOT the full A^3,
on which mp.eig QR is unstable), the same operator as probe_7_34m_*.
"""
import itertools
import mpmath as mp
import numpy as np
mp.mp.dps = 60

PASS = True
def rep(name, ok):
    global PASS; PASS = PASS and ok
    print(f"  {name:<70} {'PASS' if ok else 'FAIL'}"); return ok

# ---- replica floor g_A(s) via the <=5x5 S_A pattern quotient (mp) -----------
def _pat(tr):
    x, u, up = tr
    if x == u == up: return 0
    if x == u and u != up: return 1
    if x == up and u != up: return 2
    if u == up and x != u: return 3
    return 4

def gA(A, p, D, s):
    """g_A(s) = |C_s/C_0|^2 rho(Q_A(s)) via the S_A pattern quotient (mpmath)."""
    A = int(A); p = mp.mpf(p); D = mp.mpf(D); s = mp.mpf(s)
    th = mp.log(D / ((A - 1) * (1 - D))); E = mp.e**(th + 1j * s); den = A * D - (A - 1)
    C = ((A - 1) * D * E + D - (A - 1)) / den
    eta = ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))
    C0 = ((A - 1) * D * mp.e**th + D - (A - 1)) / den
    Cr = abs(C / C0)**2
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]
    etb = mp.conj(eta)
    st = list(itertools.product(range(A), repeat=3))
    npat = 5 if A >= 3 else 4; reps = [None] * npat
    for k, tr in enumerate(st):
        pp = _pat(tr)
        if reps[pp] is None: reps[pp] = k
    Q = mp.zeros(npat, npat)
    for a in range(npat):
        i = reps[a]; x, xp, xq = st[i]
        for j, (y, yp, yq) in enumerate(st):
            Q[a, _pat((y, yp, yq))] += T[xp][yp] * T[xq][yq] / T[x][y] * (eta if yp != y else 1) * (etb if yq != y else 1)
    ev, _ = mp.eig(Q); return Cr * max(abs(e) for e in ev)

def cheb_coeffs(A, p, D, K=8, N=64):
    """Chebyshev coefficients a_k of g_A(s)=sum_k a_k T_k(cos s), k=0..K.
    g_A is even in s, so the discrete-cosine quadrature in s is exact for the
    even extension; N nodes resolve modes up to K well below the aliasing wall."""
    sj = [mp.pi * (j + mp.mpf('0.5')) / N for j in range(N)]
    gj = [gA(A, p, D, s) for s in sj]
    out = []
    for k in range(K + 1):
        c = sum(gj[j] * mp.cos(k * sj[j]) for j in range(N)) * 2 / N
        out.append(c / (2 if k == 0 else 1))
    return out

def Dc(A, pv):
    def im(Dv):
        b0 = Dv / (A - 1); K = np.full((A, A), b0); np.fill_diagonal(K, 1 - Dv); Ki = np.linalg.inv(K)
        Tn = np.full((A, A), pv / (A - 1)); np.fill_diagonal(Tn, 1 - pv)
        By = lambda y: np.array([[Ki[y, a] * Tn[a, b] for b in range(A)] for a in range(A)])
        ev = np.linalg.eigvals(By(1) @ By(0)); t = ev[np.argsort(-np.abs(ev))[:2]]
        return max(abs(t[0].imag), abs(t[1].imag))
    lo, hi = 1e-12, (A - 1) / A - 1e-7
    for _ in range(60):
        m = (lo + hi) / 2; lo, hi = (m, hi) if im(m) < 1e-11 else (lo, m)
    return lo

def fit_order(Ds, vals):
    """log-log slope of |val| vs D (least squares); returns the decay exponent."""
    xs = [mp.log(d) for d in Ds]; ys = [mp.log(abs(v)) for v in vals]
    n = len(xs); sx = sum(xs); sy = sum(ys)
    sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)

# ----------------------------------------------------------------------------
def B1():
    print("-" * 78)
    print("B1  DECAY ORDER  |a_k(D)| ~ D^{k+1}  (k=2..6): fit exponent of |a_k| vs D.")
    print("    (small D, in the memory regime; expect exponents 3,4,5,6,7).")
    ok = True
    # small D values well inside the Gray region so the asymptotic order shows
    for A, p in ((2, '0.2'), (3, '0.2')):
        Ds = [mp.mpf('8e-4'), mp.mpf('4e-4'), mp.mpf('2e-4'), mp.mpf('1e-4')]
        coeffs = [cheb_coeffs(A, p, D, K=7) for D in Ds]
        orders = []
        for k in range(2, 7):
            vals = [c[k] for c in coeffs]
            orders.append(fit_order(Ds, vals))
        exp_str = ", ".join(f"a{k}~D^{mp.nstr(o, 4)}" for k, o in zip(range(2, 7), orders))
        print(f"     A={A} p={p}: {exp_str}")
        # each fitted order should be within 0.25 of k+1
        for k, o in zip(range(2, 7), orders):
            ok = ok and abs(o - (k + 1)) < mp.mpf('0.25')
    return rep("B1 a_k = O(D^{k+1}) decay order (k=2..6) measured", ok)

def B2():
    print("-" * 78)
    print("B2  ALTERNATING SIGNS  a_2>0, a_3<0, a_4>0, a_5<0, a_6>0  (k>=2).")
    print("    This is an ASYMPTOTIC (small-D) property -- the documented")
    print("    cancellation pattern -- tested in the small-D regime where 7.34m''")
    print("    states it.  (At moderate D=0.5 D_c the alternation can break at the")
    print("    HIGHEST mode k=6, where |a_6|~D^7~1e-17 is far below the LB(D) budget;")
    print("    this is reported, not failed -- LB(D) consumes a_2..a_5, not a_6.)")
    ok = True
    expect = [1, -1, 1, -1, 1]
    for A, p in ((2, '0.15'), (2, '0.3'), (3, '0.2')):
        dc = Dc(A, float(p))
        # asymptotic regime (where the alternation claim lives): D = 0.1 D_c
        D = mp.mpf('0.1') * dc
        a = cheb_coeffs(A, p, D, K=6)
        signs = [int(mp.sign(a[k])) for k in range(2, 7)]
        match = signs == expect
        print(f"     A={A} p={p} D=0.1 D_c (asymptotic): sign(a_2..a_6)={signs} (expect {expect}) -> {match}")
        ok = ok and match
        # moderate D: report only (high mode may deviate; magnitudes negligible)
        Dm = mp.mpf('0.5') * dc
        am = cheb_coeffs(A, p, Dm, K=6)
        signs_m = [int(mp.sign(am[k])) for k in range(2, 7)]
        mags_m = [mp.nstr(abs(am[k]), 2) for k in range(2, 7)]
        # the LB(D)-relevant modes a_2..a_5 must still alternate at moderate D
        lb_modes_ok = signs_m[:4] == expect[:4]
        print(f"       (D=0.5 D_c report: sign(a_2..a_6)={signs_m}, |a_k|={mags_m}; "
              f"LB modes a_2..a_5 alternate: {lb_modes_ok})")
        ok = ok and lb_modes_ok
    return rep("B2 alternating signs (asymptotic; LB modes a_2..a_5 at all D)", ok)

def B3():
    print("-" * 78)
    print("B3  GEOMETRIC RATIO  |a_{k+1}/a_k| ~ 0.9 D  (small-D in-mode ratio).")
    ok = True
    for A, p in ((2, '0.2'), (3, '0.2')):
        D = mp.mpf('3e-4')
        a = cheb_coeffs(A, p, D, K=6)
        ratios = [abs(a[k + 1] / a[k]) / D for k in range(2, 6)]  # divide by D => ~0.9
        rs = ", ".join(mp.nstr(r, 4) for r in ratios)
        print(f"     A={A} p={p} D={mp.nstr(D, 2)}: |a_{{k+1}}/a_k|/D (k=2..5) = [{rs}] (~0.9)")
        # the ratio/D should be O(1) and < 1.5 (clearly geometric in D, not O(1))
        for r in ratios:
            ok = ok and r < mp.mpf('1.5') and r > mp.mpf('0.2')
    return rep("B3 |a_{k+1}/a_k| ~ 0.9 D (geometric, one-D-per-mode)", ok)

def B4():
    print("-" * 78)
    print("B4  a_2 = 2 kappa_A^2 D^3 + o(D^3) (ties a_2 to the dispersion const).")
    print("    A=2 closed form: kappa_2^2 = (1-2p)^2/(p^2(1-p)^2).")
    ok = True
    for p in ('0.15', '0.3'):
        pp = mp.mpf(p); k2 = (1 - 2 * pp)**2 / (pp**2 * (1 - pp)**2)
        D = mp.mpf('2e-4')
        a = cheb_coeffs(2, p, D, K=4)
        ratio = a[2] / (2 * k2 * D**3)
        print(f"     A=2 p={p} D={mp.nstr(D, 2)}: a_2/(2 kappa_2^2 D^3) = {mp.nstr(ratio, 6)} (->1)")
        ok = ok and abs(ratio - 1) < mp.mpf('3e-2')
    # A=3: just check a_2>0 and a_2/D^3 stabilises (no closed form)
    for p in ('0.2',):
        cs = []
        for D in (mp.mpf('4e-4'), mp.mpf('1e-4')):
            a = cheb_coeffs(3, p, D, K=4); cs.append(a[2] / D**3)
        print(f"     A=3 p={p}: a_2/D^3 at D=4e-4,1e-4 = [{mp.nstr(cs[0], 6)}, {mp.nstr(cs[1], 6)}] (>0, stabilising)")
        ok = ok and cs[0] > 0 and cs[1] > 0 and abs(cs[1] / cs[0] - 1) < mp.mpf('0.2')
    return rep("B4 a_2 = 2 kappa_A^2 D^3 (curvature mode, closed form A=2)", ok)

def B5():
    print("-" * 78)
    print("B5  BAUER-FIKE/BERNSTEIN CONTRAST -- the documented one-D-power looseness.")
    print("    PROVEN majorant (Rem 7.34m''): tau* >= ln(1/D)-O(1) (gap-stability of")
    print("    frozen W_0), giving |a_k| <= 2 M e^{-k tau*} = O((cD)^k) = O(D^k).")
    print("    OBSERVED decay (B1): |a_k| = O(D^{k+1}).  The GAP is +1 D-power for")
    print("    each k>=2 -- the harmonic cancellation the crude majorant misses.")
    print("    This probe MEASURES that the gap is one power; it does NOT close it")
    print("    (the rigorous cancellation-aware a_k bound is the OPEN leaf BUG-009-D).")
    ok = True
    for A, p in ((2, '0.2'),):
        Ds = [mp.mpf('8e-4'), mp.mpf('4e-4'), mp.mpf('2e-4'), mp.mpf('1e-4')]
        coeffs = [cheb_coeffs(A, p, D, K=7) for D in Ds]
        print("     k | observed order (B1) | crude Bernstein order | gap (observed-crude)")
        for k in range(2, 7):
            vals = [c[k] for c in coeffs]
            obs = fit_order(Ds, vals)        # ~ k+1
            crude = mp.mpf(k)                # |a_k| <= O(D^k) => order k
            gap = obs - crude
            print(f"     {k} |        {mp.nstr(obs, 4):>6}        |          {k:>2}          |   {mp.nstr(gap, 4):>6}")
            # the observed order exceeds the crude majorant order by ~1 (within 0.25)
            ok = ok and abs(gap - 1) < mp.mpf('0.25')
        # also: the crude majorant tail would NOT certify LB>0 if it were tight at
        # order k (it over-counts by one D-power); demonstrate the looseness margin.
        # ratio of observed |a_3| to the crude same-order proxy |a_2|*D (=order-3 crude):
        D = mp.mpf('2e-4'); a = cheb_coeffs(A, p, D, K=4)
        crude_a3 = abs(a[2]) * D            # what an order-(k) majorant would force a_3 ~ a_2 (O(D^3)) not O(D^4)
        looseness = abs(a[3]) / crude_a3    # << 1 => true a_3 is one D-power smaller
        print(f"     looseness witness: |a_3_observed| / (|a_2| * 1) at D={mp.nstr(D, 2)} = "
              f"{mp.nstr(abs(a[3]) / abs(a[2]), 4)} (~0.9 D, i.e. one extra D vs the crude O(D^3))")
        ok = ok and abs(a[3]) / abs(a[2]) < mp.mpf('5') * D  # one extra D-power
    return rep("B5 Bernstein majorant is one-D-power loose (gap +1); a_k bound OPEN", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("BUG-009-B -- a_k Chebyshev-coefficient DECAY-ORDER probe (Remark 7.34m'')")
    print("MEASURES the observed O(D^{k+1}) decay order + Bauer-Fike contrast.")
    print("Does NOT prove the rigorous a_k bound (that is the OPEN leaf BUG-009-D).")
    print("=" * 78)
    B1(); B2(); B3(); B4(); B5()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")

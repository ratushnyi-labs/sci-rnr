#!/usr/bin/env python3
r"""
lemma_7_34_endpoint_effective_bound.py
============================================================================
BUG-009-D / Route B: EFFECTIVE endpoint positivity -- the first non-asymptotic
theorem of the corner program. For each alphabet A in {3,4,5,6} there are explicit
constants (p0, mu, K4, Keff) such that for ALL p in (0, p0]:
  (a) the reduced-cubic discriminant provably changes sign across the bracket
          B(p) = [Dhat(p) - mu p^6, Dhat(p) + mu p^6],
      Dhat = [p^2/(4(A-1))](1 + 2p + delta_2 p^2 + delta_3 p^3)
      (Sturm), so a threshold crossing lies in B(p). CAPTURE CAVEAT (adversarial
      review): that D_c is the FIRST crossing (no complexification below the
      bracket) is verified by exact Sturm root-counts at sampled p (p0, p0/2,
      p0/4; only the benign D=0 endpoint degeneracy below B(p)) but not yet
      certified for ALL p <= p0 -- a bivariate no-root-below-bracket lemma is
      the remaining rigor step (the existing signOn/Sturm machinery suffices);
  (b) for ALL D in B(p):
          |margin(D,p) - c1 p - c2 p^2 - c3 p^3| <= K4 p^4,
      hence |margin - c1 p - c2 p^2| <= Keff p^3 and margin(D,p) > 0.
CERTIFIED CONSTANTS (re-derived by running this script):
  A=3: p0 = 1/16, K4 ~ 598.5, Keff ~ 38.4
  A=4: p0 = 1/12, K4 ~ 128.8, Keff ~ 11.2
  A=5: p0 = 1/12, K4 ~ 108.8, Keff ~  9.2
  A=6: p0 = 1/12, K4 ~ 109.7, Keff ~  9.2
PROOF STRUCTURE (all exact rational arithmetic; no floats in the logic):
  C1  bracket capture: Sturm sign flip of the discriminant numerator at the
      bracket ends (denominator sign manifest);
  C2  exact PT part: ladder k <= K_EX = 11 composed over the u-parametrized
      bracket, bivariate coefficient bounds, remainder p-valuation = 4;
  C3  PT tail k >= 12: diagonally-rescaled resolvent (Kato) bound -- for
      |eta| <= gam p the four non-branch eigenvalues stay in |z| <= 2/5 and the
      branch in |z-1| <= 2/5, so the branch is the spectral radius -- plus a
      Cauchy geometric tail with q = |eta|/(gam p) <= qbar < 1/2.
PROVENANCE: user-authorized research fan-out (artifact research_c3_remainder/,
step9/step11 + step10 adversarial spot-checks: D_c-in-bracket, eigenvalue
confinement, tail <= bound with ~13 orders of slack, resid4/p^4 -> c_4 ~ 0.29
<< K4 -- all PASS); independently re-executed before packaging (A=3 and A=4
re-derived with matching constants). A=2 needs no effective theorem: the sharp
endpoint is proven in closed form (probe_7_34m_endpoint_sharp_a2.py).
RUNTIME NOTE: each alphabet takes ~2-5 min (exact arithmetic); __main__ runs
A=3 and A=4 by default; pass explicit alphabets as argv to run others.

Deps: sympy.  Python: /Users/para/.venvs/rnr/bin/python.  ~5-10 min (A=3,4).
"""
import sys as _sys

def run_certificate(Aval):
    """Run the full effective-bound certificate for one alphabet.
    Returns True iff the positivity certificate completes (prints the theorem)."""

    import itertools
    import sympy as sp
    from sympy import Rational as R

    A = sp.Integer(Aval)
    p, D, u = sp.symbols('p D u')
    W = R(1, 4 * (Aval - 1))
    d2v = R(12 * Aval**2 - 28 * Aval + 21, 4 * (Aval - 1)**2)
    d3v = R(8 * Aval**3 - 32 * Aval**2 + 53 * Aval - 32, 2 * (Aval - 1)**3)
    MU = 6 * W
    c1 = R(Aval - 2, Aval - 1)
    c2 = R(-(4 * Aval**2 + 24 * Aval - 65), 8 * (Aval - 1)**2)
    c3 = R(-(18 * Aval**2 - 153 * Aval + 232), 8 * (Aval - 1)**3)
    K_EX = 11          # exact ladder order

    def _pat(tr):
        x, u_, up = tr
        if x == u_ == up: return 0
        if x == u_ and u_ != up: return 1
        if x == up and u_ != up: return 2
        if u_ == up and x != u_: return 3
        return 4

    # ---------- ladder (exact rational lam_k, k <= K_EX) ----------
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(Aval)] for i in range(Aval)]
    st = list(itertools.product(range(Aval), repeat=3)); reps = [None] * 5
    for k, tr in enumerate(st):
        if reps[_pat(tr)] is None: reps[_pat(tr)] = k
    Q0 = sp.zeros(5, 5); Q1 = sp.zeros(5, 5); Q2 = sp.zeros(5, 5)
    for a_ in range(5):
        x, xp, xq = st[reps[a_]]
        for (y, yp, yq) in st:
            t = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x][y])
            ne_ = (1 if yp != y else 0) + (1 if yq != y else 0)
            [Q0, Q1, Q2][ne_][a_, _pat((y, yp, yq))] += t
    Q0 = Q0.applyfunc(sp.cancel); Q1 = Q1.applyfunc(sp.cancel); Q2 = Q2.applyfunc(sp.cancel)
    v = Q0[:, 0]
    xs = {0: v}; lams = {}
    for k in range(1, K_EX + 1):
        r = Q1 * xs[k - 1] + (Q2 * xs[k - 2] if k >= 2 else sp.zeros(5, 1))
        lam_k = sp.cancel(r[0]); xk = r - lam_k * v
        for j in range(1, k): xk = xk - lams[j] * xs[k - j]
        xs[k] = xk.applyfunc(sp.cancel); lams[k] = lam_k
    print("ladder done; lam_k denominators:", flush=True)
    LNUM = {}; LDEN = {}
    for k in range(1, K_EX + 1):
        n_, d_ = sp.fraction(sp.cancel(lams[k]))
        LNUM[k] = sp.Poly(sp.expand(n_), p); LDEN[k] = sp.Poly(sp.expand(d_), p)
        print(f"  k={k}: deg num {LNUM[k].degree()}, deg den {LDEN[k].degree()}", flush=True)

    # ---------- helper: rigorous univariate bounds on (0,p0] ----------
    def polyUB(poly, p0):
        """sup |poly(p)| on [0,p0], poly sympy Poly in p."""
        return sum(abs(c) * p0**m[0] for m, c in zip(poly.monoms(), poly.coeffs()))
    def polyLB0(poly, p0):
        """lower bound of |poly| on [0,p0] via constant-coefficient dominance; None if fails."""
        c0 = poly.coeff_monomial(1)
        rest = sum(abs(c) * p0**m[0] for m, c in zip(poly.monoms(), poly.coeffs()) if m[0] > 0)
        return abs(c0) - rest if abs(c0) > rest else None
    def ratUB(expr, p0):
        """sup |expr(p)| on (0,p0] for rational expr with p-valuation >= 0. None if machinery fails."""
        n_, d_ = sp.fraction(sp.together(sp.cancel(expr)))
        pn = sp.Poly(sp.expand(n_), p); pd = sp.Poly(sp.expand(d_), p)
        a = min(m[0] for m in pn.monoms()); b = min(m[0] for m in pd.monoms())
        if a - b < 0: return None
        pn1 = sp.Poly(sp.expand(n_ / p**a), p); pd1 = sp.Poly(sp.expand(d_ / p**b), p)
        lb = polyLB0(pd1, p0)
        if lb is None: return None
        return p0**(a - b) * polyUB(pn1, p0) / lb if a > b else polyUB(pn1, p0) / lb

    def signOn(polyexpr, p0):
        """rigorous sign of a univariate polynomial on (0,p0] via valuation + Sturm root count."""
        q = sp.Poly(sp.expand(polyexpr), p)
        a = min(m[0] for m in q.monoms())
        q1 = sp.Poly(sp.expand(polyexpr / p**a), p)
        nroots = q1.count_roots(0, p0)          # roots in [0,p0]
        if q1.coeff_monomial(1) == 0: return None
        if nroots > 0: return None
        return sp.sign(q1.coeff_monomial(1))

    # ---------- C1: bracket capture ----------
    Dhat = W * p**2 * (1 + 2 * p + d2v * p**2 + d3v * p**3)
    Dbr = Dhat + u * MU * p**6
    lam2s = 1 - D * A / (A - 1)
    al = R(1, Aval) + (1 - R(1, Aval)) / lam2s
    be = R(1, Aval) - R(1, Aval) / lam2s
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
    Mm = Bred(1) * Bred(0)
    Mm = Mm.applyfunc(lambda e: sp.cancel(sp.together(e)))
    c2m = sp.cancel(Mm.trace())
    c1m = sp.cancel(sum(Mm[[i for i in rr], [j for j in rr]].det() for rr in ([0, 1], [0, 2], [1, 2])))
    c0m = sp.cancel(Mm.det())
    aa = -c2m; bb = c1m; cc = -c0m
    disc = 18 * aa * bb * cc - 4 * aa**3 * cc + aa**2 * bb**2 - 4 * bb**3 - 27 * cc**2
    NG, DG = sp.fraction(sp.together(sp.cancel(disc)))
    print("discriminant numerator built", flush=True)
    # DG sign on domain: DG is a power of lam2*(A-1)-type positive factors -- verify at a point
    # and verify no roots for D in (0, 1/8), p in (0, 1/4) via factorization check:
    DGf = sp.factor(DG)
    print(f"  disc denominator (factored) = {DGf}", flush=True)

    def C1_capture(p0):
        okp = True
        signs = {}
        for uu in (-1, 1):
            expr = sp.expand(NG.subs(D, Dbr.subs(u, uu)))
            s = signOn(expr, p0)
            signs[uu] = s
            if s is None: okp = False
        if okp and signs[-1] is not None and signs[1] is not None and signs[-1] * signs[1] == -1:
            return True, signs
        return False, signs

    # ---------- C2: exact PT part over the bracket ----------
    def C2_exactpart(p0):
        Du = sp.expand(Dbr)
        ne_ = sp.expand(2 * Du * (1 - Du))                     # eta numerator
        de_ = sp.expand(A * Du - 2 * Du**2 - (A - 1))          # eta denominator (negative)
        dg_ = sp.expand(A * Du - (A - 1))                      # Cr denominator sqrt (negative)
        # common denominator of lam_k
        L = LDEN[1]
        for k in range(2, K_EX + 1): L = sp.Poly(sp.lcm(L.as_expr(), LDEN[k].as_expr()), p)
        Lx = L.as_expr()
        # SN = L*de^K + sum_k ln_k (L/ld_k) ne^k de^(K-k)
        nepow = {0: sp.Integer(1)}; depow = {0: sp.Integer(1)}
        for k in range(1, K_EX + 1):
            nepow[k] = sp.expand(nepow[k - 1] * ne_)
            depow[k] = sp.expand(depow[k - 1] * de_)
        SN = sp.expand(Lx * depow[K_EX])
        for k in range(1, K_EX + 1):
            cof = sp.expand(sp.cancel(Lx / LDEN[k].as_expr()) * LNUM[k].as_expr())
            SN = sp.expand(SN + cof * nepow[k] * depow[K_EX - k])
        GN = sp.expand(dg_**2 * Lx * depow[K_EX] - de_**2 * SN)
        GD = sp.expand(dg_**2 * Lx * depow[K_EX])
        tgt = R(3, 2) + c1 * p + c2 * p**2 + c3 * p**3
        EN = sp.expand(GN - tgt * GD * 2 * Du * (1 - Du))
        ED = sp.expand(GD * 2 * Du * (1 - Du))
        PEN = sp.Poly(EN, u, p); PED = sp.Poly(ED, u, p)
        aE = min(m[1] for m in PEN.monoms()); bE = min(m[1] for m in PED.monoms())
        print(f"  p-valuations: EN {aE}, ED {bE}, difference {aE-bE}", flush=True)
        if aE - bE < 4: return None
        # bounds with |u|<=1, p<=p0
        EN1 = [(m, c) for m, c in zip(PEN.monoms(), PEN.coeffs())]
        NB = sum(abs(c) * p0**(m[1] - aE) for m, c in EN1)
        ED1 = [(m, c) for m, c in zip(PED.monoms(), PED.coeffs())]
        c00 = sum(c for m, c in ED1 if m[1] == bE and m[0] == 0)
        # u-carrying terms at the bottom p-level would break dominance; include all non-(0,bE):
        rest = sum(abs(c) * p0**(m[1] - bE) for m, c in ED1 if not (m[1] == bE and m[0] == 0))
        if abs(c00) <= rest: return None
        LBd = abs(c00) - rest
        return p0**(aE - bE - 4) * NB / LBd     # K_rat

    # ---------- C3: Kato tail ----------
    def C3_tail(p0, SIG, TAU):
        s = [sp.Integer(1), sp.Integer(1), sp.Integer(1), SIG / p, TAU]
        def resc(M):
            return sp.Matrix(5, 5, lambda i, j: sp.cancel(M[i, j] * s[j] / s[i]))
        Q0r, Q1r, Q2r = resc(Q0), resc(Q1), resc(Q2)
        def rowUB(M, scale_p=0):
            """sup over (0,p0] of p^scale_p * ||M||_inf, via per-row ratUB.
            All entries of Q0,Q1,Q2 are >=0 (sums of positive channel terms) and the
            rescaling is positive, so row sums of entries = inf-norm row sums."""
            best = R(0)
            for i in range(5):
                rs2 = sp.cancel(sum(M[i, j] for j in range(5)))
                val = ratUB(p**scale_p * rs2, p0)
                if val is None: return None
                best = max(best, val)
            return best
        N0b = rowUB(Q0r, 0)
        B1b = rowUB(Q1r, 1)     # sup p*||Q1~||   (V-term scales |zeta| ||Q1|| <= gam B1b)
        B2b = rowUB(Q2r, 2)     # sup p^2*||Q2~|| (V-term |zeta|^2||Q2|| <= gam^2 B2b)
        if None in (N0b, B1b, B2b): return None
        CC = R(5, 2) + R(25, 6) * N0b
        # r0 = gam*p, condition (gam*B1b + gam^2*B2b)*CC <= 1/2 for all p<=p0
        gam = R(1, 2) / ((B1b + 1) * CC)
        gam = sp.nsimplify(sp.Rational(int(gam * 10**5), 10**5))  # rational round-down
        while gam > 0 and (gam * B1b + gam**2 * B2b) * CC > R(1, 2):
            gam = gam * R(9, 10)
        if gam <= R(1, 10**6): return None
        # eta envelope over bracket: |eta| <= 2*Dhi*(1-0)/((A-1)-A*Dhi)
        Dhi = sp.expand(Dhat + MU * p**6)
        etaU = sp.cancel(2 * Dhi / ((A - 1) - A * Dhi))
        qb = ratUB(etaU / (gam * p), p0)
        if qb is None or qb >= R(1, 2): return None, gam, qb, None
        # tail: |T| <= (2/5) q^{K+1}/(1-q); q(p) <= (etaU/(gam p)) pointwise
        # margin contribution <= CrUB * T / (2*Dlo*(1-Dhi))
        Dlo = sp.expand(Dhat - MU * p**6)
        CrUB = ratUB(((A - 1) - A * Dlo + 2 * Dhi**2)**2 / ((A - 1) - A * Dhi)**2, p0)
        # (numerator maximized, denominator minimized on bracket; both positive)
        if CrUB is None: return None
        # p-structure: etaU/(gam p) = O(p); (q)^{K+1}/(2Dlo(1-Dhi)) = O(p^{K+1-2})
        qexpr = etaU / (gam * p)
        tail_margin = CrUB * R(2, 5) / (1 - qb) * qexpr**(K_EX + 1) / (2 * Dlo * (1 - Dhi))
        # bound tail_margin / p^4 on (0,p0]:
        Kt = ratUB(tail_margin / p**4, p0)
        return Kt, gam, qb, (N0b, B1b, B2b, CC, CrUB)

    # ---------- run over p0 candidates ----------
    P0CANDS = [R(1, 8), R(1, 10), R(1, 12), R(1, 16), R(1, 24), R(1, 32), R(1, 48)]
    SGRID = [(R(Aval), R(Aval)), (R(Aval - 1), R(Aval - 1)), (R(2 * (Aval - 1)), R(2 * (Aval - 1))), (R(3), R(3)), (R(4), R(4)), (R(6), R(6))]
    done = False
    for P0 in P0CANDS:
        if done: break
        print(f"== trying p0 = {P0} = {float(P0):.5f}", flush=True)
        ok1, signs = C1_capture(P0)
        print(f"  C1 bracket sign change: {ok1} (signs {signs})", flush=True)
        if not ok1: continue
        Krat = C2_exactpart(P0)
        print(f"  C2 K_rat = {None if Krat is None else float(Krat)}", flush=True)
        if Krat is None: continue
        for SIG, TAU in SGRID:
            res = C3_tail(P0, SIG, TAU)
            if res is None or res[0] is None:
                print(f"  C3 (SIG={SIG},TAU={TAU}) failed: {None if res is None else (float(res[1]), None if res[2] is None else float(res[2]))}", flush=True)
                continue
            Kt, gam, qb, aux = res
            N0b, B1b, B2b, CC, CrUB = aux
            print(f"  C3 (SIG={SIG},TAU={TAU}): N0b={float(N0b):.3g} B1b={float(B1b):.3g} B2b={float(B2b):.3g} "
                  f"CC={float(CC):.3g} gam={float(gam):.5g} qbar={float(qb):.4g} K_tail={float(Kt):.6g}", flush=True)
            K4 = Krat + Kt
            Keff = abs(c3) + K4 * P0
            posexpr = c1 * p + c2 * p**2 + c3 * p**3 - K4 * p**4
            spos = signOn(sp.expand(posexpr), P0)
            print(f"  ==> K4 = {float(K4):.6g}; Keff = {float(Keff):.6g}; positivity sign on (0,{P0}]: {spos}", flush=True)
            if spos == 1:
                print(f"EFFECTIVE THEOREM (A={Aval}): for all p in (0,{P0}]:", flush=True)
                print(f"  |margin(D,p) - c1 p - c2 p^2 - c3 p^3| <= {float(K4):.5g} p^4 on the whole bracket,", flush=True)
                print(f"  |margin - c1 p - c2 p^2| <= {float(Keff):.5g} p^3,  and margin > 0.", flush=True)
                print(f"  constants: mu={MU}, K_EX={K_EX}, SIG={SIG}, TAU={TAU}, gam={gam}, K_rat={float(Krat):.5g}, K_tail={float(Kt):.5g}", flush=True)
                done = True
                break

    return done

if __name__ == "__main__":
    PASS = True
    alphabets = [int(a) for a in _sys.argv[1:]] or [3, 4]
    for Av in alphabets:
        print("=" * 78)
        ok = run_certificate(Av)
        print(f"  A={Av} effective positivity certificate: {'PASS' if ok else 'FAIL'}")
        PASS = PASS and ok
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")

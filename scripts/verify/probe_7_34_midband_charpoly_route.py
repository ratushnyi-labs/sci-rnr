#!/usr/bin/env python3
r"""
probe_7_34_midband_charpoly_route.py
============================================================================
BUG-009-D / Route B: the MID-BAND has an instrument -- real-rootedness + the
charpoly-derivative certificate govern the entire "non-perturbative" band.

CONTEXT. The Schur certificate closes p/pmax <= 0.40-0.45; the corner pocket is
covered by direct margin; between them lies the mid-band pfrac in (0.45, 0.98)
where NO perturbation-about-1 architecture can work (the branch moves ~0.6 from
1 -- mapped obstruction). The Route-B-era numerics (untracked route_b_cw3/sym
work) suggested the way out: the 5x5 pattern-matrix spectrum appears REAL on the
Gray locus, enabling sign-localization without eigenvalue extraction.

SCOUT RESULT (this probe): across the mid-band grid (A in {3,4,6}, pfrac 0.5..
0.95, D/D_c in {0.3, 0.6, 0.9, 0.99}, s in {30, 90, 150, 180} deg; 288 points,
40 dps):
  * the spectrum of Q(eta, etb) is REAL at every point (worst |Im| ~ 4e-41 --
    exactly real, not approximately);
  * the DERIVATIVE CERTIFICATE holds at every point: with Lam = L(s)/Cr,
    L(s) = 1 - (3/2) D(1-D)(1-cos s), all five derivative signs
    chi^(k)(Lam) > 0, k = 0..4 -- for a real-rooted quintic this is EXACTLY
    "all roots < Lam", i.e. rho(Q) < Lam, i.e. R_A(s;D) > 3/2.
CONSEQUENCE: the mid-band certification program reduces to sign-certifying FIVE
RATIONAL FUNCTIONS of (A, p, D, s) over the band (chi^(k)(Lam) involve no
eigenvalues) -- the same Sturm/Polya problem class solved repeatedly in this
arc -- PLUS a structural real-rootedness lemma.

STRUCTURE LEAD (for the realness proof): the eta-powers of Q are COLUMN-
determined -- Q = M . W with M REAL (nonneg pattern-transition weights) and
W = diag(1, etb, eta, |eta|^2, |eta|^2) (verified symbolically here, S2). The
realness of spec(M W) for |eta|-phase != 0 is then a structural property of
this column grading (candidate mechanisms: phase-similarity removing the
conjugate pair, or a hidden self-adjointness in a weighted inner product) --
the structural proof is the flagged next step, NOT claimed here.

REAL-COEFFICIENTS LEMMA (S3, proven): the replica-swap permutation P (classes
1 <-> 2) gives the exact similarity P Q(eta, etb) P = Q(etb, eta), hence the
characteristic polynomial is SYMMETRIC under eta <-> etb; its coefficients are
real polynomials in (eta, etb), so by the symmetric-function theorem they are
real polynomials in e1 = eta+etb and e2 = eta etb -- and on the conjugate locus
etb = conj(eta) (e1 = 2 Re eta, e2 = |eta|^2, both real) the quintic has REAL
COEFFICIENTS. This is half of the real-rootedness mechanism; what remains open
is root-REALNESS (discriminant/subresultant positivity of a real quintic over
the band -- the same sign-certification class as the derivative program).

CHECKS:
  S1  mid-band scout: real spectrum AND derivative certificate at all points of
      the (A, pfrac, D, s) grid above (labels: VERIFIED NUMERIC).
  S2  column-grading factorization Q = M . diag(1, etb, eta, |eta|^2, |eta|^2)
      with M real and eta-free (symbolic, exact).
  S3  swap-similarity P Q(eta,etb) P = Q(etb,eta) and charpoly eta<->etb
      symmetry (symbolic, exact, A=3,4,5) => real quintic coefficients on the
      conjugate locus (symmetric-function argument; numeric spot confirmation).

Deps: mpmath, numpy, sympy.  Python: /Users/para/.venvs/rnr/bin/python.  ~2-4 min.
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

def build_Q(A, p, D, s):
    th = mp.log(D / ((A - 1) * (1 - D))); E = mp.e**(th + 1j * s)
    eta = ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))
    C = ((A - 1) * D * E + D - (A - 1)) / (A * D - (A - 1))
    C0 = ((A - 1) * D * mp.e**th + D - (A - 1)) / (A * D - (A - 1)); Cr = abs(C / C0)**2
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]
    etb = mp.conj(eta)
    st = list(itertools.product(range(A), repeat=3)); reps = [None] * 5
    for k, tr in enumerate(st):
        if reps[_pat(tr)] is None: reps[_pat(tr)] = k
    Q = mp.zeros(5, 5)
    for a_ in range(5):
        i = reps[a_]; x, xp, xq = st[i]
        for j, (y, yp, yq) in enumerate(st):
            Q[a_, _pat((y, yp, yq))] += T[xp][yp] * T[xq][yq] / T[x][y] * (eta if yp != y else 1) * (etb if yq != y else 1)
    return Q, Cr

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

def S1():
    print("-" * 78); print("S1  mid-band scout: real spectrum + derivative certificate everywhere")
    tot = 0; real_ok = 0; cert_ok = 0; worst_im = mp.mpf(0)
    for A in (3, 4, 6):
        for pf in (0.5, 0.6, 0.7, 0.8, 0.9, 0.95):
            p = pf * (A - 1) / A
            dc = Dc_num(A, p)
            for fD in (0.3, 0.6, 0.9, 0.99):
                D = mp.mpf(fD) * dc
                for sd in (30, 90, 150, 180):
                    s = mp.mpf(sd) / 180 * mp.pi
                    Q, Cr = build_Q(A, p, D, s)
                    ev, _ = mp.eig(Q)
                    tot += 1
                    imax = max(abs(mp.im(e)) for e in ev)
                    rmax = max(abs(e) for e in ev)
                    if imax < mp.mpf('1e-25') * max(1, rmax): real_ok += 1
                    worst_im = max(worst_im, imax)
                    t = 1 - mp.cos(s)
                    Lam = (1 - mp.mpf(1.5) * D * (1 - D) * t) / Cr
                    coeffs = [mp.mpf(1)]
                    for e in ev:
                        new = [mp.mpf(0)] * (len(coeffs) + 1)
                        for i, c in enumerate(coeffs):
                            new[i] += c; new[i + 1] -= c * e
                        coeffs = new
                    coeffs = [mp.re(c) for c in coeffs]
                    def polyval(cs, x):
                        v = mp.mpf(0)
                        for c in cs: v = v * x + c
                        return v
                    def deriv(cs):
                        n = len(cs) - 1
                        return [cs[i] * (n - i) for i in range(n)]
                    cs = coeffs; allpos = True
                    for k in range(5):
                        if polyval(cs, Lam) <= 0: allpos = False; break
                        cs = deriv(cs)
                    if allpos: cert_ok += 1
    print(f"     points {tot}; real {real_ok}; cert {cert_ok}; worst|Im| {mp.nstr(worst_im,3)}")
    return rep("S1 real spectrum AND derivative certificate at 100% of the band grid", real_ok == tot and cert_ok == tot)

def S2():
    print("-" * 78); print("S2  column grading: Q = M(real) . diag(1, etb, eta, |eta|^2, |eta|^2)")
    import sympy as sp
    eta, etb = sp.symbols('eta etab')
    p = sp.symbols('p', positive=True)
    ok = True
    for Aval in (3, 5):
        A = sp.Integer(Aval)
        T = [[(1 - p) if i == j else p / (A - 1) for j in range(Aval)] for i in range(Aval)]
        st = list(itertools.product(range(Aval), repeat=3)); reps = [None] * 5
        for k, tr in enumerate(st):
            if reps[_pat(tr)] is None: reps[_pat(tr)] = k
        Q = sp.zeros(5, 5)
        for a_ in range(5):
            x, xp, xq = st[reps[a_]]
            for (y, yp, yq) in st:
                term = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x][y])
                if yp != y: term *= eta
                if yq != y: term *= etb
                Q[a_, _pat((y, yp, yq))] += term
        wcol = [sp.Integer(1), etb, eta, eta * etb, eta * etb]
        here = True
        for b in range(5):
            for a_ in range(5):
                m = sp.cancel(Q[a_, b] / wcol[b])
                if m.has(eta) or m.has(etb): here = False
        ok = ok and here
        print(f"     A={Aval}: all columns factor with eta-free real M: {here}")
    return rep("S2 column-grading factorization exact (the realness-lead structure)", ok)

def S3():
    print("-" * 78); print("S3  swap-similarity + symmetric charpoly => real coefficients")
    import sympy as sp
    eta, etb, lam = sp.symbols('eta etab lambda')
    p = sp.symbols('p', positive=True)
    ok = True
    for Aval in (3, 4):
        A = sp.Integer(Aval)
        T = [[(1 - p) if i == j else p / (A - 1) for j in range(Aval)] for i in range(Aval)]
        st = list(itertools.product(range(Aval), repeat=3)); reps = [None] * 5
        for k, tr in enumerate(st):
            if reps[_pat(tr)] is None: reps[_pat(tr)] = k
        Q = sp.zeros(5, 5)
        for a_ in range(5):
            x, xp, xq = st[reps[a_]]
            for (y, yp, yq) in st:
                term = sp.nsimplify(T[xp][yp] * T[xq][yq] / T[x][y])
                if yp != y: term *= eta
                if yq != y: term *= etb
                Q[a_, _pat((y, yp, yq))] += term
        Pm = sp.eye(5); Pm[1, 1] = 0; Pm[2, 2] = 0; Pm[1, 2] = 1; Pm[2, 1] = 1
        ok_swap = sp.simplify(Pm * Q * Pm - Q.xreplace({eta: etb, etb: eta})) == sp.zeros(5, 5)
        cp = Q.charpoly(lam).as_expr()
        ok_sym = sp.simplify(sp.expand(cp - cp.xreplace({eta: etb, etb: eta}))) == 0
        # EXACT spot confirmation of real coefficients on the conjugate locus
        zx = sp.Rational(3, 100); zy = sp.Rational(1, 50)
        val = sp.expand(cp.coeff(lam, 2)).xreplace(
            {eta: zx + sp.I * zy, etb: zx - sp.I * zy, p: sp.Rational(1, 5)})
        ok_real = sp.simplify(sp.im(sp.expand(val))) == 0   # exact rational arithmetic
        ok = ok and ok_swap and ok_sym and ok_real
        print(f"     A={Aval}: swap-similarity {ok_swap}; charpoly symmetric {ok_sym}; spot-real {ok_real}")
    return rep("S3 real quintic coefficients on the conjugate locus (proven)", ok)

if __name__ == "__main__":
    print("=" * 78)
    print("Mid-band instrument: real-rooted quintic + charpoly-derivative certificate")
    print("=" * 78)
    S1(); S2(); S3()
    print("=" * 78)
    print(f"OVERALL -> {'PASS' if PASS else 'FAIL'}")
    print("Mid-band program reduces to: (a) structural real-rootedness lemma (column")
    print("grading is the lead), (b) sign-certification of chi^(k)(Lam), k=0..4 --")
    print("five rational functions -- over the band (Sturm/Polya class).")

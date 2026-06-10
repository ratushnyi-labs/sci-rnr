#!/usr/bin/env python3
"""
probe_7_34_G_analytic_bounds.py
============================================================================
ANALYTIC PROOF of Hypothesis (G) of the 7.34b second-moment (replica) route
(tex block near line 12651; ground truth conventions:
probe_7_34_replica_secondmoment_tail.py).

CLAIM PROVED (full range): for every p in (0,1/2), 0 < D <= D_c(p), and
s in (0,pi],
    g(s) := |C_s/C_0|^2 * rho(W4(s))  <=  1 - (11/16) D(1-D)(1-cos s)  <  1.

PROOF SKELETON (each step machine-validated below):

  S1  STRUCTURE.  W4[(d1,d2),(e1,e2)] = B(d1^e1, d2^e2) * eta^{e1} etabar^{e2}
      with B(f1,f2) = sum_a T_{a^f1} T_{a^f2} / T_a depending only on the
      flip pair: B(0,0)=B(0,1)=B(1,0)=1, B(1,1)=1+q,
          q := (1-2p)^2 / (p(1-p))  >= 0.
      Hence  W4 = 1 m^T + q Pi M,   m=(1,etabar,eta,|eta|^2)^T, M=diag(m),
      Pi = reversal permutation (d1,d2)->(1-d1,1-d2), 1 = all-ones vector.

  S2  CHARACTERISTIC POLYNOMIAL (matrix determinant lemma).  Pi M is block-
      antidiagonal on the Pi-pairs {(0,0),(1,1)}, {(0,1),(1,0)} with both
      pair-products m1 m4 = m2 m3 = |eta|^2, so
        det(lam I - q Pi M) = (lam^2 - q^2 |eta|^2)^2,
        m^T (lam I - q Pi M)^{-1} 1 = [lam T + 4 q |eta|^2]/(lam^2-q^2|eta|^2),
        T := sum(m) = |1+eta|^2,
      and therefore, EXACTLY,
        det(lam I - W4) = (lam^2 - q^2|eta|^2) (lam^2 - T lam - kappa^2 |eta|^2),
        kappa^2 := q(q+4) = (1-2p)^2/(p^2(1-p)^2)   [q+4 = 1/(p(1-p))].
      Real coefficients (the (H1) consequence) and an EXACT quadratic x
      quadratic factorization.

  S3  PERRON ROOT.  Roots: +-q|eta| and lam_pm = [T +- sqrt(T^2+4kappa^2
      |eta|^2)]/2.  lam_+ >= kappa|eta| > q|eta| (p<1/2) and |lam_-| =
      kappa^2|eta|^2/lam_+ <= lam_+.  So rho(W4) = lam_+, always carried by
      the second quadratic factor.

  S4  DUAL OBJECTS IN CLOSED FORM.  With w = e^{is}, A=(1-D)^2, B=D^2:
        eta_s        = D(1-D)(w-1) / (A - B w),
        C_s/C_0      = (A - B w) / (1-2D),
        1 + eta_s    = (1-2D)(1-D+Dw) / (A - B w),
        |A - B w|^2  = (1-2D)^2 + 2 D^2 (1-D)^2 t,    t := 1-cos s in (0,2].
      Hence  F := |C_s/C_0|^2 |1+eta_s|^2 = |1-D+D w|^2 = 1 - 2D(1-D) t
      (the Jensen floor, exactly), G := 4 kappa^2 |eta_s|^2 |C_s/C_0|^4
      = 8 kappa^2 D^2(1-D)^2 t [(1-2D)^2 + 2D^2(1-D)^2 t]/(1-2D)^4,  and
        g(s) = [F + sqrt(F^2+G)]/2     (>= F: floor manifest).

  S5  EQUIVALENCE.  0 < F < 1 on t in (0,2], so
        g < 1  <=>  sqrt(F^2+G) < 2-F  <=>  G < 4(1-F) = 8D(1-D)t
               <=>  kappa^2 D(1-D)[(1-2D)^2 + 2D^2(1-D)^2 t] < (1-2D)^4.
      LHS affine increasing in t => worst case t=2:
        Psi(D) := kappa^2 D(1-D)[(1-2D)^2+4D^2(1-D)^2]/(1-2D)^4 < 1.

  S6  GRAY REGION.  Psi is strictly increasing on (0,1/2):
      d/dD [D(1-D)/(1-2D)^2] = 1/(1-2D)^3  (identity (1-2D)^2+4D(1-D)=1),
      and [D(1-D)]^3/(1-2D)^4 is a product of increasing positive factors.
      At D = D_c(p) (1-2D_c = sqrt(1-2p)/(1-p), D_c(1-D_c) = p^2/(4(1-p)^2)):
        Psi(D_c) = (1-2p)/(4(1-p)^2) + p^4/(16(1-p)^4) < 1/4 + 1/16 = 5/16.
      [(1-2p) < (1-p)^2 since p^2>0;  p/(1-p) < 1 for p<1/2.]
      So Psi(D) <= Psi(D_c) < 5/16 < 1 on the whole closed Gray region.

  S7  EXPLICIT MARGIN.  1-g = [4(1-F)-G] / (2[(2-F)+sqrt(F^2+G)]); the
      denominator < 4(2-F) <= 8, and 4(1-F)-G = 8D(1-D)t[1-Psi_t(D)] >=
      8D(1-D)t(1-Psi(D)) >= 8D(1-D)t * 11/16.  Hence
        1 - g(s) >= (11/16) D(1-D)(1-cos s)        QED.

  COROLLARY (curvature / vbar).  As s->0,
        g(s) = 1 - a(p,D) s^2 + O(s^4),
        a(p,D) = D(1-D) [1 - kappa^2 D(1-D)/(1-2D)^2] >= (11/16) D(1-D) > 0,
      and a = vbar = lim_n E_x[Var(K|x)]/n (per-site posterior variance).

CHECKS (PASS/FAIL each; tolerances machine-precision unless stated):
  V0   conventions anchor: brute 2^n E_x|Phi|^2 == |C|^{2n} u4 W4^{n-1} v4.
  V1   (H1) swap symmetry P W4 P == conj(W4); char-poly coefficients real.
  V2   (H2) Kraus/CP bookkeeping: W4 == sum_a T_a^{-1} kron(M_a, conj(M_a))
       (row-major vec); functional u4^T W4^{n-1} v4 == u^T Phi_CP^{n-1}(J2) ubar.
  V3   S1 structure: W4 == 1 m^T + q Pi M, q = (1-2p)^2/(p(1-p)).
  V4   S2+S3: symbolic (sympy) char-poly factorization in (lam,eta,etabar,q);
       numeric eigenvalues == {+-q|eta|, lam_pm}; rho(W4) == lam_+ everywhere.
  V5   S4 closed forms for eta_s, C_s/C_0, 1+eta_s, F, G vs probe numerics.
  V6   g closed form == probe g(s) on a hostile grid (incl. D=D_c, p=0.005).
  V7   S5+S6: sign(1-g) == sign((1-2D)^4 - kappa^2 D(1-D)[...t]) pointwise;
       sympy identities (d/dD, Psi(D_c) formula); Psi(D_c)<5/16 dense p-grid;
       Psi monotone in D.
  V8   S7 final bound g(s) <= 1-(11/16)D(1-D)(1-cos s) on the hostile grid.
  V9   curvature: closed-form a vs numeric fit; vs committed range
       [0.0014, 0.0894] on the 6 primary combos; vbar consistency
       E_x[Var(K|x)]/n -> a (transfer-matrix finite difference, n<=800).
  V10  auxiliary-bound coverage map (informational; validity is PASS/FAIL):
       (i) CP bound |C_s/C_0|^2 ||sum_a T_a^{-1} M_a^dag M_a||_op,
       (ii) abs bound |C_s/C_0|^2 rho(|W4|)
            (closed form: [(1+|eta|)^2 + sqrt((1+|eta|)^4+4kappa^2|eta|^2)]/2)
       -- both are valid upper bounds on g everywhere; report where each is
       below 1.  (Superseded by the exact formula, which covers 100%.)
  V11a FULL symbolic (sympy) derivation of S4 -- eta_s, C_s/C_0, 1+eta_s, F
       (Jensen floor identity), |denominator|^2, G -- as exact rational-
       function identities in (D, w), conjugation = w -> 1/w on |w|=1
       (upgrades V5 from 50-dps numerics to exact symbolic), PLUS the
       squared-form margin identities:
         2B - F = 1 + (5/8) dd t  > 0,
         (2B-F)^2 - F^2 = 4B(B-F),   B - F = (21/16) dd t,
       which reduce the final bound g <= B := 1-(11/16)dd t to
         Psi_t <= (21/32) B   [true: (21/32)B >= (21/32)^2 = 441/1024
                               > 5/16 >= Psi_t on Gray].
  V11b EXACT rational (fractions.Fraction) verification of the entire
       inequality chain at rational Gray points with RATIONAL D_c
       (p = (1-r^2)/2 => sqrt(1-2p) = r, D_c = (1-r)^2/(2(1+r^2))),
       down to p ~ 1e-4 / D ~ 2.5e-9 where double precision cannot resolve
       the slack: Psi(D_c) == closed form < 5/16; Psi_t <= Psi_2(D) <=
       Psi_2(D_c); 0 < F < 1; 2B-F > 0; (2B-F)^2 - F^2 - G > 0 (strict).
       Zero rounding anywhere.

VERDICT line at the end:
  proven-full-range  iff  V0..V9, V11a, V11b all PASS and V10 validity PASS.
"""
import math
import numpy as np

try:
    import sympy as sp
    HAVE_SYMPY = True
except ImportError:  # symbolic checks degraded to numeric-poly checks
    HAVE_SYMPY = False

try:
    import mpmath as mp
    mp.mp.dps = 50
    HAVE_MPMATH = True
except ImportError:
    HAVE_MPMATH = False

# ------------- committed conventions (verbatim port of the ground truth) -------------

def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))

def Dc_stable(p):
    """Cancellation-free form: D_c = p^2 / (2(1-p)[(1-p)+sqrt(1-2p)])
    (== Dc(p) exactly; sympy-verified in V7)."""
    return p * p / (2 * (1 - p) * ((1 - p) + math.sqrt(1 - 2 * p)))

def theta0(D):
    return math.log(D / (1 - D))

def eta_C(beta, D):
    eb = np.exp(beta)
    ze = (1 - eb) / ((1 + eb) * (1 - 2 * D))
    eta = (1 - ze) / (1 + ze)
    C = (1 + eb) * (1 + ze) / 2
    return eta, C

def eta_s(s, D):
    return eta_C(theta0(D) + 1j * s, D)[0]

def C_ratio(s, D):
    th = theta0(D)
    return abs(eta_C(th + 1j * s, D)[1] / eta_C(th, D)[1])

def replica_W4(s, p, D):
    et = eta_s(s, D)
    etb = np.conj(et)
    T = (1 - p, p)
    states = ((0, 0), (0, 1), (1, 0), (1, 1))
    W = np.zeros((4, 4), complex)
    for i, (d1, d2) in enumerate(states):
        for j, (e1, e2) in enumerate(states):
            c = 0.0
            for a in (0, 1):
                c += T[a ^ d1 ^ e1] * T[a ^ d2 ^ e2] / T[a]
            W[i, j] = c * (et if e1 else 1.0) * (etb if e2 else 1.0)
    u = np.array([(et if d1 else 1.0) * (etb if d2 else 1.0)
                  for (d1, d2) in states], complex)
    return W, u

def rho_num(W):
    return float(np.max(np.abs(np.linalg.eigvals(W))))

def g_probe(s, p, D):
    W, _ = replica_W4(s, p, D)
    return C_ratio(s, D) ** 2 * rho_num(W)

# ------------- closed forms proved in this probe -------------

def q_of(p):
    return (1 - 2 * p) ** 2 / (p * (1 - p))

def kappa2_of(p):
    return (1 - 2 * p) ** 2 / (p ** 2 * (1 - p) ** 2)

def FG_of(t, p, D):
    """F = 1-2D(1-D)t ;  G = 8 k^2 D^2(1-D)^2 t [(1-2D)^2+2D^2(1-D)^2 t]/(1-2D)^4."""
    k2 = kappa2_of(p)
    dd = D * (1 - D)
    c = (1 - 2 * D) ** 2
    F = 1 - 2 * dd * t
    G = 8 * k2 * dd ** 2 * t * (c + 2 * dd ** 2 * t) / c ** 2
    return F, G

def g_closed(s, p, D):
    t = 1 - math.cos(s)
    F, G = FG_of(t, p, D)
    return (F + math.sqrt(F * F + G)) / 2

def Psi_t(t, p, D):
    """Psi_t(D) := kappa^2 D(1-D)[(1-2D)^2 + 2D^2(1-D)^2 t]/(1-2D)^4;
    g(s)<1 <=> Psi_{1-cos s} < 1; full (G) <=> Psi_2 = Psi < 1."""
    k2 = kappa2_of(p)
    dd = D * (1 - D)
    c = (1 - 2 * D) ** 2
    return k2 * dd * (c + 2 * dd ** 2 * t) / c ** 2

def Psi(p, D):
    return Psi_t(2.0, p, D)

def a_closed(p, D):
    return D * (1 - D) * (1 - Psi_t(0.0, p, D))

def eta_closed(s, D):
    w = np.exp(1j * s)
    return D * (1 - D) * (w - 1) / ((1 - D) ** 2 - D ** 2 * w)

def Cratio_closed(s, D):
    w = np.exp(1j * s)
    return ((1 - D) ** 2 - D ** 2 * w) / (1 - 2 * D)

# ------------- grids -------------

def hostile_grid():
    ps = [0.005, 0.02, 0.05, 0.1, 0.25, 0.4, 0.45, 0.48]
    fDs = [0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]          # incl. D = D_c
    ss = np.concatenate([np.geomspace(1e-4, 0.1, 40, endpoint=False),
                         np.linspace(0.1, math.pi, 160)])
    return ps, fDs, ss

# ============================ V0 conventions anchor ============================

def popcount_arr(v):
    v = np.asarray(v, dtype=np.int64)
    pc = np.zeros_like(v)
    while np.any(v):
        pc += v & 1
        v = v >> 1
    return pc

def fwht(a):
    a = a.astype(float).copy()
    h = 1
    while h < len(a):
        for i in range(0, len(a), h * 2):
            u = a[i:i + h].copy(); w = a[i + h:i + 2 * h].copy()
            a[i:i + h] = u + w; a[i + h:i + 2 * h] = u - w
        h *= 2
    return a

def bsms_P(n, p):
    N = 1 << n
    v = np.arange(N)
    sw = popcount_arr((v ^ (v >> 1)) & ((1 << (n - 1)) - 1))
    return 0.5 * (1 - p) ** (n - 1 - sw) * p ** sw

def brute_E2(s_list, n, p, D):
    N = 1 << n
    P = bsms_P(n, p)
    Ph = fwht(P)
    pc_all = popcount_arr(np.arange(N))
    PYs = fwht(Ph * (1 - 2 * D) ** (-pc_all.astype(float))) / N
    assert PYs.min() > -1e-9, "deconvolution must be valid in-Gray"
    th = theta0(D)
    ks = np.arange(n + 1)
    phase = np.exp(1j * np.outer(ks, np.asarray(s_list)))
    E2 = np.zeros(len(s_list))
    ix = np.arange(N)
    for xx in range(N):
        d = pc_all[ix ^ xx]
        w0 = PYs * np.exp(th * d)
        A = np.bincount(d, weights=w0, minlength=n + 1)
        Phi = (A @ phase) / A.sum()
        E2 += P[xx] * np.abs(Phi) ** 2
    return E2

def E2_tm(s, p, D, n):
    W, u = replica_W4(s, p, D)
    val = u @ np.linalg.matrix_power(W, n - 1) @ np.ones(4)
    return C_ratio(s, D) ** (2 * n) * val

def check_V0():
    print("-" * 84)
    print("V0: conventions anchor -- brute 2^n E_x|Phi|^2 == |C|^{2n} u4 W4^{n-1} v4")
    ok = True; worst = 0.0
    n = 10
    s_list = [0.3, 1.2, 2.5, math.pi]
    for (p, fD) in [(0.25, 0.9), (0.1, 0.5)]:
        D = fD * Dc(p)
        E2b = brute_E2(s_list, n, p, D)
        for k, s in enumerate(s_list):
            v = E2_tm(s, p, D, n)
            rel = abs(E2b[k] - v.real) / max(abs(E2b[k]), 1e-300)
            worst = max(worst, rel, abs(v.imag) / max(abs(v), 1e-300))
            ok &= rel < 1e-9
    print(f"  worst rel.diff = {worst:.2e}   -> {'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V1 (H1) ============================

def check_V1():
    print("-" * 84)
    print("V1: (H1) P W4 P == conj(W4); char-poly coefficients real")
    P = np.zeros((4, 4))
    # states ((0,0),(0,1),(1,0),(1,1)) -> swap (d1,d2)->(d2,d1) swaps idx 1<->2
    P[0, 0] = P[3, 3] = P[1, 2] = P[2, 1] = 1.0
    ok = True; worst_sym = 0.0; worst_im = 0.0
    cases = [(p, fD, s) for p in (0.005, 0.1, 0.25, 0.48)
             for fD in (0.3, 0.9, 1.0) for s in (0.05, 1.3, 2.9, math.pi)]
    for (p, fD, s) in cases:
        D = fD * Dc(p)
        W, _ = replica_W4(s, p, D)
        d = np.max(np.abs(P @ W @ P - np.conj(W)))
        worst_sym = max(worst_sym, d)
        coeffs = np.poly(W)
        im = np.max(np.abs(np.imag(coeffs))) / max(np.max(np.abs(coeffs)), 1e-300)
        worst_im = max(worst_im, im)
        ok &= d < 1e-12 and im < 1e-12
    print(f"  swap-symmetry worst |PWP - conj(W)| = {worst_sym:.2e}; "
          f"char-poly worst rel imag coeff = {worst_im:.2e}   -> "
          f"{'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V2 (H2) ============================

def kraus_M(a, s, p, D):
    et = eta_s(s, D)
    T = (1 - p, p)
    M = np.zeros((2, 2), complex)
    for d in (0, 1):
        for e in (0, 1):
            M[d, e] = T[a ^ d ^ e] * (et if e else 1.0)
    return M

def check_V2():
    print("-" * 84)
    print("V2: (H2) W4 == sum_a T_a^{-1} kron(M_a, conj(M_a)) (row-major vec);")
    print("    boundary functional u4^T W4^{n-1} v4 == u^T Phi_CP^{n-1}(J2) conj(u)")
    ok = True; worst_k = 0.0; worst_f = 0.0
    n = 9
    cases = [(p, fD, s) for p in (0.1, 0.25, 0.4)
             for fD in (0.5, 1.0) for s in (0.4, 2.0, math.pi)]
    for (p, fD, s) in cases:
        D = fD * Dc(p)
        T = (1 - p, p)
        W, u4 = replica_W4(s, p, D)
        Wk = sum(np.kron(kraus_M(a, s, p, D), np.conj(kraus_M(a, s, p, D))) / T[a]
                 for a in (0, 1))
        dk = np.max(np.abs(W - Wk)) / max(np.max(np.abs(W)), 1e-300)
        worst_k = max(worst_k, dk)
        ok &= dk < 1e-12
        # CP-map iteration: Phi(X) = sum_a T_a^{-1} M_a X M_a^dagger,
        # X0 = J2 (all-ones 2x2) = unvec_r(v4); functional u^T X conj(u), u=(1,eta)
        X = np.ones((2, 2), complex)
        Ms = [kraus_M(a, s, p, D) for a in (0, 1)]
        for _ in range(n - 1):
            X = sum(M @ X @ M.conj().T / T[a] for a, M in enumerate(Ms))
        uvec = np.array([1.0, eta_s(s, D)], complex)
        val_cp = uvec @ X @ np.conj(uvec)
        val_tm = u4 @ np.linalg.matrix_power(W, n - 1) @ np.ones(4)
        df = abs(val_cp - val_tm) / max(abs(val_tm), 1e-300)
        worst_f = max(worst_f, df)
        ok &= df < 1e-10
    print(f"  worst Kraus rel.diff = {worst_k:.2e}; worst functional rel.diff = "
          f"{worst_f:.2e}   -> {'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V3 structure ============================

def check_V3():
    print("-" * 84)
    print("V3: S1 -- W4 == 1 m^T + q Pi M with q = (1-2p)^2/(p(1-p))")
    ok = True; worst = 0.0
    Pi = np.fliplr(np.eye(4))           # (d1,d2) -> (1-d1,1-d2) reversal
    for p in (0.005, 0.1, 0.25, 0.4, 0.48):
        q = q_of(p)
        # algebra check: B(1,1)-1 == q
        Bq = p ** 2 / (1 - p) + (1 - p) ** 2 / p - 1
        ok &= abs(Bq - q) < 1e-14 * max(q, 1)
        for fD in (0.5, 1.0):
            D = fD * Dc(p)
            for s in (0.05, 1.3, math.pi):
                W, _ = replica_W4(s, p, D)
                et = eta_s(s, D)
                m = np.array([1, np.conj(et), et, abs(et) ** 2], complex)
                W_struct = np.outer(np.ones(4), m) + q * Pi @ np.diag(m)
                d = np.max(np.abs(W - W_struct)) / max(np.max(np.abs(W)), 1e-300)
                worst = max(worst, d)
                ok &= d < 1e-12
    print(f"  worst rel.diff = {worst:.2e}   -> {'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V4 char poly + Perron ============================

def check_V4():
    print("-" * 84)
    print("V4: S2+S3 -- char poly == (lam^2-q^2|eta|^2)(lam^2-|1+eta|^2 lam"
          "-kappa^2|eta|^2);")
    print("    eigenvalues == {+-q|eta|, lam_pm}; rho(W4) == lam_+ everywhere")
    ok = True
    # ---- symbolic identity (general eta, etabar treated independent) ----
    if HAVE_SYMPY:
        lam, e1, e2, qs = sp.symbols('lambda eta etabar q')
        m = [sp.Integer(1), e2, e1, e1 * e2]
        Pi_idx = [3, 2, 1, 0]
        W = sp.zeros(4, 4)
        for i in range(4):
            for j in range(4):
                W[i, j] = m[j] + (qs * m[j] if Pi_idx[i] == j else 0)
        charpoly = sp.expand((lam * sp.eye(4) - W).det())
        Tsym = sp.expand((1 + e1) * (1 + e2))
        target = sp.expand((lam ** 2 - qs ** 2 * e1 * e2)
                           * (lam ** 2 - Tsym * lam - qs * (qs + 4) * e1 * e2))
        sym_ok = sp.simplify(charpoly - target) == 0
        # q algebra: q+4 == 1/(p(1-p)),  q(q+4) == (1-2p)^2/(p^2(1-p)^2)
        pp = sp.symbols('p', positive=True)
        qq = (1 - 2 * pp) ** 2 / (pp * (1 - pp))
        alg_ok = (sp.simplify(qq + 4 - 1 / (pp * (1 - pp))) == 0 and
                  sp.simplify(qq * (qq + 4)
                              - (1 - 2 * pp) ** 2 / (pp ** 2 * (1 - pp) ** 2)) == 0)
        print(f"  sympy: det(lam I - W4) factorization identity: "
              f"{'OK' if sym_ok else 'FAIL'};  q+4 = 1/(p(1-p)) and "
              f"q(q+4) = kappa^2: {'OK' if alg_ok else 'FAIL'}")
        ok &= sym_ok and alg_ok
    else:
        print("  sympy unavailable -- symbolic step degraded to numeric only")
    # ---- numeric roots on hostile grid ----
    worst_e = 0.0; worst_r = 0.0; perron_in_quad = True
    ps, fDs, _ = hostile_grid()
    ss = [1e-3, 0.05, 0.5, 1.3, 2.2, 2.9, math.pi]
    for p in ps:
        q = q_of(p); kap = math.sqrt(kappa2_of(p))
        for fD in fDs:
            D = fD * Dc(p)
            for s in ss:
                W, _ = replica_W4(s, p, D)
                ae = abs(eta_s(s, D))
                T = abs(1 + eta_s(s, D)) ** 2
                disc = math.sqrt(T * T + 4 * kap ** 2 * ae ** 2)
                lam_p = (T + disc) / 2
                lam_m = (T - disc) / 2
                ev = np.sort_complex(np.linalg.eigvals(W))
                pred = np.sort_complex(np.array(
                    [q * ae, -q * ae, lam_p, lam_m], dtype=complex))
                de = np.max(np.abs(ev - pred)) / max(np.max(np.abs(ev)), 1e-300)
                dr = abs(rho_num(W) - lam_p) / max(lam_p, 1e-300)
                worst_e = max(worst_e, de); worst_r = max(worst_r, dr)
                ok &= de < 1e-10 and dr < 1e-12
                perron_in_quad &= lam_p >= q * ae - 1e-15
    print(f"  numeric: worst eigenvalue-set rel.diff = {worst_e:.2e}; "
          f"worst rho-vs-lam_+ rel.diff = {worst_r:.2e}; "
          f"lam_+ >= q|eta| everywhere: {perron_in_quad}   -> "
          f"{'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V5 dual closed forms ============================

def check_V5():
    print("-" * 84)
    print("V5: S4 -- closed forms eta_s, C_s/C_0, 1+eta_s, F, G vs the committed")
    print("    eta_C route, evaluated at 50-digit precision (the double-precision")
    print("    eta_C reference cancels catastrophically at small s, small D;")
    print("    double-precision end-to-end agreement of g is V6)")
    if not HAVE_MPMATH:
        print("  mpmath unavailable -- SKIPPED (degraded); marking FAIL")
        return False
    ok = True
    w1 = w2 = w3 = w4 = w5 = mp.mpf(0)
    tol = mp.mpf('1e-35')
    ps, fDs, ss = hostile_grid()
    for p in ps:
        for fD in fDs:
            D = mp.mpf(fD) * mp.mpf(Dc(p))   # same double D fed to both sides
            pm = mp.mpf(p)
            th = mp.log(D / (1 - D))
            kap2 = ((1 - 2 * pm) / (pm * (1 - pm))) ** 2

            def eta_C_mp(beta):
                eb = mp.exp(beta)
                ze = (1 - eb) / ((1 + eb) * (1 - 2 * D))
                return (1 - ze) / (1 + ze), (1 + eb) * (1 + ze) / 2

            C0 = eta_C_mp(th)[1]
            for s in ss[::7]:
                sm = mp.mpf(s)
                w = mp.exp(1j * sm)
                et_n, Cs = eta_C_mp(th + 1j * sm)
                denom = (1 - D) ** 2 - D ** 2 * w
                et_c = D * (1 - D) * (w - 1) / denom
                cr_c = denom / (1 - 2 * D)
                one_eta = (1 - 2 * D) * (1 - D + D * w) / denom
                d1 = abs(et_n - et_c) / max(abs(et_n), mp.mpf('1e-300'))
                d2 = abs(Cs / C0 - cr_c) / abs(Cs / C0)
                d3 = abs(1 + et_n - one_eta) / abs(1 + et_n)
                t = 1 - mp.cos(sm)
                dd = D * (1 - D)
                c2 = (1 - 2 * D) ** 2
                F = 1 - 2 * dd * t
                G = 8 * kap2 * dd ** 2 * t * (c2 + 2 * dd ** 2 * t) / c2 ** 2
                F_n = abs(Cs / C0) ** 2 * abs(1 + et_n) ** 2
                G_n = 4 * kap2 * abs(et_n) ** 2 * abs(Cs / C0) ** 4
                d4 = abs(F - F_n) / F_n
                d5 = abs(G - G_n) / max(G_n, mp.mpf('1e-300'))
                # Jensen-floor identity F == |1-D+De^{is}|^2
                d4b = abs(F - abs(1 - D + D * w) ** 2)
                w1 = max(w1, d1); w2 = max(w2, d2); w3 = max(w3, d3)
                w4 = max(w4, d4, d4b); w5 = max(w5, d5)
                ok &= max(d1, d2, d3, d4, d4b, d5) < tol
    print(f"  worst rel.diffs (50 dps): eta {mp.nstr(w1, 3)}, "
          f"C_s/C_0 {mp.nstr(w2, 3)}, 1+eta {mp.nstr(w3, 3)}, "
          f"F {mp.nstr(w4, 3)}, G {mp.nstr(w5, 3)}   -> "
          f"{'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V6 g closed form ============================

def check_V6():
    print("-" * 84)
    print("V6: g(s) = [F+sqrt(F^2+G)]/2 == probe g(s) on the hostile grid "
          "(incl. D=D_c, p=0.005)")
    ok = True; worst = 0.0; npts = 0
    ps, fDs, ss = hostile_grid()
    for p in ps:
        for fD in fDs:
            D = fD * Dc(p)
            for s in ss[::4]:
                gn = g_probe(s, p, D)
                gc = g_closed(s, p, D)
                d = abs(gn - gc) / max(gn, 1e-300)
                worst = max(worst, d)
                ok &= d < 1e-11
                npts += 1
    print(f"  {npts} grid points; worst rel.diff = {worst:.2e}   -> "
          f"{'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V7 equivalence + Psi ============================

def check_V7():
    print("-" * 84)
    print("V7: S5+S6 -- g<1 <=> Psi_t<1 pointwise; Psi monotone in D; "
          "Psi(D_c) = (1-2p)/(4(1-p)^2) + (p/(1-p))^4/16 < 5/16")
    ok = True
    # (i) pointwise sign equivalence
    sign_ok = True
    ps, fDs, ss = hostile_grid()
    for p in ps:
        for fD in fDs:
            D = fD * Dc(p)
            for s in ss[::5]:
                t = 1 - math.cos(s)
                lhs = g_closed(s, p, D) < 1.0
                rhs = Psi_t(t, p, D) < 1.0
                sign_ok &= (lhs == rhs)
    print(f"  (i) pointwise g<1 <=> Psi_t<1: {'OK' if sign_ok else 'FAIL'}")
    ok &= sign_ok
    # (ii) sympy identities
    if HAVE_SYMPY:
        Ds, pp = sp.symbols('D p', positive=True)
        f1 = Ds * (1 - Ds) / (1 - 2 * Ds) ** 2
        id1 = sp.simplify(sp.diff(f1, Ds) - 1 / (1 - 2 * Ds) ** 3) == 0
        id1b = sp.simplify((1 - 2 * Ds) ** 2 + 4 * Ds * (1 - Ds) - 1) == 0
        k2 = (1 - 2 * pp) ** 2 / (pp ** 2 * (1 - pp) ** 2)
        Psi_sym = (k2 * Ds * (1 - Ds)
                   * ((1 - 2 * Ds) ** 2 + 4 * Ds ** 2 * (1 - Ds) ** 2)
                   / (1 - 2 * Ds) ** 4)
        Dc_sym = (1 - sp.sqrt(1 - 2 * pp) / (1 - pp)) / 2
        target = (1 - 2 * pp) / (4 * (1 - pp) ** 2) + pp ** 4 / (16 * (1 - pp) ** 4)
        id2 = sp.simplify(Psi_sym.subs(Ds, Dc_sym) - target) == 0
        # 1-2p < (1-p)^2  <=>  p^2 > 0  (first term < 1/4)
        id3 = sp.simplify(sp.expand((1 - pp) ** 2 - (1 - 2 * pp)) - pp ** 2) == 0
        # cancellation-free D_c form used in (iii)
        id4 = sp.simplify(
            (1 - sp.sqrt(1 - 2 * pp) / (1 - pp)) / 2
            - pp ** 2 / (2 * (1 - pp) * ((1 - pp) + sp.sqrt(1 - 2 * pp)))) == 0
        print(f"  (ii) sympy: d/dD[D(1-D)/(1-2D)^2]=1/(1-2D)^3: "
              f"{'OK' if id1 else 'FAIL'}; (1-2D)^2+4D(1-D)=1: "
              f"{'OK' if id1b else 'FAIL'}; Psi(D_c) closed form: "
              f"{'OK' if id2 else 'FAIL'}; (1-p)^2-(1-2p)=p^2: "
              f"{'OK' if id3 else 'FAIL'}; D_c stable form: "
              f"{'OK' if id4 else 'FAIL'}")
        ok &= id1 and id1b and id2 and id3 and id4
    else:
        print("  (ii) sympy unavailable -- skipped (numeric (iii)/(iv) still run)")
    # (iii) Psi(D_c) < 5/16 on a dense p grid incl. extremes
    #       (Dc_stable: the naive 0.5(1-sqrt(1-2p)/(1-p)) cancels at p->0)
    pgrid = np.concatenate([np.geomspace(1e-6, 0.01, 200),
                            np.linspace(0.01, 0.499999, 2000)])
    psi_dc = np.array([Psi(p, Dc_stable(p)) for p in pgrid])
    psi_dc_cf = ((1 - 2 * pgrid) / (4 * (1 - pgrid) ** 2)
                 + (pgrid / (1 - pgrid)) ** 4 / 16)
    cf_ok = np.max(np.abs(psi_dc - psi_dc_cf)
                   / np.maximum(psi_dc_cf, 1e-300)) < 1e-9
    bound_ok = bool(np.all(psi_dc < 5.0 / 16.0))
    print(f"  (iii) Psi(D_c) closed form numeric match: {'OK' if cf_ok else 'FAIL'}"
          f"; max Psi(D_c) = {psi_dc.max():.6f} < 5/16 = {5/16:.6f}: "
          f"{'OK' if bound_ok else 'FAIL'}")
    ok &= cf_ok and bound_ok
    # (iv) Psi strictly increasing in D on (0, 1/2) (numeric, dense)
    mono_ok = True
    for p in (0.005, 0.05, 0.25, 0.48):
        Dgrid = np.linspace(1e-6, 0.499, 4000)
        vals = np.array([Psi(p, d) for d in Dgrid])
        mono_ok &= bool(np.all(np.diff(vals) > 0))
    print(f"  (iv) Psi strictly increasing in D on (0,1/2): "
          f"{'OK' if mono_ok else 'FAIL'}")
    ok &= mono_ok
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V8 final margin ============================

def check_V8():
    print("-" * 84)
    print("V8: S7 -- g(s) <= 1 - (11/16) D(1-D)(1-cos s) on the hostile grid")
    ok = True
    worst_slack = float('inf')
    npts = 0
    ps, fDs, ss = hostile_grid()
    for p in ps:
        for fD in fDs:
            D = fD * Dc(p)
            dd = D * (1 - D)
            for s in ss:
                g = g_probe(s, p, D) if npts % 13 == 0 else g_closed(s, p, D)
                bound = 1 - (11.0 / 16.0) * dd * (1 - math.cos(s))
                slack = bound - g
                worst_slack = min(worst_slack, slack)
                ok &= g <= bound + 1e-13
                npts += 1
    print(f"  {npts} grid points (closed form, 1-in-13 cross-checked vs probe "
          f"numerics); min slack = {worst_slack:.3e} >= 0: "
          f"{'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V9 curvature / vbar ============================

def check_V9():
    print("-" * 84)
    print("V9: curvature a(p,D) = D(1-D)[1 - kappa^2 D(1-D)/(1-2D)^2] -- "
          "numeric fit, committed range, vbar")
    ok = True
    # (i) numeric curvature limit vs closed form, 50-digit evaluation of the
    #     g closed form (double precision drowns 1-g ~ 1e-14 at tiny s); plus
    #     a double-precision sanity at s=1e-2 vs the probe numerics
    worst = mp.mpf(0) if HAVE_MPMATH else 0.0
    worst_dp = 0.0
    for p in (0.005, 0.1, 0.25, 0.4, 0.48):
        for fD in (0.3, 0.9, 1.0):
            D = fD * Dc(p)
            a_cf = a_closed(p, D)
            if HAVE_MPMATH:
                pm, Dm = mp.mpf(p), mp.mpf(D)
                kap2 = ((1 - 2 * pm) / (pm * (1 - pm))) ** 2
                dd = Dm * (1 - Dm)
                c2 = (1 - 2 * Dm) ** 2
                for s in (mp.mpf('1e-6'), mp.mpf('1e-5')):
                    t = 1 - mp.cos(s)
                    F = 1 - 2 * dd * t
                    G = 8 * kap2 * dd ** 2 * t * (c2 + 2 * dd ** 2 * t) / c2 ** 2
                    g = (F + mp.sqrt(F * F + G)) / 2
                    a_num = (1 - g) / s ** 2
                    d = abs(a_num - mp.mpf(a_cf)) / mp.mpf(a_cf)
                    worst = max(worst, d)
                    ok &= d < mp.mpf('1e-6')   # O(s^2) correction at s<=1e-5
            # double-precision sanity vs probe eigensolver at moderate s
            s = 1e-2
            a_dp = (1 - g_probe(s, p, D)) / s ** 2
            worst_dp = max(worst_dp, abs(a_dp - a_cf) / max(a_cf, 1e-300))
    ok &= worst_dp < 5e-2              # O(s^2) corrections at s=1e-2
    print(f"  (i) (1-g)/s^2 vs a closed form: 50-dps limit (s<=1e-5) worst "
          f"rel.diff = {mp.nstr(worst, 3) if HAVE_MPMATH else 'n/a'}; "
          f"double-prec sanity (s=1e-2, vs probe eig) worst rel.diff = "
          f"{worst_dp:.2e}   {'OK' if ok else 'FAIL'}")
    # (ii) committed primary-grid range [0.0014, 0.0894]
    primary = [(p, fD) for p in (0.1, 0.25, 0.4) for fD in (0.5, 0.9)]
    avals = [a_closed(p, fD * Dc(p)) for (p, fD) in primary]
    rng_ok = (abs(min(avals) - 0.0014) / 0.0014 < 0.05
              and abs(max(avals) - 0.0894) / 0.0894 < 0.01)
    print(f"  (ii) primary-combo a range: [{min(avals):.6f}, {max(avals):.6f}] "
          f"vs committed [0.0014, 0.0894]: {'OK' if rng_ok else 'FAIL'}")
    ok &= rng_ok
    # (iii) a == vbar: E_x[Var(K|x)]/n -> a, gap = O(1/n)
    #       E_x[Var(K|x)] = (1 - E_x|Phi(h;x)|^2)/h^2 + O(h^2),  h tiny
    p, fD = 0.25, 0.9
    D = fD * Dc(p)
    a_cf = a_closed(p, D)
    h = 1e-4
    gaps = []
    for n in (100, 200, 400, 800):
        EV = (1 - E2_tm(h, p, D, n).real) / h ** 2
        gaps.append(EV / n - a_cf)
    decay_ok = all(abs(gaps[i + 1]) < 0.62 * abs(gaps[i]) + 1e-12
                   for i in range(len(gaps) - 1))
    print(f"  (iii) vbar: E[Var(K|x)]/n - a at n=100..800: "
          f"{['%.2e' % g for g in gaps]} (halving = O(1/n)): "
          f"{'OK' if decay_ok else 'FAIL'}")
    ok &= decay_ok
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V10 auxiliary-bound coverage ============================

def cp_bound(s, p, D):
    """|C_s/C_0|^2 ||sum_a T_a^{-1} M_a^dag M_a||_op  (valid: rho(Phi_CP) <=
    ||Phi_CP^*(I)||_op for CP maps, Russo--Dye)."""
    T = (1 - p, p)
    Ms = [kraus_M(a, s, p, D) for a in (0, 1)]
    S = sum(M.conj().T @ M / T[a] for a, M in enumerate(Ms))
    return C_ratio(s, D) ** 2 * float(np.max(np.linalg.eigvalsh(S)))

def abs_bound_closed(s, p, D):
    """|C_s/C_0|^2 rho(|W4|): same factorization with eta -> |eta| gives
    rho(|W4|) = [(1+|eta|)^2 + sqrt((1+|eta|)^4 + 4 kappa^2 |eta|^2)]/2."""
    ae = abs(eta_s(s, D))
    Tp = (1 + ae) ** 2
    kap2 = kappa2_of(p)
    return C_ratio(s, D) ** 2 * (Tp + math.sqrt(Tp ** 2 + 4 * kap2 * ae ** 2)) / 2

def check_V10():
    print("-" * 84)
    print("V10: auxiliary upper bounds -- validity (PASS/FAIL) + coverage map "
          "(informational)")
    ok = True
    worst_abs_cf = 0.0
    print(f"  {'(p, D/Dc)':<16} {'CP bound < 1 on':<22} {'|W4| bound < 1 on':<22}"
          f" coverage(CP, |W4|)")
    ss = np.linspace(1e-3, math.pi, 1200)
    for (p, fD) in [(0.05, 0.9), (0.1, 0.5), (0.25, 0.9), (0.4, 0.9),
                    (0.4, 1.0), (0.48, 0.9)]:
        D = fD * Dc(p)
        cp = np.array([cp_bound(s, p, D) for s in ss])
        ab = np.array([abs_bound_closed(s, p, D) for s in ss])
        # validity: each must dominate g everywhere
        g = np.array([g_closed(s, p, D) for s in ss])
        ok &= bool(np.all(cp >= g - 1e-12)) and bool(np.all(ab >= g - 1e-12))
        # abs-bound closed form vs numeric Perron root of |W4|
        for s in ss[::97]:
            W, _ = replica_W4(s, p, D)
            num = C_ratio(s, D) ** 2 * rho_num(np.abs(W))
            worst_abs_cf = max(worst_abs_cf,
                               abs(num - abs_bound_closed(s, p, D)) / num)
        def interval(v):
            below = v < 1.0
            if not below.any():
                return "never", 0.0
            i0 = int(np.argmax(below))
            return (f"s >= {ss[i0]:.3f}" if below[i0:].all()
                    else "non-interval"), float(below.mean())
        cp_iv, cp_cov = interval(cp)
        ab_iv, ab_cov = interval(ab)
        print(f"  ({p:<4}, {fD:<4})    {cp_iv:<22} {ab_iv:<22} "
              f"({cp_cov:.0%}, {ab_cov:.0%})")
    print(f"  abs-bound closed form vs numeric rho(|W4|): worst rel.diff = "
          f"{worst_abs_cf:.2e}")
    print(f"  validity (both bounds >= g everywhere): {'PASS' if ok else 'FAIL'}")
    print("  [coverage is informational: the EXACT closed form for g covers "
          "100% of (p, D, s); the auxiliary bounds are superseded]")
    return ok

# ============================ V11a symbolic S4 + margin identities ============================

def check_V11a():
    print("-" * 84)
    print("V11a: full symbolic derivation of S4 + squared-form margin identities")
    print("      (exact rational-function identities in (D, w); conj = w -> 1/w)")
    if not HAVE_SYMPY:
        print("  sympy unavailable -- FAIL (required for the symbolic chain)")
        return False
    D, w, k2 = sp.symbols('D w kappa2')
    eb = D / (1 - D) * w
    ze = (1 - eb) / ((1 + eb) * (1 - 2 * D))
    eta = (1 - ze) / (1 + ze)
    C = (1 + eb) * (1 + ze) / 2
    C0 = C.subs(w, 1)
    denom = (1 - D) ** 2 - D ** 2 * w
    eta_cf = D * (1 - D) * (w - 1) / denom
    Cr_cf = denom / (1 - 2 * D)
    t = 1 - (w + 1 / w) / 2
    ids = {}
    ids['eta = D(1-D)(w-1)/((1-D)^2-D^2 w)'] = sp.simplify(eta - eta_cf) == 0
    ids['C_0 = 1/(1-D)'] = sp.simplify(C0 - 1 / (1 - D)) == 0
    ids['C/C_0 = ((1-D)^2-D^2 w)/(1-2D)'] = sp.simplify(C / C0 - Cr_cf) == 0
    ids['1+eta = (1-2D)(1-D+Dw)/((1-D)^2-D^2 w)'] = sp.simplify(
        1 + eta_cf - (1 - 2 * D) * (1 - D + D * w) / denom) == 0
    Fw = sp.simplify(Cr_cf * (1 + eta_cf))          # = 1-D+Dw exactly
    ids['(C/C_0)(1+eta) = 1-D+Dw'] = sp.simplify(Fw - (1 - D + D * w)) == 0
    ids['F = |1-D+Dw|^2 = 1-2D(1-D)t'] = sp.simplify(
        Fw * Fw.subs(w, 1 / w) - (1 - 2 * D * (1 - D) * t)) == 0
    ids['|(1-D)^2-D^2 w|^2 = (1-2D)^2+2D^2(1-D)^2 t'] = sp.simplify(
        denom * denom.subs(w, 1 / w)
        - ((1 - 2 * D) ** 2 + 2 * D ** 2 * (1 - D) ** 2 * t)) == 0
    G_sym = 4 * k2 * (eta_cf * eta_cf.subs(w, 1 / w)) \
        * (Cr_cf * Cr_cf.subs(w, 1 / w)) ** 2
    G_cf = (8 * k2 * D ** 2 * (1 - D) ** 2 * t
            * ((1 - 2 * D) ** 2 + 2 * D ** 2 * (1 - D) ** 2 * t)
            / (1 - 2 * D) ** 4)
    ids['G = 8 k^2 D^2(1-D)^2 t [(1-2D)^2+2D^2(1-D)^2 t]/(1-2D)^4'] = \
        sp.simplify(G_sym - G_cf) == 0
    # squared-form margin identities, in independent symbols dd = D(1-D), ts = t
    dd, ts = sp.symbols('dd t_s')
    Fs = 1 - 2 * dd * ts
    Bs = 1 - sp.Rational(11, 16) * dd * ts
    ids['2B-F = 1+(5/8)dd t'] = sp.simplify(
        2 * Bs - Fs - (1 + sp.Rational(5, 8) * dd * ts)) == 0
    ids['B-F = (21/16)dd t'] = sp.simplify(
        Bs - Fs - sp.Rational(21, 16) * dd * ts) == 0
    ids['(2B-F)^2-F^2 = 4B(B-F)'] = sp.expand(
        (2 * Bs - Fs) ** 2 - Fs ** 2 - 4 * Bs * (Bs - Fs)) == 0
    ok = all(ids.values())
    for k, v in ids.items():
        print(f"    {k}: {'OK' if v else 'FAIL'}")
    print(f"  [final-bound reduction: g <= B <=> G <= 4B(B-F) <=> Psi_t <= (21/32)B;")
    print(f"   (21/32)B >= (21/32)^2 = 441/1024 = {441/1024:.4f} > 5/16 = "
          f"{5/16:.4f} >= Psi_t on Gray]")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# ============================ V11b exact rational chain ============================

def check_V11b():
    from fractions import Fraction as Fr
    print("-" * 84)
    print("V11b: exact rational verification of the full chain at rational Gray")
    print("      points with rational D_c (p=(1-r^2)/2, D_c=(1-r)^2/(2(1+r^2))),")
    print("      down to p ~ 1e-4 (slack unresolvable in double precision)")
    ok = True
    five16 = Fr(5, 16)
    rs = [Fr(1, 10), Fr(1, 2), Fr(9, 10), Fr(99, 100),
          Fr(999, 1000), Fr(9999, 10000)]
    npts = 0
    p_min = 1.0
    for r in rs:
        p = (1 - r * r) / 2
        p_min = min(p_min, float(p))
        Dc_r = (1 - r) ** 2 / (2 * (1 + r * r))
        # sanity: rational D_c agrees with the float route (use the
        # cancellation-free Dc_stable: the naive Dc() loses ~7 digits at
        # p ~ 1e-4 by subtractive cancellation; Dc_stable == Dc is
        # sympy-verified in V7(ii))
        ok &= (abs(float(Dc_r) - Dc_stable(float(p)))
               < 1e-12 * max(float(Dc_r), 1e-300))
        k2 = (1 - 2 * p) ** 2 / (p * p * (1 - p) * (1 - p))
        # exact closed form of Psi at D_c and the 5/16 bound
        ddc = Dc_r * (1 - Dc_r)
        c2c = (1 - 2 * Dc_r) ** 2
        psi_dc = k2 * ddc * (c2c + 4 * ddc * ddc) / c2c ** 2
        psi_cf = (1 - 2 * p) / (4 * (1 - p) ** 2) + p ** 4 / (16 * (1 - p) ** 4)
        ok &= (psi_dc == psi_cf) and (psi_dc < five16)
        for D in (Dc_r, Dc_r / 2, Dc_r / 10):
            dd = D * (1 - D)
            c2 = (1 - 2 * D) ** 2
            psi2 = k2 * dd * (c2 + 4 * dd * dd) / c2 ** 2
            ok &= psi2 <= psi_dc                      # Psi monotone consequence
            for t in (Fr(1, 1000), Fr(1, 2), Fr(1), Fr(2)):
                psi_t = k2 * dd * (c2 + 2 * dd * dd * t) / c2 ** 2
                F = 1 - 2 * dd * t
                G = 8 * k2 * dd * dd * t * (c2 + 2 * dd * dd * t) / c2 ** 2
                B = 1 - Fr(11, 16) * dd * t
                ok &= (psi_t <= psi2) and (psi2 < five16)
                ok &= F > 0 and F < 1
                ok &= (2 * B - F) > 0
                ok &= (2 * B - F) ** 2 - F * F - G > 0   # strict squared margin
                npts += 1
    print(f"  {npts} exact rational points: p down to {p_min:.3e}, D down to "
          f"D_c/10, t in [1e-3, 2] (t=2 is s=pi)")
    print("  all exact: Psi(D_c)==closed form < 5/16; Psi_t <= Psi(D) <= Psi(D_c);")
    print(f"  0<F<1; 2B-F>0; (2B-F)^2-F^2-G > 0 strictly   -> "
          f"{'PASS' if ok else 'FAIL'}")
    return ok

# ============================ main ============================

if __name__ == "__main__":
    print("=" * 84)
    print("7.34b Hypothesis (G): analytic proof, machine-validated step by step")
    print("=" * 84)
    v0 = check_V0()
    v1 = check_V1()
    v2 = check_V2()
    v3 = check_V3()
    v4 = check_V4()
    v5 = check_V5()
    v6 = check_V6()
    v7 = check_V7()
    v8 = check_V8()
    v9 = check_V9()
    v10 = check_V10()
    v11a = check_V11a()
    v11b = check_V11b()
    print("=" * 84)
    allv = [v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11a, v11b]
    names = ["V0", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10",
             "V11a", "V11b"]
    print(" | ".join(f"{nm} {'PASS' if v else 'FAIL'}"
                     for nm, v in zip(names, allv)))
    verdict = "proven-full-range" if all(allv) else "structure-partial"
    print(f"VERDICT: {verdict}")
    print("  g(s) <= 1 - (11/16) D(1-D)(1-cos s) for all p in (0,1/2), "
          "0 < D <= D_c(p), s in (0,pi].")

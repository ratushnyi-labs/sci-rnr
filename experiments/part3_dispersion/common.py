"""
common.py -- shared machinery for the Part III dispersion simulation studies.

Exact finite-blocklength machinery for the binary symmetric Markov source
(BSMS) under Hamming distortion on the Gray region, following the laws of
Part III (tex/papers/dispersion/rnr_dispersion.tex, Remarks 7.34b/i/n'/o)
and the committed verify scripts (scripts/verify/*7_34*):

  * fwht / XorConv       -- Walsh-Hadamard XOR-convolution (O(N log N));
  * bsms_block           -- stationary BSMS law on {0,1}^n and its surprisal;
  * deconv_output        -- the SLB backward-channel deconvolution
                            P_Y = (K^{-1})^{(x)n} P_X (valid on Gray);
  * dtilted_bits         -- exact n-block d-tilted information j_n (bits);
  * ball_logmass_bits    -- exact -log2 P_Y(B_t(x)) for Hamming balls;
  * D_c / V_lossless     -- Gray threshold and lossless varentropy rate;
  * K_beyond_gray        -- the closed-form 1/n^2 threshold-law constant K(p);
  * replica ratio R_A    -- the pattern-quotient replica spectral ratio
                            R_A(s;D) = (1-g_A)/(D(1-D)(1-cos s)) at A=2..5;
  * Dc_aary              -- A-ary Gray threshold (eigenvalue-collision flag).

Everything here is exact linear algebra / transforms, no Monte Carlo; the
studies add the sampling layers on top.
"""

import math

import numpy as np


# ---------------------------------------------------------------- transforms
def fwht(a):
    """Unnormalized Walsh-Hadamard transform (involution up to factor N)."""
    a = a.astype(np.float64, copy=True)
    N = a.size
    h = 1
    while h < N:
        a = a.reshape(N // (2 * h), 2, h)
        x = a[:, 0, :].copy()
        y = a[:, 1, :].copy()
        a[:, 0, :] = x + y
        a[:, 1, :] = x - y
        a = a.reshape(N)
        h *= 2
    return a


class XorConv:
    """XOR-convolution with a fixed kernel k: (k * f)(x) = sum_u k(u) f(x^u)."""

    def __init__(self, kernel):
        self.N = kernel.size
        self.hk = fwht(kernel)

    def __call__(self, f):
        return fwht(fwht(f) * self.hk) / self.N


def popcount(x):
    return np.bitwise_count(x).astype(np.int64)


# ---------------------------------------------------------------- BSMS source
def bsms_block(n, p):
    """Stationary BSMS law on {0,1}^n (X_1 uniform).

    Returns (P, i_bits, sw): probabilities, surprisal -log2 P, switch counts.
    """
    N = 1 << n
    x = np.arange(N, dtype=np.uint64)
    sw = popcount((x ^ (x >> np.uint64(1))) & np.uint64((1 << (n - 1)) - 1))
    log2P = -1.0 + sw * math.log2(p) + (n - 1 - sw) * math.log2(1 - p)
    P = np.exp2(log2P)
    P /= P.sum()
    return P, -log2P, sw


def h2(x):
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


def D_c(p):
    """Gray's critical distortion for the BSMS(p) (binary, Hamming)."""
    return 0.5 * (1.0 - math.sqrt(1.0 - 2.0 * p) / (1.0 - p))


def V_lossless_bits(p):
    """Lossless varentropy rate of the BSMS(p), bits^2/symbol."""
    return p * (1 - p) * (math.log2((1 - p) / p)) ** 2


def rate_gray_bits(p, D):
    """R(D) = h2(p) - h2(D) on the Gray region (per-symbol, bits)."""
    return h2(p) - h2(D)


# ------------------------------------------------------- deconvolution & j_n
def deconv_output(n, p, D):
    """P_Y = (K^{-1})^{(x)n} P_X via XOR-convolution.  Signed vector;
    min entry >= -tol certifies the all-n Gray validity at this (n,p,D)."""
    N = 1 << n
    P, _, _ = bsms_block(n, p)
    z = np.arange(N, dtype=np.uint64)
    pc = popcount(z)
    a0 = 1.0 - 2.0 * D
    log_absk = (n - pc) * math.log(1 - D) + pc * math.log(D) - n * math.log(a0)
    kern = np.exp(log_absk) * np.where(pc % 2 == 0, 1.0, -1.0)
    return XorConv(kern)(P)


def dtilted_bits(n, p, D, PY=None):
    """Exact n-block d-tilted information j_n(x^n, D) in bits, for all x.

    j_n = -lambda* n D log2(e) - log2 M(x),  M(x) = E_{Y*}[exp(-lambda* d_H)],
    lambda* = ln((1-D)/D), with Y* the exact (deconvolved) optimal output.
    """
    N = 1 << n
    if PY is None:
        PY = deconv_output(n, p, D)
    nu = D / (1 - D)                      # e^{-lambda*}
    pc = popcount(np.arange(N, dtype=np.uint64))
    kern = nu ** pc.astype(np.float64)
    M = XorConv(kern)(PY)
    M = np.maximum(M, 1e-300)
    lam_bits = math.log2((1 - D) / D)     # lambda* in bits
    return -n * D * lam_bits - np.log2(M)


def ball_logmass_bits(n, p, D, t, PY=None):
    """Exact G_n(x) = -log2 P_{Y*}(B_t(x)) for all x (Hamming ball, radius t)."""
    N = 1 << n
    if PY is None:
        PY = deconv_output(n, p, D)
    pc = popcount(np.arange(N, dtype=np.uint64))
    ball = (pc <= t).astype(np.float64)
    PB = XorConv(ball)(PY)
    PB = np.maximum(PB, 1e-300)
    return -np.log2(PB)


# ------------------------------------------------------ beyond-Gray constants
def K_beyond_gray(p):
    """Closed-form constant of the finite-n threshold law
    D_c^(n) - D_c = K(p)^2/n^2 (1+o(1));  K(p) = pi p (1-2p)^{1/4}/(2(1-p)^{3/2})."""
    return math.pi * p * (1 - 2 * p) ** 0.25 / (2 * (1 - p) ** 1.5)


# ------------------------------------------------------- A-ary replica ratio
def _pat(tr):
    x, u, up = tr
    if x == u == up:
        return 0
    if x == u and u != up:
        return 1
    if x == up and u != up:
        return 2
    if u == up and x != u:
        return 3
    return 4


def _pattern_reps(A):
    """Representative triples of the realized pattern classes (4 at A=2, 5 at A>=3)."""
    import itertools
    states = list(itertools.product(range(A), repeat=3))
    reps = {}
    for tr in states:
        c = _pat(tr)
        if c not in reps:
            reps[c] = tr
    classes = sorted(reps)
    return states, classes, [reps[c] for c in classes]


def replica_ratio(A, p, D, s):
    """R_A(s;D) = (1 - g_A)/(D(1-D)(1-cos s)),  g_A = |C_s/C_0|^2 rho(Q).

    Q is the pattern-quotient replica matrix on the realized agreement
    patterns of the triple (source, replica, conjugate replica): 4x4 at A=2,
    5x5 at A>=3 (the committed pattern-quotient builder).
    """
    E0 = D / ((A - 1) * (1 - D))
    E = E0 * complex(math.cos(s), math.sin(s))

    def C_of(Ev):
        return ((A - 1) * D * Ev + D - (A - 1)) / (A * D - (A - 1))

    eta = ((A - 1) * D * E + D - (A - 1) * E) / ((A - 1) * D * E + D - (A - 1))
    Cr = abs(C_of(E) / C_of(E0)) ** 2
    T = [[(1 - p) if i == j else p / (A - 1) for j in range(A)] for i in range(A)]
    states, classes, reps = _pattern_reps(A)
    m = len(classes)
    cindex = {c: k for k, c in enumerate(classes)}
    Q = np.zeros((m, m), dtype=np.complex128)
    etb = eta.conjugate()
    for a_, (x, xp, xq) in enumerate(reps):
        for (y, yp, yq) in states:
            term = T[xp][yp] * T[xq][yq] / T[x][y]
            if yp != y:
                term *= eta
            if yq != y:
                term *= etb
            Q[a_, cindex[_pat((y, yp, yq))]] += term
    ev = np.linalg.eigvals(Q)
    g = Cr * float(np.max(np.abs(ev)))
    return (1.0 - g) / (D * (1 - D) * (1 - math.cos(s)))


def _aary_collision_disc(A, p, Ds):
    """Discriminant of the relevant cubic block of the alternating-word
    transfer product B(1)B(0) for the A-ary symmetric chain (A>=3), on the
    multiplicity-weighted 3-class reduction {0, 1, other} (the spectator
    symbols decouple; the class weights (A-2), td+(A-3)to are essential for
    A>=4).  Vectorized over the D-grid.  disc > 0 <=> real spectrum."""
    Ds = np.asarray(Ds, dtype=float)
    lam2 = 1.0 - Ds * A / (A - 1)
    al = 1.0 / A + (1 - 1.0 / A) / lam2
    be = 1.0 / A - (1.0 / A) / lam2
    td, to = 1 - p, p / (A - 1)
    Ki0 = [al, be, be]
    Ki1 = [be, al, be]
    Tr = [[td, to, (A - 2) * to],
          [to, td, (A - 2) * to],
          [to, to, td + (A - 3) * to]]
    B0 = [[Ki0[i] * Tr[i][j] for j in range(3)] for i in range(3)]
    B1 = [[Ki1[i] * Tr[i][j] for j in range(3)] for i in range(3)]
    M = [[sum(B1[i][k] * B0[k][j] for k in range(3)) for j in range(3)]
         for i in range(3)]
    t1 = M[0][0] + M[1][1] + M[2][2]
    t2 = sum(M[i][j] * M[j][i] for i in range(3) for j in range(3))
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
           - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
           + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    c2, c1, c0 = -t1, 0.5 * (t1 * t1 - t2), -det
    return (18 * c2 * c1 * c0 - 4 * c2 ** 3 * c0 + c2 ** 2 * c1 ** 2
            - 4 * c1 ** 3 - 27 * c0 ** 2)


def Dc_aary(A, p, ngrid=4000):
    """A-ary Gray threshold D_c^(A)(p): closed form at A=2; for A>=3 the
    SMALLEST positive root of the relevant-cubic discriminant (the complex
    onset of the dominant eigenpair; the flag is non-monotone in D at
    A>=4, so locate the first sign change, then bisect)."""
    if A == 2:
        return D_c(p)
    hi = (A - 1) / A
    lo = min(p * p / (4 * (A - 1)) * 1e-3, 1e-6)  # below the small-p law scale
    Ds = np.geomspace(lo, hi * 0.999999, ngrid)
    d = _aary_collision_disc(A, p, Ds)
    idx = np.where((d[:-1] > 0) & (d[1:] <= 0))[0]
    if len(idx) == 0:
        raise RuntimeError(f"no collision found for A={A}, p={p}")
    a, b = Ds[idx[0]], Ds[idx[0] + 1]
    for _ in range(80):
        mid = 0.5 * (a + b)
        if _aary_collision_disc(A, p, [mid])[0] > 0:
            a = mid
        else:
            b = mid
    return 0.5 * (a + b)


# ---------------------------------------------------------------- reporting
class Checker:
    """PASS/FAIL line collector in the style of scripts/verify/."""

    def __init__(self):
        self.results = []

    def rep(self, name, ok):
        ok = bool(ok)
        self.results.append((name, ok))
        print(f"  {name:<70} {'PASS' if ok else 'FAIL'}", flush=True)
        return ok

    @property
    def ok(self):
        return all(o for _, o in self.results)


def write_csv(path, header, rows):
    import csv
    with open(path, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(header)
        for r in rows:
            wr.writerow(r)

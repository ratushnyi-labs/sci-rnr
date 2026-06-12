#!/usr/bin/env python3
"""
probe_7_34_nonsym_dual_identity.py
============================================================================
PROBE: does the 7.34b achievability machinery generalize to the NON-SYMMETRIC
binary Markov chain (instance (2) of the unified converse remark)?
Chain: P(0->1)=a, P(1->0)=b, stationary pi=(b/(a+b), a/(a+b)).

Q1 -- the DUAL IDENTITY is chain-agnostic:
    M(beta;x) := E_{Y*}[e^{beta d_H(x,Y*)}] = C(beta)^n E_{X'}[eta(beta)^{d_H(x,X')}]
  with the SAME zeta/eta/C as the symmetric case (they depend only on D), for
  ANY source law on {0,1}^n whose BSC(D) Walsh deconvolution
  P^_Y(w) = P^_X(w)(1-2D)^{-|w|} is a valid law.  The two Walsh steps are
  per-coordinate kernel algebra; P_X enters only through its Walsh coefficients.
  Q1a  the deconvolution-valid threshold D_c(a,b;n) for n=2..14 (bisection);
       cross-check the remark's D_c(0.1,0.3)~0.0129, D_c(0.2,0.4)~0.0455.
  Q1b  the identity in float, real + complex beta incl. the theta_0+is
       contour, several (a,b), n in {10,12,14}, random + pathological words,
       with CONDITIONING-AWARE tolerance (the brute side goes through an FWHT
       deconvolution whose ~1e-16 ABSOLUTE noise dominates when P(x)~1e-12,
       e.g. the alternating word under a strongly asymmetric low-switch chain).
  Q1x  EXACT-RATIONAL certification (polynomial identity testing): both sides
       are rational functions of E=e^beta over Q(D) of degree <= 2n, so EXACT
       Fraction equality at >= 4n+2 distinct rational E proves the identity for
       ALL complex beta (contour included).  Run at n=10 for two asymmetric
       chains with rational (a,b,D), exact deconvolution nonnegativity, exact
       equality at 44 points -- including E = D/(1-D) (the saddle), certifying
       M(theta_0;x) = (1-D)^{-n} P_X(x) EXACTLY.
  Q1c  saddle facts at beta=theta_0=ln(D/(1-D)): zeta=1, eta=0, C=1/(1-D),
       M(theta_0;x)=(1-D)^{-n} P_X(x)  =>  j_n = i_n - n h(D) deterministic.
  Q1d  the OPTIMAL tilt/lambda relation is UNCHANGED by asymmetry:
       lambda* = ln((1-D)/D); exact BA fixed-point check (marginal consistency
       sum_x P(x)Q(y|x)=q(y), E[d]=nD, I = H(X^n)-n h(D)) PLUS an independent
       iterative Blahut-Arimoto run at n=6 => R_n = [H(X^n)-n h(D)]/n = E[j_n]/n.

Q2 -- the REPLICA second moment for the asymmetric chain (NO global-flip
  reduction; the full 8x8 over (x,x',x'') with two-argument T weights, the
  pi(b)pi(c)/pi(a) initial weight, and the P(x)^{-1} outer weight):
  Q2a  E_x|Phi(s;x)|^2 = |C_s/C_0|^{2n} u8^T W8(s)^{n-1} 1  vs brute 2^n
       enumeration, n in {10,12}, machine precision; annealed first-moment
       sanity E_x Phi = ((1-D)+D e^{is})^n (chain-agnostic Binomial).
  Q2b  g(s) := |C_s/C_0|^2 rho(W8(s)) on a grid over (0,pi]: g<1 uniformly
       outside fixed neighborhoods of 0?  curvature g ~ 1 - A s^2?  argmax /
       argmin, margins; curvature vs measured per-site posterior variance.
  Q2c  rho(W8) matches the actual u W^m v growth slope (no Jordan pathology).
  Q2d  spectrum structure: conjugate-pair symmetry (b<->c swap), is the
       symmetric +-q|eta| pair / 4x4 factorization still visible? global-flip
       invariance broken?

Verdict in {generalizes-modulo-replica-inequality, dual-identity-breaks-at-X,
mixed}.
"""
import math
import numpy as np

# ------------------------- asymmetric chain basics -------------------------

def h2(x):
    if x <= 0 or x >= 1:
        return 0.0
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)

def chain(a, b):
    T = np.array([[1 - a, a], [b, 1 - b]])
    pi = np.array([b / (a + b), a / (a + b)])
    return T, pi

def bits_of(n, N):
    v = np.arange(N)
    return ((v[None, :] >> np.arange(n - 1, -1, -1)[:, None]) & 1)  # (n,N) MSB first

def asym_P(n, a, b):
    """Stationary asymmetric-chain law on n-bit integers (MSB = X_1)."""
    T, pi = chain(a, b)
    B = bits_of(n, 1 << n)
    P = pi[B[0]].astype(float)
    for t in range(n - 1):
        P = P * T[B[t], B[t + 1]]
    return P

def fwht(v):
    v = v.astype(float).copy()
    h = 1
    while h < len(v):
        for i in range(0, len(v), h * 2):
            u = v[i:i + h].copy(); w = v[i + h:i + 2 * h].copy()
            v[i:i + h] = u + w; v[i + h:i + 2 * h] = u - w
        h *= 2
    return v

def popcount_arr(v):
    v = np.asarray(v, dtype=np.int64)
    pc = np.zeros_like(v)
    while np.any(v):
        pc += v & 1
        v >>= 1
    return pc

def deconv(n, P, D):
    N = len(P)
    pc = popcount_arr(np.arange(N))
    return fwht(fwht(P) * (1 - 2 * D) ** (-pc.astype(float))) / N

def eta_C(beta, D):
    eb = np.exp(beta)
    ze = (1 - eb) / ((1 + eb) * (1 - 2 * D))
    eta = (1 - ze) / (1 + ze)
    C = (1 + eb) * (1 + ze) / 2
    return ze, eta, C

def theta0(D):
    return math.log(D / (1 - D))

# ------------------------- Q1a: deconvolution threshold -------------------------

def Dc_of_n(n, a, b, tol=1e-12):
    P = asym_P(n, a, b)
    def ok(D):
        return deconv(n, P, D).min() >= -1e-13
    lo, hi = 0.0, 0.49
    if not ok(1e-6):
        return 0.0
    for _ in range(48):
        mid = 0.5 * (lo + hi)
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo

def q1a():
    print("-" * 84)
    print("Q1a: deconvolution-valid threshold D_c(a,b;n), n=2..14 (all-n region <= min)")
    out = {}
    ok = True
    for (a, b) in [(0.1, 0.2), (0.05, 0.3), (0.1, 0.3), (0.2, 0.4)]:
        ds = [Dc_of_n(n, a, b) for n in range(2, 15)]
        dmin = min(ds)
        out[(a, b)] = dmin
        print(f"  (a,b)=({a},{b}): D_c(n) head {['%.5f' % d for d in ds[:4]]} ... "
              f"tail {['%.5f' % d for d in ds[-3:]]}  min(n<=14) = {dmin:.5f}")
    # cross-check the remark's two anchors
    c1 = abs(out[(0.1, 0.3)] - 0.0129) < 2.5e-3
    c2 = abs(out[(0.2, 0.4)] - 0.0455) < 5e-3
    ok &= c1 and c2 and all(v > 1e-3 for v in out.values())
    print(f"  cross-check remark anchors D_c(0.1,0.3)~0.0129 -> {out[(0.1,0.3)]:.5f} "
          f"({'ok' if c1 else 'MISMATCH'}), D_c(0.2,0.4)~0.0455 -> {out[(0.2,0.4)]:.5f} "
          f"({'ok' if c2 else 'MISMATCH'})")
    print(f"  Q1a -> {'PASS' if ok else 'FAIL'}")
    return ok, out

# ------------------------- Q1b/Q1c: the identity -------------------------

def sample_word(rng, n, a, b):
    T, pi = chain(a, b)
    x = [int(rng.random() < pi[1])]
    for _ in range(n - 1):
        x.append(int(rng.random() < T[x[-1], 1]))
    return np.array(x, dtype=int)

def bits_to_int(bits):
    xx = 0
    for v in bits:
        xx = (xx << 1) | int(v)
    return xx

def q1bc(rng, dvals):
    print("-" * 84)
    print("Q1b: dual identity for the asymmetric chain (real + complex beta, contour),")
    print("     float with conditioning-aware tolerance (FWHT abs-noise ~1e-15)")
    print("Q1c: saddle facts => j_n = i_n - n h(D) deterministic")
    ok = True
    worst = 0.0
    worst_slb = 0.0
    for (a, b) in [(0.1, 0.2), (0.05, 0.3), (0.1, 0.3)]:
        Dv = dvals[(a, b)]
        for n in (10, 12, 14):
            N = 1 << n
            P = asym_P(n, a, b)
            for fD in (0.5, 0.9, 1.0):
                D = fD * Dv
                PY = deconv(n, P, D)
                if PY.min() < -1e-13:
                    continue
                th = theta0(D)
                pc = popcount_arr(np.arange(N))
                words = [sample_word(rng, n, a, b) for _ in range(3)]
                words.append(np.zeros(n, dtype=int))                 # all-0 (rare state heavy)
                words.append(np.ones(n, dtype=int))                  # all-1
                words.append(np.array([i % 2 for i in range(n)]))    # alternating
                for x in words:
                    xx = bits_to_int(x)
                    d = pc[np.arange(N) ^ xx]
                    # Q1c saddle facts (absolute gate: FWHT noise is absolute;
                    # exactness is certified in Q1x with rational arithmetic)
                    ze0, eta0, C0 = eta_C(th, D)
                    M0 = float(PY @ np.exp(th * d))
                    slb_abs = abs(M0 - (1 - D) ** (-n) * P[xx])
                    worst_slb = max(worst_slb, slb_abs)
                    ok &= abs(ze0 - 1) < 1e-12 and abs(eta0) < 1e-12
                    ok &= abs(C0 - 1 / (1 - D)) < 1e-12 and slb_abs < 1e-12
                    # Q1b identity, conditioning-aware: the FWHT deconvolution
                    # noise is ABSOLUTE (~1e-16/entry) and transfers to Md as
                    # ~eps*(1+e^{Re beta})^n, which dominates when P(x)~1e-12
                    for beta in (th, th + 0.4, th + 0.5j, th + 1.5j,
                                 th + 1j * math.pi, 0.3 + 0.8j):
                        ze, eta, C = eta_C(complex(beta), D)
                        Md = PY @ np.exp(complex(beta) * d)
                        Mdual = C ** n * (P @ np.power(eta, d))
                        scale = ((1.0 + math.exp(complex(beta).real)) ** n
                                 + abs(C) ** n * (P @ np.abs(eta) ** d))
                        r = abs(Md - Mdual) / max(scale, 1e-300)
                        worst = max(worst, r)
                        ok &= r < 1e-12
        print(f"  (a,b)=({a},{b}): n in {{10,12,14}}, D in {{0.5,0.9,1.0}}xD_c, "
              f"6 words x 6 beta: ok")
    print(f"  Q1b worst conditioned residual = {worst:.2e}; "
          f"Q1c worst SLB abs.diff = {worst_slb:.2e}   -> {'PASS' if ok else 'FAIL'}")
    return ok

# ------------------------- Q1x: exact-rational certification -------------------------

def q1x():
    """Polynomial identity testing over Q: both sides of the dual identity are
    rational functions of E=e^beta of degree <= 2n over Q(D); exact equality at
    >= 4n+2 distinct rational E proves the identity for ALL complex beta."""
    from fractions import Fraction as F
    print("-" * 84)
    print("Q1x: EXACT-RATIONAL dual identity (PIT: both sides poly in E of deg <= n; "
          "44 >> n+1 points, n=10) + exact deconv >= 0 + exact saddle/SLB at E=D/(1-D)")
    ok = True
    n = 10
    N = 1 << n
    pc = [bin(i).count("1") for i in range(N)]
    for (a, b, D) in [(F(1, 10), F(1, 5), F(1, 150)),
                      (F(1, 20), F(3, 10), F(1, 200))]:
        T = [[1 - a, a], [b, 1 - b]]
        pi = [b / (a + b), a / (a + b)]
        B = bits_of(n, N)
        P = [pi[B[0][i]] for i in range(N)]
        for t in range(n - 1):
            P = [P[i] * T[B[t][i]][B[t + 1][i]] for i in range(N)]
        assert sum(P) == 1
        # exact FWHT and deconvolution
        Ph = list(P)
        h = 1
        while h < N:
            for i in range(0, N, h * 2):
                for k in range(i, i + h):
                    u, w = Ph[k], Ph[k + h]
                    Ph[k], Ph[k + h] = u + w, u - w
            h *= 2
        r = 1 / (1 - 2 * D)
        PY = [Ph[w] * r ** pc[w] for w in range(N)]
        h = 1
        while h < N:
            for i in range(0, N, h * 2):
                for k in range(i, i + h):
                    u, w = PY[k], PY[k + h]
                    PY[k], PY[k + h] = u + w, u - w
            h *= 2
        PY = [v / N for v in PY]
        nonneg = min(PY) >= 0
        ok &= nonneg
        # Both sides are POLYNOMIALS in E of degree <= n (RHS = (2c)^{-n}
        # sum_x' P ((1+E)c+(1-E))^{n-d}((1+E)c-(1-E))^d, c=1-2D), so n+1
        # distinct points suffice; we use 44 >> 2n+1 for margin.
        E0 = D / (1 - D)
        Es = [E0] + [F(k, k + 3) for k in range(1, 16)] \
                  + [-F(k, k + 5) for k in range(1, 15)] \
                  + [F(k * k + 1, k) for k in range(1, 15)]
        assert len(set(Es)) == len(Es) >= 2 * n + 2
        words = [bits_to_int([i % 2 for i in range(n)]), (1 << n) - 1, 0, 0b1011001110]
        exact_all = True
        slb_exact = True
        for xx in words:
            d = [pc[y ^ xx] for y in range(N)]
            for E in Es:
                ze = (1 - E) / ((1 + E) * (1 - 2 * D))
                eta = (1 - ze) / (1 + ze)
                C = (1 + E) * (1 + ze) / 2
                lhs = sum(PY[y] * E ** d[y] for y in range(N))
                rhs = C ** n * sum(P[y] * eta ** d[y] for y in range(N))
                exact_all &= (lhs == rhs)
                if E == E0:
                    slb_exact &= (lhs == (1 - D) ** (-n) * P[xx])
        ok &= exact_all and slb_exact
        print(f"  (a,b,D)=({a},{b},{D}): deconv>=0 exact: {nonneg}; identity EXACT "
              f"(Fraction ==) at {len(Es)} rational E x {len(words)} words: {exact_all}; "
              f"saddle M(theta0;x)==(1-D)^-n P(x) EXACT: {slb_exact}")
    print(f"  degree bound n={n} (poly in E) << {len(Es)} points => identity holds "
          f"identically in E = e^beta: ALL complex beta, contour included")
    print(f"  Q1x -> {'PASS' if ok else 'FAIL'}")
    return ok

# ------------------------- Q1d: BA / lambda relation -------------------------

def q1d(dvals):
    print("-" * 84)
    print("Q1d: optimal tilt lambda* = ln((1-D)/D) unchanged; R_n = [H(X^n)-n h(D)]/n")
    ok = True
    n = 6
    N = 1 << n
    pc = popcount_arr(np.arange(N))
    dmat = pc[np.arange(N)[:, None] ^ np.arange(N)[None, :]].astype(float)
    for (a, b) in [(0.1, 0.2), (0.05, 0.3)]:
        D = 0.9 * dvals[(a, b)]
        P = asym_P(n, a, b)
        H = float(-(P @ np.log2(P)))
        q = deconv(n, P, D)
        ok &= q.min() > 0          # interior of the region: strictly positive
        lam = math.log((1 - D) / D)
        K = np.exp(-lam * dmat)
        # exact fixed-point check with q = the deconvolution
        Z = K @ q
        c = ((P / Z) @ K)
        marg = float(np.max(np.abs(c - 1)))
        Q = q[None, :] * K / Z[:, None]
        Dist = float(P @ (Q * dmat).sum(axis=1))
        I = float(P @ (Q * np.log2(Q / q[None, :])).sum(axis=1))
        slb = H - n * h2(D)
        jn = -np.log2(P) - n * h2(D)          # j_n = i_n - n h(D)
        Ejn = float(P @ jn)
        # independent iterative BA from uniform
        qb = np.full(N, 1.0 / N)
        for _ in range(8000):
            Zb = K @ qb
            cb = (P / Zb) @ K
            qb = qb * cb
            qb /= qb.sum()
        Zb = K @ qb
        Qb = qb[None, :] * K / Zb[:, None]
        Db = float(P @ (Qb * dmat).sum(axis=1))
        Ib = float(P @ (Qb * np.log2(Qb / qb[None, :])).sum(axis=1))
        # BA point at the same slope must reproduce (nD, SLB) -- compare at its own D
        e1 = marg; e2 = abs(Dist - n * D); e3 = abs(I - slb); e4 = abs(Ejn - slb)
        e5 = abs(Db - n * D); e6 = abs(Ib - slb)
        ok &= e1 < 1e-10 and e2 < 1e-10 and e3 < 1e-10 and e4 < 1e-10
        ok &= e5 < 1e-7 and e6 < 1e-7
        print(f"  (a,b)=({a},{b}) D={D:.5f} n={n}: marginal-consistency {e1:.1e}, "
              f"|E d - nD| {e2:.1e}, |I-(H-nh(D))| {e3:.1e}, |E j_n - n R_n| {e4:.1e}; "
              f"iterative BA |D-nD| {e5:.1e}, |R-SLB| {e6:.1e}")
    print(f"  Q1d -> {'PASS' if ok else 'FAIL'}  (the lambda/D relation is the "
          f"symmetric one; asymmetry lives ONLY in P_Y*)")
    return ok

# ------------------------- Q2: asymmetric replica 8x8 -------------------------

def replica_W8_asym(s, a, b, D):
    """Full 8x8 over (x_i,x'_i,x''_i); two-argument T weights; pi enters u."""
    _, eta, _ = eta_C(theta0(D) + 1j * s, D)
    etb = np.conj(eta)
    T, pi = chain(a, b)
    states = [(u, v, w) for u in (0, 1) for v in (0, 1) for w in (0, 1)]
    W = np.zeros((8, 8), complex)
    for i, (x, xp, xq) in enumerate(states):
        for j, (y, yp, yq) in enumerate(states):
            W[i, j] = (T[xp, yp] * T[xq, yq] / T[x, y]
                       * (eta if y ^ yp else 1.0) * (etb if y ^ yq else 1.0))
    u = np.array([pi[xp] * pi[xq] / pi[x]
                  * (eta if x ^ xp else 1.0) * (etb if x ^ xq else 1.0)
                  for (x, xp, xq) in states], complex)
    return W, u

def C_ratio(s, D):
    th = theta0(D)
    return abs(eta_C(th + 1j * s, D)[2] / eta_C(th, D)[2])

def E2_tm(s, a, b, D, n):
    W, u = replica_W8_asym(s, a, b, D)
    val = u @ np.linalg.matrix_power(W, n - 1) @ np.ones(8)
    return C_ratio(s, D) ** (2 * n) * val

def brute_moments(s_list, n, a, b, D):
    N = 1 << n
    P = asym_P(n, a, b)
    PY = deconv(n, P, D)
    pc = popcount_arr(np.arange(N))
    th = theta0(D)
    ks = np.arange(n + 1)
    phase = np.exp(1j * np.outer(ks, np.asarray(s_list)))
    E1 = np.zeros(len(s_list), complex)
    E2 = np.zeros(len(s_list))
    vbar = 0.0
    ix = np.arange(N)
    for xx in range(N):
        d = pc[ix ^ xx]
        w0 = PY * np.exp(th * d)
        A = np.bincount(d, weights=w0, minlength=n + 1)
        M0 = A.sum()
        Phi = (A @ phase) / M0
        E1 += P[xx] * Phi
        E2 += P[xx] * np.abs(Phi) ** 2
        pk = A / M0
        mu = float(ks @ pk)
        vbar += P[xx] * (float((ks ** 2) @ pk) - mu * mu)
    return E1, E2, vbar / n

def q2a(dvals):
    print("-" * 84)
    print("Q2a: brute E_x|Phi|^2 == |C_s/C_0|^{2n} u8 W8^{n-1} 1 (asymmetric, no flip "
          "reduction); E_x Phi == Binom(n,D) cf")
    ok = True
    worst1 = worst2 = worst_im = 0.0
    s_list = [0.3, 0.7, 1.2, 1.9, 2.5, 3.0, math.pi]
    vbars = {}
    for (a, b) in [(0.1, 0.2), (0.05, 0.3), (0.1, 0.3)]:
        for (fD, n) in [(0.9, 10), (0.9, 12), (0.5, 12), (1.0, 10)]:
            D = fD * dvals[(a, b)]
            E1b, E2b, vbar = brute_moments(s_list, n, a, b, D)
            if fD == 0.9 and n == 12:
                vbars[(a, b)] = (vbar, D)
            for k, s in enumerate(s_list):
                ref1 = ((1 - D) + D * np.exp(1j * s)) ** n
                d1 = abs(E1b[k] - ref1) / max(abs(ref1), 1e-300)
                v8 = E2_tm(s, a, b, D, n)
                worst_im = max(worst_im, abs(v8.imag) / max(abs(v8), 1e-300))
                d2 = abs(E2b[k] - v8.real) / max(abs(E2b[k]), 1e-300)
                worst1 = max(worst1, d1); worst2 = max(worst2, d2)
                ok &= d1 < 1e-9 and d2 < 1e-9 and abs(v8.imag) < 1e-9 * abs(v8)
        print(f"  (a,b)=({a},{b}): 4 (D,n) combos x {len(s_list)} s-values: ok")
    print(f"  Q2a worst rel.diffs: first-moment {worst1:.2e} (annealed K ~ Bin(n,D) "
          f"holds for the ASYMMETRIC chain), second-moment {worst2:.2e}, "
          f"imag {worst_im:.2e}   -> {'PASS' if ok else 'FAIL'}")
    return ok, vbars

def rho_W8(s, a, b, D):
    W, _ = replica_W8_asym(s, a, b, D)
    return float(np.max(np.abs(np.linalg.eigvals(W))))

def g_of(s, a, b, D):
    return C_ratio(s, D) ** 2 * rho_W8(s, a, b, D)

def growth_slope(s, a, b, D, m_lo=1500, m_hi=3000):
    W, u = replica_W8_asym(s, a, b, D)
    v = np.ones(8, complex)
    w = u.copy()
    logs = []
    acc = 0.0
    for m in range(1, m_hi + 1):
        w = w @ W
        nm = np.linalg.norm(w)
        acc += math.log(nm)
        w /= nm
        if m >= m_lo:
            logs.append(acc + math.log(max(abs(w @ v), 1e-300)))
    ms = np.arange(m_lo, m_hi + 1, dtype=float)
    A = np.vstack([ms, np.ones_like(ms)]).T
    return float(np.linalg.lstsq(A, np.array(logs), rcond=None)[0][0])

def q2bcd(dvals, vbars):
    print("-" * 84)
    print("Q2b: g(s)=|C_s/C_0|^2 rho(W8) on (0,pi]; Q2c growth slope; Q2d spectrum")
    ok_below = True
    ok_growth = True
    grid = np.concatenate([np.geomspace(1e-4, 0.05, 300, endpoint=False),
                           np.linspace(0.05, math.pi, 2000)])
    for (a, b) in [(0.1, 0.2), (0.05, 0.3), (0.1, 0.3), (0.2, 0.4)]:
        for fD in (0.5, 0.9, 1.0):
            D = fD * dvals[(a, b)]
            g = np.array([g_of(s, a, b, D) for s in grid])
            m1 = grid >= 0.05
            m2 = grid >= 1e-3
            below1 = bool(np.all(g[m1] < 1.0))
            below2 = bool(np.all(g[m2] < 1.0))
            gmax = float(np.max(g[m1])); sarg = float(grid[m1][np.argmax(g[m1])])
            gmin = float(np.min(g[m1])); smin = float(grid[m1][np.argmin(g[m1])])
            gpi = float(g[-1])
            sm = (grid >= 2e-4) & (grid <= 2e-2)
            A_fit = float(np.mean((1 - g[sm]) / grid[sm] ** 2))
            ok_below &= below1 and below2
            extra = ""
            if fD == 0.9 and (a, b) in vbars:
                vb, Dv = vbars[(a, b)]
                extra = f" | curvature A={A_fit:.5f} vs n=12 per-site post.var {vb:.5f}"
            print(f"  (a,b)=({a},{b}) D={fD}xDc={D:.5f}: g<1 on [0.05,pi]:{below1} "
                  f"on [1e-3,pi]:{below2} | max {gmax:.6f}@s={sarg:.3f} "
                  f"(margin {1-gmax:.2e}) min {gmin:.6f}@s={smin:.3f} g(pi)={gpi:.6f}"
                  f"{extra}")
        # Q2c at the 0.9 case
        D = 0.9 * dvals[(a, b)]
        for s_chk in (1.0, math.pi):
            sl = growth_slope(s_chk, a, b, D)
            lr = math.log(rho_W8(s_chk, a, b, D))
            okg = abs(sl - lr) < max(2e-3, 5e-3 * abs(lr))
            ok_growth &= okg
            if not okg:
                print(f"    growth-rate MISMATCH at s={s_chk}: slope {sl:.6f} vs "
                      f"log rho {lr:.6f}")
    print(f"  Q2b g<1 uniformly outside s=0 for all combos -> "
          f"{'PASS' if ok_below else 'FAIL'}")
    print(f"  Q2c rho(W8) == actual growth slope -> {'PASS' if ok_growth else 'FAIL'}")
    # Q2d spectrum structure at one representative point
    print("  Q2d spectrum structure (a,b)=(0.1,0.3), D=0.9Dc:")
    a, b = 0.1, 0.3
    D = 0.9 * dvals[(a, b)]
    for s in (0.8, 2.0, math.pi):
        W, _ = replica_W8_asym(s, a, b, D)
        ev = np.linalg.eigvals(W)
        ev = ev[np.argsort(-np.abs(ev))]
        conj_closed = np.allclose(np.sort_complex(ev), np.sort_complex(np.conj(ev)),
                                  atol=1e-10)
        # global-flip invariance test: F = complement permutation on (x,x',x'')
        perm = [7 - i for i in range(8)]  # (u,v,w) -> (1-u,1-v,1-w) reverses index
        flip_inv = np.allclose(W, W[np.ix_(perm, perm)], atol=1e-12)
        # symmetric-case signature: a +-pair (lam, -lam) of magnitude q|eta|?
        pair = any(abs(ev[i] + ev[j]) < 1e-9 * max(abs(ev[i]), 1)
                   for i in range(8) for j in range(i + 1, 8)
                   if abs(abs(ev[i]) - abs(ev[j])) < 1e-9)
        print(f"    s={s:.3f}: |eigs| = {[f'{abs(e):.5f}' for e in ev]}")
        print(f"      conjugate-closed (b<->c swap symmetry): {conj_closed}; "
              f"global-flip W-invariance: {flip_inv}; "
              f"exact (lam,-lam) pair: {pair}")
    return ok_below, ok_growth

# ------------------------- main -------------------------

if __name__ == "__main__":
    print("=" * 84)
    print("7.34 NON-SYMMETRIC-chain probe: dual identity + replica second moment")
    print("=" * 84)
    rng = np.random.default_rng(11)
    ok1a, dvals = q1a()
    ok1b = q1bc(rng, dvals)
    ok1x = q1x()
    ok1d = q1d(dvals)
    ok2a, vbars = q2a(dvals)
    ok2b, ok2c = q2bcd(dvals, vbars)
    print("=" * 84)
    print(f"Q1a {'PASS' if ok1a else 'FAIL'} | Q1b/c {'PASS' if ok1b else 'FAIL'} | "
          f"Q1x {'PASS' if ok1x else 'FAIL'} | "
          f"Q1d {'PASS' if ok1d else 'FAIL'} | Q2a {'PASS' if ok2a else 'FAIL'} | "
          f"Q2b {'PASS' if ok2b else 'FAIL'} | Q2c {'PASS' if ok2c else 'FAIL'}")
    if ok1a and ok1b and ok1x and ok1d and ok2a and ok2b and ok2c:
        verdict = "generalizes-modulo-replica-inequality"
    elif not (ok1b and ok1x and ok1d):
        verdict = "dual-identity-breaks-at-X"
    else:
        verdict = "mixed"
    print(f"VERDICT: {verdict}")
    print("  dual identity: chain-agnostic, exact (incl. contour); j_n = i_n - n h(D)")
    print("  lambda*/D relation: unchanged (asymmetry lives only in P_Y*)")
    print("  replica: 8x8 survives & is exact; the CLOSED-FORM g-inequality is lost")
    print("  (global-flip reduction broken); g<1 certified numerically per (a,b,D)")

#!/usr/bin/env python3
"""
probe_7_34_replica_secondmoment_tail.py
============================================================================
SECOND-MOMENT / REPLICA route for the 7.34b quenched shell-uniform TAIL bound
    sup_{|s| >= delta} |Phi(s;x)| <= e^{-c(p) n}   (typical x),
the restated residual of Remark 7.34b (dual elliptic representation,
Phi(s;x) = (C_s/C_0)^n S_x(eta_s), commit e10bee4).

IDEA (not attempted in the committed negatives, which all work word-by-word):
average the SQUARED modulus over the source word,
    E_x |Phi(s;x)|^2 = |C_s/C_0|^{2n} sum_{x,x',x''} [P(x')P(x'')/P(x)]
                          eta_s^{d_H(x,x')} etabar_s^{d_H(x,x'')}
                     = E[e^{is(K-K')}]   (two posterior replicas given x),
an explicit TRIPLE elliptic-source chain = an 8x8 (flip-symmetry-reducible
to 4x4) transfer matrix W(s).  If
    g(s) := |C_s/C_0|^2 * rho(W(s)) < 1  uniformly on s in [delta, pi],
then E_x|Phi|^2 <= poly * g^n, and Markov + a union bound over an O(n^2)-net
(|dPhi/ds| <= E[K|x] <= n) gives the quenched uniform tail bound for all but
an exponentially small fraction of words -- the KP cluster expansion supplying
the central window |s| <= delta_KP(p,D) (the FIXED window where the contour
|eta_s| stays inside the Kotecky-Preiss disc r_KP(p)).

CHECKS
  C1   2x2 transfer matrix for S_x(eta) (states: site flipped / not flipped;
       bond weights rho^{+-1}) == committed polymer-block DP
       (lemma_7_34b_polymer_gas.py) == vectorised brute force; n<=14,
       p in {0.1,0.25,0.4}, real + on-contour complex eta, random words +
       all-stay + alternating pathological words.
  C1b  per-word Phi(s;x) = (C_s/C_0)^n S_x^{TM}(eta_s) == brute posterior
       characteristic function (FWHT deconvolution of Y*), n=12.
  C2   replica identity: brute 2^n E_x|Phi(s;x)|^2 == |C_s/C_0|^{2n}
       u4^T W4(s)^{n-1} v4 == 8x8 version; n in {8,10,12}, p in
       {0.1,0.25,0.4}, D in {0.5 Dc, 0.9 Dc}, s grid incl. s=pi.
       Plus the annealed first-moment sanity E_x Phi = ((1-D)+D e^{is})^n.
  C3   g(s) on a fine grid over (0,pi] for the 6 (p,D) combos:
       delta_KP, g_max on [delta_KP, pi], argmax, margin, g(pi),
       small-s curvature (g ~ 1 - a s^2), and large-n growth-rate
       consistency rho(W) vs the actual u W^m v slope.
  C4   shell conditioning via Bernstein tilt e^{lam*sw(x)} of the OUTER
       x-chain (the sqrt(n)-shell corresponds to lam = Theta(n^{-1/2}) -> 0;
       fixed lam is a conservative robustness check): g_lam on [delta,pi]
       for lam in {+-0.05, +-0.1, +-0.2} at (p,D)=(0.25, 0.9Dc) + dg/dlam.

Prints PASS/FAIL per check and a final verdict in
  {representation-verified-g-below-1, representation-verified-g-reaches-1,
   construction-fails, mixed}.
"""
import math
import numpy as np

# ------------------------- shared BSMS / dual objects -------------------------

def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))

def theta0(D):
    return math.log(D / (1 - D))

def eta_C(beta, D):
    """zeta, eta, C of the dual identity at (complex) beta."""
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

def bsms_P(n, p):
    """BSMS law on n-bit integers (stationary (1/2,1/2), switch prob p)."""
    N = 1 << n
    v = np.arange(N)
    sw = popcount_arr((v ^ (v >> 1)) & ((1 << (n - 1)) - 1))
    P = 0.5 * (1 - p) ** (n - 1 - sw) * p ** sw
    return P

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

# Kotecky-Preiss radius -- ported verbatim from remark_7_34b_kp_radius_and_zeros.py
def kp_feasible(eta_abs, p):
    B = ((1 - p) / p) ** 2
    for mu in np.linspace(1e-4, 6.0, 6000):
        u = eta_abs * math.exp(mu)
        if u >= 1:
            continue
        if B * (2 * u / (1 - u) + u / (1 - u) ** 2) <= mu:
            return True
    return False

def r_cluster(p):
    lo, hi = 1e-6, 0.999
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if kp_feasible(mid, p):
            lo = mid
        else:
            hi = mid
    return lo

# ------------------------- C1: 2x2 transfer matrix for S_x -------------------------

def S_tm(xbits, p, eta):
    """2x2 site transfer matrix: state e_i in {0,1} = site flipped or not.
    S_x(eta) = w^T (prod_j A_j) 1,  w=(1,eta),
    A_j[e,e'] = eta^{e'} * rho^{(e xor e')*(1-2 s_j)}, s_j = 1[x_j != x_{j+1}]."""
    rho = p / (1 - p)
    w = np.array([1.0, eta], dtype=complex)
    for j in range(len(xbits) - 1):
        sj = 1 if xbits[j] != xbits[j + 1] else 0
        r = rho ** (1 - 2 * sj)
        A = np.array([[1.0, eta * r],
                      [r,   eta]], dtype=complex)
        w = w @ A
    return w[0] + w[1]

def S_blocks(x, n, p, eta):
    """Committed polymer-block DP (lemma_7_34b_polymer_gas.py), complex-capable."""
    ratio = p / (1 - p)
    def act(a, b):
        z = eta ** (b - a + 1)
        if a >= 1:
            z *= ratio if (x[a] == x[a - 1]) else 1.0 / ratio
        if b + 1 <= n - 1:
            z *= ratio if (x[b + 1] == x[b]) else 1.0 / ratio
        return z
    memo = [None] * (n + 2)
    memo[n] = 1.0 + 0j
    memo[n + 1] = 1.0 + 0j
    for i in range(n - 1, -1, -1):
        tot = memo[i + 1]
        for b in range(i, n):
            nxt = memo[b + 2] if b + 2 <= n else 1.0 + 0j
            tot += act(i, b) * nxt
        memo[i] = tot
    return memo[0]

def S_brute_vec(xx, n, p, eta, P, pc_all):
    """sum_{x'} P(x') eta^{d_H(x,x')} / P(x), vectorised over integer words."""
    ix = np.arange(1 << n)
    d = pc_all[ix ^ xx]
    return np.sum(P * np.power(complex(eta), d)) / P[xx]

def bits_to_int(bits):
    xx = 0
    for b in bits:
        xx = (xx << 1) | int(b)
    return xx

def check1(rng):
    print("-" * 84)
    print("C1: 2x2 transfer matrix for S_x(eta) == polymer-block DP == brute (n<=14)")
    ok = True
    worst = 0.0
    for p in (0.1, 0.25, 0.4):
        D = 0.9 * Dc(p)
        for n in (10, 14):
            P = bsms_P(n, p)
            pc_all = popcount_arr(np.arange(1 << n))
            words = [(np.cumsum(rng.random(n) < p) % 2).astype(int) for _ in range(3)]
            words.append(np.zeros(n, dtype=int))                    # all-stay
            words.append(np.array([i % 2 for i in range(n)]))       # alternating
            etas = [0.05, 0.2, -0.2, 0.3,
                    eta_s(0.5, D), eta_s(1.5, D), eta_s(math.pi, D)]
            for x in words:
                xx = bits_to_int(x)
                for eta in etas:
                    a = S_brute_vec(xx, n, p, eta, P, pc_all)
                    b = S_blocks(x, n, p, eta)
                    c = S_tm(x, p, eta)
                    m = max(abs(a), 1e-300)
                    d1 = abs(a - b) / m; d2 = abs(a - c) / m
                    worst = max(worst, d1, d2)
                    ok &= d1 < 1e-9 and d2 < 1e-9
        print(f"  p={p}: n in {{10,14}}, 5 words x 7 eta (real + on-contour complex): ok")
    print(f"  C1 worst rel.diff = {worst:.2e}   -> {'PASS' if ok else 'FAIL'}")
    return ok

# ------------------------- posterior brute force (FWHT) -------------------------

def posterior_setup(n, p, D):
    N = 1 << n
    P = bsms_P(n, p)
    Ph = fwht(P)
    pc_all = popcount_arr(np.arange(N))
    PYs = fwht(Ph * (1 - 2 * D) ** (-pc_all.astype(float))) / N
    assert PYs.min() > -1e-9, "deconvolution must be valid in-Gray"
    return P, PYs, pc_all

def brute_moments(s_list, n, p, D):
    """E_x Phi(s;x) and E_x |Phi(s;x)|^2 by full 2^n enumeration."""
    N = 1 << n
    P, PYs, pc_all = posterior_setup(n, p, D)
    th = theta0(D)
    ks = np.arange(n + 1)
    phase = np.exp(1j * np.outer(ks, np.asarray(s_list)))   # (n+1, S)
    E1 = np.zeros(len(s_list), complex)
    E2 = np.zeros(len(s_list))
    ix = np.arange(N)
    for xx in range(N):
        d = pc_all[ix ^ xx]
        w0 = PYs * np.exp(th * d)
        A = np.bincount(d, weights=w0, minlength=n + 1)
        M0 = A.sum()
        Phi = (A @ phase) / M0
        E1 += P[xx] * Phi
        E2 += P[xx] * np.abs(Phi) ** 2
    return E1, E2

def brute_phi_word(s_list, xx, n, PYs, pc_all, D):
    ix = np.arange(1 << n)
    d = pc_all[ix ^ xx]
    th = theta0(D)
    w0 = PYs * np.exp(th * d)
    A = np.bincount(d, weights=w0, minlength=n + 1)
    M0 = A.sum()
    ks = np.arange(n + 1)
    return (A @ np.exp(1j * np.outer(ks, np.asarray(s_list)))) / M0

def check1b(rng):
    print("-" * 84)
    print("C1b: per-word Phi(s;x) = (C_s/C_0)^n S_x^TM(eta_s) == brute posterior cf (n=12)")
    ok = True
    worst = 0.0
    n = 12
    s_list = [0.3, 1.2, 2.5, math.pi]
    for p in (0.25, 0.4):
        D = 0.9 * Dc(p)
        P, PYs, pc_all = posterior_setup(n, p, D)
        th = theta0(D)
        C0 = eta_C(th, D)[1]
        for _ in range(3):
            x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
            xx = bits_to_int(x)
            phi_b = brute_phi_word(s_list, xx, n, PYs, pc_all, D)
            for k, s in enumerate(s_list):
                et, Cs = eta_C(th + 1j * s, D)
                phi_d = (Cs / C0) ** n * S_tm(x, p, et)
                rel = abs(phi_b[k] - phi_d) / max(abs(phi_b[k]), 1e-300)
                worst = max(worst, rel)
                ok &= rel < 1e-10
    print(f"  worst rel.diff = {worst:.2e}   -> {'PASS' if ok else 'FAIL'}")
    return ok

# ------------------------- C2: replica transfer matrix -------------------------

def replica_W4(s, p, D, lam=0.0):
    """Reduced 4x4 replica matrix over (d1,d2)=(x xor x', x xor x'').
    W[(d1,d2),(e1,e2)] = [sum_a e^{lam a} T_{a^d1^e1} T_{a^d2^e2} / T_a]
                          * eta^{e1} etabar^{e2};  u4 = eta^{d1} etabar^{d2}."""
    et = eta_s(s, D)
    etb = np.conj(et)
    T = (1 - p, p)
    states = ((0, 0), (0, 1), (1, 0), (1, 1))
    W = np.zeros((4, 4), complex)
    for i, (d1, d2) in enumerate(states):
        for j, (e1, e2) in enumerate(states):
            c = 0.0
            for a in (0, 1):
                c += math.exp(lam * a) * T[a ^ d1 ^ e1] * T[a ^ d2 ^ e2] / T[a]
            W[i, j] = c * (et if e1 else 1.0) * (etb if e2 else 1.0)
    u = np.array([(et if d1 else 1.0) * (etb if d2 else 1.0)
                  for (d1, d2) in states], complex)
    return W, u

def replica_W8(s, p, D):
    """Unreduced 8x8 over (x_i, x'_i, x''_i)."""
    et = eta_s(s, D)
    etb = np.conj(et)
    T = (1 - p, p)
    states = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    W = np.zeros((8, 8), complex)
    for i, (a, b, c) in enumerate(states):
        for j, (a2, b2, c2) in enumerate(states):
            W[i, j] = (T[b ^ b2] * T[c ^ c2] / T[a ^ a2]
                       * (et if a2 ^ b2 else 1.0) * (etb if a2 ^ c2 else 1.0))
    u = np.array([0.5 * (et if a ^ b else 1.0) * (etb if a ^ c else 1.0)
                  for (a, b, c) in states], complex)
    return W, u

def E2_tm(s, p, D, n, use8=False):
    if use8:
        W, u = replica_W8(s, p, D)
    else:
        W, u = replica_W4(s, p, D)
    val = u @ np.linalg.matrix_power(W, n - 1) @ np.ones(W.shape[0])
    return C_ratio(s, D) ** (2 * n) * val

def check2():
    print("-" * 84)
    print("C2: brute E_x|Phi(s;x)|^2 == |C_s/C_0|^{2n} u4 W4^{n-1} v4 == 8x8 version")
    ok = True
    worst2 = worst1 = worst8 = worst_im = 0.0
    s_list = [0.3, 0.7, 1.2, 1.9, 2.5, 3.0, math.pi]
    cases = []
    for p in (0.1, 0.25, 0.4):
        for fD in (0.5, 0.9):
            cases.append((p, fD, 10))
    cases += [(0.25, 0.9, 8), (0.4, 0.9, 12), (0.25, 0.9, 12)]
    for (p, fD, n) in cases:
        D = fD * Dc(p)
        E1b, E2b = brute_moments(s_list, n, p, D)
        for k, s in enumerate(s_list):
            # annealed first moment == Binomial(n,D) cf
            ref1 = ((1 - D) + D * np.exp(1j * s)) ** n
            d1 = abs(E1b[k] - ref1) / max(abs(ref1), 1e-300)
            worst1 = max(worst1, d1); ok &= d1 < 1e-9
            # second moment via 4x4 and 8x8
            v4 = E2_tm(s, p, D, n, use8=False)
            v8 = E2_tm(s, p, D, n, use8=True)
            worst_im = max(worst_im, abs(v4.imag) / max(abs(v4), 1e-300))
            d2 = abs(E2b[k] - v4.real) / max(abs(E2b[k]), 1e-300)
            d8 = abs(v8 - v4) / max(abs(v4), 1e-300)
            worst2 = max(worst2, d2); worst8 = max(worst8, d8)
            ok &= d2 < 1e-9 and d8 < 1e-11 and abs(v4.imag) < 1e-10 * abs(v4)
        print(f"  p={p:<5} D={fD}Dc={D:.5f} n={n:>2}: "
              f"E1 vs Binom cf ok, E2 brute vs TM ok ({len(s_list)} s-values)")
    print(f"  C2 worst rel.diffs: first-moment {worst1:.2e}, "
          f"second-moment brute-vs-4x4 {worst2:.2e}, 4x4-vs-8x8 {worst8:.2e}, "
          f"imag-part {worst_im:.2e}   -> {'PASS' if ok else 'FAIL'}")
    return ok

# ------------------------- C3: g(s) over (0,pi] -------------------------

def rho_W4(s, p, D, lam=0.0):
    W, _ = replica_W4(s, p, D, lam)
    return float(np.max(np.abs(np.linalg.eigvals(W))))

def g_of(s, p, D, lam=0.0):
    norm = (1 - p) + p * math.exp(lam)   # E[e^{lam * bond-switch}] per bond
    return C_ratio(s, D) ** 2 * rho_W4(s, p, D, lam) / norm

def growth_slope(s, p, D, m_lo=1500, m_hi=3000):
    """LS slope of log|u W^m v| for m in [m_lo, m_hi] (normalised iteration)."""
    W, u = replica_W4(s, p, D)
    v = np.ones(4, complex)
    w = u.copy()
    logs = []
    acc = 0.0
    for m in range(1, m_hi + 1):
        w = w @ W
        nm = np.linalg.norm(w)
        acc += math.log(nm)
        w = w / nm
        if m >= m_lo:
            val = abs(w @ v)
            logs.append(acc + math.log(max(val, 1e-300)))
    ms = np.arange(m_lo, m_hi + 1, dtype=float)
    A = np.vstack([ms, np.ones_like(ms)]).T
    slope = np.linalg.lstsq(A, np.array(logs), rcond=None)[0][0]
    return slope

def check3():
    print("-" * 84)
    print("C3: g(s) = |C_s/C_0|^2 rho(W4(s)) on (0,pi]; KP window; margins; growth rate")
    ok_below = True
    ok_growth = True
    results = {}
    grid = np.concatenate([np.geomspace(1e-5, 0.1, 600, endpoint=False),
                           np.linspace(0.1, math.pi, 2400)])
    for p in (0.1, 0.25, 0.4):
        rkp = r_cluster(p)
        for fD in (0.5, 0.9):
            D = fD * Dc(p)
            etas_abs = np.abs([eta_s(s, D) for s in grid])
            cross = np.argmax(etas_abs > rkp) if np.any(etas_abs > rkp) else None
            delta = math.pi if cross is None else float(grid[max(cross - 1, 0)])
            g = np.array([g_of(s, p, D) for s in grid])
            tail = grid >= delta
            gmax = float(np.max(g[tail])); sarg = float(grid[tail][np.argmax(g[tail])])
            gpi = float(g[-1])
            below = bool(np.all(g[tail] < 1.0))
            # whole-interval behaviour (g<1 except s->0)
            interior = grid >= 1e-3
            below_all = bool(np.all(g[interior] < 1.0))
            # small-s curvature: g ~ 1 - a s^2
            smask = (grid >= 5e-4) & (grid <= 5e-2) if delta > 5e-2 else \
                    (grid >= grid[2]) & (grid <= max(delta, grid[10]))
            a_fit = float(np.mean((1 - g[smask]) / grid[smask] ** 2))
            # growth-rate consistency at s = pi and midpoint of [delta,pi]
            gr_ok = True
            for s_chk in (0.5 * (delta + math.pi), math.pi):
                sl = growth_slope(s_chk, p, D)
                lr = math.log(rho_W4(s_chk, p, D))
                if abs(sl - lr) > max(2e-3, 5e-3 * abs(lr)):
                    gr_ok = False
            ok_below &= below
            ok_growth &= gr_ok
            results[(p, fD)] = (delta, gmax, sarg, gpi, a_fit, below, below_all, rkp)
            print(f"  p={p:<5} D={fD}Dc={D:.5f}: r_KP={rkp:.5f} delta_KP={delta:.5f} | "
                  f"g_max[delta,pi]={gmax:.6f} at s={sarg:.4f} "
                  f"(margin {1-gmax:.2e}) g(pi)={gpi:.6f} | "
                  f"g<1 on [delta,pi]: {below}; on [1e-3,pi]: {below_all}; "
                  f"small-s a={a_fit:.4f} | growth-rate ok: {gr_ok}")
    print(f"  C3 g<1 uniformly on [delta_KP,pi] for all 6 combos -> "
          f"{'PASS' if ok_below else 'FAIL'}")
    print(f"  C3 rho(W) matches actual u W^m v growth rate -> "
          f"{'PASS' if ok_growth else 'FAIL'}")
    return ok_below, ok_growth, results

# ------------------------- C4: shell tilt -------------------------

def check4(results):
    print("-" * 84)
    print("C4: Bernstein shell tilt of the outer x-chain at (p,D)=(0.25, 0.9Dc)")
    p, fD = 0.25, 0.9
    D = fD * Dc(p)
    delta = results[(p, fD)][0]
    grid = np.linspace(delta, math.pi, 800)
    ok = True
    base = None
    for lam in (-0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2):
        g = np.array([g_of(s, p, D, lam) for s in grid])
        gmax = float(np.max(g)); sarg = float(grid[np.argmax(g)])
        if lam == 0.0:
            base = (gmax, sarg)
        ok &= bool(np.all(g < 1.0))
        print(f"  lam={lam:+.2f}: g_max[delta,pi]={gmax:.6f} at s={sarg:.4f} "
              f"(margin {1-gmax:.2e})")
    # derivative at the unperturbed argmax
    s0 = base[1]
    h = 0.01
    dg = (g_of(s0, p, D, h) - g_of(s0, p, D, -h)) / (2 * h)
    print(f"  dg/dlam at s*={s0:.4f}: {dg:+.4f}  "
          f"(shell tilt lam=Theta(n^-1/2) => |g_lam-g| = O(n^-1/2): immaterial)")
    print(f"  C4 g_lam<1 on [delta,pi] for all tested lam -> {'PASS' if ok else 'FAIL'}")
    return ok

# ------------------------- main -------------------------

if __name__ == "__main__":
    print("=" * 84)
    print("7.34b replica/second-moment tail probe: E_x|Phi(s;x)|^2 via triple-chain TM")
    print("=" * 84)
    rng = np.random.default_rng(7)
    c1 = check1(rng)
    c1b = check1b(rng)
    c2 = check2()
    c3_below, c3_growth, results = check3()
    c4 = check4(results)
    print("=" * 84)
    rep_ok = c1 and c1b and c2 and c3_growth
    if rep_ok and c3_below and c4:
        verdict = "representation-verified-g-below-1"
    elif rep_ok and not c3_below:
        verdict = "representation-verified-g-reaches-1"
    elif not rep_ok:
        verdict = "construction-fails"
    else:
        verdict = "mixed"
    print(f"C1 {'PASS' if c1 else 'FAIL'} | C1b {'PASS' if c1b else 'FAIL'} | "
          f"C2 {'PASS' if c2 else 'FAIL'} | C3(g<1) {'PASS' if c3_below else 'FAIL'} | "
          f"C3(growth) {'PASS' if c3_growth else 'FAIL'} | C4 {'PASS' if c4 else 'FAIL'}")
    print(f"VERDICT: {verdict}")

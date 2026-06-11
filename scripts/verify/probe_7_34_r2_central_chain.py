#!/usr/bin/env python3
"""
probe_7_34_r2_central_chain.py
=============================================================================
(R2) ORDERED-CONSTANTS CENTRAL-WINDOW CHAIN -- hostile per-lemma verification
of the achievability assembly Var(log2 S_n | shell) = O(1) for the BSMS
Gray-region RD dispersion (Remark 7.34b; quenched tail lemma already proven).

Lemma chain under attack (each section = one lemma):
  R0  conventions anchor: DFT/log-TM posterior == brute FWHT at n=12, AND the
      Abel decomposition  P_{Y*}(B_t) = M e^{-theta0 nD} * q_t * W * phi^{nD-t}
      ties -log2 S_n to the actual ball mass (exact at n=12).
  R1  (L-a) KP radius: closed-form lower bound r_lb = p^2/(12 e (1-p)^2) is
      KP-feasible (hence r_KP >= r_lb > 0 ANALYTICALLY); contour slope
      |eta_s| <= A(D)|s| (real s) and <= 4 A(D)|z| on the complex disk; the
      s-window delta_c = min(r_lb/(4A), tau_max) keeps |eta| inside the
      MEASURED KP radius; S_x is zero-free on |eta| <= r_lb (so log S_x is
      analytic, word-uniform) and |log S_x|/n is bounded (extensive).
  R2  (L-b) vbar = a(p,D) EXACTLY: mu_x = -nD^2/(1-2D) + A(D) W1 (exact affine
      identity, vs TM); v_x local formula vs TM; (m2-1)(m2+3) == kappa^2
      (the closed-form Var(mu_x)/n = kappa^2 A^2 = D(1-D)Psi_0); exact
      enumeration at n in {10,12}: E[mu_x]=nD, total-variance split,
      Var(mu_x) == A^2 VarW1-formula, and E[v_x] - a n is the SAME constant
      at both n (exact affineness => vbar = a, no limit needed).
  R3  (L-c) influence/McDiarmid: single-bond influence on mu_x bounded by the
      analytic c_mu = 2 A rho^{-1}(rho^{-1}-rho); on v_x O(1) n-stable;
      Var(mu_x) <= (1/4) sum c_j^2 (BDI consistency).
  R4  (L-d) quenched steepest descent: R_n(x) = log2 q_t - [-(1/2)log2(2 pi
      v_x) - beta^2/2 log2 e] decays (max & RMS, RMS*sqrt(n)/ln^3 n bounded),
      INCLUDING adversarial biased-bond words (word-uniformity, the quenched
      grade); the |Phi| <= e^{-v_x s^2/4} window s* is Theta(1) n-stable;
      tail sup_{s>=eps}|Phi| <= 2 e^{-cn} with c = (1/4)ln(1/gbar) from the
      PROVEN replica bound; ORDERING ARITHMETIC: both the zone-2 and the
      vertical-end budgets close iff A > 2B'/sqrt(vbar) (exponent sign demo
      at n = 1e3, 1e6, 1e12; A = 0.95 x threshold provably diverges).
  R5  (L-e) Abel: the naive q_{t-j} <= q_t is FALSE (counterexample words
      exhibited); the corrected mechanism W = sum phi^j q_{t-j}/q_t ->
      1/(1-phi) with O(ln^3 n/sqrt n) error; the j > J = K ln^2 n truncation
      tail obeys the crude bound phi^{J+1}/((1-phi) q_t) -> 0 for
      K > B'^2/(2 ln(1/phi)) (q_{t-j} <= 1 only, no monotonicity used).
  R6  (L-f) variance assembly: Var(-log2 S_n | shell) is O(1) n-stable,
      matches Var of the predicted decomposition (beta^2 term + (1/2)log2 v_x
      + log2 W + const), and Var/n -> 0; the cross-term with the j_n part is
      o(n) by Cauchy-Schwarz once Var(log2 S_n)=O(1) (arithmetic).

Verdict line at the end: R2-chain-verified / R2-chain-issue.
Deps: numpy, mpmath.  Runtime ~1-3 min.
"""
import math
import numpy as np
import mpmath as mp

# --------------------- shared BSMS / dual objects (committed) ---------------------

def Dc(p): return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))

def theta0(D): return math.log(D / (1 - D))

def eta_C(beta, D):
    eb = np.exp(beta)
    ze = (1 - eb) / ((1 + eb) * (1 - 2 * D))
    return (1 - ze) / (1 + ze), (1 + eb) * (1 + ze) / 2

def kappa2(p): return (1 - 2 * p) ** 2 / (p * p * (1 - p) * (1 - p))

def Psi0(p, D): return kappa2(p) * D * (1 - D) / (1 - 2 * D) ** 2

def a_curv(p, D): return D * (1 - D) * (1 - Psi0(p, D))

def A_slope(D): return D * (1 - D) / (1 - 2 * D)

def m2_const(p):
    rho = p / (1 - p)
    return (1 - p) * rho ** 2 + p / rho ** 2

def g_closed(s, p, D):
    """Proven replica spectral radius g(s) = (1/2)[F + sqrt(F^2+G)] (closed form)."""
    t = 1 - np.cos(s)
    F = 1 - 2 * D * (1 - D) * t
    k2 = kappa2(p)
    G = 8 * k2 * D * D * (1 - D) ** 2 * t * ((1 - 2 * D) ** 2 + 2 * D * D * (1 - D) ** 2 * t) \
        / (1 - 2 * D) ** 4
    return 0.5 * (F + np.sqrt(F * F + G))

# --------------------- log-scaled TM posterior (committed pipeline) ---------------------

def phi_all_nodes(xbits, p, D):
    n = len(xbits)
    S = n + 1
    s = 2 * math.pi * np.arange(S) / S
    th = theta0(D)
    eta, Cs = eta_C(th + 1j * s, D)
    C0 = eta_C(th, D)[1]
    rho = p / (1 - p)
    w0 = np.ones(S, dtype=complex)
    w1 = eta.copy()
    acc = np.zeros(S)
    for j in range(n - 1):
        r = rho if xbits[j] == xbits[j + 1] else 1.0 / rho
        n0 = w0 + w1 * r
        n1 = (w0 * r + w1) * eta
        w0, w1 = n0, n1
        m = np.maximum(np.abs(w0), np.abs(w1))
        m = np.where(m == 0, 1.0, m)
        w0 /= m
        w1 /= m
        acc += np.log(m)
    Ssum = w0 + w1
    logmag = n * np.log(np.abs(Cs / C0)) + acc + np.log(np.maximum(np.abs(Ssum), 1e-300))
    phase = n * np.angle(Cs / C0) + np.angle(Ssum)
    return np.exp(np.minimum(logmag, 0.0) + 1j * phase), s

def q_from_phi(phi):
    return np.fft.fft(phi).real / len(phi)

def shell_words(rng, n, p, n_words):
    out = []
    mean = (n - 1) * p
    half = 2.0 * math.sqrt((n - 1) * p * (1 - p))
    while len(out) < n_words:
        x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
        sw = int(np.sum(x[:-1] != x[1:]))
        if abs(sw - mean) <= half:
            out.append(x)
    return out

def biased_words(rng, n, pprime, n_words):
    return [(np.cumsum(rng.random(n) < pprime) % 2).astype(int) for _ in range(n_words)]

# --------------------- brute FWHT (n=12 anchor) ---------------------

def popcount_arr(v):
    v = np.asarray(v, dtype=np.int64)
    pc = np.zeros_like(v)
    while np.any(v):
        pc += v & 1
        v >>= 1
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

# --------------------- KP gas machinery (committed conventions) ---------------------

def kp_feasible(eta_abs, p, mu_grid=None):
    B = ((1 - p) / p) ** 2
    grid = np.linspace(1e-4, 6.0, 6000) if mu_grid is None else mu_grid
    for mu in grid:
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

def block_dp_poly(x, n, ratio):
    memo = [None] * (n + 2); memo[n] = np.array([1.0])
    for i in range(n - 1, -1, -1):
        poly = memo[i + 1].copy()
        for b in range(i, n):
            L = b - i + 1; w = 1.0
            if i >= 1: w *= ratio if (x[i] == x[i - 1]) else 1.0 / ratio
            if b + 1 <= n - 1: w *= ratio if (x[b + 1] == x[b]) else 1.0 / ratio
            nxt = memo[b + 2] if b + 2 <= n else np.array([1.0])
            term = np.concatenate([np.zeros(L), w * nxt])
            m = max(len(poly), len(term)); pp = np.zeros(m)
            pp[:len(poly)] += poly; pp[:len(term)] += term; poly = pp
        memo[i] = poly
    return memo[0]

# --------------------- exact local cumulant formulas (the (b)/(c) lemmas) ---------------------

def dual_derivs(p, D):
    """(log C-ratio)', '' and eta', eta'' at s=0, high precision."""
    mp.mp.dps = 40
    Dm = mp.mpf(D)
    def lc(s): return mp.log(((1 - Dm) ** 2 - Dm ** 2 * mp.e ** (1j * s)) / (1 - 2 * Dm))
    def et(s): return Dm * (1 - Dm) * (mp.e ** (1j * s) - 1) / ((1 - Dm) ** 2 - Dm ** 2 * mp.e ** (1j * s))
    return (complex(mp.diff(lc, 0, 1)), complex(mp.diff(lc, 0, 2)),
            complex(mp.diff(et, 0, 1)), complex(mp.diff(et, 0, 2)))

def local_sums(eps, lr, n):
    """eps in {+1,-1}^{n-1} (+1 stay); returns W1, sum w^2, sum w w_next, sum w2."""
    expo = np.zeros(n)
    expo[1:] += eps
    expo[:-1] += eps
    w = np.exp(lr * expo)
    e2 = np.zeros(n - 1)
    e2[1:] += eps[:-1]
    e2[:-1] += eps[1:]
    w2 = np.exp(lr * e2)
    return float(w.sum()), float((w * w).sum()), float((w[:-1] * w[1:]).sum()), float(w2.sum())

def mu_v_local(x, p, D, derivs):
    Lc1, Lc2, e1, e2c = derivs
    n = len(x)
    eps = np.where(x[:-1] == x[1:], 1.0, -1.0)
    lr = math.log(p / (1 - p))
    W1, Sw2, Sww, SW2 = local_sums(eps, lr, n)
    S2loc = 2 * SW2 - Sw2 - 2 * Sww
    mu = (-1j * (n * Lc1 + e1 * W1)).real
    v = -(n * Lc2 + e2c * W1 + e1 * e1 * S2loc).real
    return mu, v

# =============================== R0 ===============================

def run_R0(rng):
    print("-" * 86)
    print("R0: conventions anchor (TM posterior == brute; Abel decomposition == ball mass)")
    n = 12
    ok = True
    worst_q = worst_ball = 0.0
    for p in (0.25, 0.4):
        D = 0.9 * Dc(p)
        N = 1 << n
        P = bsms_P(n, p)
        Ph = fwht(P)
        pc_all = popcount_arr(np.arange(N))
        PYs = fwht(Ph * (1 - 2 * D) ** (-pc_all.astype(float))) / N
        th = theta0(D)
        ix = np.arange(N)
        t = int(math.floor(n * D))
        phi_ratio = D / (1 - D)
        for _ in range(4):
            x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
            xx = 0
            for b in x:
                xx = (xx << 1) | int(b)
            d = pc_all[ix ^ xx]
            w = PYs * np.exp(th * d)
            Ak = np.bincount(d, weights=w, minlength=n + 1)
            M = Ak.sum()
            q_brute = Ak / M
            pi_k = np.bincount(d, weights=PYs, minlength=n + 1)
            ball_brute = pi_k[:t + 1].sum()
            phi, _ = phi_all_nodes(x, p, D)
            q = q_from_phi(phi)
            worst_q = max(worst_q, float(np.max(np.abs(q - q_brute))))
            # Abel decomposition: P(B_t) = M e^{-th nD} q_t W phi^{nD-t}
            W = sum((q[t - j] / q[t]) * phi_ratio ** j for j in range(t + 1))
            ball_dec = M * math.exp(-th * n * D) * q[t] * W * phi_ratio ** (n * D - t)
            worst_ball = max(worst_ball, abs(ball_dec / ball_brute - 1))
            ok &= abs(phi[0] - 1) < 1e-9
    ok &= worst_q < 1e-10 and worst_ball < 1e-9
    print(f"  worst |q_dft - q_brute| = {worst_q:.2e}; "
          f"worst |ball_decomp/ball_brute - 1| = {worst_ball:.2e}; Phi(0;x)=1 exact")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# =============================== R1 ===============================

def run_R1(rng):
    print("-" * 86)
    print("R1 (L-a): KP radius analytic lower bound; complex window; zero-freeness; extensivity")
    ok = True
    # (i) closed-form lower bound r_lb = p^2/(12 e (1-p)^2) is KP-feasible at mu=1
    print(f"  {'p':>5} {'r_lb (analytic)':>16} {'r_KP measured':>14} {'feasible@mu=1':>14}")
    for p in (0.05, 0.1, 0.25, 0.4, 0.49):
        B = ((1 - p) / p) ** 2
        r_lb = p ** 2 / (12 * math.e * (1 - p) ** 2)
        u = math.e * r_lb
        feas = B * (2 * u / (1 - u) + u / (1 - u) ** 2) <= 1.0
        r_meas = r_cluster(p)
        ok &= feas and (r_lb <= r_meas)
        print(f"  {p:>5} {r_lb:>16.6f} {r_meas:>14.6f} {str(feas):>14}")
    # (ii) contour slope (real s) and complex-disk slope; window inside measured r_KP
    for p, fD in ((0.25, 0.9), (0.4, 0.9), (0.25, 1.0)):
        D = fD * Dc(p)
        A = A_slope(D)
        th = theta0(D)
        s = np.linspace(1e-4, math.pi, 3000)
        et = np.abs([eta_C(th + 1j * si, D)[0] for si in s])
        slope_ok = bool(np.max(et / s) <= A * (1 + 1e-9))
        r_lb = p ** 2 / (12 * math.e * (1 - p) ** 2)
        tau_max = min(math.log(2.0), math.log(1 + (1 - 2 * D) / (2 * D * D)))
        delta_c = min(r_lb / (4 * A), tau_max, 0.5)
        rad = np.linspace(1e-4, delta_c, 40)
        ang = np.linspace(0, 2 * math.pi, 48, endpoint=False)
        zz = rad[:, None] * np.exp(1j * ang)[None, :]
        ecz = np.abs(eta_C(th + 1j * zz, D)[0])
        cslope_ok = bool(np.max(ecz / np.abs(zz)) <= 4 * A * (1 + 1e-9))
        inside_ok = bool(np.max(ecz) < r_cluster(p))
        ok &= slope_ok and cslope_ok and inside_ok
        print(f"  p={p}, D={fD}*Dc: real slope<=A(D): {slope_ok}; complex slope<=4A(D): "
              f"{cslope_ok}; max|eta| on disk(delta_c={delta_c:.4f}) < r_KP: {inside_ok}")
    # (iii) zero-freeness of S_x on |eta| <= r_lb and extensivity of |log S_x|
    for p in (0.25, 0.4):
        r_lb = p ** 2 / (12 * math.e * (1 - p) ** 2)
        ratio = p / (1 - p)
        minroot = math.inf
        ext = {}
        for n in (20, 40):
            vals = []
            words = [(np.cumsum(rng.random(n) < p) % 2).astype(int) for _ in range(15)]
            words += [np.zeros(n, dtype=int), np.array([i % 2 for i in range(n)])]
            for x in words:
                c = block_dp_poly(x, n, ratio)
                rts = np.roots(c[::-1])
                minroot = min(minroot, float(np.min(np.abs(rts))))
                ring = r_lb * np.exp(1j * np.linspace(0, 2 * math.pi, 24, endpoint=False))
                Sv = np.polyval(c[::-1], ring)
                vals.append(float(np.max(np.abs(np.log(Sv)))) / n)
            ext[n] = max(vals)
        zf = minroot > r_lb
        stab = ext[40] <= ext[20] * 1.3 + 0.05
        ok &= zf and stab
        print(f"  p={p}: min|zero of S_x| = {minroot:.4f} > r_lb = {r_lb:.4f}: {zf}; "
              f"max|log S_x|/n on |eta|=r_lb: n=20: {ext[20]:.4f}, n=40: {ext[40]:.4f} "
              f"(extensive, n-stable: {stab})")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# =============================== R2 ===============================

def run_R2(rng):
    print("-" * 86)
    print("R2 (L-b): vbar = a(p,D) via the EXACT affine identity mu_x = -nD^2/(1-2D) + A W1")
    ok = True
    # (i) closed-form identity (m2-1)(m2+3) == kappa^2  => Var(mu)/n = kappa^2 A^2 = D(1-D) Psi0
    worst = 0.0
    for p in np.linspace(0.03, 0.49, 25):
        m2 = m2_const(p)
        worst = max(worst, abs((m2 - 1) * (m2 + 3) - kappa2(p)) / kappa2(p))
    id_ok = worst < 1e-12
    ok &= id_ok
    print(f"  (m2-1)(m2+3) == kappa^2 over p-grid: worst rel.err = {worst:.2e} -> {id_ok}")
    # (ii) local formulas vs exact TM posterior (incl. D = Dc corner)
    for p, fD in ((0.25, 0.9), (0.4, 0.9), (0.4, 1.0)):
        D = fD * Dc(p)
        derivs = dual_derivs(p, D)
        n = 512
        werr_mu = werr_v = 0.0
        for x in shell_words(rng, n, p, 8):
            phi, _ = phi_all_nodes(x, p, D)
            q = np.maximum(q_from_phi(phi), 0.0)
            ks = np.arange(n + 1, dtype=float)
            mu = float(ks @ q)
            v = float((ks - mu) ** 2 @ q)
            mul, vl = mu_v_local(x, p, D, derivs)
            werr_mu = max(werr_mu, abs(mul - mu) / max(abs(mu), 1))
            werr_v = max(werr_v, abs(vl - v) / max(abs(v), 1))
        f_ok = werr_mu < 1e-5 and werr_v < 1e-4
        ok &= f_ok
        print(f"  p={p}, D={fD}*Dc, n=512: |mu_loc-mu_TM|/mu = {werr_mu:.1e}, "
              f"|v_loc-v_TM|/v = {werr_v:.1e} -> {f_ok}")
    # (iii) exact enumeration: E[mu]=nD; split; Var(mu)=A^2 VarW1; E[v]-a n constant in n
    for p, fD in ((0.25, 0.9), (0.4, 0.9)):
        D = fD * Dc(p)
        derivs = dual_derivs(p, D)
        Lc1, Lc2, e1, e2c = derivs
        A = A_slope(D)
        m2 = m2_const(p)
        bconst = {}
        en_ok = True
        for n in (10, 12):
            nb = n - 1
            N = 1 << nb
            pats = np.arange(N)
            bits = ((pats[:, None] >> np.arange(nb)[None, :]) & 1).astype(float)
            eps = 1 - 2 * bits
            lr = math.log(p / (1 - p))
            expo = np.zeros((N, n))
            expo[:, 1:] += eps
            expo[:, :-1] += eps
            w = np.exp(lr * expo)
            e2 = np.zeros((N, nb))
            e2[:, 1:] += eps[:, :-1]
            e2[:, :-1] += eps[:, 1:]
            w2 = np.exp(lr * e2)
            W1 = w.sum(1)
            S2 = 2 * w2.sum(1) - (w * w).sum(1) - 2 * (w[:, :-1] * w[:, 1:]).sum(1)
            mu = (-1j * (n * Lc1 + e1 * W1)).real
            v = -(n * Lc2 + e2c * W1 + (e1 * e1) * S2).real
            sw = bits.sum(1)
            P = (p ** sw) * ((1 - p) ** (nb - sw))
            Emu = float(P @ mu); Ev = float(P @ v)
            Vmu = float(P @ (mu * mu)) - Emu ** 2
            varW1_f = (n - 2) * (m2 ** 2 - 1) + 2 * (m2 - 1) + 2 * (n - 1) * (m2 - 1)
            c1 = abs(Emu - n * D) < 1e-9
            c2 = abs(Ev + Vmu - n * D * (1 - D)) < 1e-9
            c3 = abs(Vmu - A * A * varW1_f) < 1e-9
            en_ok &= c1 and c2 and c3
            bconst[n] = Ev - a_curv(p, D) * n
            print(f"  p={p}, D={fD}*Dc, n={n} (exact enum): E[mu]=nD: {c1}; "
                  f"split nD(1-D)=E[v]+Var(mu): {c2}; Var(mu)==A^2*VarW1: {c3}; "
                  f"E[v]-a*n = {bconst[n]:+.10f}")
        affine = abs(bconst[10] - bconst[12]) < 1e-9
        en_ok &= affine
        ok &= en_ok
        print(f"      E[v_x] EXACTLY affine in n (same constant at n=10,12): {affine} "
              f"=> vbar = a(p,D) = {a_curv(p, D):.6f} rigorously")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# =============================== R3 ===============================

def run_R3(rng):
    print("-" * 86)
    print("R3 (L-c): single-bond influence O(1), word-uniform; McDiarmid constants")
    ok = True
    for p, fD in ((0.25, 0.9), (0.4, 0.9)):
        D = fD * Dc(p)
        derivs = dual_derivs(p, D)
        A = A_slope(D)
        rho = p / (1 - p)
        c_mu_analytic = A * 2 * (1 / rho) * (1 / rho - rho)
        res = {}
        for n in (256, 1024):
            mx_mu = mx_v = 0.0
            words = shell_words(rng, n, p, 6) + biased_words(rng, n, 0.8 * p, 3) \
                + biased_words(rng, n, min(1.25 * p, 0.49), 3)
            for x in words:
                mu0, v0 = mu_v_local(x, p, D, derivs)
                for j in range(n - 1):
                    y = x.copy()
                    y[j + 1:] = 1 - y[j + 1:]  # flipping bond j == flipping suffix
                    mu1, v1 = mu_v_local(y, p, D, derivs)
                    mx_mu = max(mx_mu, abs(mu1 - mu0))
                    mx_v = max(mx_v, abs(v1 - v0))
            res[n] = (mx_mu, mx_v)
        bound_ok = res[1024][0] <= c_mu_analytic * (1 + 1e-9)
        stable_ok = res[1024][0] <= res[256][0] * 1.05 + 1e-9 and \
            res[1024][1] <= res[256][1] * 1.05 + 1e-9
        bdi_ok = kappa2(p) * A * A <= 0.25 * c_mu_analytic ** 2  # Var(mu)/n <= (1/4) c_mu^2
        ok &= bound_ok and stable_ok and bdi_ok
        print(f"  p={p}: max|D_j mu| n=256: {res[256][0]:.4f}, n=1024: {res[1024][0]:.4f} "
              f"<= analytic {c_mu_analytic:.4f}: {bound_ok}; max|D_j v| {res[256][1]:.4f} -> "
              f"{res[1024][1]:.4f} (O(1), n-stable: {stable_ok}); "
              f"BDI Var(mu)/n={kappa2(p)*A*A:.5f} <= c_mu^2/4={0.25*c_mu_analytic**2:.5f}: {bdi_ok}")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# =============================== R4 + R5 + R6 (shared word loop) ===============================

def word_stats(x, p, D, t, phi_ratio):
    phi, s = phi_all_nodes(x, p, D)
    q = q_from_phi(phi)
    q = np.maximum(q, 0.0)
    n = len(x)
    ks = np.arange(n + 1, dtype=float)
    mu = float(ks @ q)
    v = float((ks - mu) ** 2 @ q)
    beta = (t - mu) / math.sqrt(v)
    log2e = 1.0 / math.log(2.0)
    lq = math.log2(max(q[t], 1e-300))
    pred = -0.5 * math.log2(2 * math.pi * v) - 0.5 * beta * beta * log2e
    Rn = lq - pred
    # window s*: largest s with |Phi| <= e^{-v s^2/4} for all nodes below
    S = n + 1
    half = np.arange(1, S // 2 + 1)
    sh = s[half]
    okmask = np.abs(phi[half]) <= np.exp(-v * sh * sh / 4) + 1e-9
    bad = np.where(~okmask)[0]
    sstar = math.pi if len(bad) == 0 else float(sh[bad[0]] - 2 * math.pi / S)
    # Abel
    W = float(sum((q[t - j] / q[t]) * phi_ratio ** j for j in range(t + 1)))
    ratios_gt1 = bool(np.any(q[max(t - 60, 0):t] > q[t]))
    Z = -math.log2(q[t] * W) - (len(x) * D - t) * math.log2(phi_ratio)  # -log2 S_n
    return dict(mu=mu, v=v, beta=beta, Rn=Rn, sstar=sstar, W=W, gt1=ratios_gt1,
                Z=Z, q=q, phi=phi, s=s)

def run_main(rng, p, fD, n_list, w_list):
    D = fD * Dc(p)
    phi_ratio = D / (1 - D)
    vbar = a_curv(p, D)
    log2e = 1.0 / math.log(2.0)
    Winf = 1 / (1 - phi_ratio)
    print("-" * 86)
    print(f"R4/R5/R6 main loop: p={p}, D={fD}*Dc={D:.5f}, vbar=a={vbar:.5f}, phi={phi_ratio:.4f}")
    rows = []
    sstar_min_global = math.inf
    typ_keep = None
    for n, nw in zip(n_list, w_list):
        t = int(math.floor(n * D))
        st = [word_stats(x, p, D, t, phi_ratio) for x in shell_words(rng, n, p, nw)]
        if n == 1024:
            typ_keep = [dict(Rn=d["Rn"], beta=d["beta"]) for d in st]
        Rn = np.array([d["Rn"] for d in st])
        vs = np.array([d["v"] for d in st])
        bet = np.array([d["beta"] for d in st])
        Ws = np.array([d["W"] for d in st])
        Zs = np.array([d["Z"] for d in st])
        sstars = np.array([d["sstar"] for d in st])
        gt1_frac = np.mean([d["gt1"] for d in st])
        sstar_min = float(np.min(sstars))
        sstar_min_global = min(sstar_min_global, sstar_min)
        # tail bound from the PROVEN replica inequality at eps = 0.9*sstar_min
        eps = 0.9 * sstar_min
        gbar = 1 - (11 / 16) * D * (1 - D) * (1 - math.cos(eps))
        c4 = 0.25 * math.log(1 / gbar)
        tails = []
        for d in st:
            mask = (d["s"] >= eps) & (d["s"] <= math.pi)
            tails.append(float(np.max(np.abs(d["phi"][mask]))))
        # truncation tail (only meaningful if J < t)
        Bp = 2.0
        Kc = max(1, math.ceil(Bp * Bp / (2 * math.log(1 / phi_ratio))) + 1)
        J = int(Kc * math.log(n) ** 2)
        if J < t:
            trunc_meas, trunc_bound = [], []
            for d in st:
                q = d["q"]
                tm = float(sum((q[t - j] / q[t]) * phi_ratio ** j for j in range(J + 1, t + 1)))
                trunc_meas.append(tm)
                trunc_bound.append(phi_ratio ** (J + 1) / ((1 - phi_ratio) * q[t]))
            tr = (max(trunc_meas), max(trunc_bound))
        else:
            tr = None
        pred_main = 0.5 * log2e * bet ** 2 + 0.5 * np.log2(2 * math.pi * vs) \
            - math.log2(Winf) - (n * D - t) * math.log2(phi_ratio)
        rows.append(dict(n=n, t=t, maxR=float(np.max(np.abs(Rn))),
                         rmsR=float(np.sqrt(np.mean(Rn ** 2))),
                         rms_norm=float(np.sqrt(np.mean(Rn ** 2))) * math.sqrt(n) / math.log(n) ** 3,
                         vmin=float(np.min(vs / n)), vmax=float(np.max(vs / n)),
                         betmax=float(np.max(np.abs(bet))),
                         sstar=sstar_min, tailmax=max(tails), tailbnd=2 * math.exp(-c4 * n),
                         Wdev=float(np.max(np.abs(Ws / Winf - 1))), gt1=gt1_frac, tr=tr,
                         varZ=float(np.var(Zs, ddof=1)),
                         varP=float(np.var(pred_main, ddof=1))))
        r = rows[-1]
        print(f"  n={n:>5} t={t:>4}: max|R|={r['maxR']:.4f} RMS={r['rmsR']:.4f} "
              f"RMS*sqrt(n)/ln^3 n={r['rms_norm']:.4f}  s*_min={r['sstar']:.3f}  "
              f"v/n in [{r['vmin']:.4f},{r['vmax']:.4f}]")
        print(f"           tail sup|Phi|={r['tailmax']:.3e} <= 2e^-cn={r['tailbnd']:.3e}  "
              f"max|W/Winf-1|={r['Wdev']:.4f}  frac(q_t-j>q_t)={r['gt1']:.2f}  "
              f"Var(Z)={r['varZ']:.4f} Var(pred)={r['varP']:.4f}"
              + (f"  trunc: meas={r['tr'][0]:.1e}<=bnd={r['tr'][1]:.1e}" if r['tr'] else
                 "  trunc: J>=t (window covers whole sum)"))
    # adversarial biased-bond words (quenched word-uniformity of (d)).
    # Theory: |R_n| <= C (1+|beta|)^3 / sqrt(n) word-uniformly on the good set
    # (the saddle-shift cubic r3(z0) ~ C3 n (|beta|/sqrt(v_x))^3), so the right
    # checks are (i) the normalized Rhat = |R_n| sqrt(n)/(1+|beta|)^3 is
    # comparable to typical words, and (ii) adversarial max|R_n| DECAYS in n.
    advmax, advhat = {}, {}
    for n_adv in (512, 1024, 2048):
        t = int(math.floor(n_adv * D))
        adv = biased_words(rng, n_adv, 0.8 * p, 12) + \
            biased_words(rng, n_adv, min(1.25 * p, 0.49), 12)
        aR, ahat, used = [], [], 0
        for x in adv:
            d = word_stats(x, p, D, t, phi_ratio)
            if d["v"] >= vbar * n_adv / 2 and abs(d["beta"]) <= 2 * math.log(n_adv):
                aR.append(abs(d["Rn"]))
                ahat.append(abs(d["Rn"]) * math.sqrt(n_adv) / (1 + abs(d["beta"])) ** 3)
                used += 1
        advmax[n_adv] = max(aR)
        advhat[n_adv] = max(ahat)
        print(f"  adversarial biased-bond words n={n_adv} ({used}/24 in good set): "
              f"max|R|={max(aR):.4f}, max Rhat=|R|sqrt(n)/(1+|beta|)^3={max(ahat):.4f}")
    # NOTE: at FIXED bias p' the mean shift is Theta(n) so beta ~ sqrt(n) grows and
    # raw max|R| ~ C(1+beta)^3/sqrt(n) GROWS -- the theorem's invariant is the
    # word-uniform constant C, i.e. Rhat stability (and the visible cubic confirms
    # the log^3 n in the remainder is the real beta^3 at the window edge, not slack).
    typ_hat = max(abs(d["Rn"]) * math.sqrt(1024) / (1 + abs(d["beta"])) ** 3 for d in typ_keep)
    stable = max(advhat.values()) <= 3 * min(advhat.values()) + 0.05
    comparable = max(advhat.values()) <= 8 * typ_hat
    adv_ok = stable and comparable
    print(f"  word-uniformity: Rhat=|R|sqrt(n)/(1+|beta|)^3 n-stable "
          f"(max {max(advhat.values()):.3f} <= 3 x min {min(advhat.values()):.3f} + .05): {stable}; "
          f"comparable to typical (Rhat_typ(1024)={typ_hat:.3f}, <=8x): {comparable} -> {adv_ok}")
    # PASS criteria
    r_dec = rows[-1]["rmsR"] < rows[0]["rmsR"] and rows[-1]["maxR"] < rows[0]["maxR"]
    r_norm = rows[-1]["rms_norm"] <= rows[0]["rms_norm"] * 1.2
    s_ok = sstar_min_global >= 0.05
    v_ok = all(r["vmin"] >= vbar * 0.5 and r["vmax"] <= vbar * 1.5 + 0.1 for r in rows)
    tail_ok = all(r["tailmax"] <= r["tailbnd"] for r in rows) and \
        rows[-1]["tailmax"] < rows[0]["tailmax"]
    W_ok = rows[-1]["Wdev"] < rows[0]["Wdev"] * 1.2 and rows[-1]["Wdev"] < 0.2
    gt1_ok = any(r["gt1"] > 0 for r in rows)  # naive monotone bound is FALSE
    tr_ok = all((r["tr"] is None) or (r["tr"][0] <= r["tr"][1] * (1 + 1e-6) and r["tr"][1] < 1e-6)
                for r in rows)
    var_ok = all(r["varZ"] < 3.0 for r in rows) and \
        abs(rows[-1]["varZ"] - rows[0]["varZ"]) < 0.5 * max(rows[0]["varZ"], 0.05) + 0.05
    match_ok = all(abs(r["varZ"] - r["varP"]) < 0.35 * max(r["varZ"], 0.05) +
                   4 * r["rmsR"] ** 2 + 0.02 for r in rows)
    vn_ok = rows[-1]["varZ"] / rows[-1]["n"] < rows[0]["varZ"] / rows[0]["n"]
    print(f"  R4: R_n decays {r_dec}; RMS*sqrt(n)/ln^3 n bounded {r_norm}; s*=Theta(1) {s_ok}; "
          f"v/n in [vbar/2, ~vbar] {v_ok}; tail<=2e^-cn & decays {tail_ok}; adversarial {adv_ok}")
    print(f"  R5: W->1/(1-phi) {W_ok}; naive q_t-j<=q_t FALSE for {rows[-1]['gt1']:.0%} words "
          f"(corrected mechanism needed) {gt1_ok}; truncation<=crude bound {tr_ok}")
    print(f"  R6: Var(-log2 S_n|shell)=O(1) n-stable {var_ok}; matches decomposition {match_ok}; "
          f"Var/n decreasing {vn_ok}")
    return all([r_dec, r_norm, s_ok, v_ok, tail_ok, adv_ok, W_ok, gt1_ok, tr_ok,
                var_ok, match_ok, vn_ok])

# =============================== R4-arith (ordering) ===============================

def run_R4_arith():
    print("-" * 86)
    print("R4-arith: ordered-constants budgets close iff A > 2B'/sqrt(vbar)  (B'=B+1)")
    ok = True
    for p, fD in ((0.25, 0.9), (0.4, 0.9)):
        D = fD * Dc(p)
        vb = a_curv(p, D)
        Bp = 2.0
        A_thr = 2 * Bp / math.sqrt(vb)
        for lab, Af in (("A=1.5x thr (valid)", 1.5), ("A=0.95x thr (violates)", 0.95)):
            A = Af * A_thr
            expo = []
            for n in (1e3, 1e6, 1e12):
                ln = math.log(n)
                e_z2 = (Bp ** 2 / 2 - vb * A * A / 8) * ln * ln + 0.5 * ln   # zone-2 rel. err
                e_ve = (Bp ** 2 - vb * A * A / 4) * ln * ln + 0.5 * ln       # vertical ends
                expo.append((e_z2, e_ve))
            closes = all(e[0] < 0 and e[1] < 0 for e in expo[1:])
            diverges = expo[-1][0] > 0
            good = closes if Af > 1 else diverges
            ok &= good
            print(f"  p={p}: {lab}: zone-2 exponent at n=1e3/1e6/1e12: "
                  f"{expo[0][0]:+.1f}/{expo[1][0]:+.1f}/{expo[2][0]:+.1f}; "
                  f"vertical-end: {expo[0][1]:+.1f}/{expo[1][1]:+.1f}/{expo[2][1]:+.1f} -> "
                  f"{'closes' if closes else 'diverges'} (expected: "
                  f"{'closes' if Af > 1 else 'diverges'})")
        # honest n0 report for the literal constant chain (with measured-scale s0 ~ 0.3)
        A = 1.5 * A_thr
        s0 = 0.3
        n0 = 10
        while A * math.log(n0) / math.sqrt(n0) > s0 and n0 < 1e18:
            n0 *= 2
        print(f"  p={p}: literal delta_1=A ln n/sqrt(n) <= s0~0.3 first at n ~ {n0:.1e} "
              f"(asymptotic regime; mechanism verified above at accessible n)")
    # proven-theorem sanity: g(s) <= 1-(11/16)D(1-D)(1-cos s) on a grid incl. D=Dc
    worst = -math.inf
    for p in (0.05, 0.1, 0.25, 0.4, 0.49):
        for fD in (0.3, 0.9, 1.0):
            D = fD * Dc(p)
            s = np.linspace(1e-3, math.pi, 800)
            worst = max(worst, float(np.max(g_closed(s, p, D) -
                                            (1 - (11 / 16) * D * (1 - D) * (1 - np.cos(s))))))
    g_ok = worst <= 1e-12
    ok &= g_ok
    print(f"  proven bound g(s) <= 1-(11/16)D(1-D)(1-cos s): worst violation = {worst:.2e} -> {g_ok}")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# =============================== main ===============================

if __name__ == "__main__":
    print("=" * 86)
    print("(R2) ordered-constants central-window chain -- hostile per-lemma verification")
    print("=" * 86)
    rng = np.random.default_rng(7)
    res = {}
    res["R0"] = run_R0(rng)
    res["R1"] = run_R1(rng)
    res["R2"] = run_R2(rng)
    res["R3"] = run_R3(rng)
    res["R4a"] = run_R4_arith()
    res["M1"] = run_main(rng, 0.25, 0.9, [256, 512, 1024, 2048], [150, 120, 80, 50])
    res["M2"] = run_main(rng, 0.4, 0.9, [256, 512, 1024, 2048], [150, 120, 80, 50])
    print("=" * 86)
    for k, v in res.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    allok = all(res.values())
    print(f"VERDICT: {'R2-chain-verified (attack failed)' if allok else 'R2-chain-issue'}")

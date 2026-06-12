#!/usr/bin/env python3
"""
probe_7_34_gap1_writeout_v2.py
=============================================================================
Verification of EVERY numbered inequality of the GAP-1 write-out v2
(/tmp/gap1_lemma_v2.tex): the per-word, shell-uniform steepest-descent
remainder  q_t = (2 pi v_x)^{-1/2} e^{-beta^2/2} (1 + R_n),
|R_n| <= C (1+|beta|)^3 / sqrt(n)  for the BSMS Gray-region posterior
(Remark 7.34b achievability chain).

Architecture under test (two-scale, post-eps-conflation repair):
  eps_strip = delta_c/2 (analyticity scale)  !=  eps_split (tail handoff)
  zone-1+2 [0, eps2] on the SHIFTED contour Im z = y0 (exact saddle
  cancellation E2, weighted-cubic E5, in-strip domination E6/E7),
  vertical ends at Re z = +-eps2 (favorable sign, E8),
  zone-2' [eps2, eps_split] real line (quenched dyadic net, E9),
  zone-3 [eps_split, pi] the committed tail lemma (E11).

Checks (PASS/FAIL each):
  V0  machinery anchor: numpy complex-z log Phi == mpmath; Phi(0)=1
  V1  (E2) exact saddle-cancellation algebra (random instances, 1e-12)
  V2  (I2)/(I3): measured C3; cubic Taylor remainder bound on the window
  V3  (E1) contour deformation: q_t(FFT) == central + verticals + real
  V4  (E3) T_G, (E5)+(E5') inner weighted cubic, (E8) vertical ends
      incl. favorable-sign profile; (E7)/(N2) outer-zone arithmetic
  V5  zone-2': (E9) per-word band bound; (I5) annealed prefactor C1_eff;
      (E9a/E9b) net constants/margins; Jordan inequality
  V6  (E11) zone-3 committed tail bound at eps_split
  V7  final: |R_n| <= C (1+|beta|)^3/sqrt(n) with the explicit C;
      implied (tight) C reported; RMS decay
  V8  (N0)-(N5), n_a, n_b: explicit n0(p,D) computed (log-domain scan)
  V9  pitfall regressions: (i) two-scale (no tail bound at eps2),
      (ii) t|y0| = Theta(sqrt n) grows yet never needed as O(1),
      (iii) all constant orderings at both physical points

Points: (p, D) in {(0.25, 0.9 Dc), (0.4, 0.5 Dc)}; n up to 2048.
Deps: numpy, mpmath.  Runtime ~4-8 min.
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 40
LOG2PI = math.log(2 * math.pi)

# ----------------------------- closed forms -----------------------------

def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))

def A_slope(D):
    return D * (1 - D) / (1 - 2 * D)

def kappa2(p):
    return (1 - 2 * p) ** 2 / (p * p * (1 - p) * (1 - p))

def Psi0(p, D):
    return kappa2(p) * D * (1 - D) / (1 - 2 * D) ** 2

def vbar_closed(p, D):
    return D * (1 - D) * (1 - Psi0(p, D))

def g_closed(s, p, D):
    """Proven replica spectral radius g(s) (closed form, committed)."""
    t = 1 - np.cos(s)
    F = 1 - 2 * D * (1 - D) * t
    k2 = kappa2(p)
    G = 8 * k2 * D * D * (1 - D) ** 2 * t * ((1 - 2 * D) ** 2
        + 2 * D * D * (1 - D) ** 2 * t) / (1 - 2 * D) ** 4
    return 0.5 * (F + np.sqrt(F * F + G))

def eta_closed(p, D, z):
    """eta_z, complex z (w = e^{iz}); committed closed form."""
    w = np.exp(1j * np.asarray(z, dtype=complex))
    return D * (1 - D) * (w - 1) / ((1 - D) ** 2 - D ** 2 * w)

def C_ratio(p, D, z):
    w = np.exp(1j * np.asarray(z, dtype=complex))
    return ((1 - D) ** 2 - D ** 2 * w) / (1 - 2 * D)

# ------------------- log Phi on complex z (numpy, vectorized) -------------------

def bond_arrays(x, p):
    n = len(x)
    rho = p / (1.0 - p)
    l = np.ones(n)
    r = np.ones(n)
    for i in range(1, n):
        l[i] = rho if x[i] == x[i - 1] else 1.0 / rho
    for b in range(0, n - 1):
        r[b] = rho if x[b + 1] == x[b] else 1.0 / rho
    return l, r

def logPhi_grid(x, p, D, zarr):
    """Returns (logmag, phase) of Phi(z;x) on an array of complex z.
    phase is principal-valued (use exp for the function value)."""
    n = len(x)
    l, r = bond_arrays(x, p)
    eta = eta_closed(p, D, zarr)
    Cr = C_ratio(p, D, zarr)
    m = len(np.atleast_1d(eta))
    F1 = np.ones(m, dtype=complex)
    F2 = np.ones(m, dtype=complex)
    G1 = np.zeros(m, dtype=complex)
    logsc = np.zeros(m)
    for i in range(n - 1, -1, -1):
        Gi = eta * (r[i] * F2 + G1)
        Fi = F1 + l[i] * Gi
        F2 = F1
        F1 = Fi
        G1 = Gi
        if (i % 8) == 0:
            a = np.abs(F1)
            a = np.where(a > 0, a, 1.0)
            F1 = F1 / a
            F2 = F2 / a
            G1 = G1 / a
            logsc += np.log(a)
    logmag = n * np.log(np.abs(Cr)) + logsc + np.log(np.maximum(np.abs(F1), 1e-300))
    phase = n * np.angle(Cr) + np.angle(F1)
    return logmag, phase

def Phi_vals(x, p, D, zarr):
    lm, ph = logPhi_grid(x, p, D, zarr)
    return np.exp(lm + 1j * ph)

def logPhi_ray(x, p, D, z, nsteps=64):
    """Branch-tracked complex log Phi(z) along the ray 0 -> z."""
    ts = np.linspace(0, 1, nsteps + 1)[1:]
    lm, ph = logPhi_grid(x, p, D, ts * z)
    ph = np.unwrap(np.concatenate([[0.0], ph]))[1:]
    return lm[-1] + 1j * ph[-1]

# mpmath reference (from the committed saddle_strip probe conventions)

def logPhi_mp(x, p, D, z):
    n = len(x)
    w = mp.e ** (1j * mp.mpc(z))
    eta = D * (1 - D) * (w - 1) / ((1 - D) ** 2 - D ** 2 * w)
    Cr = ((1 - D) ** 2 - D ** 2 * w) / (1 - 2 * D)
    rho = mp.mpf(p) / (1 - p)
    l = [mp.mpf(1)] * n
    r = [mp.mpf(1)] * n
    for i in range(1, n):
        l[i] = rho if x[i] == x[i - 1] else 1 / rho
    for b in range(0, n - 1):
        r[b] = rho if x[b + 1] == x[b] else 1 / rho
    F1 = mp.mpc(1)
    F2 = mp.mpc(1)
    G1 = mp.mpc(0)
    logsc = mp.mpf(0)
    for i in range(n - 1, -1, -1):
        Gi = eta * (r[i] * F2 + G1)
        Fi = F1 + l[i] * Gi
        F2 = F1
        F1 = Fi
        G1 = Gi
        a = abs(F1)
        if a > 0:
            F1 /= a
            F2 /= a
            G1 /= a
            logsc += mp.log(a)
    return n * mp.log(Cr) + mp.log(F1) + logsc

# ------------------------- committed q_t pipeline (FFT) -------------------------

def phi_fft_nodes(x, p, D):
    n = len(x)
    S = n + 1
    s = 2 * math.pi * np.arange(S) / S
    lm, ph = logPhi_grid(x, p, D, s)
    return np.exp(np.minimum(lm, 0.0) + 1j * ph), s

def word_qt_mu_v(x, p, D, t):
    phi, s = phi_fft_nodes(x, p, D)
    q = np.fft.fft(phi).real / len(phi)
    qc = np.maximum(q, 0.0)
    ks = np.arange(len(q), dtype=float)
    mu = float(ks @ qc)
    v = float((ks - mu) ** 2 @ qc)
    return float(q[t]), mu, v, phi, s

def shell_words(rng, n, p, n_words):
    out = []
    mean = (n - 1) * p
    half = 2.0 * math.sqrt((n - 1) * p * (1 - p))
    while len(out) < n_words:
        x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
        sw = int(np.sum(x[:-1] != x[1:]))
        if abs(sw - mean) <= half:
            out.append([int(b) for b in x])
    return out

def biased_words(rng, n, pp, n_words):
    return [[int(b) for b in (np.cumsum(rng.random(n) < pp) % 2).astype(int)]
            for _ in range(n_words)]

# ------------------------- per-point constant pack -------------------------

class Pack:
    pass

def make_pack(p, fD, C3_meas=None):
    P = Pack()
    P.p = p
    P.D = fD * Dc(p)
    D = P.D
    P.A = A_slope(D)
    P.r_lb = p ** 2 / (12 * math.e * (1 - p) ** 2)
    P.tau_max = min(math.log(2.0), math.log(1 + (1 - 2 * D) / (2 * D * D)))
    P.delta_c = min(P.r_lb / (4 * P.A), P.tau_max, 0.5)
    P.eps_strip = 0.5 * P.delta_c
    P.vbar = vbar_closed(p, D)
    P.c_g = (11.0 / 16.0) * D * (1 - D)
    P.delta_KP_lb = P.r_lb / P.A
    rho = p / (1 - p)
    P.c_mu = 2 * P.A * (1 / rho) * (1 / rho - rho)
    P.B3 = 3 * P.c_mu / math.sqrt(P.vbar)
    P.c2 = P.c_g / (4 * math.pi ** 2)
    P.c4 = P.c_g / (16 * math.pi ** 2)
    P.c5 = P.c_g / (4 * math.pi ** 2)
    if C3_meas is not None:
        P.C3 = C3_meas
        P.eps2 = min(P.eps_strip, 3 * P.vbar / (32 * P.C3))
        P.eps_split = min(P.delta_KP_lb, P.vbar / (4 * P.C3))
        gbar = 1 - (11 / 16) * D * (1 - D) * (1 - math.cos(P.eps_split))
        P.c_tail = 0.25 * math.log(1 / gbar)
        f = (2.0 / P.vbar) ** 1.5
        P.Ca = (8 * math.e * P.C3 / 3) / math.sqrt(2 * math.pi) * f
        P.Cb = (2 * math.e * P.C3 / 3) * f
        P.C = P.Ca + P.Cb + 5.0
        P.J = max(1, math.ceil(math.log2(P.eps_split / P.eps2)))
    return P

def sigma_n(P, n):
    return (3.0 / (4 * P.C3 * n)) ** (1.0 / 3.0)

# =============================== V0 ===============================

def run_V0(rng, packs):
    print("-" * 88)
    print("V0: machinery anchor (numpy complex-z logPhi == mpmath; Phi(0;x)=1)")
    ok = True
    for P in packs:
        n = 96
        x = shell_words(rng, n, P.p, 1)[0]
        worst = 0.0
        for z in (0.01, 0.01 + 0.005j, -0.02 + 0.008j, 0.3, 2.5):
            lm, ph = logPhi_grid(x, P.p, P.D, np.array([z]))
            ref = logPhi_mp(x, P.p, P.D, z)
            v_np = complex(np.exp(lm[0] + 1j * ph[0]))
            v_mp = complex(mp.e ** ref)
            worst = max(worst, abs(v_np - v_mp) / max(abs(v_mp), 1e-30))
        lm0, _ = logPhi_grid(x, P.p, P.D, np.array([0.0]))
        ok &= worst < 1e-8 and abs(lm0[0]) < 1e-9
        print(f"  p={P.p}, D={P.D:.5f}: worst rel err vs mpmath = {worst:.2e}; "
              f"|log Phi(0)| = {abs(lm0[0]):.1e}")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# =============================== V1 ===============================

def run_V1(rng):
    print("-" * 88)
    print("V1 (E2): exact saddle-cancellation algebra "
          "i(mu-t)z - v z^2/2 |_{z=s+iy0} = -beta^2/2 - v s^2/2")
    worst = 0.0
    for _ in range(2000):
        mu = rng.uniform(1, 300)
        v = rng.uniform(0.5, 200)
        t = mu + rng.uniform(-3, 3) * math.sqrt(v)
        s = rng.uniform(-0.5, 0.5)
        y0 = (mu - t) / v
        beta = (t - mu) / math.sqrt(v)
        z = s + 1j * y0
        lhs = 1j * (mu - t) * z - v * z * z / 2
        rhs = -beta * beta / 2 - v * s * s / 2
        worst = max(worst, abs(lhs - rhs))
    ok = worst < 1e-10
    print(f"  2000 random instances: worst |lhs - rhs| = {worst:.2e} -> "
          f"{'PASS' if ok else 'FAIL'}")
    return ok

# =============================== V2 ===============================

def measure_C3(rng, P, verbose=True):
    """C3_meas = sup over words/z of |d^3_z log Phi|/n on |z|<=(4/3)eps_strip."""
    h = P.eps_strip / 24.0
    worst = 0.0
    for n in (256, 512):
        words = shell_words(rng, n, P.p, 4) + biased_words(rng, n, 0.8 * P.p, 2)
        ss = np.linspace(-1.15 * P.eps_strip, 1.15 * P.eps_strip, 7)
        ys = (-0.5 * P.eps_strip, 0.0, 0.5 * P.eps_strip)
        for x in words:
            for y in ys:
                for s0 in ss:
                    z = s0 + 1j * y
                    if abs(z) + 2 * h > (4.0 / 3.0) * P.eps_strip:
                        continue
                    zs = z + h * np.array([-2, -1, 0, 1, 2])
                    lm, ph = logPhi_grid(x, P.p, P.D, zs)
                    L = lm + 1j * np.unwrap(ph)
                    d3 = (-L[0] + 2 * L[1] - 2 * L[3] + L[4]) / (2 * h ** 3)
                    worst = max(worst, abs(d3) / n)
    return worst

def run_V2(rng, packs_raw):
    print("-" * 88)
    print("V2 (I2)/(I3): measured C3 on |z| <= (4/3) eps_strip; "
          "cubic Taylor remainder |r3(z)| <= (C3 n/6)|z|^3")
    ok = True
    packs = []
    for P0 in packs_raw:
        C3m = measure_C3(rng, P0)
        C3 = 1.05 * C3m  # 5% measurement headroom -> the admissible constant
        P = make_pack_from(P0, C3)
        # (I3) check: r3 along rays to sampled z in the rectangle region
        n = 512
        words = shell_words(rng, n, P.p, 3)
        worst_ratio = 0.0
        for x in words:
            t = int(math.floor(n * P.D))
            _, mu, v, _, _ = word_qt_mu_v(x, P.p, P.D, t)
            for z in (P.eps2 * 0.7 + 0.3j * P.eps_strip,
                      -P.eps2 + 0.45j * P.eps_strip,
                      P.eps2 + 0.0j,
                      0.5 * P.eps2 - 0.4j * P.eps_strip,
                      1.25 * P.eps_strip + 0.0j):
                L = logPhi_ray(x, P.p, P.D, z)
                r3 = L - (1j * mu * z - v * z * z / 2)
                bound = (P.C3 * n / 6) * abs(z) ** 3
                worst_ratio = max(worst_ratio, abs(r3) / bound)
        good = worst_ratio <= 1.0
        ok &= good
        print(f"  p={P.p}, D={P.D:.5f}: C3_meas={C3m:.4f} (admissible C3={P.C3:.4f}); "
              f"max |r3|/[(C3 n/6)|z|^3] = {worst_ratio:.3f} <= 1: {good}")
        packs.append(P)
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok, packs

def make_pack_from(P0, C3):
    P = Pack()
    P.__dict__.update(P0.__dict__)
    P.C3 = C3
    P.eps2 = min(P.eps_strip, 3 * P.vbar / (32 * C3))
    P.eps_split = min(P.delta_KP_lb, P.vbar / (4 * C3))
    gbar = 1 - (11 / 16) * P.D * (1 - P.D) * (1 - math.cos(P.eps_split))
    P.c_tail = 0.25 * math.log(1 / gbar)
    f = (2.0 / P.vbar) ** 1.5
    P.Ca = (8 * math.e * C3 / 3) / math.sqrt(2 * math.pi) * f
    P.Cb = (2 * math.e * C3 / 3) * f
    P.C = P.Ca + P.Cb + 5.0
    P.J = max(1, math.ceil(math.log2(P.eps_split / P.eps2)))
    return P

# =============================== V3 + V4 ===============================

def small_beta_words(rng, P, n, count, cap, tries=3000):
    out = []
    t = int(math.floor(n * P.D))
    k = 0
    while len(out) < count and k < tries:
        k += 1
        x = shell_words(rng, n, P.p, 1)[0]
        qt, mu, v, phi, s = word_qt_mu_v(x, P.p, P.D, t)
        beta = (t - mu) / math.sqrt(v)
        y0 = (mu - t) / v
        if abs(y0) <= cap and qt > 0:
            out.append((x, qt, mu, v, beta, y0, phi, s))
    return out, t

def deformed_pieces(x, P, n, t, mu, v, y0, n_real=24000):
    """Central, verticals, real-zone integrals of Phi(z) e^{-izt}."""
    # central (shifted)
    sg = np.linspace(-P.eps2, P.eps2, 2001)
    z = sg + 1j * y0
    f = Phi_vals(x, P.p, P.D, z) * np.exp(-1j * z * t)
    central = np.trapezoid(f, sg)
    # verticals
    yg = np.linspace(0.0, y0, 201)
    zm = -P.eps2 + 1j * yg
    zp = P.eps2 + 1j * yg
    fm = Phi_vals(x, P.p, P.D, zm) * np.exp(-1j * zm * t)
    fp = Phi_vals(x, P.p, P.D, zp) * np.exp(-1j * zp * t)
    Vm = 1j * np.trapezoid(fm, yg)         # -eps2: 0 -> y0
    Vp = -1j * np.trapezoid(fp, yg)        # +eps2: y0 -> 0
    # real zones [eps2, pi] both sides
    sr = np.linspace(P.eps2, math.pi, n_real // 2)
    fr_p = Phi_vals(x, P.p, P.D, sr) * np.exp(-1j * sr * t)
    fr_m = Phi_vals(x, P.p, P.D, -sr) * np.exp(1j * sr * t)
    realz = np.trapezoid(fr_p + fr_m, sr)
    return central, Vm, Vp, realz, (sg, f)

def run_V3_V4(rng, packs):
    print("-" * 88)
    print("V3 (E1): contour deformation identity;  V4 (E3)/(E5)/(E5')/(E8) "
          "+ (E7)/(N2) arithmetic")
    okV3 = True
    okV4 = True
    for P in packs:
        for n in (512, 1024):
            cap = min(P.eps2, P.eps_strip / 2)
            words, t = small_beta_words(rng, P, n, 3 if n == 512 else 2, cap)
            if not words:
                okV3 = False
                print(f"  p={P.p} n={n}: NO small-|y0| words found (cap={cap:.4f}) -> FAIL")
                continue
            for (x, qt, mu, v, beta, y0, phi, s) in words:
                central, Vm, Vp, realz, (sg, fc) = deformed_pieces(
                    x, P, n, t, mu, v, y0)
                total = (central + Vm + Vp + realz).real / (2 * math.pi)
                rel = abs(total - qt) / abs(qt)
                c_ok = rel < 2e-3
                okV3 &= c_ok
                # (E3)
                gauss = np.trapezoid(np.exp(-v * sg ** 2 / 2), sg)
                TG = math.sqrt(2 * math.pi / v) - gauss
                TG_bound = (2 / (v * P.eps2)) * math.exp(-v * P.eps2 ** 2 / 2)
                e3_ok = -1e-12 <= TG <= TG_bound * (1 + 1e-9)
                # (E5): E_in over |s|<=sigma_tilde (= eps2 at accessible n)
                sig_t = min(sigma_n(P, n), P.eps2)
                mask = np.abs(sg) <= sig_t + 1e-15
                Ein = (np.trapezoid(fc[mask], sg[mask]) * math.exp(beta * beta / 2)
                       - np.trapezoid(np.exp(-v * sg[mask] ** 2 / 2), sg[mask]))
                rhs5 = (2 * math.e * P.C3 * n / 3) * (
                    4 / v ** 2 + abs(y0) ** 3 * math.sqrt(2 * math.pi / v))
                e5_ok = abs(Ein) <= rhs5
                # |r3| <= 1 on the inner zone (lemma uses <=1 there)
                r3max = float(np.max(np.abs(
                    np.log(fc[mask] * np.exp(beta * beta / 2 + v * sg[mask] ** 2 / 2)))))
                r3_ok = r3max <= 1.0
                # (E5') arithmetic given v >= vbar n/2
                e5p_ok = (v >= P.vbar * n / 2) and (
                    rhs5 / math.sqrt(2 * math.pi / v)
                    <= P.Ca / math.sqrt(n) + P.Cb * abs(beta) ** 3 / math.sqrt(n)
                    + 1e-12)
                # (E8): vertical sup + favorable-sign profile
                yg = np.linspace(0.0, y0, 201)
                zp = P.eps2 + 1j * yg
                lmv, phv = logPhi_grid(x, P.p, P.D, zp)
                log_int = lmv + (t * yg)          # log |Phi e^{-izt}|
                bound8 = -v * P.eps2 ** 2 / 4
                e8_ok = float(np.max(log_int)) <= bound8 + 1e-9
                pred = (-v * y0 * yg + v * yg ** 2 / 2 - v * P.eps2 ** 2 / 2)
                tol = (P.C3 * n / 6) * (P.eps2 + np.abs(yg)) ** 3
                prof_ok = bool(np.all(np.abs(log_int - pred) <= tol + 1e-9))
                okV4 &= e3_ok and e5_ok and r3_ok and e5p_ok and e8_ok and prof_ok
                print(f"  p={P.p} n={n}: beta={beta:+.3f} y0={y0:+.5f} | "
                      f"E1 rel={rel:.1e} {c_ok} | TG<=bnd {e3_ok} | "
                      f"|Ein|={abs(Ein):.2e}<={rhs5:.2e} {e5_ok} | r3max={r3max:.1e}<=1 "
                      f"{r3_ok} | (E5') {e5p_ok} | vert sup {e8_ok} profile {prof_ok}")
        # (E7)/(N2): outer zone empty at accessible n; arithmetic at large n
        n_cross = 3.0 / (4 * P.C3 * P.eps2 ** 3)
        n_chk = 16 * n_cross
        sn = sigma_n(P, n_chk)
        lhsN2 = (math.log(8.0 / (sn * math.sqrt(math.pi * P.vbar * n_chk)))
                 - (P.vbar / 8) * n_chk * sn ** 2)
        n2_ok = lhsN2 <= -0.5 * math.log(n_chk)
        okV4 &= n2_ok
        print(f"  p={P.p}: outer zone empty for n < {n_cross:.2e} (sigma_n>eps2); "
              f"(N2) at n={n_chk:.1e}: log-margin = "
              f"{lhsN2 + 0.5 * math.log(n_chk):+.1f} <= 0: {n2_ok}")
    print(f"  -> V3 {'PASS' if okV3 else 'FAIL'}; V4 {'PASS' if okV4 else 'FAIL'}")
    return okV3, okV4

# =============================== V5 + V6 + V7 ===============================

def run_V5_V6_V7(rng, packs):
    print("-" * 88)
    print("V5 (E9/E9a/E9b/I5): zone-2' band bound + net constants;  "
          "V6 (E11): zone-3 tail;  V7: final |R_n| bound")
    okV5 = True
    okV6 = True
    okV7 = True
    for P in packs:
        sband = np.linspace(P.eps2, P.eps_split, 160)
        stail = np.linspace(P.eps_split, math.pi, 320)
        C1_eff = 0.0
        impliedC = 0.0
        rhat_by_n = {}
        shell_rms = {}
        for n, nw in ((256, 60), (512, 40), (1024, 16), (2048, 10)):
            t = int(math.floor(n * P.D))
            n_shell = nw
            words = shell_words(rng, n, P.p, nw) + \
                biased_words(rng, n, 0.8 * P.p, max(2, nw // 8))
            band_ok = True
            tail_ok = True
            Rs = []
            m2_band = np.zeros(len(sband))
            cnt = 0
            for iw, x in enumerate(words):
                qt, mu, v, phi, s = word_qt_mu_v(x, P.p, P.D, t)
                lm_b, _ = logPhi_grid(x, P.p, P.D, sband)
                lm_t, _ = logPhi_grid(x, P.p, P.D, stail)
                # (E9) per-word band bound
                band_ok &= bool(np.all(lm_b <= -P.c4 * n * sband ** 2 + 1e-9))
                # (E11) zone-3
                tail_ok &= bool(np.max(lm_t) <= math.log(2.0) - P.c_tail * n + 1e-9)
                m2_band += np.exp(2 * np.minimum(lm_b, 0.0))
                cnt += 1
                if qt > 0 and v > 0:
                    beta = (t - mu) / math.sqrt(v)
                    if abs(beta) <= P.B3 * math.sqrt(math.log(n)):
                        Rn = qt * math.sqrt(2 * math.pi * v) * \
                            math.exp(beta * beta / 2) - 1.0
                        Rs.append((Rn, beta, iw < n_shell))
            # (I5) annealed prefactor on the band
            m2_band /= cnt
            gb = g_closed(sband, P.p, P.D)
            C1_eff = max(C1_eff, float(np.max(m2_band / gb ** (n - 1))))
            okV5 &= band_ok
            okV6 &= tail_ok
            # V7: per-word bound + the n-invariant Rhat = |R| sqrt(n)/(1+|b|)^3
            Rarr = np.array([r for r, _, _ in Rs])
            barr = np.array([b for _, b, _ in Rs])
            sh = np.array([s for _, _, s in Rs])
            bnd = P.C * (1 + np.abs(barr)) ** 3 / math.sqrt(n)
            v7 = bool(np.all(np.abs(Rarr) <= bnd))
            okV7 &= v7
            rhats = np.abs(Rarr) * math.sqrt(n) / (1 + np.abs(barr)) ** 3
            rhat_by_n[n] = float(np.max(rhats))
            impliedC = max(impliedC, rhat_by_n[n])
            shell_rms[n] = float(np.sqrt(np.mean(Rarr[sh] ** 2)))
            print(f"  p={P.p} n={n:>4}: (E9) band {band_ok} | (E11) tail {tail_ok} | "
                  f"V7 |R|<=C(1+|b|)^3/sqrt(n) [C={P.C:.1f}] {v7}  "
                  f"max|R|={np.max(np.abs(Rarr)):.4f} maxRhat={rhat_by_n[n]:.3f} "
                  f"shellRMS={shell_rms[n]:.4f}")
        # net constants / margins / Jordan
        m1 = P.c_g / math.pi ** 2 - 2 * P.c2      # per-point margin (E9a)
        m2m = P.c_g / math.pi ** 2 - 3 * P.c2     # blockwise margin (E9b)
        sj = np.linspace(1e-6, math.pi, 4000)
        jordan = bool(np.all(1 - np.cos(sj) >= 2 * sj ** 2 / math.pi ** 2 - 1e-12))
        prefac_ok = C1_eff <= 2.0
        margins_ok = m1 > 0 and m2m > 0 and P.c4 == P.c2 / 4 and P.eps2 <= P.eps_split
        okV5 &= jordan and prefac_ok and margins_ok
        # P(E_2') crossing
        nn = 2.0 ** np.arange(8, 200)
        pe = (math.log(12 * 1.0 * P.J * P.eps_split) + np.log(nn)
              - P.c5 * nn * P.eps2 ** 2)   # with C1=1 measured grade
        n_e2p = float(nn[np.argmax(pe < math.log(1e-3))]) if np.any(
            pe < math.log(1e-3)) else float('inf')
        # n-invariance of Rhat (the committed R4 criterion: NOT raw RMS,
        # which mixes the beta-distribution across n) + shell RMS decay
        rhat_stab = max(rhat_by_n.values()) <= \
            3 * min(rhat_by_n.values()) + 0.05
        shell_dec = shell_rms[2048] < shell_rms[256]
        okV7 &= rhat_stab and shell_dec
        print(f"  p={P.p}: C1_eff={C1_eff:.3f}<=2 {prefac_ok} | margins "
              f"c_g/pi^2-2c2={m1:.2e}>0, -3c2={m2m:.2e}>0 {margins_ok} | "
              f"Jordan {jordan} | P(E_2')<1e-3 from n~{n_e2p:.1e} | "
              f"impliedC={impliedC:.3f} (explicit C={P.C:.1f}) | "
              f"Rhat n-stable {rhat_stab} | shell-RMS decay {shell_dec}")
    print(f"  -> V5 {'PASS' if okV5 else 'FAIL'}; V6 {'PASS' if okV6 else 'FAIL'}; "
          f"V7 {'PASS' if okV7 else 'FAIL'}")
    return okV5, okV6, okV7

# =============================== V8 ===============================

def n0_conditions(P, n):
    """Each (N*) as 'log-LHS <= log-RHS' margins; True iff condition holds."""
    ln = math.log(n)
    sn = sigma_n(P, n)
    y0max = P.B3 * math.sqrt(ln) * math.sqrt(2.0 / (P.vbar * n))
    N0 = y0max <= min(sn, P.eps2, P.eps_strip / 2)
    N1 = (math.log(2 / (P.eps2 * math.sqrt(math.pi * P.vbar * n)))
          - (P.vbar / 4) * n * P.eps2 ** 2) <= -0.5 * ln
    if sn < P.eps2:
        N2 = (math.log(8 / (sn * math.sqrt(math.pi * P.vbar * n)))
              - (P.vbar / 8) * n * sn ** 2) <= -0.5 * ln
    else:
        N2 = True  # outer zone empty
    N3 = (math.log(2 * P.B3 * math.sqrt(ln) / math.sqrt(2 * math.pi))
          + (P.B3 ** 2 / 2) * ln - (P.vbar / 8) * n * P.eps2 ** 2) <= -0.5 * ln
    N4 = (math.log(1 / (P.c4 * P.eps2 * 2 * math.sqrt(2 * math.pi)))
          + (P.B3 ** 2 / 2) * ln - P.c4 * n * P.eps2 ** 2) <= -0.5 * ln
    N5 = (math.log(4 * math.pi / (2 * math.sqrt(2 * math.pi))) + ln
          + (P.B3 ** 2 / 2) * ln - P.c_tail * n) <= -0.5 * ln
    return dict(N0=N0, N1=N1, N2=N2, N3=N3, N4=N4, N5=N5)

def run_V8(packs):
    print("-" * 88)
    print("V8: explicit n0(p,D) from (N0)-(N5) + n_b (log-domain scan)")
    ok = True
    for P in packs:
        cross = {}
        for key in ("N0", "N1", "N2", "N3", "N4", "N5"):
            nc = None
            for k in range(8, 300):
                n = 2.0 ** k
                if n0_conditions(P, n)[key]:
                    # require persistence at 4x and 16x
                    if n0_conditions(P, 4 * n)[key] and \
                       n0_conditions(P, 16 * n)[key]:
                        nc = n
                        break
            cross[key] = nc
        # n_b: McDiarmid E3 exponent >= 3 ln n
        nb = None
        for k in range(4, 300):
            n = 2.0 ** k
            lam = P.B3 * math.sqrt(P.vbar / 2) * math.sqrt(n * math.log(n)) - 1
            if lam > 0 and 2 * lam ** 2 / ((n - 1) * P.c_mu ** 2) >= \
                    3 * math.log(n):
                nb = n
                break
        finite = all(v is not None for v in cross.values()) and nb is not None
        ok &= finite
        n0 = max([v for v in cross.values()] + [nb, 5.0 / P.D, 2.0])
        # all conditions simultaneously at n0 and 4 n0
        all_at = all(n0_conditions(P, n0).values()) and \
            all(n0_conditions(P, 4 * n0).values())
        ok &= all_at
        cs = ", ".join(f"{k}~{v:.1e}" for k, v in cross.items())
        print(f"  p={P.p}, D={P.D:.5f}: crossings {cs}, n_b~{nb:.1e}")
        print(f"    => n0(p,D) ~ {n0:.2e} (finite: {finite}; all (N*) hold at "
              f"n0 and 4n0: {all_at})  [honest: asymptotic regime]")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# =============================== V9 ===============================

def run_V9(rng, packs):
    print("-" * 88)
    print("V9: pitfall regressions (i) two scales (ii) t|y0|=Theta(sqrt n) "
          "(iii) constant orderings")
    ok = True
    for P in packs:
        # (i) per-site rate at eps2 vs at eps_split (conflation killer)
        n = 1024
        t = int(math.floor(n * P.D))
        rates2 = []
        ratesS = []
        ty0 = {}
        for nn in (256, 1024):
            tt = int(math.floor(nn * P.D))
            vals = []
            for x in shell_words(rng, nn, P.p, 6):
                _, mu, v, _, _ = word_qt_mu_v(x, P.p, P.D, tt)
                vals.append(tt * abs(mu - tt) / v)
                if nn == 1024:
                    lm, _ = logPhi_grid(x, P.p, P.D,
                                        np.array([P.eps2, P.eps_split]))
                    rates2.append(-lm[0] / nn)
                    ratesS.append(-lm[1] / nn)
            ty0[nn] = float(np.mean(vals))
        r2 = float(np.mean(rates2))
        rS = float(np.mean(ratesS))
        sep = rS / max(r2, 1e-300)
        i_ok = sep > 10.0
        # (ii) t|y0| grows ~ sqrt(n) (ratio over 4x n should be ~2, >1.2)
        growth = ty0[1024] / max(ty0[256], 1e-30)
        ii_ok = growth > 1.2
        # (iii) orderings
        iii_ok = (P.eps2 <= P.eps_strip <= P.delta_KP_lb / 8 + 1e-15
                  and P.eps2 <= P.eps_split <= math.pi
                  and 4 * P.A * 2 * P.eps_strip <= P.r_lb + 1e-15
                  and math.sqrt(P.eps2 ** 2 + (P.eps_strip / 2) ** 2)
                  <= (4 / 3) * P.eps_strip
                  and P.vbar >= 0.75 * P.D * (1 - P.D) - 1e-12
                  and P.c_tail > 0 and P.c4 > 0 and P.C3 > 0)
        ok &= i_ok and ii_ok and iii_ok
        print(f"  p={P.p}: rate(eps_split)/rate(eps2) = {sep:.1f} > 10 {i_ok} "
              f"(tail bound NEVER invoked at eps2) | t|y0|: {ty0[256]:.2f} -> "
              f"{ty0[1024]:.2f} grows {ii_ok} | orderings {iii_ok}")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok

# =============================== main ===============================

if __name__ == "__main__":
    import sys
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 42
    print("=" * 88)
    print("GAP-1 write-out v2 -- per-inequality verification "
          "(/tmp/gap1_lemma_v2.tex)")
    print(f"points: (p, D) in {{(0.25, 0.9 Dc), (0.4, 0.5 Dc)}}; "
          f"n up to 2048; seed={seed}")
    print("=" * 88)
    rng = np.random.default_rng(seed)
    packs_raw = [make_pack(0.25, 0.9), make_pack(0.4, 0.5)]
    res = {}
    res["V0"] = run_V0(rng, packs_raw)
    res["V1"] = run_V1(rng)
    res["V2"], packs = run_V2(rng, packs_raw)
    for P in packs:
        print(f"  [constants] p={P.p} D={P.D:.5f}: eps_strip={P.eps_strip:.5f} "
              f"eps2={P.eps2:.5f} eps_split={P.eps_split:.5f} C3={P.C3:.4f} "
              f"vbar={P.vbar:.5f} c_tail={P.c_tail:.2e} c4={P.c4:.2e} "
              f"B3={P.B3:.2f} Ca={P.Ca:.1f} Cb={P.Cb:.1f} C={P.C:.1f}")
    res["V3"], res["V4"] = run_V3_V4(rng, packs)
    res["V5"], res["V6"], res["V7"] = run_V5_V6_V7(rng, packs)
    res["V8"] = run_V8(packs)
    res["V9"] = run_V9(rng, packs)
    print("=" * 88)
    for k, v in res.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    allok = all(res.values())
    print(f"VERDICT: {'gap1-writeout-v2-verified' if allok else 'gap1-writeout-v2-ISSUE'}")

#!/usr/bin/env python3
"""
probe_7_34_replica_assembly_central.py
============================================================================
ASSEMBLY probe for the 7.34b replica/second-moment route (companion of
probe_7_34_replica_secondmoment_tail.py): tests the CENTRAL-WINDOW half of
the assembled quenched argument

  log2 q_t(x) = -1/2 log2(2 pi v_x) - beta^2/2 * log2(e) + R_n(x),
      beta := (t - mu_x)/sqrt(v_x),  t = floor(nD),
      mu_x = E_Q[K|x],  v_x = Var_Q(K|x)   (posterior cumulants),

with the claim: R_n(x) -> 0 quenched-uniformly over shell words (the
steepest-descent error, polylog/sqrt(n)), Var(log2 q_t) = O(1) dominated by
the moving-center beta^2 term, and the Theta(n)-cancellation NEVER reappears
(it is performed exactly by the normalisation Phi(0;x)=1).

Also reports, per word, sup_{s in [delta_KP, pi]} |Phi(s;x)| (the quantity the
replica tail lemma bounds) for consistency with e^{-cn}.

CHECKS
  A0  DFT pipeline (log-scaled 2x2 transfer matrix at the n+1 Fourier nodes)
      == brute FWHT posterior at n=12 (q_k exact to ~1e-12).
  A1  per-word Gaussian/steepest-descent residual R_n(x): max and RMS over
      shell samples, at n in {256, 512, 1024}; expect RMS ~ n^{-1/2} polylog
      and max -> 0.
  A2  Var(log2 q_t) over shell samples is O(1) (n-stable), and matches
      Var(-beta^2/2 log2 e - 1/2 log2 v) to o(1); per-site posterior variance
      vbar = E[v_x]/n > 0 stable; sigma_mu^2 = Var(mu_x)/n > 0 (so the beta^2
      term genuinely fluctuates at Theta(1) -- the allowed O(1) budget).
  A3  tail consistency: max over shell words of sup_{[delta_KP,pi]}|Phi| at
      each n, vs the lemma's e^{-cn} with c = (1/4) ln(1/g_max).
"""
import math
import numpy as np

# ---------------- shared BSMS / dual objects (as in the verified probes) ----------------

def Dc(p):
    return 0.5 * (1 - math.sqrt(1 - 2 * p) / (1 - p))

def theta0(D):
    return math.log(D / (1 - D))

def eta_C(beta, D):
    eb = np.exp(beta)
    ze = (1 - eb) / ((1 + eb) * (1 - 2 * D))
    eta = (1 - ze) / (1 + ze)
    C = (1 + eb) * (1 + ze) / 2
    return eta, C

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

def delta_KP(p, D):
    rkp = r_cluster(p)
    grid = np.concatenate([np.geomspace(1e-5, 0.1, 600, endpoint=False),
                           np.linspace(0.1, math.pi, 2400)])
    th = theta0(D)
    etas_abs = np.abs([eta_C(th + 1j * s, D)[0] for s in grid])
    cross = np.argmax(etas_abs > rkp) if np.any(etas_abs > rkp) else None
    return math.pi if cross is None else float(grid[max(cross - 1, 0)])

# ---------------- log-scaled transfer evaluation of Phi at the DFT nodes ----------------

def phi_all_nodes(xbits, p, D):
    """Phi(s_j; x) at s_j = 2 pi j/(n+1), j=0..n, via the verified 2x2 TM,
    log-magnitude-scaled (|C_s/C_0|^n can overflow; |Phi|<=1 cannot)."""
    n = len(xbits)
    S = n + 1
    s = 2 * math.pi * np.arange(S) / S
    th = theta0(D)
    eta, Cs = eta_C(th + 1j * s, D)
    C0 = eta_C(th, D)[1]
    rho = p / (1 - p)
    # w = (1, eta) row vector; w <- w @ A_j, A_j = [[1, eta r],[r, eta]]
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
    """q_k = (1/S) sum_j Phi(s_j) e^{-i s_j k}; exact since K in {0..n}."""
    q = np.fft.fft(phi).real / len(phi)   # fft computes sum_j a_j e^{-2pi i jk/S}
    return q

# ---------------- brute cross-check at n=12 ----------------

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

def check_A0(rng):
    print("-" * 84)
    print("A0: DFT/log-TM pipeline == brute FWHT posterior q_k (n=12)")
    n = 12
    ok = True
    worst = 0.0
    for p in (0.25, 0.4):
        D = 0.9 * Dc(p)
        N = 1 << n
        P = bsms_P(n, p)
        Ph = fwht(P)
        pc_all = popcount_arr(np.arange(N))
        PYs = fwht(Ph * (1 - 2 * D) ** (-pc_all.astype(float))) / N
        th = theta0(D)
        ix = np.arange(N)
        for _ in range(4):
            x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
            xx = 0
            for b in x:
                xx = (xx << 1) | int(b)
            d = pc_all[ix ^ xx]
            w = PYs * np.exp(th * d)
            A = np.bincount(d, weights=w, minlength=n + 1)
            q_brute = A / A.sum()
            phi, _ = phi_all_nodes(x, p, D)
            q_dft = q_from_phi(phi)
            err = np.max(np.abs(q_dft - q_brute))
            worst = max(worst, err)
            ok &= err < 1e-10
    print(f"  worst abs.diff over q_k = {worst:.2e}   -> {'PASS' if ok else 'FAIL'}")
    return ok

# ---------------- main assembly experiment ----------------

def shell_words(rng, n, p, n_words):
    """Sample BSMS words conditioned on the sqrt(n) switch-count shell."""
    out = []
    mean = (n - 1) * p
    half = 2.0 * math.sqrt((n - 1) * p * (1 - p))
    while len(out) < n_words:
        x = (np.cumsum(rng.random(n) < p) % 2).astype(int)
        sw = int(np.sum(x[:-1] != x[1:]))
        if abs(sw - mean) <= half:
            out.append(x)
    return out

def run_combo(rng, p, fD, n_list, n_words, gmax_lookup):
    D = fD * Dc(p)
    dkp = delta_KP(p, D)
    log2e = 1.0 / math.log(2.0)
    print("-" * 84)
    print(f"combo p={p}, D={fD}*Dc={D:.5f}, delta_KP={dkp:.4f}")
    rows = []
    for n in n_list:
        t = int(math.floor(n * D))
        L = []      # log2 q_t
        Lpred = []  # -1/2 log2(2 pi v) - beta^2/2 log2 e
        Rres = []
        betas = []
        mus = []
        vs = []
        tails = []
        for x in shell_words(rng, n, p, n_words):
            phi, s = phi_all_nodes(x, p, D)
            q = q_from_phi(phi)
            ksum = float(np.sum(q))
            q = np.maximum(q, 0.0)
            ks = np.arange(n + 1, dtype=float)
            mu = float(ks @ q[:n + 1])
            v = float((ks - mu) ** 2 @ q[:n + 1])
            beta = (t - mu) / math.sqrt(v)
            lq = math.log2(max(q[t], 1e-300))
            pred = -0.5 * math.log2(2 * math.pi * v) - 0.5 * beta * beta * log2e
            mask = (s >= dkp) & (s <= math.pi)
            tails.append(float(np.max(np.abs(phi[mask]))))
            L.append(lq); Lpred.append(pred); Rres.append(lq - pred)
            betas.append(beta); mus.append(mu); vs.append(v)
            assert abs(ksum - 1.0) < 1e-8
        L = np.array(L); Lpred = np.array(Lpred); Rres = np.array(Rres)
        betas = np.array(betas); mus = np.array(mus); vs = np.array(vs)
        gmax = gmax_lookup
        c4 = 0.25 * math.log(1.0 / gmax)
        rows.append((n, np.var(L), np.var(Lpred), np.max(np.abs(Rres)),
                     np.sqrt(np.mean(Rres ** 2)), np.mean(vs) / n,
                     np.var(mus) / n, np.var(betas), max(tails),
                     math.exp(-c4 * n)))
        print(f"  n={n:>5} t={t:>4}: Var(log2 q_t)={rows[-1][1]:.4f}  "
              f"Var(pred)={rows[-1][2]:.4f}  max|R|={rows[-1][3]:.4f}  "
              f"RMS(R)={rows[-1][4]:.4f}")
        print(f"          vbar=E[v]/n={rows[-1][5]:.4f}  sigma_mu^2=Var(mu)/n={rows[-1][6]:.4f}  "
              f"Var(beta)={rows[-1][7]:.4f}  max tail sup|Phi|={rows[-1][8]:.3e} "
              f"(lemma e^(-n/4 ln 1/gmax)={rows[-1][9]:.3e})")
    # verdicts
    r_shrinks = rows[-1][4] < rows[0][4] and rows[-1][3] < rows[0][3]
    var_O1 = all(r[1] < 3.0 for r in rows) and \
        abs(rows[-1][1] - rows[0][1]) < 0.5 * max(rows[0][1], 0.05) + 0.05
    pred_match = all(abs(r[1] - r[2]) < 0.15 * max(r[1], 0.05) + 4 * r[4] ** 2 + 0.02
                     for r in rows)
    vbar_pos = all(r[5] > 1e-3 for r in rows)
    smu_pos = all(r[6] > 1e-4 for r in rows)
    tail_dec = rows[-1][8] < rows[0][8]
    print(f"  A1 residual shrinks with n: {'PASS' if r_shrinks else 'FAIL'} | "
          f"A2 Var(log2 q_t)=O(1) & matches Gaussian model: "
          f"{'PASS' if (var_O1 and pred_match) else 'FAIL'} | "
          f"vbar>0,sigma_mu^2>0: {'PASS' if (vbar_pos and smu_pos) else 'FAIL'} | "
          f"A3 tail decays: {'PASS' if tail_dec else 'FAIL'}")
    return r_shrinks and var_O1 and pred_match and vbar_pos and smu_pos and tail_dec

if __name__ == "__main__":
    print("=" * 84)
    print("7.34b assembly probe: central-window Gaussian form + O(1) variance budget")
    print("=" * 84)
    rng = np.random.default_rng(11)
    a0 = check_A0(rng)
    # g_max values measured in probe_7_34_replica_secondmoment_tail.py C3
    ok1 = run_combo(rng, 0.25, 0.9, [256, 512, 1024], 150, 1.0 - 4.86e-3)
    ok2 = run_combo(rng, 0.10, 0.5, [256, 512, 1024], 150, 1.0 - 1.28e-3)
    print("=" * 84)
    allok = a0 and ok1 and ok2
    print(f"A0 {'PASS' if a0 else 'FAIL'} | combos {'PASS' if (ok1 and ok2) else 'FAIL'}")
    print(f"VERDICT: {'central-window-assembly-verified' if allok else 'assembly-issue'}")

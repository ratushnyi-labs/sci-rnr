#!/usr/bin/env python3
"""
probe_7_34_cocycle_aperiodicity.py
============================================================================
APERIODICITY / COHOMOLOGICAL NON-ARITHMETICITY of the integer Hamming cocycle
    K_n(x,Y*) = sum_i 1[Y*_i != x_i]
of the (tilted) BSMS forward-filter process -- the EXACT spectral condition that
Aaronson-Denker 2001 (Gibbs-Markov LLT), Guivarc'h-Hardy, Herve-Pene 2010
(Nagaev-Guivarc'h via Keller-Liverani, lattice local limit under aperiodicity),
and Breuillard 2005 need IN PLACE OF a density floor / uniform ellipticity for a
NEAR-MEAN lattice local-limit theorem.

WHY THIS IS THE DECISIVE TEST (task A/B)
----------------------------------------
The 7.34b achievability residual reduces (saddlepoint, probe_..._saddlepoint_logS_llt.py,
GREEN gate) to a NEAR-MEAN density-floor-free LATTICE local-limit theorem for K_n,
uniform over the sqrt(n)-typical shell of schedules x. The literature dichotomy
(step iv-b) blocks the FIXED-FRACTION third-order (1/2 log n) lattice LLT on
ellipticity. The proposed escape: a NEAR-MEAN lattice LLT needs only
  (i) a Holder spectral gap (filter stability gives it; Kloeckner-Lopes-Stadlbauer
      arXiv:1412.0848), AND
  (ii) COHOMOLOGICAL NON-ARITHMETICITY of the lattice cocycle (aperiodicity),
which REPLACES the density floor (Aaronson-Denker "big images" / Herve-Pene).
The ONLY thing that could still kill this near-mean route is if condition (ii)
FAILS: if the integer cocycle is ARITHMETIC (periodic) for the BSMS filter, the
twisted transfer operator L_{is} has a peripheral eigenvalue of modulus 1 at some
s in (0,2pi), the lattice LLT carries a nontrivial periodic phase, and r_t would
NOT converge to a single constant. (The saddlepoint probe shows r_t -> 0.94
STABLY, so we expect aperiodicity -- this script CHECKS it from the operator.)

THE EXACT OBJECT (derived from the SLB-tight Bayes identity, verified to 1e-14)
------------------------------------------------------------------------------
With theta0 = ln(D/(1-D)) < 0 and the Gray-region deconvolution X = Y* (+) iid-Bern(D):
  posterior weight on y proportional to  P_{Y*}(y) e^{theta0 K(x,y)}.
The (posterior) characteristic function of K is the TWISTED-TRANSFER RATIO
  Phi_n(s; x) = E_Q[e^{is K}] = M(theta0 + i s; x) / M(theta0; x),
  M(beta; x) = E_{Y*}[ e^{beta K(x,Y*)} ] = sum_y P_{Y*}(y) prod_i e^{beta 1[y_i != x_i]}.
M(beta;x) is EXACTLY a 2x2 TRANSFER-MATRIX product over the output chain Y* (the
BSC(D) deconvolution of the BSMS source). The twisted transfer operator family is
  L_beta : per-letter twist e^{beta} on the flipped symbol; its top eigenvalue
  lambda(beta) governs M_n(beta) ~ C lambda(beta)^n.
APERIODICITY (the Aaronson-Denker/Herve-Pene condition) =
  |lambda(theta0 + i s)| < |lambda(theta0)| for all s in (0, 2pi),
  with equality ONLY at s in 2pi Z (lattice span exactly 1).
Equivalently the per-schedule, per-n decay rate
  rho_n(s; x) := |Phi_n(s; x)|^{1/n}  ->  |lambda(theta0+is)|/|lambda(theta0)| < 1.

THE OUTPUT CHAIN Y* (exact, no 2^n)
-----------------------------------
Y* is NOT finite-state Markov (rank-2 OOM, growing-distinct conditionals), so there
is no exact finite 2x2 twisted matrix. We compute M_n(beta;x) EXACTLY for complex
beta via the Walsh/Krawtchouk transfer machinery already validated in
probe_7_34_sign_var_logS.py:
  M(beta;x) = sum_y P_{Y*}(y) e^{beta K(x,y)}
            = sum_y P_{Y*}(y) prod_i (1 + (e^beta - 1) 1[y_i != x_i]).
  Let u = e^beta - 1. Expanding, group by Hamming weight pattern:
  M = sum_{k} u^k * (sum over k-subsets S: P_{Y*}( y = x except flipped on S )).
  Using the Walsh representation P_{Y*}(y) = 2^{-n} sum_w Phat_X(w)(1-2D)^{-|w|}(-1)^{<w,y>},
  and 1[y_i != x_i] selection, the exact closed form is the SAME generating identity:
  M(beta;x) = E_{Y*}[ prod_i g_i ],  g_i = e^beta if Y*_i != x_i else 1,
  which equals the z=... evaluation of the SAME Q_x machinery. Concretely
  (see _M_beta below) M(beta;x) = 2^{-n} sum_j (1-2D)^{-j} H_j(x) * c_j(beta),
  c_j(beta) = sum over outputs at Walsh-weight j of the per-letter twist =
            = the Krawtchouk-twist coefficient; we instead compute M DIRECTLY by a
  complex 2x2 transfer-matrix product for the OUTPUT law P_{Y*} (deconvolution
  realized as a complex transfer with the (1-2D)^{-1} eigen-tilt), which is exact.

We use the DIRECT exact route: M(beta;x) via the output-side DP over (Y*-state, twist),
where the output law P_{Y*} is the BSC(D)-deconvolution; the deconvolution is applied
in the Walsh domain so the DP is over j (Walsh weight). This is identical machinery to
ball_prob_formula but with the twist weight c_j(beta) replacing the Krawtchouk ball sum.

Run:   /tmp/rnr_venv/bin/python scripts/verify/probe_7_34_cocycle_aperiodicity.py [p]
Deps:  numpy, mpmath
"""
import math
import sys
import cmath
import numpy as np
import mpmath as mp


def h2(x):
    if x <= 0 or x >= 1:
        return 0.0
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


def D_c(p):
    return 0.5 * (1.0 - math.sqrt(1.0 - 2.0 * p) / (1.0 - p))


def V_lossless(p):
    return p * (1 - p) * (math.log2((1 - p) / p)) ** 2


# --------------------------------------------------------------------------
# H_j(x): coefficients of Q_x(z) (real), from probe_7_34_sign_var_logS.py
# --------------------------------------------------------------------------
def _mul_deg1(c, s):
    m = len(c)
    o = [None] * (m + 1)
    o[0] = c[0]
    for k in range(1, m):
        o[k] = c[k] + s * c[k - 1]
    o[m] = s * c[m - 1]
    return o


def H_coeffs(x, p):
    """Real coeffs [H_0..H_n] of Q_x(z)=sum_{x'} P_X(x') prod_i ((1+z) if x'_i==x_i else (1-z))."""
    n = len(x)
    half = mp.mpf('0.5')
    x0 = int(x[0])
    v = [_mul_deg1([half], 1 if x0 == 0 else -1),
         _mul_deg1([half], 1 if x0 == 1 else -1)]
    P = mp.mpf(p); Qq = 1 - P
    for i in range(1, n):
        xi = int(x[i])
        vb, va = v[0], v[1]
        L = len(vb)
        a0 = [Qq * vb[k] + P * va[k] for k in range(L)]
        a1 = [P * vb[k] + Qq * va[k] for k in range(L)]
        v = [_mul_deg1(a0, 1 if xi == 0 else -1),
             _mul_deg1(a1, 1 if xi == 1 else -1)]
    Q = [v[0][k] + v[1][k] for k in range(len(v[0]))]
    while len(Q) < n + 1:
        Q.append(mp.mpf(0))
    return Q[:n + 1]


# --------------------------------------------------------------------------
# EXACT posterior law q_k(x) = P(K = k | X = x)  (real), from saddlepoint probe.
#   pi_k(x) = P_{Y*}(d_H(x,Y*) = k)  [untilted output-side weight enumerator]
#   q_k     propto pi_k e^{theta0 k}
# Then M(beta;x) = sum_k pi_k e^{beta k}  (the twisted transfer = output-side MGF of K),
#   so Phi_n(s;x) = M(theta0+is)/M(theta0) = sum_k q_k e^{isk}  (posterior char. fn).
# This is EXACT and avoids any 2x2 finite-state approximation of the non-Markov Y*.
# --------------------------------------------------------------------------
def kraw_full(n):
    from math import comb
    K = [[0] * (n + 1) for _ in range(n + 1)]
    for j in range(n + 1):
        for k in range(n + 1):
            s = 0
            lo = max(0, k - (n - j)); hi = min(k, j)
            for l in range(lo, hi + 1):
                s += ((-1) ** l) * comb(j, l) * comb(n - j, k - l)
            K[k][j] = s
    return K


def pi_k(x, p, D, Km):
    """Exact pi_k = P_{Y*}(d_H(x,Y*)=k), k=0..n.  (mpmath real)"""
    n = len(x)
    H = H_coeffs(x, p)
    inv = mp.mpf(1) / (1 - 2 * mp.mpf(D))
    ip = [mp.mpf(1)] * (n + 1)
    for j in range(1, n + 1):
        ip[j] = ip[j - 1] * inv
    out = [mp.mpf(0)] * (n + 1)
    for k in range(n + 1):
        s = mp.mpf(0)
        for j in range(n + 1):
            s += ip[j] * Km[k][j] * H[j]
        out[k] = s / (mp.mpf(2) ** n)
    return out


def sample(n, p, rng):
    st = (rng.random(n) < p).astype(np.int8)
    st[0] = (rng.random() < 0.5)
    return np.cumsum(st) % 2


# --------------------------------------------------------------------------
# CORE: per-schedule twisted-transfer ratio and posterior characteristic function.
#   M(beta;x) = sum_k pi_k e^{beta k}   (beta complex)
#   Phi_n(s;x) = M(theta0 + i s)/M(theta0) = E_Q[e^{i s K}]
# --------------------------------------------------------------------------
def char_fn_and_ratemod(pis, theta0, s_vals):
    """Return per-s complex Phi_n(s)=M(theta0+is)/M(theta0) and modulus, exact-from-pi_k.
    pis: list of mpmath reals length n+1. Returns numpy complex array over s_vals."""
    n = len(pis) - 1
    ks = np.arange(n + 1)
    pivec = np.array([float(v) for v in pis], dtype=np.float64)  # pi_k can be tiny but >=0
    # M(theta0) real
    w0 = pivec * np.exp(theta0 * ks)
    M0 = w0.sum()
    out = np.empty(len(s_vals), dtype=np.complex128)
    for idx, s in enumerate(s_vals):
        # M(theta0+is) = sum_k pi_k e^{theta0 k} e^{i s k} = sum_k w0_k e^{i s k}
        Ms = np.sum(w0 * np.exp(1j * s * ks))
        out[idx] = Ms / M0
    return out


# --------------------------------------------------------------------------
# Per-schedule SPECTRAL-RADIUS proxy (the Aaronson-Denker / Herve-Pene quantity):
#   rho_n(s;x) = |Phi_n(s;x)|^{1/n_eff}, n_eff = posterior std support ~ effective length.
# We report the modulus |Phi_n(s;x)| (the twisted-cocycle char. fn) directly AND its
# per-letter decay |Phi_n|^{1/n}. Aperiodic <=> |Phi_n(s)| -> 0 for s in (0,2pi),
# = 1 only at s = 2pi (lattice span 1). Peripheral peak structure is the lattice span.
# --------------------------------------------------------------------------
def analyze(p, D, n, R, seed, s_grid, dps=80):
    rng = np.random.default_rng(seed)
    Km = kraw_full(n)
    th0 = float(mp.log(mp.mpf(D) / (1 - mp.mpf(D))))
    t = int(math.floor(n * D))
    nD = n * D
    # accumulate |Phi_n(s)| per schedule; also the ANNEALED M ratio (sum over schedules
    # of M(theta0+is) weighted by the schedule's own M(theta0), = E_shell with the proper
    # near-mean operative weight).
    mod_list = []           # per-schedule |Phi_n(s)| arrays
    Mnum = np.zeros(len(s_grid), dtype=np.complex128)  # sum_x M(th0+is;x) (annealed numerator)
    Mden = 0.0
    vqs = []
    kept = 0
    for _ in range(R):
        x = sample(n, p, rng)
        with mp.workdps(dps):
            pis = pi_k(x, p, D, Km)
            # operative tilt weights w0_k = pi_k e^{th0 k}
            ks = np.arange(n + 1)
            w0 = np.array([float(pis[k]) * math.exp(th0 * k) for k in range(n + 1)])
        M0 = w0.sum()
        if M0 <= 0:
            continue
        q = w0 / M0
        mq = float((ks * q).sum())
        vq = float(((ks - mq) ** 2 * q).sum())
        if vq <= 0:
            continue
        # per-schedule char fn
        phi = np.array([np.sum(w0 * np.exp(1j * s * ks)) / M0 for s in s_grid])
        mod_list.append(np.abs(phi))
        # annealed numerator/denominator (weight each schedule by M0 so it's E over the
        # operative posterior-normalizer measure)
        Mnum += np.array([np.sum(w0 * np.exp(1j * s * ks)) for s in s_grid])
        Mden += M0
        vqs.append(vq)
        kept += 1
    if kept < 10:
        return None
    mods = np.array(mod_list)               # (kept, len(s_grid))
    mean_mod = mods.mean(axis=0)
    max_mod = mods.max(axis=0)              # worst-case (uniform-over-shell) modulus
    annealed = np.abs(Mnum / Mden)         # |E_shell M(th0+is)|/|E_shell M(th0)|
    vq = float(np.mean(vqs))
    return dict(n=n, t=t, nD=nD, kept=kept, vq=vq, th0=th0,
                mean_mod=mean_mod, max_mod=max_mod, annealed=annealed)


def main():
    p = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4
    D = 0.9 * D_c(p)
    th0 = math.log(D / (1 - D))
    print("#" * 100)
    print(f"# BSMS p={p}  D=0.9 D_c={D:.5f}  V_lossless={V_lossless(p):.5f}  theta0={th0:.3f}")
    print(f"# APERIODICITY of the integer Hamming cocycle K_n -- the Aaronson-Denker / Herve-Pene")
    print(f"# lattice-LLT condition that REPLACES the density floor in the NEAR-MEAN regime.")
    print(f"# Object: Phi_n(s;x)=E_Q[e^{{isK}}]=M(theta0+is;x)/M(theta0;x), exact from pi_k(x).")
    print(f"# Aperiodic  <=>  |Phi_n(s)| -> 0 for s in (0,2pi),  =1 only at s in 2pi Z (span 1).")
    print("#" * 100)

    # dual-lattice phase grid: dense near 0 (Gaussian zone), across (0,2pi), peak at 2pi.
    s_grid = np.concatenate([
        np.linspace(0.0, 0.6, 7),          # Gaussian zone near 0
        np.array([math.pi / 2, math.pi, 1.5 * math.pi]),  # interior
        np.linspace(2 * math.pi - 0.6, 2 * math.pi + 0.6, 13),  # peripheral peak window
    ])
    s_grid = np.unique(np.round(s_grid, 6))
    labels = {0.0: "0", round(math.pi / 2, 6): "pi/2", round(math.pi, 6): "pi",
              round(1.5 * math.pi, 6): "3pi/2", round(2 * math.pi, 6): "2pi"}

    sched = [(64, 300), (96, 250), (128, 200), (160, 160),
             (220, 120), (300, 90), (400, 60), (520, 45)]
    rows = []
    # pick a few diagnostic s values to print per n (interior + peripheral)
    s_show = [round(math.pi / 2, 6), round(math.pi, 6), round(1.5 * math.pi, 6),
              round(2 * math.pi, 6)]
    hdr = f"{'n':>4} {'kept':>4} {'vq':>7} | " + " ".join(
        f"|Phi|@{labels.get(s, f'{s:.2f}'):>5}" for s in s_show)
    print(hdr)
    print("-" * len(hdr))
    for n, R in sched:
        t = int(math.floor(n * D))
        if t < 8:
            continue
        r = analyze(p, D, n, R, seed=7000 + n, s_grid=s_grid)
        if r is None:
            print(f"{n:>4}  few")
            continue
        rows.append(r)
        idxs = [int(np.argmin(np.abs(s_grid - s))) for s in s_show]
        vals = " ".join(f"{r['max_mod'][i]:>10.5f}" for i in idxs)
        print(f"{n:>4} {r['kept']:>4} {r['vq']:>7.2f} | {vals}", flush=True)

    print("-" * 100)
    if len(rows) < 3:
        print("INSUFFICIENT rows.")
        return
    last = rows[-1]
    ns = np.array([r['n'] for r in rows])

    # ----- (1) INTERIOR DECAY: |Phi_n(s)| -> 0 for s in (0,2pi) ? rate in n -----
    print("APERIODICITY DIAGNOSTICS:")
    print("(1) INTERIOR decay |Phi_n(s)| -> 0 for s in (0,2pi) (worst-case over shell = max_mod):")
    interior_s = [s for s in s_grid if 0.3 < s < 2 * math.pi - 0.3]
    spans = []
    for s in interior_s:
        i = int(np.argmin(np.abs(s_grid - s)))
        seq = np.array([r['max_mod'][i] for r in rows])
        # fit log|Phi| ~ -rate * n  (geometric decay => aperiodic, rate>0)
        good = seq > 1e-14
        if good.sum() >= 3:
            rate = -np.polyfit(ns[good], np.log(seq[good]), 1)[0]
        else:
            rate = float('inf')
        spans.append(rate)
        lbl = labels.get(round(s, 6), f"{s:.2f}")
        print(f"    s={lbl:>5}: |Phi| n=[{rows[0]['n']}->{rows[-1]['n']}] "
              f"{seq[0]:.4f}->{seq[-1]:.2e}   per-letter decay rate={rate:+.4f}/n")
    interior_decays = all(r > 1e-3 for r in spans)
    print(f"    ==> ALL interior s decay geometrically (aperiodic, no periodic phase): "
          f"{'YES' if interior_decays else 'NO'}")

    # ----- (2) PERIPHERAL PEAK at s=2pi only (lattice span exactly 1) -----
    print("(2) PERIPHERAL structure: |Phi_n(s)| near s=2pi (lattice span check):")
    peri = [s for s in s_grid if abs(s - 2 * math.pi) <= 0.6]
    i2pi = int(np.argmin(np.abs(s_grid - 2 * math.pi)))
    seq2pi = np.array([r['max_mod'][i2pi] for r in rows])
    # at exactly 2pi: e^{i 2pi k}=1 for integer k => |Phi|=1 EXACTLY (span 1 lattice). Check.
    print(f"    |Phi_n(2pi)| (should be ~1, the lattice span-1 peak): "
          f"{seq2pi[0]:.6f} -> {seq2pi[-1]:.6f}")
    # immediate neighbourhood: does it dip below 1 as soon as s != 2pi ? (no sub-lattice)
    near = []
    for s in peri:
        if abs(s - 2 * math.pi) < 1e-6:
            continue
        i = int(np.argmin(np.abs(s_grid - s)))
        near.append((s, last['max_mod'][i]))
    near_max = max(v for _, v in near) if near else 0.0
    print(f"    max |Phi_n(s)| for s in [2pi-0.6,2pi+0.6]\\{{2pi}} at n={last['n']}: {near_max:.5f} "
          f"(<1 => isolated peak, span exactly 1)")
    span1 = abs(seq2pi[-1] - 1.0) < 1e-3 and near_max < 0.99
    print(f"    ==> lattice span exactly 1, isolated peripheral peak at 2pi: {'YES' if span1 else 'NO'}")

    # ----- (3) NO sub-peak in (0,2pi) (no arithmetic sub-lattice d|1, trivial, but verify min) -----
    print("(3) NO interior peripheral eigenvalue (no nontrivial periodic phase):")
    midmax = max(last['max_mod'][int(np.argmin(np.abs(s_grid - s)))]
                 for s in s_grid if 0.3 < s < 2 * math.pi - 0.3)
    print(f"    max |Phi_n(s)| over interior (0.3, 2pi-0.3) at n={last['n']}: {midmax:.5f} "
          f"(<<1 => no interior peripheral eigenvalue)")
    no_interior_peak = midmax < 0.5
    print(f"    ==> no interior modulus-1 phase: {'YES' if no_interior_peak else 'NO'}")

    # ----- (4) GAUSSIAN ZONE near 0 (the near-mean regime sanity): |Phi| ~ exp(-vq s^2/2) -----
    print("(4) GAUSSIAN zone near s=0 (near-mean): |Phi_n(s)| ~ exp(-vq s^2 / 2)?")
    for s in [s for s in s_grid if 0 < s <= 0.6][:4]:
        i = int(np.argmin(np.abs(s_grid - s)))
        pred = math.exp(-last['vq'] * s * s / 2)
        print(f"    s={s:.3f}: |Phi|={last['max_mod'][i]:.5f}  exp(-vq s^2/2)={pred:.5f}")

    print("-" * 100)
    verdict = interior_decays and span1 and no_interior_peak
    print("VERDICT (aperiodicity / cohomological non-arithmeticity of the integer Hamming cocycle):")
    if verdict:
        print("  APERIODIC CONFIRMED. |Phi_n(s)| decays geometrically for ALL s in (0,2pi), with an")
        print("  isolated peripheral peak |Phi_n(2pi)|=1 (lattice span exactly 1) and NO interior")
        print("  modulus-1 phase. This is EXACTLY the Aaronson-Denker / Herve-Pene non-arithmeticity")
        print("  hypothesis. => the lattice cocycle is NOT arithmetic; barrier (B-i) is REFUTED;")
        print("  the near-mean lattice LLT has NO periodic-phase obstruction. LINK-HOLDS for the")
        print("  aperiodicity input. (Whether the near-mean relaxation truly removes the density-floor")
        print("  NEED -- barrier (B-ii) -- is the SEPARATE analytic question argued in the report.)")
    else:
        print("  BARRIER FOUND: aperiodicity FAILS at the printed s -- the cocycle carries a periodic")
        print("  phase => lattice LLT obstructed; r_t cannot converge to a single constant.")
        print("  This would BLOCK the near-mean route with a precise arithmetic barrier.")


if __name__ == "__main__":
    main()

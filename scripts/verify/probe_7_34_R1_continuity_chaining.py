#!/usr/bin/env python3
"""
probe_7_34_R1_continuity_chaining.py
============================================================================
TRACK (R1): is r_t(x)=q_t(x) sqrt(2 pi v_q) BOUNDED above/below and LOCALLY
LIPSCHITZ, UNIFORMLY over the source-typical sqrt(n)-shell of schedules x --
and does that quenched/uniform statement REDUCE to (annealed boundedness) +
(R2 single-flip stability) + a CONTINUITY/CHAINING argument, or does it need a
genuinely new QUENCHED lattice LLT?

The annealed statement (E_x[r_t]~0.94, r_t -> bounded const in distribution) is
the Aaronson-Denker 2001 / Gouezel 2010 finite-variance lattice LLT (Gibbs-Markov
+ aperiodicity, density-floor-free) -- ALREADY established (aperiodicity verified
to n=520, twisted radius<1 on (0,2pi)).

The QUENCHED/uniform upgrade route (a) [CONTINUITY/CHAINING] needs the map
x |-> log2 r_t(x) to be Holder/Lipschitz with a GEOMETRICALLY-LOCALIZED modulus,
i.e. single-flip influence D_i(log r_t) DECAYING geometrically with the distance
of coordinate i from the chain's "active region" (filter stability => two-sided
geometric coordinate forgetting). IF that localization holds, a covering/chaining
over the shell (or a McDiarmid/Azuma sup-control with summable Lipschitz weights
c_i, sum c_i^2 = O(1)) upgrades annealed boundedness to SHELL-UNIFORM boundedness.

THIS PROBE tests the LOAD-BEARING localization in the GENUINE NEAR-MEAN REGIME
(large n, t=floor(nD) LARGE, NOT the frozen-radius t in {1,2} artifact that the
small-n influence probe is trapped in):

 (T-A) RIGHT-DISTANCE influence profile c_d = mean|D_i log2 r_t| at right-distance
       d=n-1-i, in the near-mean regime (n up to ~140, t>=8). Does c_d -> 0
       geometrically in d?  (Filter stability would give this.)
 (T-B) BULK vs the LOCALIZED-active influence: is the influence concentrated near
       the SWITCH positions of x (the chain's active structure), or spread
       uniformly across all coordinates (=> global non-additive object, no
       chaining)?  Report Gini/participation-ratio of the influence vector.
 (T-C) The DECISIVE chaining quantity: the SUM of squared single-flip influences
       sum_i c_i^2 where c_i is the WORST-CASE (over the shell, not averaged)
       single-flip change. McDiarmid/Azuma needs sum_i c_i^2 = O(1). Report it.
 (T-D) Direct CONTINUITY test: take two shell schedules x, x' at Hamming distance
       h and measure |log2 r_t(x)-log2 r_t(x')| vs h. Holder/Lipschitz <=> bounded
       slope; report the empirical Lipschitz constant and whether it is n-stable.

VERDICT LOGIC:
  * If (T-A) c_d decays geometrically AND (T-C) sum_i c_i^2 = O(1) n-stably AND
    (T-D) Lipschitz constant n-stable: (R1) REDUCES to annealed+(R2)+continuity;
    the quenched upgrade is a chaining/McDiarmid corollary -- VERDICT PROVES-route.
  * If (T-A) is FLAT (no geometric decay) or (T-C) sum_i c_i^2 = Theta(n): the
    worst-case single-flip influence is NOT summable; chaining/McDiarmid only gives
    Var=O(n); the quenched upgrade needs a genuinely new quenched lattice LLT
    (sequential Nagaev-Guivarch, route (b)) -- VERDICT OPEN.

Exact O(n^2) local mass via Walsh/Krawtchouk (same machinery as
probe_7_34_saddlepoint_logS_llt.py / probe_7_34_localmass_influence_regularity.py).

Run:  /tmp/rnr_venv/bin/python scripts/verify/probe_7_34_R1_continuity_chaining.py [p]
Deps: numpy, mpmath
"""
import math, sys
import numpy as np
import mpmath as mp
from math import comb


def h2(x):
    if x <= 0 or x >= 1:
        return 0.0
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


def D_c(p):
    return 0.5 * (1.0 - math.sqrt(1.0 - 2.0 * p) / (1.0 - p))


def V_lossless(p):
    return p * (1 - p) * (math.log2((1 - p) / p)) ** 2


def _mul_deg1(c, s):
    m = len(c)
    o = [None] * (m + 1)
    o[0] = c[0]
    for k in range(1, m):
        o[k] = c[k] + s * c[k - 1]
    o[m] = s * c[m - 1]
    return o


def H_coeffs(x, p):
    n = len(x)
    half = mp.mpf('0.5')
    x0 = int(x[0])
    v = [_mul_deg1([half], 1 if x0 == 0 else -1),
         _mul_deg1([half], 1 if x0 == 1 else -1)]
    P = mp.mpf(p)
    Qq = 1 - P
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


_KCACHE = {}
def kraw_full(n):
    if n in _KCACHE:
        return _KCACHE[n]
    K = [[0] * (n + 1) for _ in range(n + 1)]
    for j in range(n + 1):
        for k in range(n + 1):
            s = 0
            lo = max(0, k - (n - j)); hi = min(k, j)
            for l in range(lo, hi + 1):
                s += ((-1) ** l) * comb(j, l) * comb(n - j, k - l)
            K[k][j] = s
    _KCACHE[n] = K
    return K


def pi_k_vec(x, p, D, Km):
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


def logr_t(x, p, D, Km, dps=70):
    """Return log2 r_t = log2(q_t sqrt(2 pi v_q)) at the boundary t=floor(nD),
    tilted q_k propto pi_k e^{theta0 k}.  None if degenerate."""
    n = len(x)
    t = int(math.floor(n * D))
    th0 = float(mp.log(mp.mpf(D) / (1 - mp.mpf(D))))
    with mp.workdps(dps):
        pis = pi_k_vec(x, p, D, Km)
        M = sum(pis[k] * mp.e ** (mp.mpf(th0) * k) for k in range(n + 1))
        if M <= 0:
            return None
        q = [float(pis[k] * mp.e ** (mp.mpf(th0) * k) / M) for k in range(n + 1)]
    if q[t] <= 0:
        return None
    mq = sum(k * q[k] for k in range(n + 1))
    vq = sum((k - mq) ** 2 * q[k] for k in range(n + 1))
    if vq <= 0:
        return None
    rt = q[t] * math.sqrt(2 * math.pi * vq)
    if rt <= 0:
        return None
    return math.log2(rt)


def sample(n, p, rng):
    st = (rng.random(n) < p).astype(np.int8)
    st[0] = (rng.random() < 0.5)
    return np.cumsum(st) % 2


def near_mean_ok(x, p, D, Km, dps=70):
    """Keep only schedules whose posterior mean mq ~ nD (the near-mean regime),
    so we are NOT measuring the frozen-radius artifact."""
    n = len(x); t = int(math.floor(n * D)); nD = n * D
    th0 = float(mp.log(mp.mpf(D) / (1 - mp.mpf(D))))
    with mp.workdps(dps):
        pis = pi_k_vec(x, p, D, Km)
        M = sum(pis[k] * mp.e ** (mp.mpf(th0) * k) for k in range(n + 1))
        if M <= 0:
            return False
        q = [float(pis[k] * mp.e ** (mp.mpf(th0) * k) / M) for k in range(n + 1)]
    mq = sum(k * q[k] for k in range(n + 1))
    return abs(mq - nD) < 1.0 + 0.1 * nD


# ---------------------------------------------------------------- (T-A)+(T-B)+(T-C)
def influence_analysis(p, n, R, seed, Dfrac=0.9, dps=70):
    D = Dfrac * D_c(p)
    Km = kraw_full(n)
    rng = np.random.default_rng(seed)
    prof_sum = np.zeros(n); prof_cnt = np.zeros(n)   # mean|D_i| by right-distance
    energies = []                                    # sum_i (D_i)^2 per word (averaged-DIRICHLET)
    worst_c = np.zeros(n)                            # WORST-CASE |D_i| over the shell, per right-distance
    part_ratios = []                                 # participation ratio of the |D_i| vector
    switch_corr = []                                 # corr(|D_i|, near-switch indicator)
    kept = 0
    for _ in range(R):
        x = sample(n, p, rng)
        if not near_mean_ok(x, p, D, Km, dps):
            continue
        base = logr_t(x, p, D, Km, dps)
        if base is None:
            continue
        Di = np.zeros(n)
        ok = True
        for i in range(n):
            xf = x.copy(); xf[i] ^= 1
            lf = logr_t(xf, p, D, Km, dps)
            if lf is None:
                ok = False; break
            Di[i] = base - lf
        if not ok:
            continue
        kept += 1
        a = np.abs(Di)
        energies.append(float(np.sum(Di ** 2)))
        for i in range(n):
            d = n - 1 - i
            prof_sum[d] += a[i]; prof_cnt[d] += 1
            worst_c[d] = max(worst_c[d], a[i])
        # participation ratio: (sum a)^2 / (n * sum a^2) in [1/n,1]; ~1 => spread, ~0 => localized
        s1 = a.sum(); s2 = (a ** 2).sum()
        if s2 > 0:
            part_ratios.append((s1 ** 2) / (n * s2))
        # switch correlation: is influence near switch positions x_i!=x_{i-1}?
        sw = np.zeros(n); sw[1:] = (x[1:] != x[:-1]).astype(float)
        # smear switch indicator by +-1 (active region)
        act = sw.copy()
        act[:-1] = np.maximum(act[:-1], sw[1:]); act[1:] = np.maximum(act[1:], sw[:-1])
        if a.std() > 0 and act.std() > 0:
            switch_corr.append(np.corrcoef(a, act)[0, 1])
    if kept < 5:
        return None
    prof = np.divide(prof_sum, np.maximum(prof_cnt, 1))
    energies = np.array(energies)
    return dict(n=n, kept=kept, t=int(math.floor(n * D)),
                prof=prof, worst_c=worst_c,
                E=energies.mean(), E_se=energies.std() / math.sqrt(kept),
                sumworst2=float(np.sum(worst_c ** 2)),
                part_ratio=float(np.mean(part_ratios)) if part_ratios else float('nan'),
                switch_corr=float(np.mean(switch_corr)) if switch_corr else float('nan'))


# ---------------------------------------------------------------- (T-D) continuity slope
def continuity_slope(p, n, R, seed, Dfrac=0.9, dps=70, max_h=6):
    """Pairs of shell schedules at Hamming distance h: |log2 r_t(x)-log2 r_t(x')| vs h.
    Empirical Lipschitz constant = mean(|delta|/h). n-stable & bounded => Holder/Lipschitz."""
    D = Dfrac * D_c(p)
    Km = kraw_full(n)
    rng = np.random.default_rng(seed)
    by_h = {h: [] for h in range(1, max_h + 1)}
    tries = 0
    while tries < R:
        x = sample(n, p, rng)
        if not near_mean_ok(x, p, D, Km, dps):
            tries += 1; continue
        lx = logr_t(x, p, D, Km, dps)
        if lx is None:
            tries += 1; continue
        # flip h random coordinates
        h = rng.integers(1, max_h + 1)
        idx = rng.choice(n, size=h, replace=False)
        xp = x.copy(); xp[idx] ^= 1
        # keep x' near-mean too (so we stay in regime)
        if not near_mean_ok(xp, p, D, Km, dps):
            tries += 1; continue
        lxp = logr_t(xp, p, D, Km, dps)
        if lxp is None:
            tries += 1; continue
        by_h[int(h)].append(abs(lx - lxp))
        tries += 1
    rows = []
    for h in range(1, max_h + 1):
        v = np.array(by_h[h])
        if len(v) >= 5:
            rows.append((h, len(v), v.mean(), v.mean() / h))
    return rows


def main():
    p = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4
    D = 0.9 * D_c(p)
    print("#" * 100)
    print(f"# BSMS p={p}  D=0.9 D_c={D:.5f}  V_lossless={V_lossless(p):.5f}")
    print(f"# (R1) CONTINUITY/CHAINING route: does shell-uniform boundedness of r_t reduce to")
    print(f"#   annealed + (R2 single-flip stability) + a Holder-continuity/chaining argument?")
    print(f"#   Tested in the GENUINE NEAR-MEAN regime (t>=8, NOT the frozen-radius t in {{1,2}} artifact).")
    print("#" * 100)

    print("\n(T-A)+(T-B)+(T-C) single-flip influence of log2 r_t, near-mean schedules only:")
    print(f"{'n':>4} {'t':>3} {'kept':>5} | {'E=avg sum Di^2':>14} | {'sum worst_c^2':>13} | "
          f"{'partRatio':>9} {'swCorr':>7} | right-distance worst_c[d] (d=0..6)")
    rows = []
    for n in (60, 80, 100, 120, 140):
        t = int(math.floor(n * D))
        if t < 6:
            continue
        R = 220 if n <= 100 else 140
        r = influence_analysis(p, n, R, seed=4100 + n)
        if r is None:
            print(f"{n:>4} (insufficient near-mean)"); continue
        rows.append(r)
        wc = r['worst_c']
        wcstr = " ".join(f"{wc[d]:.3f}" for d in range(0, min(7, n)))
        print(f"{n:>4} {r['t']:>3} {r['kept']:>5} | {r['E']:>8.4f}+-{r['E_se']:.3f} | "
              f"{r['sumworst2']:>13.4f} | {r['part_ratio']:>9.3f} {r['switch_corr']:>+7.3f} | {wcstr}",
              flush=True)

    if len(rows) >= 3:
        ns = np.array([r['n'] for r in rows], float)
        E = np.array([r['E'] for r in rows])
        sw2 = np.array([r['sumworst2'] for r in rows])
        pr = np.array([r['part_ratio'] for r in rows])
        sE = np.polyfit(ns, E, 1)[0]
        ssw = np.polyfit(ns, sw2, 1)[0]
        # geometric decay of the worst-case profile for the largest n
        wc = rows[-1]['worst_c']
        ratios = [wc[d + 1] / wc[d] for d in range(0, 6) if wc[d] > 1e-9]
        print("-" * 100)
        print(f"(T-A) worst-case right-distance profile decay ratios (largest n): "
              f"{', '.join(f'{x:.2f}' for x in ratios)}")
        decays = all(x < 0.85 for x in ratios) if len(ratios) >= 3 else False
        print(f"      >>> worst-case influence GEOMETRICALLY LOCALIZED in distance: "
              f"{'YES' if decays else 'NO (flat => global object)'}")
        print(f"(T-B) participation ratio (1=spread, 0=localized): "
              f"mean over n = {pr.mean():.3f}  => influence is "
              f"{'LOCALIZED' if pr.mean() < 0.4 else 'SPREAD across all coordinates'}")
        print(f"(T-C) averaged Dirichlet energy E[sum Di^2] slope vs n = {sE:+.5f}  "
              f"(BOUNDED if ~0): {'BOUNDED' if abs(sE) < 0.01 else 'GROWS Theta(n)'}")
        print(f"      WORST-CASE chaining sum_i c_i^2 slope vs n = {ssw:+.5f}  "
              f"(McDiarmid/Azuma needs O(1)): {'O(1) -- summable' if abs(ssw) < 0.02 else 'Theta(n) -- NOT summable'}")

    print("\n(T-D) CONTINUITY: |log2 r_t(x)-log2 r_t(x')| vs Hamming distance h (near-mean pairs):")
    for n in (80, 120):
        t = int(math.floor(n * D))
        if t < 6:
            continue
        print(f"  n={n} (t={t}):")
        rows_h = continuity_slope(p, n, 600, seed=5100 + n)
        for h, cnt, md, slope in rows_h:
            print(f"    h={h} cnt={cnt:>4} mean|delta|={md:.4f} mean|delta|/h={slope:.4f}", flush=True)

    print()
    print("VERDICT LOGIC (R1):")
    print("  PROVES-route (R1 reduces to annealed+(R2)+continuity) IFF:")
    print("    (T-A) worst-case influence geometrically localized in distance, AND")
    print("    (T-C) worst-case chaining sum_i c_i^2 = O(1) n-stable, AND")
    print("    (T-D) Lipschitz constant mean|delta|/h bounded & n-stable.")
    print("  OPEN (needs new quenched lattice LLT, route (b) sequential Nagaev-Guivarch) IF:")
    print("    (T-A) flat / (T-C) sum c_i^2 = Theta(n): worst-case single-flip influence NOT")
    print("    summable => McDiarmid/chaining only gives Var=O(n), not the O(1) (R1) needs.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
probe_7_34_localmass_nearmean_influence.py
============================================================================
DECISIVE measurement for the QUENCHED-UNIFORM route (Remark 7.34b achievability):
the coordinate-influence ENERGY of log2 q_t in the GENUINE NEAR-MEAN regime
(t = floor(nD) >= 8, escaping the frozen-radius artifact), where the LLT zone
actually lives.

The route (annealed Aaronson-Denker local CLT + concentration of q_t(X)) CLOSES
only if log2 q_t(x) is a bounded-influence / geometrically-forgetting Holder
functional of the schedule x -- i.e. (i) influence energy E=sum_i (D_i log q_t)^2
is O(1) (n-independent) and (ii) per-coordinate influence DECAYS with distance.
Then McDiarmid/Kontorovich-Ramanan gives Var(log q_t)=O(1)=o(n), and annealed +
concentration => the uniform-over-shell LLT.

The route FAILS (same wall as the prior ball-CDF Efron-Stein route) if E=Theta(n)
or the influence does NOT decay: then the small observed Var(log q_t) is a
RADIUS-ASYMPTOTIC (t->infty) cancellation that NO concentration inequality sees.

This probe pairs, at the SAME (n, near-mean t):
  E    = E[ sum_i (D_i log2 q_t)^2 ]            (the McDiarmid/Efron-Stein energy)
  Vlq  = Var_x( log2 q_t(X) )                   (the quantity to bound)
and the Efron-Stein INEQUALITY says Vlq <= (1/2) E. If E=Theta(n) while Vlq=O(1),
the inequality is SLACK by Theta(n) -- proof that concentration over-counts and
the true smallness is a cancellation, not concentration. Also reports the
per-coordinate influence profile by distance (geometric decay test).

Exact O(n^2) tilted local mass via Walsh/Krawtchouk transfer matrix.
Run:  python scripts/verify/probe_7_34_localmass_nearmean_influence.py [p]
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


_KC = {}
def kraw_full(n):
    if n in _KC:
        return _KC[n]
    K = [[0] * (n + 1) for _ in range(n + 1)]
    for j in range(n + 1):
        for k in range(n + 1):
            s = 0
            lo = max(0, k - (n - j)); hi = min(k, j)
            for l in range(lo, hi + 1):
                s += ((-1) ** l) * comb(j, l) * comb(n - j, k - l)
            K[k][j] = s
    _KC[n] = K
    return K


def logq_t(x, p, D, Km, th0, t, dps=60):
    n = len(x)
    H = H_coeffs(x, p)
    inv = mp.mpf(1) / (1 - 2 * mp.mpf(D))
    ip = [mp.mpf(1)] * (n + 1)
    for j in range(1, n + 1):
        ip[j] = ip[j - 1] * inv
    # only need pi_k; compute all (O(n^2))
    pis = [mp.mpf(0)] * (n + 1)
    for k in range(n + 1):
        s = mp.mpf(0)
        row = Km[k]
        for j in range(n + 1):
            s += ip[j] * row[j] * H[j]
        pis[k] = s / (mp.mpf(2) ** n)
    M = mp.mpf(0)
    for k in range(n + 1):
        M += pis[k] * mp.e ** (mp.mpf(th0) * k)
    if M <= 0:
        return None
    qt = pis[t] * mp.e ** (mp.mpf(th0) * t) / M
    if qt <= 0:
        return None
    return float(mp.log(qt, 2))


def sample(n, p, rng):
    st = (rng.random(n) < p).astype(np.int8)
    st[0] = (rng.random() < 0.5)
    return np.cumsum(st) % 2


def run(p, n, R_E, R_V, seed, Dfrac=0.9, dps=60):
    D = Dfrac * D_c(p)
    Km = kraw_full(n)
    th0 = float(mp.log(mp.mpf(D) / (1 - mp.mpf(D))))
    t = int(math.floor(n * D))
    rng = np.random.default_rng(seed)
    # ---- Vlq over many samples (cheap: one logq per sample) ----
    lqs = []
    for _ in range(R_V):
        x = sample(n, p, rng)
        with mp.workdps(dps):
            v = logq_t(x, p, D, Km, th0, t, dps)
        if v is not None:
            lqs.append(v)
    lqs = np.array(lqs)
    Vlq = lqs.var(ddof=1) if len(lqs) > 5 else float('nan')
    # ---- E over fewer samples (n flips each: expensive) ----
    energies = []
    prof_sum = np.zeros(n); prof_cnt = np.zeros(n)
    for _ in range(R_E):
        x = sample(n, p, rng)
        with mp.workdps(dps):
            lq0 = logq_t(x, p, D, Km, th0, t, dps)
        if lq0 is None:
            continue
        Di = np.zeros(n); ok = True
        for i in range(n):
            xf = x.copy(); xf[i] ^= 1
            with mp.workdps(dps):
                lf = logq_t(xf, p, D, Km, th0, t, dps)
            if lf is None:
                ok = False; break
            Di[i] = lq0 - lf
        if not ok:
            continue
        energies.append(float(np.sum(Di ** 2)))
        for i in range(n):
            d = n - 1 - i
            prof_sum[d] += abs(Di[i]); prof_cnt[d] += 1
    energies = np.array(energies)
    prof = np.divide(prof_sum, np.maximum(prof_cnt, 1))
    return dict(n=n, t=t, nD=n * D, kE=len(energies), kV=len(lqs),
                E=energies.mean() if len(energies) else float('nan'),
                E_se=energies.std() / math.sqrt(max(len(energies), 1)),
                Vlq=Vlq, prof=prof)


def main():
    p = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4
    D = 0.9 * D_c(p)
    print("#" * 100)
    print(f"# BSMS p={p}  D=0.9 D_c={D:.5f}  V_lossless={V_lossless(p):.5f}")
    print(f"# NEAR-MEAN influence energy of log2 q_t (t=floor(nD)>=8) vs direct Var_x(log2 q_t).")
    print(f"# Efron-Stein: Var(log2 q_t) <= (1/2) E.  Route needs E=O(1) (=> concentration).")
    print(f"# If E=Theta(n) while Var=O(1): inequality slack by Theta(n) => smallness is a")
    print(f"#   radius-asymptotic CANCELLATION, NOT concentration (route FAILS, same wall).")
    print("#" * 100)
    print(f"{'n':>4} {'t':>3} {'nD':>6} {'kE':>4} {'kV':>5} | {'E=sumDi^2':>11} {'E/n':>7} | "
          f"{'Var(lq)':>9} {'ES gap E/2-Var':>13} | per-coord |D_i| by dist d (0..6)")
    rows = []
    # near-mean schedule: pick n so floor(nD)>=8; keep n<=128 for tractability of the flip loop
    for n in (72, 88, 104, 120, 136):
        t = int(math.floor(n * D))
        if t < 8:
            continue
        R_E = 60 if n <= 104 else 36
        R_V = 600 if n <= 104 else 300
        r = run(p, n, R_E, R_V, seed=5000 + n)
        if r is None or not (r['kE'] > 3):
            print(f"{n:>4} (insufficient)"); continue
        rows.append(r)
        prof = r['prof']
        ds = " ".join(f"{prof[d]:.3f}" for d in range(0, min(7, n)))
        es_gap = 0.5 * r['E'] - r['Vlq']
        print(f"{r['n']:>4} {r['t']:>3} {r['nD']:>6.2f} {r['kE']:>4} {r['kV']:>5} | "
              f"{r['E']:>7.3f}+-{r['E_se']:.2f} {r['E']/r['n']:>7.4f} | "
              f"{r['Vlq']:>9.4f} {es_gap:>13.3f} | {ds}", flush=True)
    if len(rows) >= 3:
        ns = np.array([r['n'] for r in rows], float)
        E = np.array([r['E'] for r in rows])
        Vlq = np.array([r['Vlq'] for r in rows])
        sE = np.polyfit(ns, E, 1)[0]
        sV = np.polyfit(ns, Vlq, 1)[0]
        # influence-decay over the largest n
        prof = rows[-1]['prof']
        ratios = [prof[d + 1] / prof[d] for d in range(0, 5) if prof[d] > 1e-9]
        print("-" * 100)
        print(f"E slope vs n = {sE:+.4f}  (O(1) if ~0, Theta(n) if >0)")
        print(f"Var(log2 q_t) slope vs n = {sV:+.5f}  (O(1) if ~0)")
        print(f"per-coord influence ratios prof[d+1]/prof[d] (largest n): "
              f"{', '.join(f'{x:.2f}' for x in ratios)}")
        E_bounded = abs(sE) < 0.02
        decays = all(x < 0.85 for x in ratios) if ratios else False
        print(f">>> influence energy E BOUNDED: {'YES' if E_bounded else 'NO -> Theta(n)'}")
        print(f">>> per-coordinate influence GEOMETRICALLY DECAYING: {'YES' if decays else 'NO (flat)'}")
        print()
        if (not E_bounded or not decays) and abs(sV) < 0.02:
            print("VERDICT: Var(log2 q_t)=O(1) (small) BUT influence energy Theta(n) / non-decaying.")
            print("  The Efron-Stein/McDiarmid inequality is SLACK by Theta(n): concentration CANNOT")
            print("  reproduce the observed smallness. The route (annealed AD + concentration) FAILS")
            print("  on ingredient (b): log2 q_t is NOT a bounded-influence Holder functional; its")
            print("  smallness is the radius-asymptotic lattice cancellation the LLT itself supplies,")
            print("  not a concentration property. Same wall as the prior ball-CDF Efron-Stein route.")
        elif E_bounded and decays:
            print("VERDICT: influence energy BOUNDED and decaying => concentration route VIABLE.")
        else:
            print("VERDICT: AMBIGUOUS at these n; extend.")


if __name__ == "__main__":
    main()

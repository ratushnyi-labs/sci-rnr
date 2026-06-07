#!/usr/bin/env python3
"""
probe_7_34_R1_chaining_fast.py
============================================================================
FAST decisive test for TRACK (R1): does shell-uniform boundedness of
r_t(x)=q_t sqrt(2 pi v_q) REDUCE to annealed + (R2 single-flip stability) +
a McDiarmid/chaining argument, or does it need a NEW quenched lattice LLT?

The McDiarmid/Azuma route can upgrade annealed boundedness to SHELL-UNIFORM
(quenched) boundedness with Var(log r_t)=O(1) ONLY IF the WORST-CASE single-flip
bounded-difference coefficients c_i satisfy sum_i c_i^2 = O(1) (n-independent).
Equivalently, a continuity/chaining/covering argument over the shell needs the
map x|->log r_t to be Holder with a GEOMETRICALLY-LOCALIZED modulus (single-flip
influence decaying with distance from a localized active region).

This probe (smaller R, near-mean schedules only, larger n) reports:
  sumworst2 = sum_i (worst-case_x |D_i log2 r_t|)^2   -- the McDiarmid budget
  partRatio = participation ratio of |D_i| (1=spread, 0=localized)
  worst_c profile by right-distance d -- geometric decay <=> filter-localized
  Lipschitz slope mean|delta|/h on near-mean pairs at Hamming distance h.

VERDICT: sumworst2 O(1) & profile decaying & slope bounded => PROVES-route.
         sumworst2 grows / profile flat / influence spread => OPEN (new LLT).
"""
import math, sys
import numpy as np
import mpmath as mp
from math import comb


def D_c(p): return 0.5 * (1.0 - math.sqrt(1.0 - 2.0 * p) / (1.0 - p))
def V_lossless(p): return p * (1 - p) * (math.log2((1 - p) / p)) ** 2


def _mul_deg1(c, s):
    m = len(c); o = [None] * (m + 1); o[0] = c[0]
    for k in range(1, m): o[k] = c[k] + s * c[k - 1]
    o[m] = s * c[m - 1]; return o


def H_coeffs(x, p):
    n = len(x); half = mp.mpf('0.5'); x0 = int(x[0])
    v = [_mul_deg1([half], 1 if x0 == 0 else -1), _mul_deg1([half], 1 if x0 == 1 else -1)]
    P = mp.mpf(p); Qq = 1 - P
    for i in range(1, n):
        xi = int(x[i]); vb, va = v[0], v[1]; L = len(vb)
        a0 = [Qq * vb[k] + P * va[k] for k in range(L)]
        a1 = [P * vb[k] + Qq * va[k] for k in range(L)]
        v = [_mul_deg1(a0, 1 if xi == 0 else -1), _mul_deg1(a1, 1 if xi == 1 else -1)]
    Q = [v[0][k] + v[1][k] for k in range(len(v[0]))]
    while len(Q) < n + 1: Q.append(mp.mpf(0))
    return Q[:n + 1]


_KC = {}
def kraw_full(n):
    if n in _KC: return _KC[n]
    K = [[0] * (n + 1) for _ in range(n + 1)]
    for j in range(n + 1):
        for k in range(n + 1):
            s = 0; lo = max(0, k - (n - j)); hi = min(k, j)
            for l in range(lo, hi + 1):
                s += ((-1) ** l) * comb(j, l) * comb(n - j, k - l)
            K[k][j] = s
    _KC[n] = K; return K


def pi_k_vec(x, p, D, Km):
    n = len(x); H = H_coeffs(x, p); inv = mp.mpf(1) / (1 - 2 * mp.mpf(D))
    ip = [mp.mpf(1)] * (n + 1)
    for j in range(1, n + 1): ip[j] = ip[j - 1] * inv
    out = [mp.mpf(0)] * (n + 1)
    for k in range(n + 1):
        s = mp.mpf(0)
        for j in range(n + 1): s += ip[j] * Km[k][j] * H[j]
        out[k] = s / (mp.mpf(2) ** n)
    return out


def rt_and_mean(x, p, D, Km, dps=60):
    """Return (log2 r_t, mq) ; None if degenerate."""
    n = len(x); t = int(math.floor(n * D))
    th0 = float(mp.log(mp.mpf(D) / (1 - mp.mpf(D))))
    with mp.workdps(dps):
        pis = pi_k_vec(x, p, D, Km)
        M = sum(pis[k] * mp.e ** (mp.mpf(th0) * k) for k in range(n + 1))
        if M <= 0: return None
        q = [float(pis[k] * mp.e ** (mp.mpf(th0) * k) / M) for k in range(n + 1)]
    if q[t] <= 0: return None
    mq = sum(k * q[k] for k in range(n + 1))
    vq = sum((k - mq) ** 2 * q[k] for k in range(n + 1))
    if vq <= 0: return None
    rt = q[t] * math.sqrt(2 * math.pi * vq)
    if rt <= 0: return None
    return math.log2(rt), mq


def sample(n, p, rng):
    st = (rng.random(n) < p).astype(np.int8); st[0] = (rng.random() < 0.5)
    return np.cumsum(st) % 2


def analyze(p, n, R, seed, Dfrac=0.9, dps=60):
    D = Dfrac * D_c(p); Km = kraw_full(n); rng = np.random.default_rng(seed)
    nD = n * D; t = int(math.floor(nD))
    worst_c = np.zeros(n); part = []; ens = []
    pair_slopes = {h: [] for h in (1, 2, 3, 4)}
    kept = 0
    cache = []  # store (x, lr0) of near-mean words for pair test
    while kept < R:
        x = sample(n, p, rng)
        b = rt_and_mean(x, p, D, Km, dps)
        if b is None: continue
        lr0, mq = b
        if abs(mq - nD) > 1.0 + 0.1 * nD: continue
        Di = np.zeros(n); ok = True
        for i in range(n):
            xf = x.copy(); xf[i] ^= 1
            bf = rt_and_mean(xf, p, D, Km, dps)
            if bf is None: ok = False; break
            Di[i] = lr0 - bf[0]
        if not ok: continue
        kept += 1
        a = np.abs(Di); ens.append(float(np.sum(Di ** 2)))
        for i in range(n):
            d = n - 1 - i
            worst_c[d] = max(worst_c[d], a[i])
        s1 = a.sum(); s2 = (a ** 2).sum()
        if s2 > 0: part.append((s1 ** 2) / (n * s2))
        cache.append((x, lr0))
    # pair-continuity from cache
    for _ in range(4 * R):
        if len(cache) < 2: break
        i, j = rng.integers(0, len(cache), size=2)
        if i == j: continue
        x1, l1 = cache[i]; x2, l2 = cache[j]
        h = int(np.sum(x1 != x2))
        if h in pair_slopes: pair_slopes[h].append(abs(l1 - l2) / h)
    slopes = {h: (np.mean(v) if v else float('nan')) for h, v in pair_slopes.items()}
    return dict(n=n, t=t, kept=kept, sumworst2=float(np.sum(worst_c ** 2)),
                worst_c=worst_c, part=float(np.mean(part)) if part else float('nan'),
                E=float(np.mean(ens)), slopes=slopes)


def main():
    p = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4
    D = 0.9 * D_c(p)
    print(f"# BSMS p={p} D=0.9Dc={D:.5f} V_lossless={V_lossless(p):.5f} -- (R1) chaining budget")
    print(f"{'n':>4} {'t':>3} {'kept':>4} | {'sumWorst_c^2':>12} {'E[sumDi^2]':>11} {'partRatio':>9} | "
          f"worst_c[d] d=0..5 | Lip slope h=1,2,3,4")
    rows = []
    for n in (60, 80, 100, 120):
        t = int(math.floor(n * D))
        if t < 6: continue
        R = 80 if n <= 80 else (50 if n == 100 else 35)
        r = analyze(p, n, R, seed=6100 + n)
        rows.append(r)
        wc = " ".join(f"{r['worst_c'][d]:.3f}" for d in range(0, 6))
        sl = " ".join(f"{r['slopes'][h]:.3f}" for h in (1, 2, 3, 4))
        print(f"{n:>4} {r['t']:>3} {r['kept']:>4} | {r['sumworst2']:>12.4f} {r['E']:>11.4f} "
              f"{r['part']:>9.3f} | {wc} | {sl}", flush=True)
    if len(rows) >= 3:
        ns = np.array([r['n'] for r in rows], float)
        sw = np.array([r['sumworst2'] for r in rows])
        E = np.array([r['E'] for r in rows])
        ssw = np.polyfit(ns, sw, 1)[0]; sE = np.polyfit(ns, E, 1)[0]
        wc = rows[-1]['worst_c']
        rat = [wc[d + 1] / wc[d] for d in range(0, 5) if wc[d] > 1e-9]
        print("-" * 90)
        print(f"sum worst_c^2 slope/n = {ssw:+.4f}  -> "
              f"{'O(1) summable (McDiarmid budget bounded)' if abs(ssw) < 0.02 else 'GROWS Theta(n): worst-case influence NOT summable'}")
        print(f"E[sum Di^2] slope/n   = {sE:+.4f}  -> "
              f"{'BOUNDED (averaged Dirichlet energy O(1))' if abs(sE) < 0.01 else 'grows'}")
        print(f"worst_c profile decay ratios = {', '.join(f'{x:.2f}' for x in rat)} -> "
              f"{'LOCALIZED' if all(x<0.85 for x in rat) else 'FLAT (no geometric localization)'}")
        print()
        print("READING: the DISCRIMINATOR is the WORST-CASE chaining budget sum_i c_i^2.")
        print(" - If it grows Theta(n) while the AVERAGED energy E[sum Di^2] stays O(1), then the")
        print("   averaged (Efron-Stein/Glauber-Poincare) variance is O(1) [annealed], but the")
        print("   WORST-CASE bounded-difference budget is Theta(n): McDiarmid/Azuma sup-control")
        print("   over the shell gives only Var=O(n), NOT the O(1) the QUENCHED (R1) needs.")
        print(" - The averaged O(1) is exactly the ANNEALED Aaronson-Denker/Gouezel statement;")
        print("   the worst-case Theta(n) is exactly why the quenched upgrade does NOT follow from")
        print("   continuity/McDiarmid alone and needs a genuinely new quenched lattice LLT.")


if __name__ == "__main__":
    main()

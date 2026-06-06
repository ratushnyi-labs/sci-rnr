#!/usr/bin/env python3
"""
probe_7_34_localmass_influence_regularity.py
============================================================================
THE QUENCHED-UNIFORM ROUTE for the BSMS Gray-region RD-dispersion achievability
(sole open gap of Remark 7.34b).

The achievability needs the LOCAL lattice mass to be near 1/sqrt(2 pi v_q)
UNIFORMLY over the typical sqrt(n)-shell of schedules x:
    q_t(x) := P_{Y*}( d_H(., x) = t  |  tilt s* ) ~ 1/sqrt(2 pi v_q),
equivalently the LLT ratio r_t(x) := q_t(x) sqrt(2 pi v_q) -> const, PER schedule.
The annealed (averaged-over-x) local CLT (Aaronson-Denker for Gibbs-Markov maps)
gives only E_x[q_t(X)]. The proposed route upgrades annealed -> quenched-uniform
by CONCENTRATION of the random variable r_t(X) (equiv. log q_t(X)) as a regular
(Holder / bounded-influence) functional of the schedule x.

THIS PROBE tests the LOAD-BEARING REGULARITY claim of that route -- the thing a
McDiarmid / Kontorovich-Ramanan concentration inequality needs:

  (A) COORDINATE INFLUENCE of log q_t: D_i := log2 q_t(x) - log2 q_t(x^{flip i}).
      Does |D_i| DECAY GEOMETRICALLY as coordinate i moves away from the chain's
      "active" region (filter stability => geometric coordinate forgetting)?
  (B) INFLUENCE ENERGY  E[ sum_i (D_i)^2 ].  This is THE discriminator:
        - If O(1) (bounded, n-independent): bounded-difference / Efron-Stein gives
          Var(log q_t) = O(1) and McDiarmid gives concentration of r_t to o(1)
          => the annealed local CLT + concentration CLOSES the quenched-uniform LLT.
        - If Theta(n): the local mass has the SAME global-shell-object pathology as
          the ball CDF G_n (prior Efron-Stein route, BLOCKED) => route fails the
          same way; concentration only gives Var=O(n), not o(1).
  (C) DIRECT  Var_x( log2 q_t(X) )  over the shell vs n: bounded or growing?
      And the LLT ratio  r_t(x)  spread (cross-check the paper's 0.095 figure).

CONTRAST WITH THE PRIOR (BLOCKED) ROUTE (Session 2026-06-06): that route tried to
concentrate the BALL CDF G_n = -log2 P(d_H <= t) (or log2 S_n = G_n - j_n), which
has Theta(n) influence energy (it is a GLOBAL shell-slab CDF, the leading additive
influence of j_n only cancels at leading order; O(1) remainder per coord does not
decay). The OPEN QUESTION here is whether the SINGLE local mass q_t (one lattice
shell, the LLT object, NOT a CDF/cumulative slab) is MORE regular -- this is the
whole bet of the prompt's route (1).

Exact O(n^2) local mass via the Walsh/Krawtchouk transfer-matrix machinery
(pi_k = untilted shell mass; q_k propto pi_k e^{theta0 k}; the tilted local mass).
Run:    python scripts/verify/probe_7_34_localmass_influence_regularity.py [p]
Deps:   numpy, mpmath
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
    """Weight enumerator coeffs H_j(x) = [z^j] sum_{x'} P_X(x') prod((1+z)|x'_i=x_i,(1-z)|else)."""
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
    """Untilted output-side shell masses pi_k = P_{Y*}(d_H(.,x)=k), k=0..n."""
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


def local_mass_logq_t(x, p, D, Km, dps=70):
    """Return (log2 q_t, t, vq, mq, rt) for the TILTED local mass at the saddle.
    We tilt q_k propto pi_k e^{theta0 k} (theta0=ln(D/(1-D))) so the lower-tail
    boundary t=floor(nD) sits near the posterior mean (the near-mean LLT zone),
    matching the saddlepoint reduction. Returns log2 of the boundary local mass q_t."""
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
    logq = math.log2(q[t])
    rt = q[t] * math.sqrt(2 * math.pi * vq)
    return dict(logq=logq, t=t, vq=vq, mq=mq, rt=rt)


def sample(n, p, rng):
    st = (rng.random(n) < p).astype(np.int8)
    st[0] = (rng.random() < 0.5)
    return np.cumsum(st) % 2


# ---------------------------------------------------------------- (A)+(B) influence
def influence_profile(p, n, R, seed, Dfrac=0.9, dps=70):
    """Measure D_i = log2 q_t(x) - log2 q_t(x^{flip i}) for every coordinate i,
    averaged |D_i| profile (by distance from right end) + influence energy sum_i D_i^2."""
    D = Dfrac * D_c(p)
    Km = kraw_full(n)
    rng = np.random.default_rng(seed)
    energies = []
    # profile of mean |D_i| by RIGHT-distance d = n-1-i (filter forgets the past;
    # the boundary integer t is read at the end of the schedule scan)
    prof_sum = np.zeros(n); prof_cnt = np.zeros(n)
    kept = 0
    for _ in range(R):
        x = sample(n, p, rng)
        base = local_mass_logq_t(x, p, D, Km, dps)
        if base is None:
            continue
        lq0 = base['logq']
        Di = np.zeros(n)
        ok = True
        for i in range(n):
            xf = x.copy(); xf[i] ^= 1
            bf = local_mass_logq_t(xf, p, D, Km, dps)
            if bf is None:
                ok = False; break
            Di[i] = lq0 - bf['logq']
        if not ok:
            continue
        kept += 1
        energies.append(float(np.sum(Di ** 2)))
        for i in range(n):
            d = n - 1 - i
            prof_sum[d] += abs(Di[i]); prof_cnt[d] += 1
    if kept < 5:
        return None
    energies = np.array(energies)
    prof = np.divide(prof_sum, np.maximum(prof_cnt, 1))
    return dict(n=n, kept=kept, t=int(math.floor(n * D)),
                E=energies.mean(), E_se=energies.std() / math.sqrt(kept),
                prof=prof)


# ---------------------------------------------------------------- (C) direct variance
def direct_variance(p, n, R, seed, Dfrac=0.9, dps=70):
    """Var_x(log2 q_t(X)) and r_t spread over the source-typical shell, vs n."""
    D = Dfrac * D_c(p)
    Km = kraw_full(n)
    rng = np.random.default_rng(seed)
    lqs = []; rts = []; vqs = []
    for _ in range(R):
        x = sample(n, p, rng)
        b = local_mass_logq_t(x, p, D, Km, dps)
        if b is None:
            continue
        lqs.append(b['logq']); rts.append(b['rt']); vqs.append(b['vq'])
    if len(lqs) < 20:
        return None
    lqs = np.array(lqs); rts = np.array(rts); vqs = np.array(vqs)
    return dict(n=n, kept=len(lqs), t=int(math.floor(n * D)),
                Vlq=lqs.var(ddof=1), meanlq=lqs.mean(),
                rt_mean=rts.mean(), rt_sd=rts.std(), rt_relsd=rts.std() / rts.mean(),
                vq_mean=vqs.mean())


def main():
    p = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4
    D = 0.9 * D_c(p)
    print("#" * 100)
    print(f"# BSMS p={p}  D=0.9 D_c={D:.5f}  V_lossless={V_lossless(p):.5f}")
    print(f"# QUENCHED-UNIFORM ROUTE: is log2 q_t(x) (the LOCAL mass) a bounded-influence")
    print(f"# functional of the schedule x?  (A) influence decay, (B) energy O(1) vs Theta(n),")
    print(f"# (C) direct Var_x(log2 q_t) over the shell.")
    print("#" * 100)

    # -------- (B) THE DISCRIMINATOR: influence energy E[sum_i (D_i log q_t)^2] vs n
    print("\n(A)+(B) INFLUENCE of log2 q_t under single-bit flips of the schedule x.")
    print(f"{'n':>4} {'t':>3} {'kept':>5} | {'E=sum Di^2':>12} {'E/n':>8} | tail of mean|D_i| by right-distance d")
    Brows = []
    for n in (12, 14, 16, 18, 20, 22):
        R = 400 if n <= 16 else (200 if n <= 20 else 120)
        r = influence_profile(p, n, R, seed=7000 + n)
        if r is None:
            print(f"{n:>4} (insufficient)"); continue
        Brows.append(r)
        prof = r['prof']
        tailstr = " ".join(f"d{d}:{prof[d]:.3f}" for d in range(0, min(8, n)))
        print(f"{n:>4} {r['t']:>3} {r['kept']:>5} | {r['E']:>8.4f}+-{r['E_se']:.3f} "
              f"{r['E']/n:>8.4f} | {tailstr}", flush=True)

    if len(Brows) >= 4:
        ns = np.array([r['n'] for r in Brows], float)
        E = np.array([r['E'] for r in Brows])
        # fit E ~ a*n + b  and  E ~ const
        A = np.vstack([ns, np.ones_like(ns)]).T
        a, b = np.linalg.lstsq(A, E, rcond=None)[0]
        # geometric decay check: ratio prof[d+1]/prof[d] over small d for largest n
        prof = Brows[-1]['prof']
        ratios = [prof[d + 1] / prof[d] for d in range(0, 5) if prof[d] > 1e-9]
        print("-" * 100)
        print(f"(B) influence energy fit  E ~ {a:+.4f}*n + {b:+.3f}")
        print(f"    geometric-decay ratios prof[d+1]/prof[d] (largest n): "
              f"{', '.join(f'{x:.2f}' for x in ratios)}")
        bounded = (abs(a) < 0.02) and (b < 5.0)
        decays = all(x < 0.85 for x in ratios) if ratios else False
        print(f"    >>> influence energy BOUNDED (slope~0): {'YES' if bounded else 'NO (Theta(n))'}")
        print(f"    >>> per-coordinate influence GEOMETRICALLY DECAYING: {'YES' if decays else 'NO'}")

    # -------- (C) direct Var_x(log2 q_t) over the shell, larger n
    print("\n(C) DIRECT  Var_x(log2 q_t(X))  and r_t spread over the source shell, vs n.")
    print(f"{'n':>4} {'t':>3} {'kept':>5} | {'Var(log2 q_t)':>14} | {'r_t mean':>9} {'r_t sd':>7} "
          f"{'r_t relsd':>9} | {'vq':>7}")
    Crows = []
    for n in (64, 96, 128, 160, 220, 300, 400):
        t = int(math.floor(n * D))
        if t < 8:
            continue
        R = 400 if n <= 128 else (250 if n <= 220 else 130)
        r = direct_variance(p, n, R, seed=8000 + n)
        if r is None:
            print(f"{n:>4} (insufficient)"); continue
        Crows.append(r)
        print(f"{n:>4} {r['t']:>3} {r['kept']:>5} | {r['Vlq']:>14.5f} | "
              f"{r['rt_mean']:>9.4f} {r['rt_sd']:>7.4f} {r['rt_relsd']:>9.4f} | {r['vq_mean']:>7.2f}",
              flush=True)
    if len(Crows) >= 3:
        ns = np.array([r['n'] for r in Crows], float)
        Vlq = np.array([r['Vlq'] for r in Crows])
        relsd = np.array([r['rt_relsd'] for r in Crows])
        sv = np.polyfit(ns, Vlq, 1)[0]
        srel = np.polyfit(ns, relsd, 1)[0]
        print("-" * 100)
        print(f"(C) Var(log2 q_t) slope vs n = {sv:+.5f}  (BOUNDED if ~0)")
        print(f"    r_t relative spread slope vs n = {srel:+.2e}  (CONCENTRATING if <=0)")
        Cbounded = abs(sv) < 0.01
        print(f"    >>> Var_x(log2 q_t) BOUNDED over the shell: {'YES' if Cbounded else 'NO (grows)'}")

    print()
    print("VERDICT LOGIC:")
    print("  Route CLOSES-modulo-concentration IFF (B) energy is O(1)/bounded AND geometrically")
    print("  decaying (=> McDiarmid/Kontorovich-Ramanan gives Var(log q_t)=O(1), concentration")
    print("  of r_t to o(1), and annealed local CLT + concentration => quenched-uniform LLT).")
    print("  Route FAILS the SAME way as the ball-CDF Efron-Stein route IFF (B) energy is Theta(n)")
    print("  (the local mass is then ALSO a global non-additive shell object, not a regular")
    print("  Holder/bounded-influence functional).")


if __name__ == "__main__":
    main()

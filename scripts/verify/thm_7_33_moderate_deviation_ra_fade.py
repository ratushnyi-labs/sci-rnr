#!/usr/bin/env python3
"""
Verification for Theorem 7.33 (Moderate-deviation rate of the RNR archive, and the
fading of the random-access penalty across the CLT->LDP bridge).

Setting (Theorems 7.29/7.31): per-position cost l_i = -log2 P(X_i | past), realised
length L^rep = sum l_i, source entropy rate h, varentropy V.  Moderate-deviation
scaling: budget = N h + a_N * sqrt(N V) with a_N -> oo, a_N = o(sqrt N).

Claims:
  V1  MDP: -(1/a_N^2) log2 P(L^rep >= Nh + a_N sqrt(NV)) -> 1/(2 ln2)
      (base-2; equivalently P = exp(-(1+o(1)) a_N^2 / 2) in nats), uniformly over the
      moderate range.  [Standard MDP for bounded weakly-dependent sums; Petrov 1975.]
  V2  RA penalty = constant standardized SHIFT.  With sub-block overhead (N/K)delta
      at K=c sqrt(N) (= sqrt(N) delta/c), the RNR overflow at the same budget has
      effective deviation (a_N - beta), beta := delta/(c sqrt V):
        P(L^rep + (N/K)delta >= Nh + a_N sqrt(NV)) = exp(-(1+o(1))(a_N-beta)^2/2).
  V3  The FADE (three regimes unified):
       (a) exponent ratio (a_N-beta)^2 / a_N^2 -> 1 as a_N -> oo  (RA invisible at the
           exponent level -> recovers Thm 7.31 LDP-invisibility);
       (b) at a_N = O(1) the shift beta is comparable to a_N (RA visible -> Thm 7.29);
       (c) the overflow PROBABILITY still carries the growing factor
           exp(beta a_N - beta^2/2) -> oo: the RA penalty fades from the EXPONENT but
           multiplies the PROBABILITY.
  V4  Cross-check vs the exact Cramer/LDP exponent E(R) of Thm 7.31: near the mean,
      E(h + a_N sqrt(V/N)) * N  ->  (a_N-... )  matches the MDP quadratic (the LDP
      rate's quadratic expansion equals the MDP rate).
"""

import numpy as np
from math import lgamma, log, log2

np.random.seed(0)
LN2 = np.log(2.0)
ok_all = True


def report(tag, ok, msg):
    global ok_all
    ok_all = ok_all and ok
    print(f"[{'PASS' if ok else 'FAIL'}] {tag}: {msg}")


# i.i.d. source; ideal predictor l = -log2 P0(x).  Exact tail via multinomial.
P0 = np.array([0.6, 0.3, 0.1]); A = len(P0)
l = -np.log2(P0)
h = float((P0 * l).sum())
V = float((P0 * (l - h) ** 2).sum())            # varentropy


def exact_tail_log2p(N, thresh):
    """log2 P(sum l_i >= thresh) exactly via the multinomial law (small alphabet)."""
    terms = []
    logfacN = lgamma(N + 1); logP = np.log(P0)

    def rec(idx, rem, counts):
        if idx == A - 1:
            c = counts + [rem]
            s = sum(ci * li for ci, li in zip(c, l))
            if s >= thresh - 1e-9:
                lp = logfacN + sum(-lgamma(ci + 1) + ci * lpx for ci, lpx in zip(c, logP))
                terms.append(lp)
            return
        for n in range(rem + 1):
            rec(idx + 1, rem - n, counts + [n])

    rec(0, N, [])
    if not terms:
        return -np.inf
    m = max(terms)
    return (m + log(sum(np.exp(np.array(terms) - m)))) / LN2


# ----------------------------------------------------------------------------
# V1 : MDP  -(1/a_N^2) log2 P -> 1/(2 ln2)
# ----------------------------------------------------------------------------
print("=" * 70)
print("V1  moderate-deviation rate -(1/a_N^2) log2 P(L>=Nh+a_N sqrt(NV)) -> 1/(2 ln2)")
target = 1.0 / (2 * LN2)                          # base-2 MDP constant
# choose (N, a_N) along the moderate path a_N ~ N^{1/4} (-> oo, = o(sqrt N))
rows = []
for N in [400, 1600, 6400]:
    aN = N ** 0.25                               # moderate scaling
    thr = N * h + aN * np.sqrt(N * V)
    log2p = exact_tail_log2p(N, thr)
    rate = -log2p / aN ** 2
    rows.append((N, aN, rate))
    print(f"     N={N:5d}  a_N={aN:5.2f}  -(1/a_N^2)log2 P = {rate:.4f}   target 1/(2ln2)={target:.4f}")
# converging toward target from above (finite-N + prefactor corrections)
conv = abs(rows[-1][2] - target) < abs(rows[0][2] - target) and abs(rows[-1][2] - target) < 0.15
report("V1", conv and rows[-1][2] > 0,
       f"MDP rate -> 1/(2ln2)={target:.4f} along a_N~N^1/4 (rate {rows[0][2]:.3f}->{rows[-1][2]:.3f})")


# ----------------------------------------------------------------------------
# V2 : RA penalty = constant standardized shift beta
# ----------------------------------------------------------------------------
print("=" * 70)
print("V2  RA penalty (N/K)delta at K=c sqrt N  ==  constant standardized shift beta")
delta = 0.25; c = 1.0
beta = delta / (c * np.sqrt(V))
ok2 = True
for N in [400, 1600, 6400]:
    aN = N ** 0.25
    ra = (N / (c * np.sqrt(N))) * delta          # = sqrt(N) delta/c
    # RNR overflow at budget Nh + a_N sqrt(NV): P(sum l >= budget - ra)
    thr_rnr = N * h + aN * np.sqrt(N * V) - ra
    # this equals threshold with effective deviation (a_N - beta):
    thr_eff = N * h + (aN - beta) * np.sqrt(N * V)
    ok2 = ok2 and abs(thr_rnr - thr_eff) < 1e-6
print(f"     beta = delta/(c sqrt V) = {beta:.4f}  (delta={delta}, V={V:.4f})")
report("V2", ok2,
       f"P(L_RNR>=Nh+a_N sqrt(NV)) = P(L>=Nh+(a_N-beta)sqrt(NV)); RA = constant shift "
       f"beta={beta:.4f} in standardized units")


# ----------------------------------------------------------------------------
# V3 : the FADE -- exponent ratio ->1, but probability factor -> oo
# ----------------------------------------------------------------------------
print("=" * 70)
print("V3  RA fades from the EXPONENT (ratio->1) but multiplies the PROBABILITY (->oo)")
exp_ratios, logratio_over_aN, prob_factors = [], [], []
for aN in [1.0, 3.0, 10.0, 30.0, 100.0]:
    exp_ratio = (aN - beta) ** 2 / aN ** 2       # ->1 as aN->oo (RA invisible in exponent)
    # rigorous: (1/a_N) log(ratio) = (a_N^2/2 - (a_N-beta)^2/2)/a_N = beta - beta^2/(2 a_N) -> beta
    lr = (aN ** 2 / 2 - (aN - beta) ** 2 / 2) / aN
    prob_factor = np.exp(beta * aN - beta ** 2 / 2)   # CLT-scale value; order exp(beta a_N) -> oo
    exp_ratios.append(exp_ratio); logratio_over_aN.append(lr); prob_factors.append(prob_factor)
    print(f"     a_N={aN:6.1f}  (a_N-beta)^2/a_N^2={exp_ratio:.4f}  (1/a_N)log(ratio)={lr:.4f}(->beta={beta:.4f})"
          f"  prob factor~{prob_factor:11.2f}")
fade_ok = (abs(exp_ratios[-1] - 1) < 1e-2          # exponent ratio -> 1 (RA invisible in exponent)
           and exp_ratios[0] < 0.9                  # at a_N=O(1): RA visible (ratio << 1)
           and abs(logratio_over_aN[-1] - beta) < 1e-2   # (1/a_N)log(ratio) -> beta (rigorous)
           and prob_factors[-1] > prob_factors[0] * 1e3)  # probability factor diverges
report("V3", fade_ok,
       f"exponent ratio {exp_ratios[0]:.3f}(a_N=1,RA visible,~Thm7.29)->{exp_ratios[-1]:.3f}"
       f"(a_N=100,RA invisible,~Thm7.31); rigorously (1/a_N)log(ratio)->beta={beta:.3f}; "
       f"prob factor {prob_factors[0]:.1f}->{prob_factors[-1]:.2e} diverges")
# V3b: EXACT-tail confirmation of the probability-ratio claim (not just algebra) --
# (1/a_N) log[P(L>=B-ra)/P(L>=B)] rises toward beta along the moderate path a_N~N^1/4,
# with the slow O(N^-1/4) Cramer-Petrov convergence the proof predicts.
print("-" * 70)
print("V3b  EXACT-tail ratio (1/a_N)log P(L>=B-ra)/P(L>=B) trends UP toward beta")
exact_lr = []
for N in [400, 900, 1600]:
    aN = N ** 0.25
    B = N * h + aN * np.sqrt(N * V); ra = np.sqrt(N) * delta / c
    lpP = exact_tail_log2p(N, B); lpR = exact_tail_log2p(N, B - ra)
    lr = (lpR - lpP) * LN2 / aN                  # (1/a_N) log_e(ratio)
    exact_lr.append(lr)
    print(f"     N={N:5d} a_N={aN:5.2f}  exact (1/a_N)log(ratio)={lr:.4f}  (-> beta={beta:.4f})")
trend_ok = (exact_lr[0] < exact_lr[1] < exact_lr[2] < beta
            and (beta - exact_lr[2]) < (beta - exact_lr[0]))   # monotone up toward beta, gap shrinks
report("V3b", trend_ok,
       f"exact (1/a_N)log(ratio) {exact_lr[0]:.3f}->{exact_lr[2]:.3f} rises toward beta={beta:.3f} "
       f"(O(N^-1/4) Cramer-Petrov convergence) -- the probability inflation is REAL, not algebra")


# ----------------------------------------------------------------------------
# V4 : the LDP rate's quadratic-at-mean equals the MDP rate (consistency)
# ----------------------------------------------------------------------------
print("=" * 70)
print("V4  LDP rate E(R) quadratic expansion at the mean == MDP rate")
def E_of_R(R, t_hi=200.0):                       # base-2 LDP rate, iid matched
    ts = np.linspace(1e-6, t_hi, 200000)
    Lam = np.log2((P0[None, :] * P0[None, :] ** (-ts[:, None])).sum(axis=1))
    return float(np.max(ts * R - Lam))
# near the mean: E(h+x) ~ x^2/(2 ln2 V); set x = a_N sqrt(V/N) (moderate), check
# N*E(h+x) ~ a_N^2/(2 ln2) = the MDP exponent
xs = [0.02, 0.01, 0.005]
ok4 = True
for x in xs:
    Ex = E_of_R(h + x)
    quad = x ** 2 / (2 * LN2 * V)
    ok4 = ok4 and abs(Ex - quad) / quad < 0.05
    print(f"     x={x}: E(h+x)={Ex:.3e}  x^2/(2 ln2 V)={quad:.3e}  (MDP quadratic)")
report("V4", ok4, "LDP exponent's quadratic-at-mean = MDP rate (the two regimes agree at the bridge)")


print("=" * 70)
print("ALL PASS" if ok_all else "SOME FAILED")
import sys
sys.exit(0 if ok_all else 1)

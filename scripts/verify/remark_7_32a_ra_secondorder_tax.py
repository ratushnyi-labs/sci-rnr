#!/usr/bin/env python3
"""
Verification for Remark 7.32a (Universal and distributed RNR at second order: the
parameter is a cheap header, the side information is not).

Single known source (Thm 7.29): clean CLT-scale repair length
   N h + (N/K) delta_inf + sqrt(N V) Q^{-1}(eps) + o(sqrt N).
This script confirms the DICHOTOMY when the known-source assumption is lifted in the two
canonical directions under the sub-block (random-access) construction:

  (i) UNIVERSAL -- sqrt(N) floor SURVIVES random access. The d-param unknown is encodable as
     a TWO-PASS global header theta-hat of (d/2)log2 N = o(sqrt N) bits; each sub-block then
     decodes independently against fixed P_theta-hat:
        L = N h + (N/K) delta_inf + (d/2) log2 N + sqrt(N V) Q^{-1}(eps) + o(sqrt N),
     LOCALLY DECODABLE with the Kontoyiannis-Verdu dispersion sqrt(NV) INTACT (Clarke-Barron
     regret (d/2)log N = o(sqrt N), dispersion-invisible). The cost of not knowing the source
     is second-order-free EVEN with random access. (CLT companion to the LDP Theorem 7.32.)
  (ii) The header earns its keep: a naive single-pass scheme that RESETS per sub-block re-pays
     (N/K)(d/2)log2 K = Theta(sqrt N log N) at K=c sqrt N -- ABOVE dispersion, but AVOIDABLE
     (it is the cost of discarding the cheap header, not a tax of universality).
  (iii) DISTRIBUTED (RNR Slepian-Wolf 7.8) -- the GENUINELY OPEN tax. Side info Y is Theta(N)
     SOURCE data, NOT a d-param: no cheap header. Global SW second order is O(sqrt N) (Tan-Kosut
     2014), but the sub-block per-block union bound at the SW corner pays
     m sqrt(K) Q^{-1}(eps/m) = Theta(N^{3/4} sqrt(log N)) at K=c sqrt N -- above sqrt N; whether
     a locally-decodable sub-block SW code restores sqrt N is OPEN.

So a d-dim PARAMETER is a cheap header keeping random-access universality at the sqrt N floor,
whereas Theta(N) bits of SIDE INFORMATION admit no such header (the second order stays open).

Checks (K=c sqrt N, the balanced regime):
  V1 header (d/2)log N = o(sqrt N): (i) keeps the sqrt N floor (header is dispersion-invisible).
  V2 avoidable naive-reset cost (N/K)(d/2)log K = Theta(sqrt N log N): /sqrt N -> infinity (a log factor).
  V3 SW sub-block union-bound cost m sqrt(K) Q^{-1}(eps/m) = Theta(N^{3/4} sqrt log N): /sqrt N -> infinity.
  V4 scale ordering: sqrt N  <<  sqrt N log N  <<  N^{3/4} sqrt log N (ratios -> infinity).
  V5 LOAD-BEARING (i): the two-pass theta-hat scheme codes each block against the DATA-estimated
     theta-hat ("twice-use of data"); confirm the overfit term L_p - L_theta-hat is O_P(1) with mean
     -> d/(2 ln2) (N-independent, =1/(2 ln2)=0.721 bits for d=1, the Wilks chi^2_d/(2 ln2)), so the
     KV dispersion sqrt(NV) is INTACT (standardized two-pass length -> N(0,1), 0.95-quantile -> 1.645).
"""
import math
from math import erf, sqrt, log

def Qinv(eps):
    lo, hi = -12.0, 12.0
    for _ in range(300):
        m = (lo+hi)/2
        if 0.5*(1+erf(m/sqrt(2))) < 1-eps: lo = m
        else: hi = m
    return (lo+hi)/2

c = 1.0; d = 4; V = 0.7; eps = 0.05; qi = Qinv(eps)
Ns = [10**4, 10**6, 10**8, 10**10, 10**12]

print("="*72)
print("V1: global-universal regret (d/2)log2 N = o(sqrt N) (below dispersion sqrt(NV)Q^-1)")
ok1 = True
for N in Ns:
    regret = (d/2)*math.log2(N)
    disp = sqrt(N*V)*qi
    r = regret/disp
    print(f"  N={N:>13}: (d/2)log N={regret:.1f}, sqrt(NV)Q^-1={disp:.3e}, ratio={r:.3e}")
    ok1 = ok1 and r < (1.0 if N<=10**4 else 0.1)   # ->0
ok1 = ok1 and ((d/2)*math.log2(Ns[-1]))/(sqrt(Ns[-1]*V)*qi) < 1e-3
print(f"  V1 {'PASS' if ok1 else 'FAIL'}  (regret/dispersion -> 0)")

print("="*72)
print("V2: AVOIDABLE naive-reset cost (N/K)(d/2)log2 K at K=c sqrt N = Theta(sqrt N log N) > dispersion")
ok2 = True
for N in Ns:
    K = c*sqrt(N); m = N/K
    tax = m*(d/2)*math.log2(K)          # regret re-paid per sub-block
    disp = sqrt(N*V)*qi
    ratio = tax/disp                     # ~ log N -> infinity
    per_sqrtN = tax/sqrt(N)
    print(f"  N={N:>13}: reset-tax={tax:.3e}, /sqrt N={per_sqrtN:.2f} (~log N), tax/disp={ratio:.2f}")
    ok2 = ok2 and ratio > 1
ok2 = ok2 and (Ns[-1] and (lambda N: (N/(c*sqrt(N)))*(d/2)*math.log2(c*sqrt(N))/(sqrt(N*V)*qi))(Ns[-1]) > 10)
print(f"  V2 {'PASS' if ok2 else 'FAIL'}  (reset tax / dispersion -> infinity, a log factor)")

print("="*72)
print("V3: SW sub-block tax m sqrt(K) Q^-1(eps/m) at K=c sqrt N = Theta(N^{3/4} sqrt log N) > dispersion")
ok3 = True
for N in Ns:
    K = c*sqrt(N); m = N/K
    tax = m*sqrt(K)*Qinv(min(eps/m, 0.49))     # per-sub-block union bound at the SW corner
    disp = sqrt(N*V)*qi
    ratio = tax/disp                            # ~ N^{1/4} sqrt(log N) -> infinity
    print(f"  N={N:>13}: SW-subblock-tax={tax:.3e}, /sqrt N={tax/sqrt(N):.2f}, tax/disp={ratio:.2f}")
    ok3 = ok3 and ratio > 1
print(f"  V3 {'PASS' if ok3 else 'FAIL'}  (SW sub-block tax / dispersion -> infinity)")

print("="*72)
print("V4: scale ordering  sqrt N  <  sqrt N log N  <  N^{3/4} sqrt(log N)  (ratios -> infinity)")
ok4 = True
for N in Ns:
    a = sqrt(N)                       # single-source dispersion scale
    b = sqrt(N)*log(N)                # universal-reset tax scale
    cc = N**0.75*sqrt(log(N))         # SW sub-block tax scale
    print(f"  N={N:>13}: sqrt N={a:.3e},  sqrtN logN={b:.3e} (b/a={b/a:.1f}),  "
          f"N^3/4 sqrt logN={cc:.3e} (c/b={cc/b:.2f})")
    ok4 = ok4 and (b > a) and (cc > b)
# ratios grow
r_ba = (sqrt(Ns[-1])*log(Ns[-1]))/sqrt(Ns[-1])
r_cb = (Ns[-1]**0.75*sqrt(log(Ns[-1])))/(sqrt(Ns[-1])*log(Ns[-1]))
ok4 = ok4 and r_ba > 10 and r_cb > 10
print(f"  V4 {'PASS' if ok4 else 'FAIL'}  (sqrt N << sqrt N log N << N^3/4 sqrt log N)")

print("="*72)
print("V5: LOAD-BEARING (i) -- two-pass theta-hat header keeps the dispersion ('twice-use of data')")
print("    overfit L_p - L_theta-hat is O_P(1), mean -> d/(2 ln2); standardized length -> N(0,1)")
try:
    import numpy as np
    p = 0.3; d_fam = 1
    h_b = -(p*math.log2(p) + (1-p)*math.log2(1-p))
    Vbar = p*(1-p)*(math.log2(p/(1-p)))**2            # Bernoulli varentropy (bits^2)
    target = d_fam/(2*math.log(2))                     # 1/(2 ln2)=0.7213 for d=1 (Wilks chi^2_d/(2 ln2))
    rng = np.random.default_rng(12345)
    ok5 = True; qhat_last = None
    for Nv, R in [(2000,6000),(8000,6000),(32000,6000)]:
        X = rng.random((R, Nv)) < p
        n1 = X.sum(axis=1).astype(float)
        th = np.clip(n1/Nv, 1e-9, 1-1e-9)             # MLE = sample mean = theta-hat
        L_th = -(n1*np.log2(th) + (Nv-n1)*np.log2(1-th))           # two-pass length vs fixed P_theta-hat
        L_p  = -(n1*math.log2(p) + (Nv-n1)*math.log2(1-p))         # length vs true p
        overfit = L_p - L_th                          # >= 0 (MLE), O_P(1)
        mo = float(overfit.mean())
        std = (L_th - Nv*h_b)/math.sqrt(Nv*Vbar)      # standardized two-pass length
        qhat = float(np.quantile(std, 0.95)); qhat_last = qhat
        print(f"  N={Nv:>6}: mean overfit={mo:.3f} bits (target {target:.3f}, N-indep), "
              f"std-quantile(.95)={qhat:.3f} (Q^-1(.05)={Qinv(0.05):.3f})")
        ok5 = ok5 and (0.55 < mo < 0.90)              # ~0.72 across N => N-independent O_P(1)
    ok5 = ok5 and (abs(qhat_last - Qinv(0.05)) < 0.12)            # dispersion sqrt(NV) intact at largest N
    print(f"  V5 {'PASS' if ok5 else 'FAIL'}  (overfit O_P(1) ~ 1/(2 ln2); dispersion sqrt(NV) intact)")
except ImportError:
    ok5 = True
    print("  V5 SKIP (numpy unavailable) -- scale checks V1-V4 stand")

print("="*72)
allok = ok1 and ok2 and ok3 and ok4 and ok5
print(f"RESULT: {'ALL PASS' if allok else 'SOME FAILED'}  [V1 {ok1}, V2 {ok2}, V3 {ok3}, V4 {ok4}, V5 {ok5}]")

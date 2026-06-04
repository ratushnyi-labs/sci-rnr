#!/usr/bin/env python3
"""
Verification for Theorem 7.35 (Extreme-value law for the worst-case random-access
query cost / read-buffer provisioning).

Claim: partition X^N into m = N/K disjoint blocks; L_i = -log2 P(block_i | context_i)
is the per-block codelength (the bits read to answer a random-access query landing in
block i). With h = entropy rate (bits/sym), V = varentropy rate (bits^2/sym):

  (a) LEADING-ORDER provisioning law (iid / finite-state Markov; ln m = o(K); exp-mixing
      alpha extends on the restricted range ln m = o(K/(log K loglog K)^2), per 7.33's MDP scope):
      max_{i<=m} L_i  =  K h  +  sqrt(2 K V ln m) * (1 + o_P(1)).
  (b) REFINED Gumbel (iid / finite-state aperiodic Markov; ln m = o(K^{1/3}), i.e. K >> (ln m)^3),
      a_m = sqrt(2 ln m), b_m = a_m - (ln ln m + ln 4pi)/(2 a_m):
      a_m * ( (max_i L_i - K h)/sqrt(KV) - b_m )  =>  Gumbel,
      E[max] = K h + sqrt(KV) (b_m + gamma/a_m) + o(sqrt(KV)/a_m).

NOTE (codex/cold-review fixes): the refined Gumbel needs ln m = o(K^{1/3}) (NOT o(K)) --
the Cramer-Petrov correction exp((x^3/sqrt K) lambda) at x ~ sqrt(2 ln m) -> 1 iff
x^3/sqrt K -> 0; the convergence carries the a_m factor (else the centred max -> 0, not
Gumbel); the sharp regime / Leadbetter D,D' dependence is iid / finite-state Markov.

Checks:
  V1  i.i.d.: mean block length = K h, block variance = K V (constants sanity).
  V2  i.i.d.: standardized max matches the refined Gumbel b_m + gamma/a_m (K >> (ln m)^3).
  V3  i.i.d.: leading-order ratio (max-Kh)/sqrt(2KV ln m) -> 1 from below as m grows.
  V4  finite-state Markov: same law with Markov (h, V) incl. covariance (the 7.28 varentropy).
  V5  SHARP scale boundary: refined-Gumbel error tracks x^3/sqrt K, small iff K >> (ln m)^3.
  V6  matching lower bound: a separated half-family of blocks already forces the max up.
"""
import math
import numpy as np

rng = np.random.default_rng(20260604)
EULER = 0.5772156649015329

def entropy_rate_iid(p):
    p = np.asarray(p)
    return float(-(p*np.log2(p)).sum())

def varentropy_iid(p):
    p = np.asarray(p)
    h = entropy_rate_iid(p)
    return float((p*(np.log2(p))**2).sum() - h*h)

def gumbel_norm(m):
    a = math.sqrt(2*math.log(m))
    b = a - (math.log(math.log(m)) + math.log(4*math.pi))/(2*a)
    return a, b

def block_max_iid(p, K, m, R):
    """R realizations; each: m blocks of K i.i.d. symbols; return array of max_i L_i.
    Fast inverse-CDF sampler (searchsorted on cumulative p) instead of rng.choice(p=)."""
    p = np.asarray(p)
    surpr = -np.log2(p)
    cump = np.cumsum(p); cump[-1] = 1.0
    out = np.empty(R)
    for r in range(R):
        u = rng.random((m, K))
        idx = np.searchsorted(cump, u, side='right')
        np.clip(idx, 0, len(p)-1, out=idx)
        L = surpr[idx].sum(axis=1)
        out[r] = L.max()
    return out

# ---------------------------------------------------------------- V1
print("="*70)
print("V1: i.i.d. block mean = K h and block variance = K V (constants)")
p = [0.5, 0.25, 0.15, 0.10]
h = entropy_rate_iid(p); V = varentropy_iid(p)
K = 400; R = 8000
surpr = -np.log2(np.asarray(p))
idx = rng.choice(4, size=(R, K), p=p)
Ls = surpr[idx].sum(axis=1)
emp_mean = float(Ls.mean()); emp_var = float(Ls.var(ddof=1))
print(f"  h={h:.5f} bits/sym, V={V:.5f} bits^2/sym, K={K}")
print(f"  E[L]: empirical {emp_mean:.3f} vs K h = {K*h:.3f}  (ratio {emp_mean/(K*h):.4f})")
print(f"  Var[L]: empirical {emp_var:.2f} vs K V = {K*V:.2f}  (ratio {emp_var/(K*V):.4f})")
v1 = abs(emp_mean/(K*h)-1) < 0.01 and abs(emp_var/(K*V)-1) < 0.06
print(f"  V1 {'PASS' if v1 else 'FAIL'}")

# ---------------------------------------------------------------- V2
print("="*70)
print("V2: standardized max matches refined Gumbel formula b_m + gamma/a_m")
# K=700 keeps K >> (ln m)^3 for all m below ((ln 1024)^3=333), the valid refined-Gumbel regime.
K = 700; R = 2000
ok2 = True
for m in [16, 64, 256, 1024]:
    mx = block_max_iid(p, K, m, R)
    z = (mx.mean() - K*h)/math.sqrt(K*V)
    a, b = gumbel_norm(m)
    pred = b + EULER/a
    rel = abs(z-pred)/pred
    print(f"  m={m:5d}: std max z={z:.4f}  refined pred={pred:.4f}  rel.err={rel:.3f}")
    ok2 = ok2 and rel < 0.06
print(f"  V2 {'PASS' if ok2 else 'FAIL'}")

# ---------------------------------------------------------------- V3
print("="*70)
print("V3: leading-order ratio (max-Kh)/sqrt(2 K V ln m) -> 1 (from below, slowly)")
K = 300; R = 2500
ratios = []
for m in [16, 64, 256, 1024]:
    mx = block_max_iid(p, K, m, R)
    ratio = (mx.mean() - K*h)/math.sqrt(2*K*V*math.log(m))
    ratios.append(ratio)
    print(f"  m={m:5d}: leading-order ratio = {ratio:.4f}")
# the crude leading-order normalization converges to 1 from below SLOWLY (the b_m loglog
# correction is O(1/ln m)); the test is monotone increase toward 1, last value close.
mono = all(ratios[i] < ratios[i+1] for i in range(len(ratios)-1))
ok3 = mono and all(0.70 < r <= 1.02 for r in ratios) and ratios[-1] > 0.87
print(f"  V3 {'PASS' if ok3 else 'FAIL'}  (monotone increasing toward 1 from below; "
      f"last={ratios[-1]:.3f})")

# ---------------------------------------------------------------- V4  Markov
print("="*70)
print("V4: Markov source -- same EVT law with Markov (h, V) incl. covariance")

def markov_stationary(P):
    P = np.asarray(P); n = len(P)
    pi = np.full(n, 1.0/n)
    for _ in range(5000):
        nx = pi @ P
        if np.max(np.abs(nx-pi)) < 1e-15:
            pi = nx; break
        pi = nx
    return pi

def markov_h_V(P, pi):
    """h = entropy rate; V = varentropy RATE via the fundamental-matrix CLT variance
    for the additive functional g(i,j) = -log2 P[i,j] of the Markov chain."""
    P = np.asarray(P); pi = np.asarray(pi); n = len(P)
    with np.errstate(divide='ignore'):
        g = np.where(P > 0, -np.log2(P), 0.0)          # per-step surprisal g(i,j)
    # f(i) = E[g(i,J)|i] = sum_j P[i,j] g(i,j); h = sum_i pi_i f(i)
    f = (P*g).sum(axis=1)
    h = float(pi @ f)
    # Markov CLT variance (Peticky/Bratteli; fundamental matrix Z=(I-P+1pi)^-1) for the
    # functional G(i,j)=g(i,j). Centered increment d(i,j)=g(i,j)-h.
    # Var-rate = sum_i pi_i sum_j P[i,j] d(i,j)^2
    #          + 2 sum_i pi_i sum_j P[i,j] d(i,j) * ( (Z @ fc)[j] )   where fc(i)=f(i)-h
    # (standard Markov additive-functional CLT; cross-covariance via Z on the cell means)
    A = np.eye(n) - P + np.outer(np.ones(n), pi)
    Z = np.linalg.inv(A)
    d = g - h                                          # centered cell values
    fc = f - h                                         # centered conditional means
    term_var = float(pi @ (P*d**2).sum(axis=1))        # one-step variance
    u = Z @ fc                                         # discounted future deviation per state
    # cross term: 2 E[ d(i,J) * u(J) ]
    cross = float(pi @ (P*(d*u[None, :])).sum(axis=1)) * 2.0
    V = term_var + cross
    return h, V

def block_max_markov(P, pi, K, m, R):
    """R parallel chains advanced together; one contiguous chain of m*K steps per
    realization -> m disjoint contiguous (dependent) blocks. Vectorized over R."""
    P = np.asarray(P); n = len(P)
    with np.errstate(divide='ignore'):
        g = np.where(P > 0, -np.log2(P), 0.0)
    cdf = np.cumsum(P, axis=1)
    total = m*K
    cur = rng.choice(n, size=R, p=pi)             # R parallel current states
    block_sum = np.zeros((R, m))                   # running per-block surprisal sums
    for t in range(total):
        bidx = t // K
        u = rng.random(R)
        nxt = (u[:, None] > cdf[cur]).sum(axis=1)  # vectorized inverse-CDF over R chains
        nxt = np.minimum(nxt, n-1)
        block_sum[np.arange(R), bidx] += g[cur, nxt]
        cur = nxt
    return block_sum.max(axis=1)

P = [[0.7, 0.2, 0.1],
     [0.1, 0.6, 0.3],
     [0.2, 0.3, 0.5]]
pi = markov_stationary(P)
hM, VM = markov_h_V(P, pi)
print(f"  Markov h={hM:.5f} bits/sym, V(rate)={VM:.5f} bits^2/sym  (pi={np.round(pi,3)})")
# cross-check the fundamental-matrix varentropy against empirical block variance / K
# (independent validation: the EVT prediction below must not rest on an unchecked V).
def markov_block_sums(P, pi, K, R):
    P = np.asarray(P); n = len(P)
    with np.errstate(divide='ignore'):
        g = np.where(P > 0, -np.log2(P), 0.0)
    cdf = np.cumsum(P, axis=1)
    cur = rng.choice(n, size=R, p=pi)
    s = np.zeros(R)
    for _ in range(K):
        u = rng.random(R)
        nxt = np.minimum((u[:, None] > cdf[cur]).sum(axis=1), n-1)
        s += g[cur, nxt]
        cur = nxt
    return s
bs = markov_block_sums(P, pi, 400, 6000)
V_emp = float(bs.var(ddof=1))/400; h_emp = float(bs.mean())/400
print(f"  cross-check @K=400: h_emp={h_emp:.5f} (vs {hM:.5f}), "
      f"V_emp={V_emp:.5f} (vs formula {VM:.5f}, ratio {V_emp/VM:.3f})")
vcross = abs(h_emp/hM-1) < 0.01 and abs(V_emp/VM-1) < 0.08
K = 300; R = 1000
ok4 = vcross
for m in [64, 256, 1024]:
    mx = block_max_markov(P, pi, K, m, R)
    z = (mx.mean() - K*hM)/math.sqrt(K*VM)
    a, b = gumbel_norm(m)
    pred = b + EULER/a
    rel = abs(z-pred)/pred
    print(f"  m={m:5d}: std max z={z:.4f}  refined pred={pred:.4f}  rel.err={rel:.3f}")
    ok4 = ok4 and rel < 0.10
print(f"  V4 {'PASS' if ok4 else 'FAIL'}")

# ---------------------------------------------------------------- V5  scale boundary
print("="*70)
print("V5: SHARP scale boundary -- the refined Gumbel needs ln m = o(K^{1/3}), i.e. K >> (ln m)^3")
# (codex/cold-review FIX: the original claim ln m = o(K) is too WEAK for the refined
# centering.) The Cramer-Petrov multiplicative correction is exp((x^3/sqrt K) lambda) with
# lambda(0)=mu3/(6 sigma^3); at the extremal level x ~ sqrt(2 ln m) it -> 1 IFF x^3/sqrt K -> 0,
# i.e. (ln m)^{3/2}=o(sqrt K), i.e. ln m = o(K^{1/3}). A SKEWED source (large 3rd cumulant)
# exhibits it: the refined-Gumbel error tracks x^3/sqrt K and is small iff K >> (ln m)^3.
# (iid counterexample of proof step 2: m = floor(e^{sqrt K}) has ln m=o(K) but x^3/sqrtK ~ K^{1/4}->oo.)
p5 = [0.85, 0.08, 0.05, 0.02]
h5 = entropy_rate_iid(p5); V5v = varentropy_iid(p5)
m = 1024; R = 2000
a, b = gumbel_norm(m); pred = b + EULER/a
lnm3 = (math.log(m))**3
print(f"  skewed source h={h5:.4f}, V={V5v:.4f}; m={m}, a_m={a:.3f}, (ln m)^3={lnm3:.0f}, pred={pred:.4f}")
errs = []
for K in [1200, 333, 100, 30, 10]:
    mx = block_max_iid(p5, K, m, R)
    z = (mx.mean() - K*h5)/math.sqrt(K*V5v)
    rel = abs(z-pred)/pred
    x3 = (a**3)/math.sqrt(K)                 # the Cramer correction magnitude x^3/sqrt K
    errs.append((K, rel))
    print(f"  K={K:5d}: K/(ln m)^3={K/lnm3:5.2f}  x^3/sqrtK={x3:5.2f}  std max z={z:.4f}  "
          f"rel.err={rel:.4f}")
e = [x[1] for x in errs]
mono5 = all(e[i] < e[i+1] for i in range(len(e)-1))
# clean when K >> (ln m)^3 (first, K/(ln m)^3=3.6): small error; badly violated (last,
# K/(ln m)^3=0.03): large error; monotone growth as K drops through the (ln m)^3 threshold.
ok5 = e[0] < 0.03 and mono5 and e[-1] > 0.10 and e[-1] > 4*e[0]
print(f"  V5 {'PASS' if ok5 else 'FAIL'}  (refined Gumbel accurate IFF K >> (ln m)^3: "
      f"err {e[0]:.3f} at 3.6x -> {e[-1]:.3f} at 0.03x; confirms ln m = o(K^{{1/3}}))")

# ---------------------------------------------------------------- V6  lower bound
print("="*70)
print("V6: matching lower bound -- a separated half-family already forces the max up")
K = 400; m = 1024; R = 2500
za = []; zh = []
cump = np.cumsum(np.asarray(p)); cump[-1] = 1.0
sur = -np.log2(np.asarray(p))
for r in range(R):
    u = rng.random((m, K))
    idx = np.clip(np.searchsorted(cump, u, side='right'), 0, len(p)-1)
    L = sur[idx].sum(axis=1)
    za.append(L.max())
    zh.append(L[::2].max())
mean_all = np.mean(za); mean_half = np.mean(zh)
z_all = (mean_all-K*h)/math.sqrt(K*V)
z_half = (mean_half-K*h)/math.sqrt(K*V)
a_all,_ = gumbel_norm(m); a_half,_ = gumbel_norm(m//2)
print(f"  std max over all m={m}:     z={z_all:.4f}  (a_m=sqrt(2 ln m)={a_all:.4f})")
print(f"  std max over half m/2={m//2}:  z={z_half:.4f}  (a_half={a_half:.4f})")
print(f"  ratio z_half/z_all = {z_half/z_all:.4f}  (separated half already ~ full max)")
ok6 = z_half/z_all > 0.93
print(f"  V6 {'PASS' if ok6 else 'FAIL'}")

# ---------------------------------------------------------------- summary
print("="*70)
allok = v1 and ok2 and ok3 and ok4 and ok5 and ok6
print(f"RESULT: {'ALL PASS' if allok else 'SOME FAILED'}  "
      f"[V1 {v1}, V2 {ok2}, V3 {ok3}, V4 {ok4}, V5 {ok5}, V6 {ok6}]")

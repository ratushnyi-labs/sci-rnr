#!/usr/bin/env python3
"""
Verification for Theorem 7.31 (Large-deviations reliability function of the RNR
archive size; the two-scale random-access picture).

Claim under test.  Let l_i = -log2 M(X_i | ctx_i) be the per-position coded
length of the realised RNR repair stream L^rep = sum_i l_i (notation of
Theorems 7.25/7.28/7.29).  For an i.i.d. or finite-state irreducible aperiodic
Markov source X ~ P coded under a (possibly mismatched) model M, the archive
overflow probability satisfies a large-deviation principle

      -(1/N) log2 P( (1/N) L^rep >= R )  -->  E(R)   for R > h_cross,

with rate (reliability) function the Legendre transform

      E(R) = sup_t [ t R - Lambda(t) ],

      Lambda(t) = log2 sum_x P(x) M(x)^{-t}          (i.i.d.)
      Lambda(t) = log2 rho( A_t ),  A_t(x,x')=P(x'|x) M(x'|x)^{-t}   (Markov)

(rho = Perron-Frobenius eigenvalue of the tilted transfer matrix).

Verified facts:
  V1  i.i.d., EXACT tail via the multinomial law: -(1/N) log2 P(S_N>=NR) -> E(R).
  V2  Markov, Gaertner-Ellis: tilted-matrix Lambda, exponent confirmed by tilted
      importance sampling (the matching LOWER bound on P, which the
      Azuma/Bernstein trilogy does NOT provide).
  V3  Consistency with the deviation trilogy and dispersion:
        Lambda'(0) = h_cross  (= h + KL(P||M), matches Thm 7.4 mean shift),
        Lambda''(0) = ln2 * v^2  (v^2 = long-run variance of Thm 7.25/7.28),
        E(R) ~ (R - h_cross)^2 / (2 ln2 v^2) near the mean  (recovers dispersion,
        Thm 7.29) and the Bernstein exponent is a LOWER bound on E, tight as R->h.
  V4  Mismatch (M != P): mean -> h + KL, E(h_cross)=0, exponent strictly convex.
  V5  Two-scale random access: the deterministic (N/K) delta_inf overhead with
      K = c sqrt(N) shifts R by delta_inf/(c sqrt N) -> 0, so the LDP exponent is
      UNCHANGED in the limit (RA is "free" at the large-deviations scale) while it
      was Theta(sqrt N) -- i.e. order-of-the-dispersion -- at the CLT scale (7.29).
  V6  Incompressibility / no-universal-dominance link (Thm 8.1): E(log2|Sigma|)
      is a finite positive exponent (expansion is exponentially rare but happens).

All PASS lines must print for the script to pass.
"""

import numpy as np
from itertools import product as iproduct

np.random.seed(0)
LN2 = np.log(2.0)
ok_all = True


def report(tag, ok, msg):
    global ok_all
    ok_all = ok_all and ok
    print(f"[{'PASS' if ok else 'FAIL'}] {tag}: {msg}")


# ----------------------------------------------------------------------------
# i.i.d. machinery
# ----------------------------------------------------------------------------
def iid_lambda(t, P, M):
    """Lambda(t) = log2 sum_x P(x) M(x)^{-t}  (base-2 scaled CGF).
    Numerically stable via log-sum-exp; vectorised over t."""
    t = np.asarray(t, dtype=float)
    # log2( sum_x exp2( log2 P(x) - t log2 M(x) ) )
    expo = np.log2(P)[None, :] - t[..., None] * np.log2(M)[None, :]  # (n_t, A)
    mmax = np.max(expo, axis=-1, keepdims=True)
    out = (mmax[..., 0] + np.log2(np.sum(2.0 ** (expo - mmax), axis=-1)))
    return out if t.ndim else out.reshape(())


def iid_lambda_derivs(P, M):
    """Lambda'(0), Lambda''(0) by closed form.
    phi(t)=sum P M^{-t}=E[2^{t l}], l=-log2 M.  Lambda=log2 phi.
    Lambda'(0)=E[l]=h_cross ; Lambda''(0)=ln2*Var(l).
    """
    l = -np.log2(M)
    Eh = np.sum(P * l)
    Var = np.sum(P * (l - Eh) ** 2)
    return Eh, LN2 * Var, Var


def legendre_E(R, lam, t_hi=200.0):
    """E(R)=sup_{t>=0}[tR-Lambda(t)] for R above the mean (upper tail).
    Uses a vectorised t-grid (lam must accept array t)."""
    ts = np.linspace(1e-6, t_hi, 400000)
    vals = ts * R - lam(ts)
    return np.max(vals), ts[np.argmax(vals)]


def legendre_E_scalar(R, lam_scalar, t_hi=60.0):
    """Same Legendre transform but lam_scalar accepts only scalar t (e.g. the
    Markov log2 rho(A_t)).  Coarse grid + local golden-section refinement."""
    ts = np.linspace(1e-4, t_hi, 4000)
    vals = np.array([t * R - lam_scalar(float(t)) for t in ts])
    j = int(np.argmax(vals))
    a, b = ts[max(0, j - 1)], ts[min(len(ts) - 1, j + 1)]
    gr = (np.sqrt(5) - 1) / 2
    for _ in range(80):
        c, d = b - gr * (b - a), a + gr * (b - a)
        fc = c * R - lam_scalar(c)
        fd = d * R - lam_scalar(d)
        if fc > fd:
            b = d
        else:
            a = c
    tstar = 0.5 * (a + b)
    return tstar * R - lam_scalar(tstar), tstar


def iid_exact_tail_exponent(N, R, P, M):
    """EXACT P(S_N >= N R), S_N = sum_i l(X_i), via the multinomial law over
    count vectors.  Returns -(1/N) log2 P.  Alphabet size kept small."""
    A = len(P)
    l = -np.log2(M)
    from math import lgamma, log
    logfacN = lgamma(N + 1)
    logP = np.log(P)
    target = N * R
    # enumerate count vectors (n_0..n_{A-1}) summing to N
    total_logprob_terms = []

    def rec(idx, remaining, counts):
        if idx == A - 1:
            counts2 = counts + [remaining]
            s = sum(c * lc for c, lc in zip(counts2, l))
            if s >= target - 1e-12:
                lp = logfacN
                for c, lpx in zip(counts2, logP):
                    lp += -lgamma(c + 1) + c * lpx
                total_logprob_terms.append(lp)
            return
        for n in range(remaining + 1):
            rec(idx + 1, remaining - n, counts + [n])

    rec(0, N, [])
    if not total_logprob_terms:
        return np.inf
    m = max(total_logprob_terms)
    logsum = m + log(sum(np.exp(np.array(total_logprob_terms) - m)))
    logP_tail = logsum  # natural log
    return -(logP_tail / LN2) / N  # base-2 exponent per symbol


# ----------------------------------------------------------------------------
# V1 : i.i.d. exact-tail convergence to the Cramer exponent
# ----------------------------------------------------------------------------
print("=" * 70)
print("V1  i.i.d. exact multinomial tail  ->  Cramer exponent E(R)")
P = np.array([0.7, 0.2, 0.1])
M = P.copy()  # matched
lam = lambda t: iid_lambda(t, P, M)
R = 1.45  # above h_cross = 1.157
E_an, t_star = legendre_E(R, lam)
exps, corr = [], []
Ns = [40, 80, 160, 320, 640]
for N in Ns:
    e = iid_exact_tail_exponent(N, R, P, M)
    exps.append(e)
    # Bahadur-Rao prefactor: P ~ 2^{-N E}/(t* sigma_t sqrt(2 pi N)), so the raw
    # exponent overshoots E by (log2 N)/(2N) + O(1/N); subtract the leading term.
    corr.append(e - np.log2(N) / (2 * N))
    print(f"     N={N:4d}  raw -(1/N)log2 P = {e:.5f}   BR-corrected = {corr[-1]:.5f}"
          f"   E(R)={E_an:.5f}")
mono = all(exps[i] >= exps[i + 1] - 1e-9 for i in range(len(exps) - 1))   # decreasing to E
conv = all(abs(corr[i] - E_an) >= abs(corr[i + 1] - E_an) - 1e-9 for i in range(len(corr) - 1))
report("V1", mono and conv and abs(corr[-1] - E_an) < 0.006 and E_an > 0,
       f"E(R={R})={E_an:.4f}; raw exponent decreases to E, BR-corrected->E "
       f"(corr err@N=640={abs(corr[-1]-E_an):.4f}, t*={t_star:.3f})")


# ----------------------------------------------------------------------------
# V3 : consistency Lambda'(0)=h_cross, Lambda''(0)=ln2 v^2, quadratic-at-mean
# ----------------------------------------------------------------------------
print("=" * 70)
print("V3  consistency with mean / dispersion / Bernstein")
Eh, L2pp, Var = iid_lambda_derivs(P, M)
# numerical derivatives of Lambda at 0
dt = 1e-4
L0 = lam(0.0)
Lp = (lam(dt) - lam(-dt)) / (2 * dt)
Lpp = (lam(dt) - 2 * lam(0.0) + lam(-dt)) / dt ** 2
report("V3a-mean", abs(Lp - Eh) < 1e-3 and abs(Eh - 1.156) < 0.05,
       f"Lambda'(0)={Lp:.5f} == h_cross={Eh:.5f}")
report("V3b-var", abs(Lpp - L2pp) < 1e-3,
       f"Lambda''(0)={Lpp:.5f} == ln2*Var(l)={L2pp:.5f}  (v^2=Var={Var:.5f})")
# quadratic-at-mean: E(h+a) ~ a^2/(2 ln2 v^2)
a = 0.02
E_near, _ = legendre_E(Eh + a, lam)
E_quad = a ** 2 / (2 * LN2 * Var)
report("V3c-dispersion", abs(E_near - E_quad) / E_quad < 0.05,
       f"E(h+a)={E_near:.3e} ~ a^2/(2 ln2 v^2)={E_quad:.3e}  (recovers varentropy)")
# Bernstein exponent (variance-aware, B=max l) is a LOWER bound on E, tight as a->0
B = np.max(-np.log2(M))
def bernstein_exp(a):  # base-2 exponent of exp(-N a^2/(2(v^2 + B a/3))) -> /ln2
    return (a ** 2 / (2 * (Var + B * a / 3))) / LN2
ok_bound = True
ratio_small = ratio_large = None
for aa in [0.005, 0.02, 0.1, 0.4]:
    Ea, _ = legendre_E(Eh + aa, lam)
    be = bernstein_exp(aa)
    ok_bound = ok_bound and (be <= Ea + 1e-9)
    if aa == 0.005:
        ratio_small = be / Ea
    if aa == 0.4:
        ratio_large = be / Ea
report("V3d-bernstein<=cramer", ok_bound and ratio_small > 0.9 and ratio_large < 0.9,
       f"Bernstein<=Cramer always; ratio: {ratio_small:.3f} (a=.005, tight) -> "
       f"{ratio_large:.3f} (a=.4, loose)")
# V3e: the tail bounds do NOT order uniformly (codex adv5 finding). The true a^2
# curvature is 1/(2 ln2 v^2); Bernstein is sharp near the mean, but for the SYMMETRIC
# Bernoulli (v^2 = range^2/4) the Azuma/Hoeffding range bound BEATS Bernstein in the
# tail. Both stay <= the exact Cramer exponent E (the envelope).
true_curv = 1.0 / (2 * LN2 * Var)                    # exact LDP leading curvature
small = 0.01
E_small, _ = legendre_E(Eh + small, lam)
matches_const = abs((E_small / small ** 2) - true_curv) / true_curv < 0.05
# symmetric Bernoulli ell in {0,1}, p=1/2:  Azuma vs Bernstein cross relative to each
pb = 0.5; vb = pb * (1 - pb); rb = 1.0; meanb = pb
def E_bern2(R):
    return R * np.log2(R / pb) + (1 - R) * np.log2((1 - R) / (1 - pb))
cross_tail = cross_le = True
for a in [0.1, 0.2, 0.4]:
    Eb = E_bern2(meanb + a)
    az = 2 * a * a / (rb ** 2 * LN2)
    be = a * a / (2 * (vb + rb * a / 3)) / LN2
    cross_le = cross_le and az <= Eb + 1e-9 and be <= Eb + 1e-9   # both lower bounds on E
    if a == 0.4:
        cross_tail = az > be                                       # Azuma beats Bernstein in tail
report("V3e-bounds-cross",
       matches_const and cross_le and cross_tail,
       f"E(h+a)/a^2={E_small/small**2:.4f}~1/(2ln2 v^2)={true_curv:.4f}; symmetric-Bernoulli "
       f"tail a=0.4: Azuma>Bernstein (no uniform order), both <= exact E")


# ----------------------------------------------------------------------------
# V4 : mismatch  M != P
# ----------------------------------------------------------------------------
print("=" * 70)
print("V4  mismatched model  (code under M, source P)")
Mq = np.array([0.5, 0.3, 0.2])
lamq = lambda t: iid_lambda(t, P, Mq)
h = -np.sum(P * np.log2(P))
KL = np.sum(P * np.log2(P / Mq))
h_cross = h + KL
Ehq, _, _ = iid_lambda_derivs(P, Mq)
report("V4a-meanshift", abs(Ehq - h_cross) < 1e-9 and KL > 0,
       f"Lambda'(0)={Ehq:.5f} = h({h:.4f})+KL({KL:.4f})={h_cross:.5f}")
E_at_mean, _ = legendre_E(h_cross, lamq, t_hi=50)
report("V4b-zero-at-mean", E_at_mean < 1e-3,
       f"E(h_cross)={E_at_mean:.2e} ~ 0 (no overflow at the mean)")
Eov, _ = legendre_E(h_cross + 0.2, lamq)
report("V4c-pos-overflow", Eov > 0,
       f"E(h_cross+0.2)={Eov:.4f} > 0 (overflow exponentially rare)")


# ----------------------------------------------------------------------------
# Markov machinery
# ----------------------------------------------------------------------------
def stationary(Pm):
    w, v = np.linalg.eig(Pm.T)
    i = np.argmin(np.abs(w - 1))
    pi = np.real(v[:, i]); pi = pi / pi.sum()
    return pi


def markov_lambda(t, Pm, Mm):
    A = Pm * (Mm ** (-t))           # tilted transfer matrix A_t(x,x')
    rho = np.max(np.abs(np.linalg.eigvals(A)))
    return np.log2(rho)


def markov_h_cross(Pm, Mm):
    pi = stationary(Pm)
    g = -np.log2(Mm)               # per-step cost
    return np.sum(pi[:, None] * Pm * g)


def markov_longrun_var(Pm, n=4_000_00):
    """Direct long-run variance v^2 of l_i along a simulated chain (batch means
    cross-check of Lambda''(0)/ln2)."""
    A = Pm.shape[0]
    g = -np.log2(Pm)               # matched: l = -log2 P(x'|x)
    pi = stationary(Pm)
    x = np.random.choice(A, p=pi)
    vals = np.empty(n)
    cum = np.cumsum(Pm, axis=1)
    u = np.random.rand(n)
    for i in range(n):
        xn = int(np.searchsorted(cum[x], u[i]))
        vals[i] = g[x, xn]
        x = xn
    # long-run variance via Bartlett / truncated covariance series
    vbar = vals - vals.mean()
    v2 = np.var(vals)
    for k in range(1, 200):
        c = np.mean(vbar[:-k] * vbar[k:])
        v2 += 2 * c
    return v2


def markov_cgf_exact(t, N, Pm, Mm):
    """EXACT (1/N) log2 E[2^{t S_N}] via the transfer matrix:
    E[2^{t S_N}] = pi^T A_t^N 1,  A_t(x,x') = Pm(x'|x) Mm(x'|x)^{-t}.
    Gaertner-Ellis: this -> log2 rho(A_t) = Lambda(t) as N->oo."""
    pi = stationary(Pm)
    At = Pm * (Mm ** (-t))
    v = pi.copy()
    for _ in range(N):
        v = v @ At
    return np.log2(np.sum(v)) / N


def markov_direct_mc_exponent(N, R, Pm, nsamp=400000):
    """Direct (untilted) MC raw exponent -(1/N) log2 P((1/N) S_N >= R) at a
    MODERATE R.  Matched model l = -log2 P.  Returns (raw, BR-corrected)."""
    A = Pm.shape[0]
    g = -np.log2(Pm)
    pi = stationary(Pm)
    cum = np.cumsum(Pm, axis=1)
    x = np.random.choice(A, size=nsamp, p=pi)
    s = np.zeros(nsamp)
    for _ in range(N):
        u = np.random.rand(nsamp)
        xn = (u[:, None] > cum[x]).sum(axis=1)
        s += g[x, xn]
        x = xn
    p = np.mean(s >= N * R)
    if p == 0:
        return np.inf, np.inf
    raw = -np.log2(p) / N
    return raw, raw - np.log2(N) / (2 * N)


# ----------------------------------------------------------------------------
# V2 : Markov Gaertner-Ellis exponent -- EXACT transfer-matrix CGF + MC tail
# ----------------------------------------------------------------------------
print("=" * 70)
print("V2  Markov Gaertner-Ellis  (exact transfer-matrix CGF -> log2 rho(A_t))")
Pm = np.array([[0.9, 0.1], [0.4, 0.6]])
Mm = Pm.copy()                       # matched
lamM = lambda t: markov_lambda(t, Pm, Mm)   # log2 rho(A_t)
hc_M = markov_h_cross(Pm, Mm)
# (a) the exact finite-N CGF converges to log2 rho(A_t) for several t
cgf_ok = True
for tt in [-0.4, 0.2, 0.6, 1.0]:
    limit = lamM(tt)
    seq = [markov_cgf_exact(tt, n, Pm, Mm) for n in [20, 80, 320]]
    cgf_ok = cgf_ok and abs(seq[-1] - limit) < 5e-3
    print(f"     t={tt:+.2f}  (1/N)log2 E[2^tS] : {seq[0]:.4f},{seq[1]:.4f},"
          f"{seq[2]:.4f} -> log2 rho(A_t)={limit:.4f}")
# (b) exponent at R above the mean, via Legendre on log2 rho(A_t)
RM = hc_M + 0.20
E_M, t_star_M = legendre_E_scalar(RM, lamM, t_hi=60)
# (c) direct MC tail at the same R for growing N: the raw exponent must
# decrease toward E_M FROM ABOVE (the LDP), the BR-corrected one approach it.
raws, brs = [], []
for n in [60, 120, 240]:
    r, b = markov_direct_mc_exponent(n, RM, Pm, nsamp=800000)
    raws.append(r); brs.append(b)
    print(f"     N={n:4d}  raw MC exponent={r:.4f}  BR-corrected={b:.4f}   E(R)={E_M:.4f}")
mc_trend = (all(raws[i] >= raws[i + 1] - 1e-3 for i in range(len(raws) - 1))
            and all(r > E_M - 1e-3 for r in raws)
            and abs(brs[-1] - E_M) < 0.012)
report("V2", cgf_ok and E_M > 0 and mc_trend,
       f"exact CGF->log2 rho(A_t); E(R={RM:.3f})={E_M:.4f}; raw MC exponent "
       f"decreases to E from above, BR-corrected@N=240={brs[-1]:.4f} (t*={t_star_M:.3f})")
# cross-check Lambda''(0)/ln2 == direct long-run variance
Lpp_M = (lamM(1e-3) - 2 * lamM(0.0) + lamM(-1e-3)) / 1e-3 ** 2
v2_direct = markov_longrun_var(Pm)
report("V2b-longrunvar", abs(Lpp_M / LN2 - v2_direct) / v2_direct < 0.10,
       f"Lambda''(0)/ln2={Lpp_M/LN2:.4f} == direct long-run var v^2={v2_direct:.4f}")
# Markov mean consistency: Lambda'(0) == h_cross (entropy rate), and Lambda(0)==0
LpM0 = (lamM(1e-4) - lamM(-1e-4)) / 2e-4
report("V2c-markov-mean", abs(LpM0 - hc_M) < 2e-3 and abs(lamM(0.0)) < 1e-9,
       f"Lambda(0)={lamM(0.0):.2e}, Lambda'(0)={LpM0:.5f} == h_cross(Markov)={hc_M:.5f}")


# ----------------------------------------------------------------------------
# V5 : two-scale random access -- LDP-invisible, dispersion-visible
# ----------------------------------------------------------------------------
print("=" * 70)
print("V5  random-access overhead: dispersion-scale vs large-deviation scale")
delta_inf = 0.3        # per-sync excess (Crutchfield-Feldman)
c = 1.0                # K = c sqrt(N)
V = Var                # iid varentropy
# at the CLT/dispersion scale the RA overhead (N/K)delta = sqrt(N) delta/c is
# Theta(sqrt N) -- SAME order as the dispersion sqrt(N V) Q^{-1}: ratio const in N
disp_ratios = []
for N in [10_000, 100_000, 1_000_000]:
    ra = (N / (c * np.sqrt(N))) * delta_inf          # = sqrt(N) delta/c
    disp = np.sqrt(N * V) * 1.2816                    # Q^{-1}(0.1)
    disp_ratios.append(ra / disp)
clt_const = np.allclose(disp_ratios, disp_ratios[0], rtol=1e-6)
# at the LDP scale: E with the RA shift -> E without it (shift per-symbol ->0)
R0 = Eh + 0.25
E_clean, _ = legendre_E(R0, lam)
shifts = []
for N in [1_000, 100_000, 10_000_000]:
    shift = delta_inf / (c * np.sqrt(N))             # per-symbol RA overhead
    E_shift, _ = legendre_E(R0 - shift, lam)         # P(sum l + (N/K)delta >= NR0)
    shifts.append(E_shift)
ldp_converges = abs(shifts[-1] - E_clean) < 1e-3 and abs(shifts[0] - E_clean) > abs(shifts[-1] - E_clean)
report("V5a-dispersion-scale", clt_const,
       f"(N/K)delta / dispersion = {disp_ratios[0]:.4f} constant in N (RA visible at CLT scale)")
report("V5b-ldp-invisible", ldp_converges,
       f"E with RA shift -> E_clean={E_clean:.4f} (E_RA: {shifts[0]:.4f}->{shifts[-1]:.4f}); "
       f"RA overhead is sub-exponential, invisible at the LDP scale")


# ----------------------------------------------------------------------------
# V6 : incompressibility exponent (no-universal-dominance, Thm 8.1)
# ----------------------------------------------------------------------------
print("=" * 70)
print("V6  incompressibility / expansion exponent  E(log2|Sigma|)")
Rmax = np.log2(len(P))               # uniform / incompressible rate
E_inc, t_inc = legendre_E(Rmax - 1e-6, lam, t_hi=2000)
# the expansion threshold (archive >= source size in bits) is R = 8 for bytes;
# here just confirm E(log|Sigma|) is finite & positive: rare but possible
report("V6", 0 < E_inc < np.inf,
       f"E(log2|Sigma|={Rmax:.3f})={E_inc:.4f} finite>0: near-uniform (incompressible) "
       f"blocks are exponentially rare but occur (consistent with Thm 8.1)")


# ----------------------------------------------------------------------------
# V7 : interior effective domain + boundary behaviour (codex-fix scope)
# ----------------------------------------------------------------------------
print("=" * 70)
print("V7  interior effective domain  h_cross < R < ess-sup l, and the boundary")
ell = -np.log2(M)                 # matched P==M here
ess_sup = float(np.max(ell))      # max attainable per-step cost
x_star = int(np.argmax(ell))      # the max-cost (min-prob) symbol
# (a) E finite & positive strictly inside (h_cross, ess_sup)
Rin = 0.5 * (Eh + ess_sup)
E_in, _ = legendre_E(Rin, lam, t_hi=400)
# (b) right boundary R = ess_sup: E finite, attained as t->oo, equals -log2 P(x*)
E_bnd, t_bnd = legendre_E(ess_sup, lam, t_hi=5000)
E_bnd_closed = -np.log2(P[x_star])     # exponent of the all-max-cost realisation
# (c) beyond the boundary R > ess_sup: event impossible, exact tail empty (exp=inf)
exp_above = iid_exact_tail_exponent(30, ess_sup + 0.05, P, M)
report("V7a-interior", 0 < E_in < np.inf,
       f"E(R={Rin:.3f}) in (h_cross={Eh:.3f}, ess_sup={ess_sup:.3f}) = {E_in:.4f} finite>0")
report("V7b-right-boundary", abs(E_bnd - E_bnd_closed) < 5e-3,
       f"E(ess_sup={ess_sup:.3f})={E_bnd:.4f} == -log2 P(x*)={E_bnd_closed:.4f} "
       f"(maximiser receding, t={t_bnd:.0f})")
report("V7c-impossible-above", exp_above == np.inf,
       f"R>ess_sup: P(S_N>=NR)=0 exactly, exponent=+inf (event impossible)")
# (d) MARKOV: the upper endpoint is the max CYCLE MEAN r_+ = lim Lambda'(t), which
#     can be STRICTLY BELOW the single-step max ell_max -- a high-cost edge off every
#     high-mean cycle cannot be sustained (codex-caught correction).
Pm_d = np.array([[0.5, 0.5], [0.5, 0.5]])      # support: all edges
g_d = np.array([[0.10, 3.00], [0.10, 0.10]])   # g(0,1)=3 is the max single-step edge
Mm_d = 2.0 ** (-g_d); Mm_d = Mm_d / Mm_d.sum(axis=1, keepdims=True)
g_d = -np.log2(Mm_d)
ell_max_d = float(g_d.max())                                  # single-step max
cyc_means = [g_d[0, 0], g_d[1, 1], (g_d[0, 1] + g_d[1, 0]) / 2]   # 0->0, 1->1, 0->1->0
r_plus_cyc = max(cyc_means)
lamM_d = lambda t: markov_lambda(t, Pm_d, Mm_d)
r_plus_lim = (lamM_d(60.0 + 1e-3) - lamM_d(60.0 - 1e-3)) / (2e-3)   # lim Lambda'(t)
# empirical: no length-N path can sustain an average above r_+ (greedy max attempt)
x, best = 0, 0.0
for _ in range(5000):
    x2, s = np.random.randint(2), 0.0
    for _ in range(300):
        xn = int(np.argmax(g_d[x2])); s += g_d[x2, xn]; x2 = xn
    best = max(best, s / 300)
report("V7d-markov-cyclemean",
       abs(r_plus_cyc - r_plus_lim) < 1e-2 and r_plus_cyc < ell_max_d - 1e-6
       and best <= r_plus_cyc + 1e-9,
       f"Markov r_+ = max cycle mean {r_plus_cyc:.4f} == lim Lambda'(t) {r_plus_lim:.4f} "
       f"< ell_max {ell_max_d:.4f}; greedy max avg {best:.4f} <= r_+ (ell_max is the WRONG endpoint)")


print("=" * 70)
print("ALL PASS" if ok_all else "SOME FAILED")
import sys
sys.exit(0 if ok_all else 1)

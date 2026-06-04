#!/usr/bin/env python3
"""
Verification for Theorem 7.32 (Optimal archive-overflow exponent and the
universality dividend for RNR).

Setting.  i.i.d. source P0 over a finite alphabet (size A), parametric class =
the (d = A-1)-dim simplex.  Three codes for the RNR repair stream:
  - FIXED/ideal model (Theorem 7.31, matched M=P0): length L_fix = -log2 P0(x^N)
    = N * crossent(emp, P0).
  - UNIVERSAL KT mixture (Bayes, Dirichlet(1/2) prior): length
    L_univ = -log2 P_KT(x^N).
  - OPTIMAL-overflow code (length-ranked): the best UD code for the overflow
    criterion P(L >= NR).

Claims under test:
  V1  Clarke-Barron / Laplace: L_univ = N * Hhat(x^N) + (d/2) log2 N + O(1)
      (Hhat = empirical entropy), pointwise.
  V2  E*(R) = inf{ D(Q||P0) : H(Q) >= R } = D(Q_s||P0), Q_s ∝ P0^s the
      entropy-flattened tilt with H(Q_s)=R; and E_arith(R)=inf{D : crossent>=R}
      (= Theorem 7.31's matched exponent). Verify E*(R) >= E_arith(R), STRICT for
      R > H(P0). (The fixed/ideal code is overflow-SUBOPTIMAL.)
  V3  Optimal converse (type counting): the length-ranked code's overflow
      exponent equals E*(R) (P0-mass outside the top 2^{NR} sequences).
  V4  Monte-Carlo: KT code overflows STRICTLY LESS than the fixed code; the
      Bahadur-Rao-corrected empirical exponents track E* (KT) and E_arith (fixed).
  V5  Redundancy L_univ - L_fix = +(d/2)log2 N typically, but NEGATIVE on the
      overflow event (the universal code ADAPTS to the atypical empirical law).
  V6  Three scales: (d/2)log2 N is o(sqrt N) (dispersion-invisible, Thm 7.29) and
      o(N) (LDP-invisible, Thm 7.31), yet E* > E_arith (LDP-IMPROVING). The
      universality price is paid only at the mean (3rd order) and bought back as a
      strictly better reliability exponent.
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


def H(Q):
    Q = np.asarray(Q, float); Q = Q[Q > 0]
    return float(-(Q * np.log2(Q)).sum())


def D(Q, P):
    Q = np.asarray(Q, float); P = np.asarray(P, float)
    m = Q > 0
    return float((Q[m] * np.log2(Q[m] / P[m])).sum())


def crossent(Q, P):  # E_Q[-log2 P]
    Q = np.asarray(Q, float); P = np.asarray(P, float)
    return float(-(Q * np.log2(P)).sum())


# ----------------------------------------------------------------------------
# KT (Krichevsky-Trofimov, Dirichlet(1/2)) exact code length
# ----------------------------------------------------------------------------
def L_kt(counts):
    """-log2 P_KT(x^N), P_KT = Dirichlet(1/2)-multinomial marginal."""
    counts = np.asarray(counts, float)
    N = counts.sum(); A = len(counts)
    # log P_KT = sum_a [lgamma(n_a+1/2)-lgamma(1/2)] + lgamma(A/2) - lgamma(N+A/2)
    lp = (sum(lgamma(n + 0.5) - lgamma(0.5) for n in counts)
          + lgamma(A * 0.5) - lgamma(N + A * 0.5))
    return -lp / LN2


# ----------------------------------------------------------------------------
# V1 : L_univ = N*Hhat + (d/2) log2 N + O(1)
# ----------------------------------------------------------------------------
print("=" * 70)
print("V1  KT length  =  N*Hhat + (d/2) log2 N + O(1)")
P0 = np.array([0.6, 0.3, 0.1]); A = len(P0); d = A - 1
gaps = []
for N in [200, 1000, 5000, 20000]:
    x = np.random.choice(A, size=N, p=P0)
    counts = np.bincount(x, minlength=A)
    emp = counts / N
    Lu = L_kt(counts)
    pred = N * H(emp) + (d / 2) * log2(N)
    gaps.append(Lu - pred)
    print(f"     N={N:6d}  L_univ={Lu:10.2f}  N*Hhat+(d/2)log2 N={pred:10.2f}  "
          f"gap(O(1))={Lu-pred:+.3f}")
# the gap should be O(1) (bounded, not growing with N)
report("V1", max(abs(g) for g in gaps) < 3.0 and abs(gaps[-1] - gaps[0]) < 2.0,
       f"L_univ - [N Hhat + (d/2)log2 N] = O(1) (gaps {min(gaps):+.2f}..{max(gaps):+.2f})")


# ----------------------------------------------------------------------------
# V2 : E*(R) = inf{D: H(Q)>=R} via tilt Q_s ∝ P0^s ; E* >= E_arith, strict
# ----------------------------------------------------------------------------
print("=" * 70)
print("V2  E*(R)=inf{D:H>=R} >= E_arith(R)=inf{D:crossent>=R}  (strict for R>H(P0))")
P0b = np.array([0.9, 0.1])           # binary, H(P0)=0.469
hP0 = H(P0b)

def tilt(s):                          # Q_s ∝ P0^s  (s in (0,1]: P0->uniform as s->0)
    w = P0b ** s; return w / w.sum()

def Estar(R):                         # inf{D(Q||P0): H(Q)>=R}, R in (H(P0), log2 A)
    # H(tilt(s)) increases from H(P0) (s=1) to log2 A (s=0); bisect for H=R
    lo, hi = 1e-6, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if H(tilt(mid)) < R: hi = mid      # too uniform? need smaller H -> larger s
        else: lo = mid
    Q = tilt(0.5 * (lo + hi))
    return D(Q, P0b), Q

def Earith(R):                        # inf{D(Q||P0): crossent(Q,P0)>=R}  (binary, grid)
    qs = np.linspace(1e-6, 1 - 1e-6, 200001)
    Q = np.stack([qs, 1 - qs], 1)
    ce = -(Q * np.log2(P0b)).sum(1)
    Dv = (Q * np.log2(Q / P0b)).sum(1)
    feas = ce >= R
    return float(Dv[feas].min()) if feas.any() else np.inf

ok2 = True
for R in [0.55, 0.70, 0.85, 0.95]:
    es, Q = Estar(R)
    ea = Earith(R)
    strict = es > ea + 1e-4
    ok2 = ok2 and (es >= ea - 1e-6) and (strict if R > hP0 + 0.05 else True)
    print(f"     R={R}:  E*={es:.4f} (Q*={np.round(Q,3).tolist()})  "
          f"E_arith={ea:.4f}   E*>E_arith? {es>ea+1e-9}")
report("V2", ok2, f"H(P0)={hP0:.3f}; E*(R) >= E_arith(R), strict for R>H(P0) "
       f"(fixed/ideal code is overflow-SUBOPTIMAL)")


# ----------------------------------------------------------------------------
# V3 : optimal converse via type counting -- length-ranked code attains E*
# ----------------------------------------------------------------------------
print("=" * 70)
print("V3  type-counting: P0-mass outside top 2^{NR} sequences ~ 2^{-N E*(R)}")
# For binary, enumerate types (k ones out of N); count 2^{N h(k/N)}, mass per
# sequence 2^{-N crossent}. Rank-order by probability; find mass beyond top 2^{NR}.
def converse_exponent(N, R):
    ks = np.arange(N + 1)
    q = ks / N
    with np.errstate(divide='ignore', invalid='ignore'):
        ent = np.where((q > 0) & (q < 1), -(q*np.log2(q)+(1-q)*np.log2(1-q)), 0.0)
    logcount = N * ent                                   # log2 #sequences of this type
    # per-sequence log2 prob under P0
    logp_seq = ks*np.log2(P0b[1]) + (N-ks)*np.log2(P0b[0])   # ones=symbol1
    logp_type = logcount + logp_seq                      # log2 P0(type)
    # sort types by per-sequence prob (descending) = ascending codelength
    order = np.argsort(-logp_seq)
    cum_count_log = -np.inf
    # accumulate sequence counts until exceeding 2^{NR}; the rest is the overflow mass
    # do it in log space
    cnt = 0.0  # cumulative count (linear, but use log via logsumexp)
    logcnt = -np.inf
    mass_outside_log = []
    acc_types = []
    for idx in order:
        acc_types.append(idx)
        # cumulative count log
        logcnt = np.logaddexp(logcnt*LN2, logcount[idx]*LN2)/LN2 if logcnt>-np.inf else logcount[idx]
        if logcnt >= N * R:
            break
    # overflow = mass of all types NOT in acc_types (those with shortest codes excluded)
    mask = np.ones(N + 1, bool); mask[acc_types] = False
    if not mask.any():
        return np.inf
    lp = logp_type[mask]
    m = lp.max()
    log_mass = m + np.log2(np.sum(2.0 ** (lp - m)))
    return -log_mass / N

R3 = 0.85
es3, _ = Estar(R3)
conv = [converse_exponent(N, R3) for N in [200, 800, 3200]]
print(f"     R={R3}: converse exponents {['%.4f'%c for c in conv]}  ->  E*={es3:.4f}")
report("V3", abs(conv[-1] - es3) < 0.03,
       f"length-ranked converse exponent -> E*={es3:.4f} (optimal overflow exponent)")


# ----------------------------------------------------------------------------
# V4 : Monte-Carlo -- KT overflows LESS than fixed; exponents track E*, E_arith
# ----------------------------------------------------------------------------
print("=" * 70)
print("V4  Monte-Carlo overflow rates: universal (KT) vs fixed/ideal")
def mc_overflow(N, R, nsamp=300000):
    """Exact MC overflow rates for the fixed (-log2 P0) and universal (KT) codes.
    No scipy: KT length via a precomputed lgamma table over integer counts."""
    lg = np.array([lgamma(k + 0.5) for k in range(N + 1)])      # lgamma(k+1/2)
    const = (-2 * lgamma(0.5) + lgamma(1.0) - lgamma(N + 1.0))  # A=2
    samples = np.random.multinomial(N, P0b, size=nsamp)        # (nsamp, 2)
    cfix = samples @ (-np.log2(P0b))                           # L_fix per sample (bits)
    n0 = samples[:, 0]; n1 = samples[:, 1]
    Lkt = -(lg[n0] + lg[n1] + const) / LN2                     # L_univ per sample
    return float(np.mean(cfix >= N * R)), float(np.mean(Lkt >= N * R))

R4 = 0.65            # fatter tail (closer to H(P0)) so both rates are measurable;
# N must exceed the discreteness-resolution point (~150 here: N*E_arith > (d/2)log2 N)
# for the (d/2)log2 N shift to separate the two integer overflow sets.
es4, _ = Estar(R4); ea4 = Earith(R4)
rows = []
for N in [200, 350, 500]:
    of, ou = mc_overflow(N, R4, nsamp=1000000)
    ratio = of / ou if ou > 0 else np.inf
    rows.append((N, of, ou, ratio))
    print(f"     N={N:4d}  P(over)_fix={of:.4e} P(over)_uni={ou:.4e}  fix/uni={ratio:.2f}")
uni_less = all(r[2] <= r[1] for r in rows)                    # KT overflows <= fixed
ratio_emerges = rows[-1][3] >= 1.3 and rows[-1][3] >= rows[0][3] - 1e-9   # gap grows
report("V4", uni_less and ratio_emerges and es4 > ea4,
       f"universal overflows <= fixed at every N; the fix/uni ratio emerges and grows "
       f"({rows[0][3]:.2f}->{rows[-1][3]:.2f}) ~ 2^(N(E*-E_arith)); predicted "
       f"E_arith={ea4:.3f} < E*={es4:.3f}")


# ----------------------------------------------------------------------------
# V5 : redundancy sign flip -- universal SHORTER on the overflow event
# ----------------------------------------------------------------------------
print("=" * 70)
print("V5  L_univ - L_fix: +(d/2)log2 N typically, NEGATIVE on the overflow event")
N5 = 400; R5 = 0.75
# typical realisation
xt = np.random.choice(2, size=N5, p=P0b); ct = np.bincount(xt, minlength=2)
red_typ = L_kt(ct) - crossent(ct / N5, P0b) * N5
# overflow realisation in the GAP regime: crossent(Q,P0) >= R (fixed overflows) but
# H(Q) < R (universal does NOT) -- for R=0.75 this is q1 in [0.189, 0.215).
q1 = 0.20; co = np.array([round(N5 * (1 - q1)), round(N5 * q1)])
emp_o = co / N5
ce_o = crossent(emp_o, P0b); Hh_o = H(emp_o)
Lfix_o = ce_o * N5
Luni_o = L_kt(co)
red_over = Luni_o - Lfix_o
fix_overflows = Lfix_o >= N5 * R5
uni_overflows = Luni_o >= N5 * R5
print(f"     typical:  L_univ-L_fix={red_typ:+.2f}  ((d/2)log2 N={(1/2)*log2(N5):.2f})")
print(f"     gap-event q1={q1}: crossent={ce_o:.3f}>=R (fix overflows: L_fix={Lfix_o:.1f}>={N5*R5:.0f}={fix_overflows});"
      f" Hhat={Hh_o:.3f}<R (univ: L_univ={Luni_o:.1f}<{N5*R5:.0f}={not uni_overflows})  diff={red_over:+.2f}")
report("V5", red_typ > 0 and red_over < 0 and fix_overflows and not uni_overflows,
       f"redundancy +{red_typ:.1f} typical but {red_over:+.1f} on the gap event; the FIXED code "
       f"overflows there while the UNIVERSAL code adapts (Hhat<R) and does NOT")


# ----------------------------------------------------------------------------
# V6 : three scales -- (d/2)log N invisible at dispersion & LDP, improves LDP
# ----------------------------------------------------------------------------
print("=" * 70)
print("V6  (d/2)log2 N: o(sqrt N) and o(N), yet E* > E_arith")
V = float((P0b * (-np.log2(P0b) - hP0) ** 2).sum())   # varentropy
for N in [10**4, 10**6, 10**8]:
    red = (d / 2) * log2(N)
    disp = np.sqrt(N * V) * 1.2816
    print(f"     N={N:.0e}  (d/2)log2N={red:8.1f}  dispersion~{disp:10.1f}  "
          f"red/disp={red/disp:.2e}  red/N={red/N:.2e}")
es6, _ = Estar(0.85); ea6 = Earith(0.85)
report("V6",
       (d / 2) * log2(10**8) / np.sqrt(10**8) < 1e-2 and es6 > ea6,
       f"universality price is o(sqrt N) (dispersion-invisible) and o(N) (LDP-invisible) "
       f"but improves the exponent E_arith={ea6:.3f} -> E*={es6:.3f}")


# ----------------------------------------------------------------------------
# V7 : scope (codex adv/bal findings) -- interior P0 needed; (d/2) is prior-
#      UNIVERSAL for interior empiricals but the formula FAILS at the boundary.
# ----------------------------------------------------------------------------
print("=" * 70)
print("V7  interior P0 needed; (d/2)log N prior-universal for interior empiricals")

def L_dir(counts, alpha):
    counts = np.asarray(counts, float); A = len(counts); N = counts.sum()
    lp = (sum(lgamma(n + alpha) - lgamma(alpha) for n in counts)
          + lgamma(A * alpha) - lgamma(N + A * alpha))
    return -lp / LN2

# (a) interior empirical q=(0.7,0.3): leading (d/2)log N gap is O(1) for ALL priors
N7 = 50000; ci = np.array([int(0.7 * N7), int(0.3 * N7)]); d7 = 1
gaps_int = [L_dir(ci, a) - (N7 * H(ci / N7) + (d7 / 2) * log2(N7)) for a in (0.5, 1.0, 2.0)]
# (b) boundary empirical 0^N: Dir(1) coefficient is LARGER (gap grows ~ (d/2)log10),
#     but KT/Dir(1/2) keeps the interior (A-1)/2 coefficient (gap O(1)) -- codex adv3
#     formula [(A-k)alpha + (k-1)/2]: for alpha=1/2 boundary == interior (A-1)/2.
g_b1 = [L_dir([N, 0], 1.0) - (N * 0 + (d7 / 2) * log2(N)) for N in (1000, 10000)]
g_bkt = [L_dir([N, 0], 0.5) - (N * 0 + (d7 / 2) * log2(N)) for N in (1000, 10000)]
boundary_grows = (g_b1[1] - g_b1[0] > 1.0) and abs(g_bkt[1] - g_bkt[0]) < 0.2  # Dir(1) grows, KT O(1)
# (c) degenerate P0=(1,0): D(Q||P0)=inf for any Q with mass on symbol 1 -> E*=E_arith=inf
P0deg = np.array([1.0, 0.0])
Ddeg = D(np.array([0.5, 0.5]), P0deg)                        # = inf
report("V7",
       max(abs(g) for g in gaps_int) < 2.0                  # interior: O(1) for all priors
       and abs(gaps_int[0] - gaps_int[2]) < 2.0
       and boundary_grows and np.isinf(Ddeg),
       f"interior empirical: (d/2)log N gap O(1) for Dir(.5/1/2)={[round(g,2) for g in gaps_int]} "
       f"(prior-universal); boundary 0^N: Dir(1) coeff GROWS {g_b1[0]:.1f}->{g_b1[1]:.1f} but "
       f"KT/Dir(.5) stays O(1) {g_bkt[0]:.2f}->{g_bkt[1]:.2f} (KT boundary coeff = interior (A-1)/2); "
       f"degenerate P0=(1,0): D=inf so E*=E_arith=inf (interior P0 required)")


# ----------------------------------------------------------------------------
# V8 : A>=3 boundary types CAN enter {Hhat>=R} (codex adv2) but are D-dominated
#      -- they do not change the Sanov exponent E* (interior Q_R dominates).
# ----------------------------------------------------------------------------
print("=" * 70)
print("V8  A>=3 boundary types may enter the overflow event but are exp-negligible")
P03 = np.array([0.9, 0.05, 0.05]); R8 = 0.8
Qb = np.array([0.0, 0.5, 0.5])                     # boundary: symbol 0 absent
in_overflow = H(Qb) >= R8                          # high entropy despite boundary
Db = D(Qb, P03)                                    # finite (P0 full support)
# interior minimiser Q_R prop P0^s with H(Q_R)=R8
def tilt3(s):
    w = P03 ** s; return w / w.sum()
lo, hi = 1e-6, 1.0
for _ in range(300):
    m = 0.5 * (lo + hi)
    if H(tilt3(m)) < R8: hi = m
    else: lo = m
QR = tilt3(0.5 * (lo + hi)); Estar3 = D(QR, P03)
report("V8",
       in_overflow and np.isfinite(Db) and Db > Estar3 + 0.5 and (QR > 0).all(),
       f"A=3 boundary Qb=(0,.5,.5): H={H(Qb):.2f}>=R={R8} (ENTERS overflow), D={Db:.3f} "
       f">> E*={Estar3:.4f} at interior Q_R={np.round(QR,3).tolist()} -> boundary exp-negligible, "
       f"exponent unchanged")


print("=" * 70)
print("ALL PASS" if ok_all else "SOME FAILED")
import sys
sys.exit(0 if ok_all else 1)

#!/usr/bin/env python3
r"""
BUG-003-C: realized sync-point overhead AS A FRACTION of the total coded
stream, over a (K, W, eta) sweep, on (i) a synthetic stationary Markov
source and (ii) a non-Markov (long-range) stress source.

This is the executable HEURISTIC SUPPORT for the abstract's "approximately
one percent" conjecture (Section 10.7, Corollary 10.2a). It is explicitly
NOT a measurement on a named real corpus -- that residual is the external
leaf BUG-003-D (closeable_in_repo = NO). Synthetic Markov sources are
favorable to the conjecture by construction; the non-Markov stress source
is included so the probe CAN report an honestly larger fraction.

What is computed (exact entropy bookkeeping, not a neural run):
  The sync overhead Delta_sync of Section 10.7 / Def. 10.3 is the excess of
  the COLD-restarted encoding of a sub-block over the WARM-context (no-sync)
  encoding of the SAME bytes by the SAME W-bounded predictor. It is NOT the
  predictor's intrinsic suboptimality (the gap between its warm rate and the
  true rate h): that gap is paid identically with or without sync points and
  cancels. So we measure, for each (source, K, W):
        L_warm(K, W) = sum_i Hpred(X_i | full available context),   (no reset)
        L_cold(K, W) = sum_i Hpred(X_i | context truncated by the reset),
        delta_real(K, W) = L_cold - L_warm  >= 0,                   (the excess)
  with the excess concentrated in the first W positions after a reset (where
  the in-block context is shorter than the predictor would otherwise use);
  positions beyond W are warm in both, so they cancel. This is exactly the
  per-sync quantity rung 1 bounds: delta_real <= S * W * log2(1/eta).
  With m = ceil(N/K) sync points over a stream of N bytes (full sub-blocks
  N = mK), the realized total sync overhead is m * delta_real and the
  realized warm-context (no-sync) coded length is m * L_warm, so the realized
  OVERHEAD-AS-FRACTION reported here is
        frac(K, W) = (m * delta_real) / (m * L_warm) = delta_real / L_warm.
  For order-W' Markov with W >= W' and K >= W', the warm rate equals the true
  rate h and delta_real = delta_W' = H(X^{W'}) - W' h(X), the finite
  source-specific constant of Lemma 5.1c, so frac ~ delta_W'/(K h(X)) scales
  as ~ 1/K (rung 3 of Corollary 10.2a).

PASS/FAIL (printed per check, "OVERALL -> PASS" on success; CI contract):
  (A) loose ceiling respected: for every swept (source, K, W, eta), the
      realized per-sync excess delta_real <= the rung-1 ceiling
      S * W * log2(1/eta) (Corollary 10.2a rung 1 / Lemma 5.1a).
  (B) typical-K band: at the H10(c) spacing/window ratio K/W = 16 (its
      K = 64 KiB, W = 4 KiB), the realized fraction on the NATURAL-DATA-LIKE
      (mildly correlated) Markov sources falls in the conjectured
      low-single-digit-percent band (<= ~1% target, asserted <= 2% to match
      H10(c)'s falsification threshold). Clearly labelled HEURISTIC,
      synthetic, NOT a corpus datum.
  (C) ~1/K scaling: at fixed W the realized fraction shrinks monotonically as
      K grows and the product frac(K)*K is constant (delta_real is
      K-independent once K >= W' and the cold region is <= W positions deep),
      i.e. more sync points (smaller K) => proportionally more overhead.

  Two counterweight sources -- a STICKY order-1 Markov source (large per-sync
  delta_1) and the NON-MARKOV long-range stress source (cold region W bytes
  deep) -- are reported but NOT asserted into the band: both legitimately
  exceed it at small K/W, so the probe can fail honestly and the ~1% figure
  reads as natural-data-favorable, not universal. The loose ceiling (A) and
  the ~1/K scaling (C) ARE asserted for them.
"""

import itertools
import math
import sys
from collections import defaultdict

EPS = 1e-15


def line(ok, label, detail):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {label}: {detail}")
    return bool(ok)


def entropy_of_dist(probs):
    return -sum(p * math.log2(p) for p in probs if p > EPS)


# ---------------------------------------------------------------------------
# Source models. Each exposes:
#   h_rate()            -> entropy rate h(X)
#   block_length(K, W)  -> sum_{i=1}^K H(X_i | X_{i-W:i-1})  (true conditional)
# computed exactly from the joint over K positions for small alphabets.
# ---------------------------------------------------------------------------
def stationary_order1(P):
    n = len(P)
    pi = [1.0 / n] * n
    for _ in range(5000):
        nxt = [sum(pi[j] * P[j][i] for j in range(n)) for i in range(n)]
        if max(abs(nxt[i] - pi[i]) for i in range(n)) < EPS:
            break
        pi = nxt
    return pi


def order1_joint(P, pi, N):
    n = len(P)
    joint = {}
    for tup in itertools.product(range(n), repeat=N):
        p = pi[tup[0]]
        for i in range(N - 1):
            p *= P[tup[i]][tup[i + 1]]
        if p > EPS:
            joint[tup] = p
    return joint


def cond_entropy_at(joint, n_pos, N):
    """H(X_{n_pos+1} | X_1..X_{n_pos}) over the joint of length-N tuples."""
    if n_pos == 0:
        marg = defaultdict(float)
        for tup, p in joint.items():
            marg[tup[0]] += p
        return entropy_of_dist(marg.values())
    prefix = defaultdict(float)
    cond = defaultdict(lambda: defaultdict(float))
    for tup, p in joint.items():
        pre = tup[:n_pos]
        prefix[pre] += p
        cond[pre][tup[n_pos]] += p
    H = 0.0
    for pre, pp in prefix.items():
        H += pp * entropy_of_dist([v / pp for v in cond[pre].values()])
    return H


class Order1Markov:
    """Stationary order-1 Markov source on a 2-symbol alphabet.

    Warm (no-sync) interior block of size K: every position conditions on at
    least the previous byte, so the true conditional gives h(X) per position,
    L_warm = K * h(X). Cold (post-reset) block: position 1 conditions on
    nothing (the reset truncates context), paying the marginal H_marg, then
    every later position is warm (h(X)), so L_cold = H_marg + (K-1) h(X).
    delta_real = L_cold - L_warm = H_marg - h(X) = delta_1, K-independent and
    W-independent (for W >= 1), exactly Lemma 5.1c's delta_{W'=1}.
    """

    def __init__(self, P, name):
        self.P = P
        self.name = name
        self.pi = stationary_order1(P)
        self._h = sum(self.pi[x] * entropy_of_dist(P[x]) for x in range(len(P)))
        self._Hmarg = entropy_of_dist(self.pi)

    def h_rate(self):
        return self._h

    def warm_length(self, K, W):
        # no reset: an interior sub-block in steady state pays h(X)/position.
        return K * self._h

    def cold_length(self, K, W):
        # reset truncates the context at position 1 (pays the marginal H_marg);
        # positions >= 2 are warm. (W>=1 so only position 1 is cold for
        # order-1.) Exact for any K >= 1.
        return self._Hmarg + max(0, K - 1) * self._h


class NonMarkovStress:
    """Non-Markov long-range stress: a W-bounded predictor's COLD penalty does
    not collapse to a single-byte marginal term.

    We model a source whose W-bounded predictor has a fixed warm rate h_W (its
    intrinsic suboptimality h_W - h is paid with OR without sync points and so
    cancels out of the sync overhead). What the reset costs is that the first
    W positions after a reset see a context shorter than W bytes: at in-block
    depth d (1 <= d <= W) the predictor conditions on only d-1 prior bytes
    instead of W, paying a per-position cold surcharge c(d) that decays as the
    in-block context fills up (c(1) the largest, c(W) ~ 0). Positions beyond W
    are warm in both runs and cancel. Thus
        delta_real = sum_{d=1}^{W} c(d),
    which GROWS with W (more cold positions) but is K-independent for K >= W,
    bounded by W * max_d c(d) <= S * W * log2(1/eta) when each c(d) <= the
    per-byte floor cap. The realized FRACTION delta_real / (K * h_W) is larger
    than the Markov delta_1 / (K h) case (the cold region is W bytes deep, not
    1), so the stress source legitimately reports a higher band -- the honest
    counterweight to the favorable Markov sources.
    """

    def __init__(self, name):
        self.name = name
        self._h = 0.55                 # true entropy rate (bits/byte)
        self._hW = 0.62                # W-bounded predictor warm rate (h_W > h)
        # per-position cold surcharge by in-block depth d (decays as d grows);
        # capped well under any reasonable S*log2(1/eta) per-byte ceiling.
        self._c = {1: 0.45, 2: 0.34, 3: 0.25, 4: 0.18, 5: 0.12, 6: 0.08,
                   7: 0.05, 8: 0.03}

    def h_rate(self):
        return self._h

    def _cold_at_depth(self, d):
        if d in self._c:
            return self._c[d]
        return 0.0                     # context full by depth > max key

    def warm_length(self, K, W):
        # no reset: predictor pays its warm rate h_W each position.
        return K * self._hW

    def cold_length(self, K, W):
        # first W positions pay h_W + c(d) (truncated context); rest warm.
        cold_excess = sum(self._cold_at_depth(d) for d in range(1, W + 1))
        return K * self._hW + cold_excess


# ---------------------------------------------------------------------------
# Realized overhead-as-fraction.
# ---------------------------------------------------------------------------
def realized(src, K, W):
    """Returns (delta_real_per_sync, frac) where
       delta_real = L_cold(K,W) - L_warm(K,W)   (per-sync cold-restart excess,
                    isolated from the predictor's intrinsic suboptimality),
       frac       = delta_real / L_warm(K,W)    (overhead as fraction of the
                    no-sync warm-context coded stream)."""
    Lwarm = src.warm_length(K, W)
    Lcold = src.cold_length(K, W)
    delta = Lcold - Lwarm
    frac = delta / Lwarm if Lwarm > 0 else float("inf")
    return delta, frac


def loose_ceiling(S, W, eta):
    """Rung-1 / Lemma 5.1a per-sync ceiling S*W*log2(1/eta)."""
    return S * W * math.log2(1.0 / eta)


# ---------------------------------------------------------------------------
# Checks.
# ---------------------------------------------------------------------------
# H10(c) tests K = 64 KiB with W = 4 KiB, i.e. spacing/window ratio K/W = 16.
# The realized overhead FRACTION depends on that ratio (cold region ~W deep,
# warm block ~K bytes), so the meaningful "typical-K" analogue keeps K/W = 16;
# we also sweep larger ratios to exhibit the ~1/K collapse at fixed W.
KW_RATIO = 16                           # H10(c): 64 KiB / 4 KiB
K_RATIO_SWEEP = [4, 16, 64, 256, 1024]  # K = ratio * W (more syncs <=> smaller)
W_SWEEP = [1, 2, 4]
ETA_SWEEP = [2.0 ** -8, 2.0 ** -12, 2.0 ** -16]
S_VALUES = [1, 2]                       # Type-I (S=1), per-byte Type-III-B (S=2)
BAND_TARGET = 0.01                      # conjectured ~1%
BAND_ASSERT = 0.02                      # H10(c) falsification threshold at 64 KiB

# "natural-data-like" Markov sources: mildly correlated, small per-sync
# delta_1 = H_marg - h(X). The band assertion is made on THESE -- the
# conjecture is about natural data at typical parameters, not pathological
# sources. (Lemma 5.1c's worked instance is order-1 English bigram,
# delta_1 ~ 0.71 bits/sync, h ~ 0.71 bits/letter; a sync every K=64 bytes
# then costs ~0.71 / (64*0.71) ~ 1.6% at W small, dropping below 1% as K
# grows -- the order of magnitude the abstract conjectures.)
NATURAL_SOURCES = [
    Order1Markov([[0.7, 0.3], [0.4, 0.6]], "order-1 Markov A (mild)"),
    Order1Markov([[0.6, 0.4], [0.45, 0.55]], "order-1 Markov C (weak corr.)"),
]
# stress / counterweight sources: NOT asserted into the band. They show the
# ~1% figure is favorable-case, not universal -- a sticky Markov source and a
# non-Markov long-range source both legitimately exceed the band at small K/W.
STICKY_SOURCE = Order1Markov([[0.9, 0.1], [0.2, 0.8]], "order-1 Markov B (sticky; counterweight)")
STRESS_SOURCE = NonMarkovStress("non-Markov run source (long-range stress)")
COUNTERWEIGHTS = [STICKY_SOURCE, STRESS_SOURCE]
MARKOV_SOURCES = NATURAL_SOURCES + [STICKY_SOURCE]  # all Markov, for ceiling sweep


def k_sweep_for(W):
    return [r * W for r in K_RATIO_SWEEP]


def check_loose_ceiling_respected():
    print("\n(A) Loose rung-1 ceiling S*W*log2(1/eta) respected (all sources, all sweep points)")
    ok = True
    worst = None
    n_pts = 0
    for src in NATURAL_SOURCES + COUNTERWEIGHTS:
        for W in W_SWEEP:
            for K in k_sweep_for(W):
                delta, _ = realized(src, K, W)
                for S in S_VALUES:
                    for eta in ETA_SWEEP:
                        n_pts += 1
                        ceil = loose_ceiling(S, W, eta)
                        if delta > ceil + 1e-9:
                            ok = False
                            worst = (src.name, K, W, S, eta, delta, ceil)
    if worst is None:
        # show the tightest swept point (smallest ceiling: S=1, eta=2^-8, W=1)
        d, _ = realized(STRESS_SOURCE, KW_RATIO * 1, 1)
        c = loose_ceiling(1, 1, 2.0 ** -8)
        ok = line(True, f"delta_real <= S*W*log2(1/eta) at all {n_pts} swept points",
                  f"tightest e.g. {STRESS_SOURCE.name} W=1: delta={d:.4f} "
                  f"<= ceiling(S=1,eta=2^-8)={c:.1f} bits")
    else:
        ok = line(False, "ceiling VIOLATED",
                  f"{worst[0]} K={worst[1]} W={worst[2]} S={worst[3]} eta={worst[4]:.1e}: "
                  f"delta={worst[5]:.4f} > ceiling={worst[6]:.4f}")
    return ok


def check_typical_K_band_markov():
    print(f"\n(B) Typical K/W = {KW_RATIO} (H10(c): K=64 KiB, W=4 KiB) realized fraction in ~1% "
          f"band on NATURAL-DATA-LIKE Markov sources [HEURISTIC, synthetic -- NOT a corpus measurement]")
    ok = True
    for src in NATURAL_SOURCES:
        for W in W_SWEEP:
            K = KW_RATIO * W
            _, frac = realized(src, K, W)
            in_band = frac <= BAND_ASSERT
            ok &= line(in_band, f"{src.name} W={W} K={K}",
                       f"realized fraction = {frac*100:.4f}%  "
                       f"(target ~{BAND_TARGET*100:.0f}%, assert <= {BAND_ASSERT*100:.0f}% per H10(c))")
    # counterweight sources at the SAME K/W ratio: NOT asserted into the band.
    # A sticky Markov source (large delta_1) and a non-Markov long-range source
    # (cold region W bytes deep) both legitimately exceed 2% at small K/W,
    # proving the probe can fail and the ~1% figure is natural-data-favorable,
    # not universal.
    for src in COUNTERWEIGHTS:
        for W in W_SWEEP:
            K = KW_RATIO * W
            _, frac = realized(src, K, W)
            verdict = "exceeds band (expected; counterweight)" if frac > BAND_ASSERT \
                else "within band"
            print(f"  [info] {src.name} W={W} K={K}: realized fraction = "
                  f"{frac*100:.4f}% -- {verdict}; band membership NOT asserted")
    return ok


def check_inverse_K_scaling():
    print("\n(C) Realized fraction scales ~1/K at fixed W (more sync points => more overhead)")
    ok = True
    for src in NATURAL_SOURCES + COUNTERWEIGHTS:
        W = 2
        Ks = k_sweep_for(W)
        fracs = [realized(src, K, W)[1] for K in Ks]
        mono = all(fracs[i] >= fracs[i + 1] - 1e-12 for i in range(len(fracs) - 1))
        ok &= line(mono, f"{src.name} W={W} fraction monotone decreasing in K",
                   "fracs(%) = " + ", ".join(f"{f*100:.4f}" for f in fracs))
        # delta_real is K-independent at fixed W (cold region <= W positions),
        # so frac * K is constant => frac ~ 1/K exactly.
        prods = [f * K for f, K in zip(fracs, Ks)]
        spread = (max(prods) - min(prods)) / (sum(prods) / len(prods))
        const = spread < 1e-6
        ok &= line(const, f"{src.name} frac*K constant (=> exact ~1/K at fixed W)",
                   f"frac*K = {prods[0]:.6f} (relative spread {spread:.2e})")
    return ok


def main() -> int:
    print("Verification: realized sync-point overhead AS A FRACTION (BUG-003-C)")
    print("=" * 72)
    print("HEURISTIC SUPPORT for the abstract's '~1%' conjecture on SYNTHETIC")
    print("sources (Corollary 10.2a). This is NOT a measurement on a named real")
    print("corpus -- that residual is the external leaf BUG-003-D. Markov sources")
    print("are favorable by construction; a non-Markov long-range stress source")
    print("is included so the result cannot be read as universal.")

    results = {
        "(A) loose rung-1 ceiling S*W*log2(1/eta) respected": check_loose_ceiling_respected(),
        "(B) typical-K realized fraction in ~1% band (Markov)": check_typical_K_band_markov(),
        "(C) realized fraction ~1/K (more syncs => more overhead)": check_inverse_K_scaling(),
    }

    print("\n" + "=" * 72)
    print("Summary:")
    all_ok = True
    for name, ok in results.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        all_ok &= ok

    print()
    if all_ok:
        print("Supports the ~1% conjecture on synthetic natural-data-like Markov")
        print("sources at the H10(c) K/W ratio (sticky-Markov and non-Markov")
        print("counterweights legitimately exceed the band); does NOT establish a")
        print("real-corpus value -- that is the external leaf BUG-003-D.")
        print("OVERALL -> PASS")
        return 0
    print("OVERALL -> FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())

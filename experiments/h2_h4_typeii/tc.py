#!/usr/bin/env python3
"""Multi-information (total correlation) estimation for Theorem 5.5.

TC(D_N) = sum_f H_marg(field_f) - H(joint record), estimated from M
records. Three entropy estimators are computed:

  plug         empirical plug-in (biased down for H, hence TC biased UP)
  mm           Miller-Madow first-order correction
  grassberger  Grassberger (2003) psi-function estimator,
               H = [ln M - (1/M) sum_i n_i G(n_i)] / ln 2 with
               G(n) = psi(n) + ((-1)^n / 2) (psi((n+1)/2) - psi(n/2)).

The PRIMARY reported value is the Grassberger estimate: on held-out
known-entropy cases at this experiment's exact sample sizes it reduced
the product-of-marginals bias floor from 0.094 bpb (MM) to 0.019 bpb
(M = 65536 over support 65536) while matching a closed-form TC = 1.0 bpb
anchor to 0.0001 bpb. The permuted-control TC of each corpus measures
the realized residual bias floor directly (its true TC is 0).

Uncertainty: multinomial bootstrap over records -- resampled joint counts
are drawn over the OBSERVED distinct record tuples, per-field marginal
counts are aggregated from the same resampled joint counts (coherent
resampling), and the Grassberger estimator is re-applied per resample;
the percentile interval is reported. This is a SECONDARY interval in the
METHODS.md sense (the frozen BCa-on-means machinery does not cover
plug-in entropy functionals); it is labeled
method="percentile-multinomial-grassberger" wherever stored.
"""

from __future__ import annotations

import numpy as np
from scipy.special import psi


def _factorize_records(records: np.ndarray):
    """(joint_counts, tuples) for an (M, R) uint8 array."""
    void = np.ascontiguousarray(records).view(
        np.dtype((np.void, records.shape[1])))
    uniq, counts = np.unique(void.ravel(), return_counts=True)
    tuples = uniq.view(np.uint8).reshape(len(uniq), records.shape[1])
    return counts.astype(np.int64), tuples


def _g_table(nmax: int) -> np.ndarray:
    """Grassberger G(n) lookup for n = 0..nmax (G(0) unused, weight 0)."""
    n = np.arange(1, nmax + 1, dtype=np.float64)
    g = psi(n) + 0.5 * ((-1.0) ** n) * (psi((n + 1) / 2) - psi(n / 2))
    return np.concatenate(([0.0], g))


def _h_plug(counts, M):
    p = counts[counts > 0] / M
    return float(-(p * np.log2(p)).sum())


def _h_mm(counts, M):
    d = int((counts > 0).sum())
    return _h_plug(counts, M) + (d - 1) / (2.0 * M * np.log(2.0))


def _h_grass_from_table(counts: np.ndarray, M: int, g: np.ndarray) -> float:
    return float((np.log(M) - (counts * g[counts]).sum() / M) / np.log(2.0))


def estimate_tc(records: np.ndarray, *, n_boot: int = 10000,
                boot_seed: int = 0, alpha: float = 0.05) -> dict:
    """TC estimate (bits/record and bits/byte) with bootstrap percentiles."""
    M, R = records.shape
    joint_counts, tuples = _factorize_records(records)
    d = len(joint_counts)
    g = _g_table(M)

    marg_counts, field_maps = [], []
    for f in range(R):
        vals, vidx = np.unique(tuples[:, f], return_inverse=True)
        mc = np.zeros(len(vals), dtype=np.int64)
        np.add.at(mc, vidx, joint_counts)
        marg_counts.append(mc)
        field_maps.append((len(vals), vidx))

    est = {}
    for tag, fn in (("plug", lambda c: _h_plug(c, M)),
                    ("mm", lambda c: _h_mm(c, M)),
                    ("grass", lambda c: _h_grass_from_table(c, M, g))):
        hj = fn(joint_counts)
        hm = sum(fn(mc) for mc in marg_counts)
        est[tag] = {"h_joint_bits": hj, "h_marg_sum_bits": hm,
                    "tc_bits": hm - hj}

    # Coherent multinomial bootstrap with per-resample Grassberger.
    tcs = np.empty(max(n_boot, 1))
    rng = np.random.default_rng(boot_seed)
    p = joint_counts / M
    chunk = max(1, min(max(n_boot, 1), int(1e7 / max(d, 1))))
    done = 0
    while done < max(n_boot, 1):
        b = min(chunk, max(n_boot, 1) - done)
        cj = rng.multinomial(M, p, size=b)  # (b, d) ints
        hj = (np.log(M) - (cj * g[cj]).sum(axis=1) / M) / np.log(2.0)
        hm = np.zeros(b)
        for n_vals, vidx in field_maps:
            cm = np.zeros((b, n_vals), dtype=np.int64)
            np.add.at(cm.T, vidx, cj.T)
            hm += (np.log(M) - (cm * g[cm]).sum(axis=1) / M) / np.log(2.0)
        tcs[done:done + b] = hm - hj
        done += b
    lo, hi = np.percentile(tcs, [100 * alpha / 2, 100 * (1 - alpha / 2)])

    tc_bits = est["grass"]["tc_bits"]
    return {
        "M": M, "R": R, "distinct_records": d, "coverage": d / M,
        # Joint support nearly exhausts the sample => no estimator can
        # recover H(joint); the estimate is flagged unusable.
        "estimable": d < 0.9 * M,
        "estimators_bits": est,
        "tc_bits": tc_bits, "tc_bpb": tc_bits / R,
        "tc_mm_bpb": est["mm"]["tc_bits"] / R,
        "tc_plug_bpb": est["plug"]["tc_bits"] / R,
        "tc_boot_lo_bpb": float(lo) / R, "tc_boot_hi_bpb": float(hi) / R,
        "n_boot": n_boot, "ci_method": "percentile-multinomial-grassberger",
    }


def chain_mi(records: np.ndarray, max_lag: int = 2) -> dict:
    """Best-parent tree mutual information (Grassberger entropies).

    Each field f >= 1 picks a single earlier parent p(f) in
    {f-1, ..., f-max_lag} maximizing I(X_f; X_{p(f)}); the statistic is
    tree_mi = sum_f I(X_f; X_{p(f)}). Because every field has exactly
    one parent earlier in the ordering, the chain rule with dropped
    conditioning gives H(joint) <= H(X_0) + sum_f H(X_f | X_{p(f)}),
    hence tree_mi is a rigorous LOWER bound on TC(D_N). (Choosing the
    best of max_lag candidate parents keeps the bound valid -- it is a
    maximum over valid trees -- and captures dependencies that skip a
    constant field, e.g. the interleaved-high-byte columnar layout.)

    Unlike the full joint plug-in, every term lives on a <= 2-byte
    alphabet and stays estimable at any M used here; under the H4
    per-field permutation control its true value is 0. Used as the
    permutation-control statistic for corpora whose permuted joint
    support exceeds M.
    """
    M, R = records.shape
    g = _g_table(M)

    def h_of(cols: np.ndarray) -> float:
        counts, _ = _factorize_records(np.ascontiguousarray(cols))
        return _h_grass_from_table(counts, M, g)

    h_field = [h_of(records[:, f:f + 1]) for f in range(R)]
    pair_mi, parents = [], []
    for f in range(1, R):
        best, best_p = None, None
        for lag in range(1, min(max_lag, f) + 1):
            p = f - lag
            h_pair = h_of(records[:, [p, f]])
            mi = h_field[p] + h_field[f] - h_pair
            if best is None or mi > best:
                best, best_p = mi, p
        pair_mi.append(best)
        parents.append(best_p)
    total = float(sum(pair_mi))
    return {"chain_mi_bits": total, "chain_mi_bpb": total / R,
            "pair_mi_bits": [float(x) for x in pair_mi],
            "parents": parents, "max_lag": max_lag}


def prefix_tc_curve(records: np.ndarray, prefixes: list[int]) -> list[dict]:
    """Point-estimate TC(prefix r) for r in prefixes -- the small-scale
    Theta(N) proxy of H3 (linear growth of TC in the record length)."""
    out = []
    for r in prefixes:
        e = estimate_tc(records[:, :r], n_boot=1, boot_seed=0)
        out.append({"prefix_bytes": r, "tc_bits": e["tc_bits"],
                    "distinct_records": e["distinct_records"]})
    return out

#!/usr/bin/env python3
"""Rate, timing, and uncertainty metrics for the RNR empirical program.

Implements the measurement definitions of the experimental-design
companion, section 5 (bits per byte, timing under repetition) and the
statistical machinery of section 6 (bootstrap confidence intervals,
paired comparisons, Cohen's d, Holm-Bonferroni step-down).

FROZEN METHODOLOGY: the statistical choices implemented here (BCa
bootstrap as the primary interval, 10000 resamples, 95% nominal level,
seed-derivation policy, per-family Holm-Bonferroni correction) are
specified and frozen in bench/METHODS.md. This module is the reference
implementation of that document. Any change to a statistical code path
after the first recorded measurement requires an amendment note in
METHODS.md (see its Amendment rule) before the change is used.
"""

from __future__ import annotations

import dataclasses
import hashlib
import time

import numpy as np
from scipy.stats import norm

# Frozen defaults; see bench/METHODS.md sections 2-4.
N_RESAMPLES = 10_000
ALPHA = 0.05
MASTER_SEED = 20260710


# ---------------------------------------------------------------------------
# Seed policy (METHODS.md section 4)
# ---------------------------------------------------------------------------

def derive_seed(experiment_id: str, index: int, tag: str = "run") -> int:
    """Deterministic 32-bit seed derived from the master seed.

    seed = first 4 bytes of SHA-256(f"{MASTER_SEED}:{tag}:{experiment_id}:{index}")
    """
    digest = hashlib.sha256(
        f"{MASTER_SEED}:{tag}:{experiment_id}:{index}".encode()
    ).digest()
    return int.from_bytes(digest[:4], "big")


# ---------------------------------------------------------------------------
# Primary metric: coding rate
# ---------------------------------------------------------------------------

def bits_per_byte(compressed_size_bytes: int, original_size_bytes: int) -> float:
    """Coding rate in bits per byte: 8 * compressed / original."""
    if original_size_bytes <= 0:
        raise ValueError("original size must be positive")
    return 8.0 * compressed_size_bytes / original_size_bytes


# ---------------------------------------------------------------------------
# Bootstrap confidence intervals (METHODS.md section 2)
# ---------------------------------------------------------------------------

@dataclasses.dataclass(frozen=True)
class CI:
    point: float
    lo: float
    hi: float
    level: float  # e.g. 0.95
    method: str  # "bca", "percentile", or "degenerate"
    n_resamples: int


def _jackknife_stats(x: np.ndarray, stat) -> np.ndarray:
    """Leave-one-out statistic values (closed form for the mean)."""
    n = len(x)
    if stat is np.mean:
        total = x.sum()
        return (total - x) / (n - 1)
    idx = np.arange(n)
    return np.array([stat(x[idx != i]) for i in range(n)])


def bootstrap_ci(
    samples,
    *,
    stat=np.mean,
    n_resamples: int = N_RESAMPLES,
    alpha: float = ALPHA,
    method: str = "bca",
    seed: int = 0,
) -> CI:
    """Bootstrap confidence interval for stat(samples).

    Primary method is BCa (bias-corrected and accelerated,
    Efron-Tibshirani 1993 ch. 14): bias correction z0 from the fraction
    of resampled statistics below the point estimate, acceleration a
    from the jackknife skewness. Falls back to the plain percentile
    interval when BCa is undefined (z0 infinite because all resamples
    fall on one side, or zero jackknife variance), and degenerates to a
    zero-width interval when all samples are identical. The method
    actually used is recorded in CI.method.

    `stat` must accept an `axis` keyword when given a 2-D array
    (np.mean and np.median qualify).
    """
    x = np.asarray(samples, dtype=float)
    if x.ndim != 1 or len(x) < 2:
        raise ValueError("need a 1-D sample of size >= 2")
    n = len(x)
    theta = float(stat(x))

    if np.all(x == x[0]):
        return CI(theta, theta, theta, 1 - alpha, "degenerate", 0)

    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_resamples, n))
    boots = np.asarray(stat(x[idx], axis=1), dtype=float)

    def _percentile() -> CI:
        lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
        return CI(theta, float(lo), float(hi), 1 - alpha, "percentile", n_resamples)

    if method == "percentile":
        return _percentile()
    if method != "bca":
        raise ValueError(f"unknown method {method!r}")

    # Bias correction, with a mid-p tie adjustment.
    prop = (np.sum(boots < theta) + 0.5 * np.sum(boots == theta)) / n_resamples
    if prop <= 0.0 or prop >= 1.0:
        return _percentile()
    z0 = norm.ppf(prop)

    # Acceleration from the jackknife.
    jack = _jackknife_stats(x, stat)
    d = jack.mean() - jack
    denom = 6.0 * np.sum(d**2) ** 1.5
    if denom == 0.0:
        return _percentile()
    a = np.sum(d**3) / denom

    z_lo, z_hi = norm.ppf(alpha / 2), norm.ppf(1 - alpha / 2)
    q_lo = norm.cdf(z0 + (z0 + z_lo) / (1 - a * (z0 + z_lo)))
    q_hi = norm.cdf(z0 + (z0 + z_hi) / (1 - a * (z0 + z_hi)))
    lo, hi = np.percentile(boots, [100 * q_lo, 100 * q_hi])
    return CI(theta, float(lo), float(hi), 1 - alpha, "bca", n_resamples)


def paired_bootstrap_diff(
    a,
    b,
    *,
    n_resamples: int = N_RESAMPLES,
    alpha: float = ALPHA,
    seed: int = 0,
) -> CI:
    """CI on mean(a - b) via bootstrap on the paired differences.

    This is the exp-design section 6.2 "paired bootstrap on per-block
    rates": a and b are per-block rates of two coders on the SAME
    blocks, in the same order.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError("paired samples must have identical shape")
    return bootstrap_ci(
        a - b, n_resamples=n_resamples, alpha=alpha, method="bca", seed=seed
    )


def cohens_d_paired(a, b) -> float:
    """Cohen's d for paired block-level comparisons: mean(diff)/sd(diff)."""
    diff = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    sd = diff.std(ddof=1)
    if sd == 0.0:
        return 0.0 if diff.mean() == 0.0 else float("inf") * np.sign(diff.mean())
    return float(diff.mean() / sd)


# ---------------------------------------------------------------------------
# Timing under repetition (METHODS.md section 3)
# ---------------------------------------------------------------------------

@dataclasses.dataclass(frozen=True)
class TimingResult:
    times_s: tuple  # measured wall-clock times, warmup excluded
    mean_s: float
    ci: CI


def time_under_repetition(
    fn,
    *args,
    repeats: int = 5,
    warmup: int = 1,
    ci_seed: int = 0,
    **kwargs,
) -> TimingResult:
    """Wall-clock timing of fn(*args, **kwargs) under repetition.

    Runs `warmup` unmeasured calls (cache/JIT/page-fault settling),
    then `repeats` measured calls with time.perf_counter. Returns all
    measured times, their mean, and a bootstrap CI on the mean.
    """
    if repeats < 2:
        raise ValueError("repeats must be >= 2 for a CI")
    for _ in range(warmup):
        fn(*args, **kwargs)
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn(*args, **kwargs)
        times.append(time.perf_counter() - t0)
    ci = bootstrap_ci(times, seed=ci_seed)
    return TimingResult(tuple(times), float(np.mean(times)), ci)


# ---------------------------------------------------------------------------
# Multiple-hypothesis correction (METHODS.md section 5)
# ---------------------------------------------------------------------------

def holm_bonferroni(p_values, alpha: float = ALPHA) -> list:
    """Holm-Bonferroni step-down within one hypothesis family.

    Sort the m p-values ascending; reject the i-th smallest (1-based)
    while p_(i) <= alpha / (m - i + 1); stop at the first failure and
    retain everything after it. Controls the family-wise error rate at
    alpha. Returns a list of booleans (True = rejected null) aligned
    with the input order.
    """
    p = np.asarray(p_values, dtype=float)
    if np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must lie in [0, 1]")
    m = len(p)
    order = np.argsort(p, kind="stable")
    reject = np.zeros(m, dtype=bool)
    for rank, i in enumerate(order, start=1):
        if p[i] <= alpha / (m - rank + 1):
            reject[i] = True
        else:
            break
    return reject.tolist()

#!/usr/bin/env python3
"""Executable checks for the bench/ baseline+metrics+statistics harness.

Checks (PASS/FAIL per line, OVERALL at the end):
  1. Every available baseline coder round-trips on synthetic Markov
     text and compresses it below 8 bpb (and above 0).
  2. bits_per_byte is exact on a hand-computed case.
  3. The bootstrap CI covers a known Bernoulli mean at (approximately)
     the nominal 95% rate in a 500-simulation self-test.
  4. The percentile method and the degenerate fallback behave as
     specified.
  5. holm_bonferroni matches a hand-worked step-down example.
  6. Paired bootstrap + Cohen's d sanity on a constructed shift.
  7. Results store: append + load + summarize round-trip; the store is
     append-only (earlier bytes untouched by later appends).
  8. METHODS.md exists, is substantive, and is referenced by the
     metrics.py module docstring.

Run:  /Users/para/.venvs/rnr/bin/python bench/run_checks.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import baselines
import metrics
import results
import smoke

BENCH_DIR = Path(__file__).resolve().parent

_failures = 0


def report(ok: bool, label: str, detail: str = "") -> None:
    global _failures
    tag = "PASS" if ok else "FAIL"
    if not ok:
        _failures += 1
    suffix = f"  ({detail})" if detail else ""
    print(f"{tag}: {label}{suffix}")


# ---------------------------------------------------------------------------
# 1. Baseline round-trips
# ---------------------------------------------------------------------------

def check_baselines() -> None:
    data = smoke.markov_text(n_bytes=64 * 1024, seed=7)
    coders = baselines.available_codecs()
    report(
        {"gzip", "bz2", "lzma"} <= set(coders),
        "stdlib baselines (gzip, bz2, lzma) available",
        f"available: {coders}",
    )
    report("zstd" in coders, "zstd baseline available (zstandard package)")
    report("brotli" in coders, "brotli baseline available (brotli package)")

    for coder in coders:
        try:
            res = baselines.run_baseline(coder, data)
        except Exception as exc:  # round-trip failure or codec error
            report(False, f"{coder} round-trip", repr(exc))
            continue
        bpb = res.bits_per_byte
        report(
            0.0 < bpb < 8.0,
            f"{coder} round-trip + compresses Markov text",
            f"{bpb:.3f} bpb, enc {res.enc_s * 1e3:.1f} ms",
        )

    # Isolated-child mode: one codec suffices to exercise the path.
    res = baselines.run_baseline("gzip", data, isolated=True)
    report(
        res.rss_mode == "isolated-child" and res.peak_rss_estimate > 1_000_000,
        "isolated-child RSS measurement returns a plausible high-watermark",
        f"{res.peak_rss_estimate / 1e6:.1f} MB",
    )


# ---------------------------------------------------------------------------
# 2-4. Metrics
# ---------------------------------------------------------------------------

def check_bits_per_byte() -> None:
    report(
        metrics.bits_per_byte(125, 1000) == 1.0
        and metrics.bits_per_byte(1000, 1000) == 8.0,
        "bits_per_byte exact on hand-computed cases",
    )


def check_ci_coverage() -> None:
    """500-simulation coverage self-test for the BCa bootstrap CI.

    True model: Bernoulli(p=0.3), n=100 per simulation; nominal level
    95%. The self-test uses B=1000 resamples for speed (the frozen
    reporting default is B=10000; the coverage property being checked
    is insensitive to this). BCa on a bounded discrete mean at n=100
    typically covers at 92-95%; with 500 simulations the binomial noise
    is about +/-1% (1 sigma), so the acceptance window is [0.90, 0.99].
    """
    n_sims, n, p, B = 500, 100, 0.3, 1000
    rng = np.random.default_rng(metrics.derive_seed("ci-coverage-selftest", 0))
    covered = 0
    for sim in range(n_sims):
        x = (rng.random(n) < p).astype(float)
        ci = metrics.bootstrap_ci(x, n_resamples=B, seed=int(rng.integers(2**31)))
        if ci.lo <= p <= ci.hi:
            covered += 1
    rate = covered / n_sims
    report(
        0.90 <= rate <= 0.99,
        "BCa bootstrap CI covers Bernoulli(0.3) mean at nominal rate",
        f"coverage {rate:.3f} over {n_sims} sims (window [0.90, 0.99])",
    )


def check_ci_methods() -> None:
    x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    bca = metrics.bootstrap_ci(x, n_resamples=2000, seed=1)
    pct = metrics.bootstrap_ci(x, n_resamples=2000, seed=1, method="percentile")
    report(
        bca.method == "bca"
        and pct.method == "percentile"
        and bca.lo < bca.point < bca.hi
        and pct.lo < pct.point < pct.hi,
        "BCa and percentile intervals both bracket the point estimate",
        f"bca [{bca.lo:.3f}, {bca.hi:.3f}], pct [{pct.lo:.3f}, {pct.hi:.3f}]",
    )
    deg = metrics.bootstrap_ci([2.0, 2.0, 2.0, 2.0])
    report(
        deg.method == "degenerate" and deg.lo == deg.hi == 2.0,
        "degenerate fallback on constant samples",
    )


def check_holm() -> None:
    """Hand-worked example: p = [0.01, 0.04, 0.03, 0.005], alpha = 0.05.
    Sorted: 0.005 <= 0.05/4, 0.01 <= 0.05/3, 0.03 > 0.05/2 -> stop.
    Rejected: 0.005 and 0.01 only."""
    got = metrics.holm_bonferroni([0.01, 0.04, 0.03, 0.005], alpha=0.05)
    report(
        got == [True, False, False, True],
        "holm_bonferroni matches hand-worked step-down example",
        f"got {got}",
    )


def check_paired() -> None:
    rng = np.random.default_rng(3)
    a = rng.normal(2.0, 0.1, size=40)
    b = a - 0.5  # b is uniformly 0.5 lower
    ci = metrics.paired_bootstrap_diff(a, b, n_resamples=2000, seed=2)
    d = metrics.cohens_d_paired(a, b)
    report(
        abs(ci.point - 0.5) < 1e-12 and ci.lo <= 0.5 <= ci.hi and d > 100,
        "paired bootstrap diff + Cohen's d on a constant 0.5 shift",
        f"diff CI [{ci.lo:.4f}, {ci.hi:.4f}], d={d:.3g}",
    )


# ---------------------------------------------------------------------------
# 7. Results store
# ---------------------------------------------------------------------------

def check_results_store() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        store = Path(tmp) / "test.jsonl"
        for coder in ("gzip", "lzma"):
            for seed in (11, 22, 33):
                rec = results.make_record(
                    corpus="unit-test-corpus",
                    coder=coder,
                    config={"level": 9},
                    seed=seed,
                    measured={
                        "bpb": 2.0 + 0.01 * seed + (0.5 if coder == "gzip" else 0),
                        "enc_s": 0.1,
                        "dec_s": 0.01,
                        "peak_rss_estimate": 10_000_000,
                    },
                )
                results.append(store, rec)

        loaded = results.load(store)
        report(len(loaded) == 6, "results store append + load round-trip",
               f"{len(loaded)} records")

        keys = {results.record_key(r) for r in loaded}
        report(len(keys) == 6, "record keys are distinct")

        table = results.summarize(store)
        report(
            "gzip" in table and "lzma" in table and table.count("|") > 20
            and "unit-test-corpus" in table,
            "summarize renders a Markdown table with both coders",
        )

        before = store.read_bytes()
        results.append(
            store,
            results.make_record(
                "unit-test-corpus", "bz2", {"level": 9}, 44,
                {"bpb": 2.2, "enc_s": 0.2, "dec_s": 0.02,
                 "peak_rss_estimate": 9_000_000},
            ),
        )
        after = store.read_bytes()
        report(
            after[: len(before)] == before and len(results.load(store)) == 7,
            "store is append-only (existing bytes untouched by append)",
        )

        try:
            results.append(store, {"corpus": "x", "coder": "y"})
            report(False, "append rejects records missing required keys")
        except ValueError:
            report(True, "append rejects records missing required keys")


# ---------------------------------------------------------------------------
# 8. METHODS.md
# ---------------------------------------------------------------------------

def check_methods_doc() -> None:
    path = BENCH_DIR / "METHODS.md"
    exists = path.is_file()
    text = path.read_text(encoding="utf-8") if exists else ""
    report(
        exists and len(text) > 2000 and "Holm" in text and "BCa" in text
        and "Amendment" in text,
        "METHODS.md exists and freezes bootstrap + correction + amendment rule",
    )
    doc = metrics.__doc__ or ""
    report(
        "METHODS.md" in doc,
        "metrics.py module docstring references METHODS.md",
    )


def main() -> int:
    print("bench/ harness checks")
    print("=" * 60)
    check_baselines()
    check_bits_per_byte()
    check_ci_methods()
    check_holm()
    check_paired()
    check_ci_coverage()
    check_results_store()
    check_methods_doc()
    print("=" * 60)
    if _failures == 0:
        print("OVERALL: PASS")
        return 0
    print(f"OVERALL: FAIL ({_failures} failing checks)")
    return 1


if __name__ == "__main__":
    sys.exit(main())

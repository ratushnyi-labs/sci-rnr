# Frozen statistical methodology for the RNR empirical program

Version 1.0 — frozen 2026-07-10, before the first recorded measurement.
Companion to the experimental-design document (`tex/rnr_experimental_design.tex`),
sections 5–6; reference implementation in `bench/metrics.py`.

Sections 1–7 below are FROZEN. They are never edited in place after the
first measurement is recorded; changes go through the Amendment rule
(section 7) and the Amendment log (section 8).

## 1. Scope

This document fixes, prior to any measurement, the statistical procedures
used to report every primary hypothesis test H1–H18 of the experimental
design: point estimates, uncertainty intervals, timing protocol, seed
derivation, and multiple-hypothesis correction. Exploratory or secondary
analyses may use other procedures but must be labeled "secondary,
uncorrected" wherever reported.

## 2. Uncertainty intervals: BCa bootstrap

- Primary interval: **BCa** (bias-corrected and accelerated) bootstrap
  (Efron–Tibshirani 1993, ch. 14), two-sided, nominal level **95%**
  (alpha = 0.05).
- Resample count: **B = 10000** (exp-design section 6.1).
- Fallback ladder, recorded in the reported interval's `method` field:
  1. **bca** — default.
  2. **percentile** — used automatically when BCa is undefined: the
     bias-correction z0 is infinite (all resamples on one side of the
     point estimate) or the jackknife variance is zero.
  3. **degenerate** — all samples identical; zero-width interval.
- Statistic resampled: the mean, unless a hypothesis explicitly names a
  quantile (H10/H13 report 95th-percentile latency; the same BCa scheme
  applies with the quantile as the plug-in statistic).
- Paired comparisons (RNR variant vs baseline on the same blocks,
  exp-design section 6.2): bootstrap on the **paired per-block
  differences** (`metrics.paired_bootstrap_diff`), never on the two
  margins separately. Cohen's d for paired block-level comparisons is
  mean(diff)/sd(diff) (`metrics.cohens_d_paired`).

## 3. Timing protocol

- Each (coder, corpus) measurement is repeated **5 times** with distinct
  deterministic run seeds (exp-design section 6.1), after **1 unmeasured
  warmup run** (cache/page-fault settling); the warmup is not one of
  the 5.
- Clock: `time.perf_counter` wall-clock, single-threaded unless the
  coder is inherently multi-threaded (then thread count is recorded in
  the config field).
- Report mean and BCa 95% CI over the 5 measured runs.
- Peak memory: resident-set high-watermark (`ru_maxrss`). Reported
  numbers use the **isolated-child** mode of `bench/baselines.py`
  (fresh process per measurement); the cheaper in-process delta mode is
  a lower bound reserved for smoke tests and is labeled as such in the
  results store (`rss_mode` field).

## 4. Seed policy

- Master seed: **20260710**.
- Every derived seed is the first 4 bytes (big-endian) of
  `SHA-256("{master}:{tag}:{experiment_id}:{index}")`
  (`metrics.derive_seed`), with tags:
  - `run` — per-run seeds, index 0..4;
  - `boot` — bootstrap RNG seed for each reported interval;
  - `corpus` — deterministic corpus subsetting (e.g. CommonCrawl
    sampling), per exp-design section 6.1.
- All seeds actually used are recorded in the results store alongside
  the measurement. No unlogged randomness in any reported number.

## 5. Multiple-hypothesis correction: Holm–Bonferroni within families

Primary hypothesis tests are corrected by **Holm–Bonferroni step-down
within each hypothesis family** at family-wise alpha = 0.05
(`metrics.holm_bonferroni`): sort the m family p-values ascending,
reject the i-th smallest while p_(i) <= alpha/(m − i + 1), stop at the
first failure.

Families (fixed here, before any measurement):

| Family | Hypotheses | Character |
|---|---|---|
| RATE   | H1, H2, H5, H7, H15 | coding-rate comparisons vs baselines |
| INFO   | H3, H4, H6          | information-theoretic condition measurements |
| ACCESS | H10, H11, H12, H13, H17 | random-access / latency / cache scaling |
| REPRO  | H8, H9, H14         | byte-identity conformance (binary) |
| SAFETY | H16, H18            | hard worst-case guarantees (binary) |

Rules:

- Each hypothesis contributes **one pre-specified primary test** to its
  family (one designated coder/corpus pair, named in the
  pre-registration record before the run). Additional coder/corpus
  pairs under the same hypothesis are secondary and reported with
  uncorrected CIs, clearly labeled (exp-design section 6.4).
- REPRO and SAFETY families are **pass/fail conformance tests**, not
  p-value tests (a single byte divergence or bound violation fails
  them outright); they are excluded from the p-value correction and
  reported as exact counts.
- Note: exp-design section 6.4 (written for the original 14 hypotheses)
  prescribes Holm–Bonferroni across all primary tests jointly. This
  document refines that to within-family correction for the extended
  H1–H18 list; the refinement predates the first measurement and is
  therefore part of the frozen protocol, not an amendment.

## 6. Results storage

- Append-only JSON-lines store (`bench/results.py`), one measurement
  per line, keyed by (corpus, coder, config, seed, timestamp).
- Records are never edited or deleted. A correction is a new record
  whose `supersedes` field names the replaced record's key; summaries
  must drop superseded records.
- Corpora are referenced by SHA-256 hash in a `corpus_sha256` field
  (exp-design section 3).

## 7. Amendment rule

After the first measurement is recorded (the first line written to any
`bench/results/*.jsonl` intended for reporting), ANY change to:

- sections 1–7 of this document, or
- a statistical code path in `bench/metrics.py` or the aggregation
  logic in `bench/results.py`

requires, **before the changed procedure is used**, an entry in the
Amendment log below stating: date, what changed, why, and which
already-recorded results are affected (and whether they are re-analyzed
or left as-is with the old procedure). Frozen sections are never edited
in place; amendments only add to the log. A change without a log entry
invalidates every result produced after it.

## 8. Amendment log

(empty — no amendments)

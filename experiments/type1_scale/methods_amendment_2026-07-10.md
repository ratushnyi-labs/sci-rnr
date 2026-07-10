# Dated amendment note — Type-I scale campaign timing protocol

Date: 2026-07-10. Scope: `experiments/type1_scale/` only. This note
implements the METHODS.md §7 amendment rule for the scale campaign
without editing the frozen METHODS.md (which stays untouched); it is
registered before the first paper-tier measurement.

## What deviates

METHODS.md §3 (frozen) prescribes, per (coder, corpus) measurement,
**5 timed runs + 1 unmeasured warmup**, with BCa 95% CIs over the 5
runs, and in-process `time.perf_counter` timing with isolated-child RSS.

The scale campaign (protocol.md §2, §5) uses instead:

1. **3 timed runs at the 500 MB tier; 1 timed run at 1500 MB and
   5000 MB; no warmup run.**
2. Wall time measured by the parent around a dedicated child process
   per phase (encode / decode), not around an in-process library call.
3. Peak RSS from `/usr/bin/time -l` on the same child
   (`rss_mode: time-l-child`), not `bench/baselines.py` isolated mode.

## Why

* At 500–5000 MB, cell wall times range from minutes to the 90-minute
  cap. Five runs plus warmup would multiply a full tier to multiple
  days on this host (8-core arm64; measured smoke throughputs down to
  ~0.25 MB/s for zstd-19/brotli-q11-class cells) for negligible
  variance information: at >= 100 s per run, run-to-run wall-time
  variation on an otherwise idle machine is far below the effect sizes
  of interest (coders differ by integer factors), and rate/RSS metrics
  are exactly reproducible for these deterministic coders.
* A warmup run is also *misleading* at this scale: a 5 GB corpus does
  not fit the page cache reliably, so the "warm" state the warmup rule
  was written for (small corpora, METHODS §3) does not exist; the cold
  first pass is the realistic regime.
* External CLI archivers cannot be timed in-process; the child-process
  clock plus `/usr/bin/time -l` high-watermark is the uniform
  instrument across the whole matrix (including the RNR CLIs), which
  is more comparable than mixing library and subprocess timings.

## Consequences for reporting

* BCa CIs (METHODS §2) are reported only where >= 3 runs exist (the
  500 MB tier); 1500/5000 MB timings are single-run point estimates
  and must be labeled as such wherever shown.
* No already-recorded result is affected: this note predates the first
  paper-tier record. Smoke-tier records (machinery validation) also use
  the 1-run policy and are never reported as benchmark data.

## What does NOT change

Seed derivation (METHODS §4), the append-only results store and
supersede rule (§6), the BCa procedure itself (§2), Holm–Bonferroni
family structure (§5), and corpus-by-hash referencing all apply
unchanged. This campaign feeds H1/H2-family secondary scale evidence;
any primary hypothesis test keeps the frozen METHODS procedures.

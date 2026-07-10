# Type-I scale benchmark — pre-registered protocol

Version 1.0 — registered 2026-07-10, before any tier measurement.
Companion to `bench/METHODS.md` (frozen statistical methodology) and
`tex/rnr_experimental_design.tex` §5–6. Deviations from METHODS.md §3
are logged in `methods_amendment_2026-07-10.md` (dated note in this
directory; METHODS.md itself is frozen and is not edited).

Smoke-tier results (16/48 MiB) validate the machinery only and are
never reported as benchmark data. The pre-registered claims concern the
500 / 1500 / 5000 MB tiers exclusively.

## 1. Design

Full-factorial campaign over

* **size** — exactly 500 MB, 1500 MB, 5000 MB (10^6-byte MB:
  500,000,000 / 1,500,000,000 / 5,000,000,000 bytes);
* **content type** — white_noise, text, video, images, archive_mix,
  precompressed (recipes in §3);
* **coder** — the matrix in §4;

executed by `run_campaign.py` (per-cell checkpointing, §5), summarized
by `make_figures.py` (§6), gated by `run_checks.py` (§7).

A **cell** is one (size, type, coder, run-index) tuple: one encode,
one size/hash measurement, one decode-verify, on a single dedicated
child process per phase.

## 2. Tiers

| tier | bytes | role | runs/cell |
|---|---|---|---|
| smoke | 16 MiB | machinery validation only | 1 |
| smoke48 | 48 MiB | machinery validation (second size point for figure paths) | 1 |
| 500 | 500,000,000 | paper tier | 3 |
| 1500 | 1,500,000,000 | paper tier | 1 |
| 5000 | 5,000,000,000 | paper tier | 1 |

The repetition policy (3 runs at 500 MB; 1 run at 1500/5000 MB)
deviates from the frozen 5-runs-plus-warmup rule of METHODS.md §3 and
is documented, with rationale, in the dated amendment note.

## 3. Content-type recipes (corpora.py)

Every corpus is materialized deterministically to the exact tier byte
count and pinned by sha256 in `corpora/<tier>/manifest.json` before any
cell runs; each results record carries the corpus sha256. Seeds follow
METHODS.md §4 (`derive_seed`, tag `corpus`, experiment id
`type1_scale:<tier>:<type>`).

### 3.1 white_noise
numpy Philox counter-based PRG (`Generator.bytes`), keyed by the derived
seed. Reproducible at any size; numpy version recorded. Negative
control: the known optimum is 8 bpb and every honest coder must stay
within its documented framing overhead of it.

### 3.2 text
Byte slice of the standard Hutter corpus files, cycled to the target
size: smoke = enwik8 prefix; paper tiers = enwik9 (deferred `data/`
entry; must be fetched before the 500 MB tier runs). Cycling period is
1 GB at paper tiers, which exceeds the largest window in the coder
matrix (zstd --ultra -22: 128 MiB), so no coder can span a full period;
the 1500/5000 MB text cells therefore measure repeated-content corpora
whose repetition is *invisible to every coder in the matrix* — recorded
as a design note, not a confound. Using enwik8 (100 MB period) at paper
tiers is forbidden: the zstd-22u window would cover a period.

### 3.3 video
Synthetic raw-video model v1: 640x360 8-bit luma frames; each frame is
an AR(1) temporal mixture (rho = 0.94) over spatially box-smoothed
(window 15) Philox innovations, quantized to uint8. Approximates the
strong spatial+temporal correlation of raw video while being seeded,
self-contained, and generable at any size. **Caveat (honest):** this is
a synthetic stand-in, not camera footage; real Y4M material can be
added later only as a *new, separately named* content type — recipe v1
stays fixed.

### 3.4 images
Smoke: the pinned `telemetry_grid.f32` payload (16 MiB, exactly the
smoke size). Paper tiers: a sequence of independent seeded correlated
2048x2048 float32 rasters (box window 49 + smooth deterministic trend
per raster; 16 MiB each), the same raster/numeric class as the data/
telemetry generator. Same synthetic-stand-in caveat as §3.3.

### 3.5 archive_mix
Cyclic concatenation calgary.tar | canterbury.tar | rnr_scripts_src.tar
| sqlite_synth.db (cycle length 21,798,912 bytes), truncated to the
tier size. **Design note:** at paper tiers the cycle repeats and every
coder window (>= 4 MiB) can span cycle boundaries, so this type
deliberately measures *long-range redundancy exploitation* — the
realistic backup-archive workload — not one-shot heterogeneous coding.
Figures and tables must carry this note.

### 3.6 precompressed
`zstd -3 -T1` (the zstd CLI default level — the most common
precompressed payload in practice) applied to the §3.2 text stream,
truncated to the tier size. Inner zstd CLI version recorded in the
corpus manifest. Negative control: near-incompressible input with
non-trivial internal structure.

## 4. Coder matrix

### 4.1 RNR lane

| coder | config | tiers | rationale |
|---|---|---|---|
| rnr1-fast | W=3, K=64 KiB | all | scale implementation from `impl_fast/` |
| rnr1-cm-fast | impl_fast CM variant | all | only if `impl_fast/` delivers it |
| rnr1-ref | W=3, K=64 KiB (`impl/rnr1.py`) | smoke only | normative reference; measured ~0.05–0.07 MB/s on this host, so a 500 MB cell would need ~2.5 h > the 90 min cap — pre-registered as infeasible at paper tiers |

`impl_fast/` integration contract: `impl_fast/BENCH_ADAPTER.json` maps
coder names to `{"pack": argv, "unpack": argv, "config": {...}}`
templates with `{in}`/`{out}` placeholders (relative argv[0] resolved
inside `impl_fast/`); fallback probe is `impl_fast/rnr1_fast.py` with an
`impl/rnr1.py`-compatible `pack`/`unpack` CLI. If neither exists the
RNR-fast lane is simply not scheduled (never silently substituted).

### 4.2 Baseline CLI archivers

Single-threaded (`-T1` where the tool has threads), pinned binaries,
versions recorded per record.

| coder | invocation | tiers |
|---|---|---|
| zstd-19 | `zstd -19 -T1` | all |
| zstd-22u | `zstd --ultra -22 -T1` (decode `--memory=2048MB`) | all |
| xz-6 | `xz -6 -T1` | all |
| xz-9e | `xz -9 -e -T1` | all |
| brotli-q11 | `brotli -q 11 --lgwin=24` | smoke, 500, 1500 |
| brotli-q9 | `brotli -q 9 --lgwin=24` | smoke, 5000 |
| gzip-9 | `gzip -9` | all |
| bzip2-9 | `bzip2 -9` | all |

brotli -q11 at 5000 MB is pre-registered as infeasible: measured q11
throughput on this host is ~0.25–1 MB/s (16 MiB noise took ~68 s under
load), which extrapolates to multiple hours per 5000 MB cell — beyond
the 90 min cap; q9 is the pre-registered 5000 MB stand-in and is also
run at smoke so the substitution is characterized. Several other cells
are *expected* censoring candidates at 5000 MB on this host (smoke
throughput bands, sensitive to background load by up to ~4x: zstd-19
~0.25–1 MB/s on text, xz-9e ~0.6–2 MB/s): cap hits are recorded as
censored cells and reported as such, not dropped. Paper-tier runs MUST
be executed on an otherwise idle machine (see §8).

**Exclusions (pre-registered):** paq8/zpaq/cmix-class context-mixing
archivers are excluded from the paper tiers because their throughput
(typically 10–500 KB/s) implies multiple days per 500 MB cell —
irreconcilable with the 90 min cap; they are rate-champions but not
scale-feasible. The RNR CM lane (rnr1-cm-fast) covers the
context-mixing comparison if impl_fast delivers it.

## 5. Execution rules

### 5.1 Measurement
Per cell (coders.py): encode child writes the compressed stream to a
scratch tempfile (size then measured and recorded); decode child streams
its output directly into sha256 (no decoded file on disk; RNR CLIs
decode to a tempfile that is hashed then deleted — recorded per record
as `decode_io`). Round-trip passes iff the decoded sha256 and length
equal the corpus manifest values. Wall time = parent `perf_counter`
around each child; peak RSS = `/usr/bin/time -l` child high-watermark
(bytes; `rss_mode = time-l-child`). Scratch files are deleted on every
code path, so disk stays bounded by corpus + one compressed file.

### 5.2 Repetition and seeds
Runs per tier per §2. Every record carries a METHODS §4 derived seed
(tag `run`, experiment id `type1_scale:<tier>:<type>:<coder>`, index =
run index) even though all coders in the matrix are deterministic.

### 5.3 Time cap and censoring
Per-cell cap: 90 minutes wall (encode + decode combined). On expiry the
child process group is SIGKILLed and the cell is recorded with
`censored: true` and the phase reached; censored cells are terminal in
the journal (not retried) and appear in tables as censored, never as
missing. The cap applies per run-cell, so a 3-run 500 MB cell may
censor some runs and complete others.

### 5.4 Checkpointing and resume
`out/<tier>/cells.journal` records one JSON line per terminal cell
(done / censored / failed) *after* its results record is appended to
`out/<tier>/results.jsonl` (bench/results.py append-only schema). A
re-run skips journaled cells. A crash between record and journal write
re-runs at most one cell; summaries deduplicate on (corpus, coder,
config, seed). Do not run two campaign processes against the same
out-dir concurrently.

### 5.5 Integrity
Corpora are re-hashed against the manifest pin at every campaign start;
drift aborts the campaign. A round-trip failure marks the cell
`failed` (terminal) and the campaign exits nonzero at the end.

## 6. Outputs and figures

`make_figures.py` regenerates deterministically from the JSONL stores:
rate-vs-size per type; encode and decode time vs size (log-log);
rate-vs-encode-throughput Pareto per (size, type); one summary CSV per
(size, type); `captions.md`. PDF+PNG, colorblind-safe validated
palette, no in-figure titles, censored cells listed in the CSVs.

## 7. Gate

`run_checks.py` must print `OVERALL: PASS` before any paper-tier run is
launched: coder availability, corpus determinism, the full smoke matrix
(every coder on >= 2 types at 16 MiB) with all metrics present and
round-trips verified, kill-and-resume checkpointing, figure/CSV
generation, and figure determinism (byte-identical on regeneration).

## 8. Host

Apple Silicon (arm64) macOS, 8 cores, 24 GB RAM; binaries:
/opt/homebrew/bin/{zstd 1.5.7, xz 5.8.3, brotli 1.2.0}, Apple gzip 479,
bzip2 1.0.8, /usr/bin/time. Python: /Users/para/.venvs/rnr/bin/python
(3.14.x, numpy 2.4.x). All throughput expectations quoted in this
protocol were measured on this host at smoke scale, some under
non-trivial background load (observed wall-time inflation up to ~4x);
they are planning bands, not results. Paper tiers require an otherwise
idle machine, and the 500 MB tier's 3 runs provide the spread check.

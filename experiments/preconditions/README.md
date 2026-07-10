# Precondition statistics — which corpora qualify (BUG-001-D / BUG-003-D / BUG-010-D, dev-slice scale)

Measures the qualifying data-class precondition of the main paper
(`tex/papers/core/rnr_core.tex` §10.2.1, inequality (10.1)) on every
non-deferred corpus of `data/MANIFEST.json`, using the reference Type-I
coder (`impl/rnr1.py`) as the side-information instrument:

* **(a) conditional-advantage ratio** (BUG-001-D): Definition-2.1 ratio
  `E[L(Enc(X))]/E[L(B(X))]` at `W ∈ {2,3,4}`, `K = 64 KiB`, paired
  per-block vs zstd/lzma at matched sync granularity (seekable-format
  proxy) and vs streaming baselines, with the full (10.1) term
  breakdown (a)–(e); plus the sync-point overhead datum at typical
  `(K, W)` (BUG-003-D / exp-design H10(c));
* **(b) stationarity / mixing indicators**: per-third entropy-ladder
  spread; byte-stream autocorrelation at structural lags
  (77 = base64 line, 4096 = SQLite page, 8192 = telemetry row);
* **(c) entropy-rate ladder**: held-out (train-half/eval-half) and
  plug-in order-0..4 conditional empirical entropy, and the
  Type-III-B `I(H;L|C)` instrument for H6 scoping.

Statistics follow the frozen methodology (`bench/METHODS.md` v1.0):
BCa bootstrap, B = 10000, derived seeds, Holm–Bonferroni within the
screening family — reported as a labeled SECONDARY analysis (the
screen is not one of the frozen H1–H18 primaries; no amendment).

## Running

```sh
# gate, reduced scale (< 10 min): unit checks + reduced campaign
/Users/para/.venvs/rnr/bin/python -u run_checks.py

# gate + the real dev-slice campaign (writes findings.md,
# qualification.csv, out/manifest_qualifies_patch.json):
/Users/para/.venvs/rnr/bin/python -u run_checks.py --full

# campaign only:
/Users/para/.venvs/rnr/bin/python -u measure_preconditions.py [--full]
```

Every gate check prints `[PASS]/[FAIL]` and the run ends with
`OVERALL: PASS/FAIL` (exit code 0/1).

## Deliverables

* `findings.md` — numbers, CIs, qualification table, per-hypothesis
  verdicts, honest caveats (written by the `--full` campaign).
* `qualification.csv` — the corpus × precondition table.
* `out/results.jsonl` — append-only bench-compatible records
  (`bench/results.py` schema, `corpus_sha256` per METHODS §6).
* `out/manifest_qualifies_patch.json` — proposed `qualifies_for`
  updates for `data/MANIFEST.json` (data/ is **not** edited directly).
* `out/*_{scale}.csv`, `out/summary_{scale}.json` — per-scale detail.

## Scale honesty

Coding runs on the first 1 MiB (full) / 256 KiB (reduced) per corpus,
indicators on the first 8 MiB / 2 MiB; the instrument is the order-W
counting predictor, not the paper's neural predictor. This is the
LOCAL dev-slice version of the measurements BUG-001-D / BUG-010-D call
external; those items remain open at full scale. enwik9 is deferred
and not measured.

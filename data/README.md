# RNR corpus suite (small-classics scale)

Corpus infrastructure for the RNR empirical program, implementing
`tex/rnr_experimental_design.tex` §3 (corpus selection) at "small classics"
scale. Payloads live in `data/payloads/` and are **not** tracked by git
(see `.gitignore`); only the code, this README, and `MANIFEST.json` are.
Every payload is pinned by sha256 in the manifest, per §3's rule that
"experimental scripts reference corpora by hash, not filename".

## Contents

| name           | §3 class                   | bytes       | provenance | qualifies for |
|----------------|----------------------------|------------:|------------|---------------|
| enwik8         | 3.1 text/source            | 100,000,000 | download (mattmahoney.net) | H1 H3 H4 H6 H8 |
| enwik9         | 3.1 text/source            | 1,000,000,000 (deferred) | download (mattmahoney.net) | H1 H3 H4 H6 H15 |
| canterbury     | 3.1 text/source            | 2,821,120   | download (corpus.canterbury.ac.nz) | H1 H3 H4 H8 |
| scripts_src    | 3.1 text/source            | 3,358,720   | derived: pinned tar of repo `scripts/` | H1 H3 H4 |
| sqlite_synth   | 3.2 structured binary      | 12,353,536  | generated: seeded synthetic SQLite DB | H2 H3 H4 H8 |
| telemetry_grid | 3.3 image/numeric          | 16,777,216  | generated: seeded correlated float32 field | H4 H6 |
| calgary        | 3.4 heterogeneous/negative | 3,265,536   | download (corpus.canterbury.ac.nz) | H1 H5 H8 |
| ctrl_urandom   | 3.4 negative control       | 16,777,216  | generated: `os.urandom` | H16 |
| ctrl_zstd      | 3.4 negative control       | 9,111,253   | derived: enwik8[:32 MiB] at zstd -19 | H16 |
| ctrl_base64    | 3.4 negative control       | 16,997,969  | generated: base64 of seeded PRNG bytes | H1 (calibration) |

Byte sizes above are those of the current pinned instances; `MANIFEST.json`
is authoritative (it also carries the full sha256, source URL, license note,
and provenance notes for each entry).

Class mapping notes:

* **enwik8** is the small-scale stand-in for §3.1's enwik9; **enwik9**
  itself is a manifest-only *deferred* entry (1 GB) fetched on demand.
* **scripts_src** (this repository's own `scripts/` tree, ~245 Python
  files, deterministic tar) stands in for §3.1's Linux-kernel-source
  tarball as the source-code sample.
* **sqlite_synth** stands in for §3.2's database pages: B-tree pages of a
  real SQLite file with a seeded sensor/readings/logs schema and an index.
* **telemetry_grid** stands in for §3.3's HDF5 scientific arrays: a raw
  headerless 2048x2048 little-endian float32 grid, Gaussian-filtered white
  noise (correlation length 24 px) plus a smooth trend, fully seeded.
* **calgary** is classed 3.4 because its 14+ file mix (text, sources,
  geophysical binary, bitmap, object files) is a heterogeneous-scheduling
  corpus in the sense of §3.4's ISO entry; it also serves 3.1 measurements.
* The three `ctrl_*` entries are §3.4's negative controls at reduced scale
  (16 MiB urandom vs the spec's 64 MB; 32 MiB-input zstd blob vs the
  spec's 128 MB). AES-CTR ciphertext from §3.4 is not materialized;
  `ctrl_urandom` covers the "appears random" control at this scale.

## Usage

```sh
# fetch/build everything non-deferred (idempotent; sha256-checked)
/Users/para/.venvs/rnr/bin/python data/fetch.py

# also fetch enwik9 (~0.3 GB download, 1 GB payload)
/Users/para/.venvs/rnr/bin/python data/fetch.py --include-deferred

# verify the suite
/Users/para/.venvs/rnr/bin/python data/run_checks.py
```

Experiment code reads corpora only through the loader:

```python
from loader import iter_blocks, corpus_names, corpus_bytes, read_range

for block in iter_blocks("enwik8", K=4096):      # streaming, O(K) memory
    ...
read_range("sqlite_synth", offset=8192, length=4096)  # random access
```

## Reproducibility and pinning

* Downloads are streamed, then hashed; `fetch.py` re-running verifies the
  existing payload against the manifest and skips it (idempotent).
* Generated entries are seeded (suite seed `20260710`; per-entry derived
  seeds recorded in `params`). Two caveats are inherent and documented in
  the manifest: `ctrl_urandom` is non-reproducible by design (the hash pins
  the drawn instance), and `sqlite_synth`'s byte layout may differ across
  SQLite library versions (generation version recorded in its entry).
* `scripts_src` is pinned at generation time; if `scripts/` changes and the
  payload is regenerated from scratch, the hash changes and the manifest is
  updated — the manifest sha256 is always the reference.
* If a download fails, `fetch.py` synthesizes a deterministic text stand-in
  and marks the entry `"standin": true`; `run_checks.py` reports stand-ins
  in its summary and skips content checks that require the real payload.

## Checks

`run_checks.py` prints PASS/FAIL per check and an OVERALL verdict
(exit code 0 iff all pass):

* C1 manifest schema (fields required by §3 and §12 open-artifacts list),
* C2 payload integrity (existence, byte size, sha256) for every
  non-deferred entry,
* C3 loader round-trip (block reassembly at two K values, block-count
  formula, `drop_last` accounting),
* C4 class coverage (at least one materialized corpus in each of
  §3.1-3.4),
* C5 deferred handling (enwik9 manifest-only; loader refuses cleanly),
* C6 content sanity (tar readability, SQLite magic and row counts,
  telemetry spatial autocorrelation > 0.9, negative controls incompressible
  under zlib, zstd blob round-trips to its exact source prefix, base64
  decodes to its seeded source, text-vs-random compressibility gap).

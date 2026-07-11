# rnr-cli — production CLI for the RNR1 Type-I coder

Statically compilable Rust port of the normative reference coder
`impl/rnr1.py` (archive format: `impl/FORMAT.md`).  The Python reference
stays normative; this binary is accepted only through **bit-identity**:
its archives must equal `rnr1.pack()` byte for byte and it must decode
exactly the archives the reference accepts.  The validated C port
(`impl_fast/`) serves as a second independent cross-decode engine.

## Layout

```
Cargo.toml            crate manifest (only dependency: sha2)
src/main.rs           CLI (pack / unpack / read / info)
src/coder.rs          E12-A map, arithmetic coder, order-W predictor,
                      sub-block encode/decode (integer-only path)
src/container.rs      archive container, seek index, store mode,
                      random-access read, multi-threaded pack/unpack
run_checks.py         NON-NEGOTIABLE acceptance gate (see below)
build.sh              cross-target release builds + dist/MANIFEST.json
BENCH_ADAPTER.json    coder lanes for the scale campaign (rnr-rust,
                      rnr-rust-mt8), same schema as impl_fast's
```

## CLI

```
rnr pack   <input> <output.rnr> [--W 3] [--K 65536] [--threads 1]
rnr unpack <archive.rnr> <output> [--no-verify] [--threads 1]
rnr read   <archive.rnr> --pos P --len k [--out FILE] [--verify]
rnr info   <archive.rnr>
```

`--threads T` parallelizes over sub-blocks (independent restart points,
Definition 10.3); output is byte-identical to `--threads 1`.  `read`
decodes only the sub-block(s) covering `[P, P+k)` (decode window
`<= K + k`, Theorem 10.3).

## Acceptance gate

```
/Users/para/.venvs/rnr/bin/python rnr-cli/run_checks.py
```

R0 build; R1 smoke-class bit-identity + four-way cross-decode
(rust/reference/C engine) at W in {2,3}, K in {64,16} KiB; R2 corpus
matrix (every non-deferred `data/` manifest entry, W in {2,3},
K = 64 KiB: rust archive == reference archive byte for byte, MT == ST,
rust fully decodes reference archives, C engine decodes rust archives);
R3 error paths (tamper/truncation/magic/model-hash rejected); R4
two-process determinism; R5 random-access read conformance incl.
warmup/decoded statistics; R6 enwik8 throughput report.

Reference corpus archives come from the `impl_fast` cache
(`impl_fast/workdir/refarc`, produced by the F1-gated parallel
reference driver); missing entries are regenerated through the same
driver.  On the corpus matrix, "reference decodes the rust archive" is
implied by byte-identity (the rust archive *is* the reference archive,
which the reference already decoded in the impl_fast run) and is
additionally exercised directly on the R1 smoke matrix.

## Cross-target builds

```
sh rnr-cli/build.sh        # RNR_SKIP_DOCKER=1 to skip the Docker smoke
```

Builds `dist/rnr-<target>` for aarch64/x86_64 apple-darwin and
aarch64/x86_64 linux-musl (fully static: `file` reports
static-pie/statically linked, no interpreter), runs a cross-OS
determinism smoke (linux binaries under Docker must produce archives
byte-identical to the macOS binary on the same 4 MiB input), and writes
`dist/MANIFEST.json` (per-target sha256, linkage, smoke status).

Windows MSVC targets (`x86_64-pc-windows-msvc`,
`aarch64-pc-windows-msvc`): the sources compile (`cargo check` passes;
no platform-specific code), but linking on this host fails for lack of
the MSVC import libraries (`kernel32.lib`, `ntdll.lib`, `userenv.lib`,
`ws2_32.lib`, `dbghelp.lib` — Windows SDK components).  `cargo-xwin`
was not installed and has no prebuilt binary available here.  On a host
with the SDK (or xwin-fetched libs), `cargo build --release --target
x86_64-pc-windows-msvc` needs no source changes.  The recorded status
is in `dist/MANIFEST.json`.

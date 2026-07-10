# impl_fast — engineering notes, measurements, and amendments

Date: 2026-07-10.  Host: macOS arm64, 8 cores, 24 GB RAM, Apple clang,
`/Users/para/.venvs/rnr/bin/python`.

## What this is

`rnr1_fast.c` is a C performance port of the **normative** Python
reference Type-I coder `impl/rnr1.py` (format per `impl/FORMAT.md`),
with a thin ctypes wrapper `rnr1fast.py`.  The reference stays
normative; the port is accepted only through **byte-identity**:

- fast `pack()` output == reference `pack()` output, byte for byte,
  on every non-deferred `data/MANIFEST.json` entry at W in {2,3},
  K = 64 KiB, plus the smoke classes of `impl/run_checks.py` C2, at
  K in {64 KiB, 16 KiB} (checks F2/F3 in `run_checks.py`);
- cross-decode in both directions with full verification (sub-block
  hashes + archive H_v): fast decodes reference archives, the
  reference decodes fast archives (check F4);
- the multi-threaded paths must reproduce the single-threaded archive
  bytes exactly (F3b/F4b, F6).

Because the corpus-scale reference encodes/decodes would take ~90 CPU
minutes per direction in pure Python, the checks drive the reference
through `ref_driver.py`: a process pool over **independent sub-blocks**
calling `rnr1.encode_subblock` / `rnr1.decode_subblock`, with the
container assembled by the reference's own code paths.  Sub-blocks are
self-contained restart points (fresh predictor + coder state), so the
driver's output equals sequential `rnr1.pack()` output by construction;
this equality is itself verified at runtime on multi-sub-block inputs
before any parallel result is trusted (check F1).

## Bit-identity-preserving performance structure

The naive C port ran at 1.15 MB/s because the reference recomputes a
256-entry stable argsort and full prefix sums per coded byte.  The port
keeps the identical integer semantics but restructures the state:

1. Per-context **top-16 cache** (count desc, byte asc — the exact
   stable-argsort prefix) maintained incrementally under +1 updates.
   Exactness argument: one key changes per update, so re-inserting that
   key into the sorted prefix reproduces the full sort's prefix; an
   outside byte enters iff it now precedes the 16th element in the
   total order.  Cap-halving can reorder equal-count ties, so it
   triggers a full stable rebuild (rare: a context must accumulate
   2^16 observations inside one sub-block).
2. Per-context **16 bucket sums** so ESC-class prefix sums are O(32)
   instead of O(256).
3. **uint16 counts** (any count <= context total < 2^16 at rest; the
   single transient 2^16 case wraps and is repaired inside the halving
   step that fires in the same update).
4. One **exact-key table pass** per position shared by dist() and
   update(): entries created before coding have total == 0 and are
   invisible to the backoff rule (total > 0), exactly like absent dict
   entries in the reference.  Tables are open-addressing with **full
   64-bit key compare** (the reference dict is exact; lossy hashing is
   not allowed), generation-tagged for O(1) per-sub-block reset,
   load-factor <= 1/2 with deterministic growth.
5. Arithmetic-coder division skips for `cum_lo == 0` / `cum_hi ==
   total` (algebraically identical results).
6. Optional **multi-threaded** pack/unpack over independent sub-blocks
   (pthread pool, atomic work index, ordered assembly) — byte-identical
   output by the same argument as (F1), and verified.

## Measured performance (enwik8, W=3, K=64 KiB)

Single process, this host, best of 2 runs (run_checks F7 gate line,
2026-07-10, idle machine):

| configuration | encode | decode |
|---------------|--------|--------|
| 1 thread      | 9.07 MB/s | 9.90 MB/s |
| 8 threads     | 45.44 MB/s | 50.47 MB/s |

The >= 10 MB/s target is met by the deliverable in its 8-thread
configuration (45 MB/s); the single-thread rate is reported honestly
at 9.1 MB/s, just below target (the remaining hot spots are the
arithmetic coder's per-symbol 64-bit divisions and renormalization
loop and the per-position count-table updates, which are
semantics-bound).  Compression ratio is identical to the reference by
construction (3.6319 bpb on enwik8 at W=3, K=64 KiB; the reference
produces the same bytes).

## Memory

- Encoder/decoder context state: 640 B per live context entry, pool
  grown geometrically, reset (not freed) per sub-block; slot tables
  grow to <= 2x live entries per order.  Measured on enwik8-class text
  at W=3, K=64 KiB: ~15 MB coder state per instance (34 MB total RSS
  for a 10 MB standalone encode including input/output buffers); the
  hard bound per instance is (W+1) x min(K, 2^16-ish contexts) x 640 B
  approx 84 MB at W=3, K=64 KiB on adversarial data.  Multi-threaded
  paths use one instance per thread (~120 MB text / ~0.7 GB adversarial
  at 8 threads, well within 24 GB).
- `pack()`/`unpack()` additionally hold the input and output buffers in
  memory (the wrapper reads whole files); at the 5000 MB benchmark
  scale budget ~2x file size + coder state.  Peak RSS on the enwik8
  gate run is recorded in the F7 line of the run_checks log.

## Amendment note (timing methodology; METHODS.md is frozen)

`bench/METHODS.md` prescribes 5 measured runs (+1 warmup) with BCa CIs
for benchmark timing.  The F7 throughput gate in this directory is an
*engineering acceptance check*, not a paper benchmark number: it uses
best-of-2 runs per configuration at the full 100 MB enwik8 to keep the
gate's wall time inside the task budget.  Paper-grade Type-I scale
timings at 500/1500/5000 MB should follow METHODS.md (or its dated
amendment for big sizes) and use `rnr1fast.py` CLI in isolated-child
mode; thread count must be recorded per METHODS.md's multi-thread
provision (line "single-threaded unless the coder is inherently
multi-threaded (then thread count is recorded)").

## Caveats

- macOS-only build as delivered (SHA-256 via CommonCrypto, which lives
  in libSystem; no extra link flags).  Portability would need a vendored
  SHA-256 (~60 lines) — isolated in `sha256_trunc8`/`model_hash`.
- `rnr1_unpack` matches the reference's accept/reject behavior on all
  well-formed archives and on the tamper classes exercised by F5.  Two
  deliberate divergences on *malformed* input: (a) archives whose
  header W > 8 with a correct model hash (unproducible by any
  conformant writer; reference would attempt them, we return
  "unsupported"), and (b) truncated raw-mode sub-blocks decoded with
  `verify=False` (reference silently returns short output; we return
  "truncated archive").  With verification on, both implementations
  reject identically.
- The stretch goal (porting the impl_cm context-mixing predictor) was
  not attempted within budget; the container/AC layers here are
  reusable for it.
- Reference archives for F3/F4 are cached (`--refdir`); a cold re-run
  regenerates them via `ref_driver.py` in ~25 min on 8 cores.

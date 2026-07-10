# Determinism & conformance suite — findings (H8 / H9 / H14, R-3.x)

Experiment: `conformance-h8-h9-h14` — Part I (tex/papers/core/rnr_core.tex)
Theorems 10.1/10.2; exp-design §2.7/§2.10/§8; engineering-spec R-3.1–R-3.8,
R-5.2.  Gate: `experiments/conformance/run_checks.py` (permanent; fast mode
< 1 min including the container leg, `--full` = the campaign reported
here).  Raw measurements: `experiments/conformance/results/conformance.jsonl`
(bench/results.py format).  Statistics per the frozen methodology:
H8/H9/H14 form the REPRO family of `bench/METHODS.md` §5 — **binary
pass/fail conformance, reported as exact counts, excluded from p-value
correction**; no bootstrap CI is applicable to an identity test (a single
byte divergence fails outright).

## Campaign setup (full mode, 2026-07-10)

- **Test corpus** (exp-design §8.1 operationalized): **107 blocks,
  827,070 bytes total** — 99 corpus blocks (11 per corpus x 9 manifest
  corpora: enwik8, calgary, canterbury, scripts_src, sqlite_synth,
  telemetry_grid, ctrl_urandom, ctrl_zstd, ctrl_base64), sizes 1–16 KiB,
  offsets derived from the frozen master seed (`metrics.derive_seed`,
  tag `corpus`); plus 7 synthetic boundary blocks (empty, 1 byte, n=K,
  n=K+1, minimum K=256, zeros, periodic) and 1 overflow-stress block
  (66 KiB single repeated byte in one sub-block).  The §8.1 100-block
  census = 99 corpus blocks + stress; the 7 boundary blocks are extras.
  Spec SHA-256 `80b3bc044116ca0b…`; per-block SHA-256 in the generated
  `workdir/blocks/spec.json`.  Configs span W in {1,2,3,4,5}, K in
  {256, 512, 1024, 2048, 4096, 131072}; **403 sub-blocks** total
  (228 coded, 175 fail-safe raw — both container modes exercised).
- **Legs** (each leg re-encodes the identical corpus in a fresh process):
  - A, A2: native **arm64 macOS**, Python 3.14.3, numpy 2.4.6,
    PYTHONHASHSEED=0 (A2 = independent process, identical env);
  - B: PYTHONHASHSEED=1; C: PYTHONHASHSEED unset (randomized hashing);
  - D: **x86-64 Linux** via `docker run --platform linux/amd64`, image
    `rnr-conformance-amd64` (Id `sha256:fe4bd34fdd96…26b4bc5a7`, built
    from `python:3.12-slim` digest `sha256:423ed6ab25b1…852199fbf` with
    numpy pinned 2.4.6), Python 3.12.13 — leg D therefore varies **ISA
    (arm64 → x86-64), OS (macOS → Linux), and Python version
    (3.14 → 3.12)** simultaneously; leg wall time 88.1 s;
  - E: native arm64 decode-only pass over leg D's archive files
    (exp-design §8.3 steps 5–6 cross-decode matrix).
- Compared objects per block: archive bytes (H8); per-sub-block SHA-256
  of the **full integer inference stream** — the count-predictor's
  freq[256] integer "logits", total, argmax byte, candidate list, class
  map, exact class-partition sums, coded class, and payload coding
  interval at every position (H14); and the per-sub-block
  **adaptive-state hash ladder** — SHA-256 of the predictor's complete
  count-table state after each sub-block, chained (H9).  The instrumented
  encoder is cross-checked byte-for-byte against `rnr1.pack()` on every
  block of every leg, so all verdicts apply to the actual reference
  implementation, not a divergent copy.

## Results (exact counts; gate output 12 PASS / 0 FAIL / 0 SKIP)

| Check | Quantity | Result |
|---|---|---|
| K2 | round-trip + H_v/sub-block-hash verify (leg A) | 107/107 |
| K3 | same-process double-encode identical | 2/2 probed |
| K4 | across-process identity A==A2 (full report projection) | identical |
| K5 | PYTHONHASHSEED invariance A==B==C | identical |
| K6 | **H8** cross-ISA archive **file bytes** A vs D | **107/107 identical** |
| K7 | **H14** cross-ISA inference-trace identity (827,070 positions) | **107/107 blocks** |
| K8 | **H9** cross-ISA adaptive-state ladder (403 rungs) | **107/107 blocks** |
| K9 | cross-decode: x86-64 decodes native archives / native decodes x86-64 archives | 107/107 and 107/107 |
| K10 | overflow headroom: max coder total on stress block | 2,097,376 = 256+32·(2^16−1), 8.0x below 2^24 |
| K11 | guard live: accept 2^24−1 / reject 2^24; fault-injected scale (2^20) trips dist() guard | all trip correctly |
| K12 | reduction-order probe: shuffled integer reductions == fixed order | 11 probe points × 3 shuffles, exact |

Global comparison chains agreed across all legs (archives chain
`d409011f92deeb14…`, traces chain `c91dba5b9572ec47…`, states chain
`1afd6bacaaa1a75b…`).  Encoder-vs-decoder adaptation agreement (the
decoder recomputes the same inference trace and the same state ladder)
held on all 228 coded sub-blocks on **both** ISAs.

**Secondary, uncorrected** (METHODS.md §1): with 0 divergences in 107
blocks, the one-sided 95% Clopper–Pearson upper bound on the per-block
divergence probability is 1−0.05^(1/107) ≈ **2.8%**; at trace-position
granularity (0 divergences in 827,070 positions) it is ≈ **3.6×10⁻⁶**
per position.  These bounds are descriptive only; the REPRO-family
verdict is the exact count.

## Per-hypothesis verdicts

- **H8 (bit-exact reproduction, Thm 10.1) — SUPPORTED on tested
  platforms.** Byte-identical archives across repeated runs, independent
  processes, three PYTHONHASHSEED regimes, and arm64-macOS vs
  x86-64-Linux (2 of the 6 platforms of exp-design §8.2; see caveats).
  Cross-decode matrix verified in both directions, 107/107.
- **H9 (online-adaptation determinism, Thm 10.2 / Def 10.1) — SUPPORTED
  for the count-based compatible-online path.** The adaptive-count state
  evolution — hashed after every sub-block, 403 rungs — is identical
  across all legs and identical between encoder and decoder.  Sub-blocks
  with up to 67,584 adaptation steps (one integer count-state update per
  byte) exceed H9's 1000-step threshold.
- **H14 (integer inference-path portability) — SUPPORTED for the two CPU
  targets tested.** The complete per-position integer inference stream is
  bit-identical arm64 vs x86-64: identity rate 100% on 827,070 positions
  (H14 asks for 100% on a 1000-input set; this is 827x that scale on the
  CPU targets).
- **R-3.2/R-5.2 (overflow guard) — CONFORMANT.** The stress input drives
  every context into the COUNT_CAP halving regime; the maximum
  arithmetic-coder total observed is exactly 256+32·(2^16−1) = 2,097,376
  (headroom 2^24/2,097,376 = 8.00x), reproduced bit-exactly on both
  ISAs.  The guard is live at the exact boundary (accepts 2^24−1,
  rejects 2^24) and trips under fault injection.
- **R-3.4 (reduction order) — CONFORMANT, with documented coverage.**
  Direct probe passed (class-partition sums, totals, prefix sums
  recomputed in randomized orders equal the fixed-order values exactly).
  Coverage argument: the coding path is integer-only, so every reduction
  is exact and its value order-independent; the implementation fixes the
  order syntactically (np.cumsum / np.add.at / left-to-right loops).  An
  order-sensitive reduction could only arise from floating point and
  would produce different bits on at least one ISA — numpy dispatches
  different SIMD kernels (NEON vs SSE/AVX) with different blocking on
  the two legs — so the per-position cross-ISA trace identity of K7 is
  an end-to-end witness for R-3.4 on real kernels, and K12 is the
  in-situ probe.  Trace records also dtype-audit the path (int64 only,
  R-3.1).

## Caveats (honest scope limits)

1. **Platform coverage is 2 of 6** (exp-design §8.2): ARM64 macOS native
   and x86-64 Linux under Docker emulation.  The x86-64 leg executes
   genuine x86-64 numpy SIMD kernels but through binary translation on
   Apple silicon, not physical x86-64 hardware; AVX-512-vs-baseline,
   RISC-V/QEMU, and Windows/MSVC legs remain untested.  Both legs are
   little-endian.
2. **H14's accelerator targets are untested**: GPU int8 tensor cores,
   ANE, TPU (exp-design targets c–f) are out of reach here; the verdict
   covers the CPU scalar/SIMD integer path only.
3. **H9 is validated for the adaptive-counts instantiation** of
   Definition 10.1 (deterministic integer count update, reset at sync
   points) — the compatible-online path of the reference coder.  The
   integer SGD/signSGD/Adam optimizer variants named in H9 are not
   implemented in `impl/rnr1.py` and are not tested.
4. **Single implementation, single language.** Cross-ISA identity of one
   Python+numpy implementation does not demonstrate cross-implementation
   conformance (two independent codebases meeting the spec).  Python's
   arbitrary-precision integers cannot exhibit the fixed-width wraparound
   a C port could; the numpy int64 path can, and is dtype-audited.
5. **The overflow guard is assertion-based** (`assert` in dist()/
   encode()); under `python -O` it is stripped — the gate detects this
   and SKIPs K11 rather than passing it.  A production implementation
   needs a hard trap per R-5.2.  The guard-trip demonstration is fault
   injection (patched count scale), since by construction the reference
   coder's totals cannot reach 2^24 (8x headroom, K10).
6. **PYTHONHASHSEED legs are a guard, not a proof**: CPython int hashing
   is seed-independent, so hash-seed variation exercises dict-order
   sensitivity only indirectly (the coder uses int-keyed dicts and never
   iterates them in hash order).
7. The predictor is the order-W integer count model of the reference
   coder, not a W-bounded transformer; conformance of a neural-predictor
   implementation to R-3.x is a separate future test.
8. Operational note: during the campaign the docker client process of
   leg D was interrupted once; the container completed normally and the
   gate was re-run with `--resume` (crash recovery that reuses completed
   leg reports; comparisons always re-run over the full leg outputs).
   All quantities above come from complete leg reports and byte-level
   file comparisons.

## Reproduction

```
/Users/para/.venvs/rnr/bin/python -u experiments/conformance/run_checks.py          # fast gate
/Users/para/.venvs/rnr/bin/python -u experiments/conformance/run_checks.py --full   # this campaign
```

Docker unavailable => the four cross-ISA checks report SKIP (never
PASS); a missing image is rebuilt from `docker/Dockerfile`, and if the
build is impossible the gate falls back to pinned `python:3.12-slim`
with a runtime `pip install numpy==2.4.6`, recording digests either way.

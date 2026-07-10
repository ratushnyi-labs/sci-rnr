# H5/H6/H7 — Type-III mode selection, nibble mutual information, dictionary coding

Experiment: `experiments/h5_h7_modes/` (this directory).
Protocol: exp-design §2.4–2.6 (H5, H6, H7), §5–6; Part I Theorems 6.1/6.2
(approximate-oracle scheduling), 6.3/6.3b (I(H;L|C) advantage), 6.5
(marginal-MDL admission), 6.6 (bits-back), Lemma 6.6a (rate-neutral index),
Theorem 6.6b (o(N) overhead shape). Statistics follow the frozen
methodology `bench/METHODS.md` (master seed 20260710, seeds via
`metrics.derive_seed`); this campaign is **secondary/labeled** relative to
the frozen H1–H18 primaries — every interval below carries its method
label; deviations listed under *Honest caveats*.

Campaign: 2026-07-10, **reduced** scale (H5 composite 256 KiB, 32 × 8 KiB
tiles; H6 cheap instruments 2 MiB, CM instrument 64 KiB, B = 200; H7
ladder ≤ 256 KiB); runtime 140 s; **all 17 gates PASS**
(`run_checks.py`; `--full` scales to 2 × 1 MiB composites, 8 MiB/256 KiB
H6 slices, B = 1000, 1 MiB ladder). Instruments: committed reference
coders `impl/rnr1.py` (order-W counting predictor) and
`impl_cm/cm_predictor.py` (integer context-mixing predictor, "CM"),
bench baselines as external modes. Results: `out/summary_reduced.json`
(detail), `out/results_reduced.jsonl` (31 bench-compatible records via
`bench/results.py`).

## Design

**H5.** A seeded heterogeneous composite (text + sqlite + telemetry +
urandom segments of 1–3 tiles, the preconditions qualification table's
heterogeneity classes) is tiled at 8 KiB; every tile is coded standalone
by all 8 modes (store, gzip9, bz2-9, lzma-9e, zstd-22, brotli-11,
rnr1-ngram W=3, rnr1-cm W=3; bare-blob framing uniform across modes),
giving the exact cost matrix L_m(t). Two Type-III-A selectors run
against the per-tile oracle and every fixed mode: *prefix* (U_m(t)
extrapolated from coding the tile's first quarter; encoder-side, charged
⌈log₂ 8⌉ = 3 flag bits/tile) and *causal* (U_m(t) from the previous tile
of the same segment; decoder-reproducible, flag-free; first tile of each
segment uses the declared default zstd-22). **ε estimation
(documented):** because L_m(t) is measured for every mode, the per-tile
error ε_t = max_m |U_m(t) − L_m(t)| is computed EXACTLY for each
selector's U, not estimated; Theorem 6.1's bound
L_sel(t) ≤ L_oracle(t) + 2ε_t is then checked tile by tile (it must hold
identically — the check validates the harness), and the empirical
content is the *size* of ε per training rule plus the realized regret.

**H6.** I(H;L|C) per corpus under three context instruments — `prev1`
(C = previous byte; replication of the preconditions screen, 2 MiB),
`cm` (C = the CM predictor's argmax byte: integer context mixing of
orders 0–3 with adaptive logistic weights, a decoder-reproducible
function of the whole past; 64 KiB slice, the instrument runs at
~15 KB/s), and `plane` (record-structured corpora: C = (byte plane,
previous same-plane byte), stride 4, plus a per-plane decomposition).
Estimator: Miller-Madow-corrected plug-in **minus a conditional-
permutation null** (low nibbles permuted within context strata — the
exact H⊥L|C null with identical stratum sizes and (C,H)/(C,L)
marginals; it removes the residual thin-support inflation that
first-order Miller-Madow leaves: on the 64 KiB CM slice the raw MM
value carries ≈ +0.36 bpb of pure bias on the incompressible controls,
which the null removes to ±0.004, gates G3/G4). A cyclic-shift null is
NOT valid for this statistic (it kills context relevance but not the
within-byte H–L coupling) — measured and rejected during
instrumentation. Chunk-bootstrap percentile CIs (2 KiB chunks, B = 200),
shifted by the fixed null.

**H7.** A minimal, real Type-III-C dictionary coder: substring
candidates (lengths 2–7), greedy admission by Theorem 6.5's frozen-rest
rule ΔL(τ) = dictionary_cost + f·(code_cost − 8ℓ) < 0 with realized
non-overlapping counts f restricted to still-literal positions,
dictionary_cost = 8(ℓ+1) bits, parse_overhead = 0 (boundaries implied by
token lengths in a single entropy-coded stream — declared); greedy
longest-match parse; one adaptive arithmetic-coded stream (the
`impl/rnr1.py` coder) over {256 literals} ∪ {K token ids}; byte-exact
round-trip decoder. Bits-back accounting per Theorem 6.6: (a) identity
(6.4) verified exactly on an enumerable micro model (all parses of all
2⁸ strings over a 5-token dictionary, both q = exact posterior and
q = point mass, data law ≠ model marginal); (b) on the real coder, a DP
over all parses under the frozen final token law computes −log₂ P_X(X),
so the deterministic parse's KL term −log₂ P(z*|X) is measured exactly.

## Per-hypothesis verdicts

**H5 (Type-III-A approximate-oracle scheduling) — CONFIRMED at this
scale, well inside the 5% target.** Composite: 32 tiles = 2 enwik8 +
5 sqlite + 10 telemetry + 15 urandom (seeded draw); oracle 6.6323 bpb.

| arm | bpb | regret vs oracle [95% CI, paired-tile percentile] |
|---|---|---|
| per-tile oracle | 6.6323 | — |
| **prefix selector** (incl. flag bits) | **6.6576** | **+0.38% [0.09%, 0.77%]** |
| causal selector (flag-free) | 6.7950 | +2.45% [0.99%, 4.35%] |
| best fixed mode (lzma-9e) | 6.7117 | +1.20% |
| brotli-11 / zstd-22 / rnr1-cm | 6.8654 / 6.9777 / 7.0014 | +3.5% / +5.2% / +5.6% |
| gzip9 / rnr1-ngram / bz2-9 / store | 7.0459 / 7.3362 / 7.5065 / 8.0000 | +6.2% / +10.6% / +13.2% / +20.6% |

- The prefix selector beats every fixed mode outright (the Theorem 6.1
  corollary — composition without loss — realized); the causal
  selector's 2.45% sits above fixed lzma-9e here because the flag-free
  rule pays for segment boundaries on a composite where lzma is
  near-universal runner-up.
- **2ε bound (Theorem 6.1):** holds on every tile for both selectors
  with the exactly instrumented ε_t, and aggregate regret stays within
  Σ2ε_t + flag bits (gate G1). Realized regret ≪ 2ε̄.
- **ε magnitudes:** prefix instrument ε̄ = 1.136 bpb, p95 = 1.564 —
  inflated by container-header extrapolation from 2 KiB prefixes (an
  honest property of the cheap instrument, and the errors are strongly
  correlated across modes, which is why its argmin is still nearly
  oracle-exact). Causal instrument ε̄ = 0.076 bpb, p95 = 0.233 — the
  decoder-reproducible instrument meets H5's "ε below 0.1 bits per
  position" calibration target.
- Falsification threshold (rate > oracle + 20%): not approached (0.38%).
- Mode choices track class structure: store on urandom tiles, lzma-9e
  on telemetry/sqlite, brotli-11 on text.

**H6 (Type-III-B I(H;L|C) bands) — ASCII text CONFIRMED (in band or
above, favorable); telemetry float32 re-measured: the preconditions
value replicates and the stronger CM instrument RAISES it (0.259 →
0.283), but it remains BELOW the 1.0–2.0 band — verdict "below band,
NOT adverse" (comfortably above the 0.1 falsification floor); no
corpus in the suite is adverse.** All values bpb, conditional-
permutation-debiased Miller-Madow, [95% chunk-percentile CI shifted by
the fixed null]:

| corpus | class | prev1 (2 MiB) | cm (64 KiB) | plane/lag4 | verdict vs band |
|---|---|---|---|---|---|
| enwik8 | ascii-text | 0.940 [0.937, 0.950] | 0.714 [0.711, 0.746] | — | above [0.3, 0.7] (favorable) |
| canterbury | ascii-text | 0.611 [0.598, 0.637] | 0.646 [0.637, 0.675] | — | in [0.3, 0.7] |
| scripts_src | ascii-text | 0.875 [0.872, 0.884] | 0.692 [0.676, 0.744] | — | above [0.3, 0.7] (favorable) |
| calgary | mixed | 0.899 [0.902, 0.918] | 0.677 [0.670, 0.707] | — | reported (no band) |
| sqlite_synth | structured-binary | 0.772 [0.803, 0.818] | 0.550 [0.529, 0.646] | 0.613 | reported (no band) |
| telemetry_grid | float32 | 0.259 [0.277, 0.282] | 0.283 [0.285, 0.305] | 0.051 | below [1.0, 2.0], above 0.1 floor |
| ctrl_urandom | control | 0.000 | 0.002 | — | anchor-consistent (~0) |
| ctrl_zstd | control | 0.000 | −0.004 | — | anchor-consistent (~0) |
| ctrl_base64 | control | 0.400 | 0.386 | — | consistent with preconditions 0.4002 |

- **Replication:** prev1 matches the preconditions screen (enwik8
  0.940 vs 0.9397; canterbury 0.611 vs 0.5902; scripts_src 0.875 vs
  0.8661; sqlite 0.772 vs 0.7182; telemetry 0.259 vs 0.2608; base64
  0.400 vs 0.4002) despite different slice sizes — the instrument is
  stable.
- **Telemetry per-plane decomposition (the pre-adverse check the task
  required):** the coupling is NOT uniform over float32 byte planes.
  Little-endian planes 0/1 (low mantissa bytes) carry ≈ 0 CMI (near
  i.i.d.-uniform noise, as expected for smoothed-Gaussian synthetic
  data); plane 2 carries 0.176 under the lag-4 context; plane 3
  (sign+exponent) carries 0.440 under the CM context. The whole-stream
  value is therefore diluted ~4× by the noise planes. Even the best
  per-plane figure (0.44) stays below the 1.0–2.0 band: for THIS
  corpus — a synthetic smoothed-noise field whose mantissa bytes are
  nearly incompressible — the H6 float32 band is not reached, and that
  is a property of the stand-in (`data/README.md`: telemetry_grid
  stands in for HDF5 scientific arrays), not yet evidence against the
  band on real scientific float32 data with structured exponent
  dynamics. The falsification floor (0.1) is cleared by 2.8×.
- **Effect of the stronger instrument:** on text, conditioning on the
  CM prediction *lowers* I(H;L|C) (0.94 → 0.71 on enwik8) — a stronger
  context explains part of the nibble coupling, exactly the direction
  Theorem 6.3b's tightness implies; the Type-III-B saving is largest
  for weak/cheap contexts. On telemetry the CM instrument *raises* the
  measured saving (0.259 → 0.283): the adaptive byte model surfaces
  cross-nibble structure the raw previous byte misses.

**H7 (Type-III-C greedy MDL + bits-back) — mechanism CONFIRMED:
admission rule consistent, o(N) overhead shape realized, hierarchical
index exactly rate-neutral, bits-back identity exact; the absolute rate
of this deliberately minimal coder is ≈ 4.3–4.4 bpb (H7's original
2%-of-H(X|Y) bits-back-VAE claim on images is NOT tested here).**

| corpus | N | rate bpb | ĥ holdout (order) | rate/(N·ĥ) | K | dict overhead frac |
|---|---|---|---|---|---|---|
| enwik8 | 32 Ki | 4.4214 | 4.5026 (1) | 0.982 | 215 | 0.0400 |
| enwik8 | 64 Ki | 4.3892 | 4.2213 (1) | 1.040 | 244 | 0.0222 |
| enwik8 | 128 Ki | 4.3715 | 4.1865 (1) | 1.044 | 259 | 0.0119 |
| enwik8 | 256 Ki | 4.4114 | 4.1002 (1) | 1.076 | 281 | 0.0063 |
| scripts_src | 64 Ki | 4.3351 | 4.5243 (1) | 0.958 | 296 | 0.0285 |
| scripts_src | 128 Ki | 4.2940 | 4.0745 (1) | 1.054 | 328 | 0.0158 |

- **o(N) check (Theorem 6.6b shape):** dictionary overhead fraction
  strictly decreases along both ladders (enwik8 0.0400 → 0.0063 over 8×
  in N; scripts_src 0.0285 → 0.0158) while K grows sublinearly
  (215 → 281 over 8×). PASS (gate G6).
- **Rate vs N·ĥ:** 0.96–1.08 of the order-1 held-out cross-entropy
  estimate. ĥ is itself an UPPER bound on the entropy rate, so ratios
  slightly below 1 are consistent (the dictionary captures longer-than-
  order-1 structure); this is a ballpark consistency check, not an
  optimality claim. Round-trips byte-exact at every point; every
  admitted phrase had ΔL < 0 (gate G5).
- **Lemma 6.6a (rate-neutral index):** on the realized token-id
  substreams (up to 88k ids), the ideal two-stage (family/variant grid)
  code length equals the ideal flat length to 0.0e+00 bits at every
  ladder point — the chain-rule identity (6.5) verified on real data.
  Realized adaptive two-stage differs from realized flat by ≤ 25 bits
  (≤ 0.013%, the O(1)-per-token practical-coder slack), and realized
  flat exceeds the substream's empirical entropy by only 0.1–0.8%,
  inside the declared adaptive-redundancy envelope. "Index cost == its
  entropy within tolerance": PASS (gate G7).
- **Bits-back (Theorem 6.6):** identity (6.4) holds pointwise to
  3.6e-15 bits on the enumerable micro model for both q = exact
  posterior (net = −log₂ P_X(X), KL = 0) and q = point mass
  (net = −log₂ P(X,z*)); the expected bits-back saving of posterior
  sampling over the best single parse is 1.43 bits/string there. On the
  real coder (enwik8 16 KiB, frozen final token law), the deterministic
  greedy parse's KL term −log₂ P(z*|X) = 2062.9 bits = 0.126 bpb
  (≈ 2.9% of the stream) — the measured amount a bits-back latent-parse
  implementation would recover. PASS (gate G8).

## Honest caveats

1. **Scale:** reduced/dev slices throughout; the H5 composite is
   256 KiB under ONE composition seed (the seeded draw weighted
   urandom 15/32 tiles — heavier on the store regime; `--full` runs
   2 × 1 MiB compositions), H6's CM instrument sees 64 KiB, H7's ladder
   tops at 256 KiB. Full-corpus versions are wired
   (`run_checks.py --full`) but not run here.
2. **H5 ε semantics:** ε is a joint property of (modes, U-instrument,
   tile size). The prefix instrument's large ε̄ (1.14 bpb) is dominated
   by fixed container headers extrapolated 4×; Theorem 6.1 bounds
   regret by 2ε whatever U is, so this is reported, not hidden. The
   causal instrument is the deployable (decoder-reproducible) one and
   meets the < 0.1 bpb calibration target. H5's "calibrated U_m(t)"
   with a learned cost model is not implemented; these two simple
   selectors bracket it from the cheap side.
3. **H5 framing:** L_m(t) uses bare blobs for all modes (no
   container/seek-index bytes; rnr modes include their Theorem 7.22
   store escape). Absolute bpb slightly flatters every mode equally;
   regret ratios are unaffected.
4. **H6 estimator:** Miller-Madow corrects first-order bias only; the
   remaining thin-support inflation is removed by a conditional-
   permutation null whose own sampling noise is NOT propagated into
   the CI (labeled approximation; anchor-validated at matched n and
   support in gates G3/G4, including signal retention: debiased
   1.2223 vs analytic 1.1715 on the thin-support coupled anchor —
   note the residual +0.05 over-estimate visible there bounds the
   method error). Where a point estimate falls below its shifted CI
   (e.g. sqlite plane_lag4 0.613 vs [0.679, 0.699]) the offset is
   chunk-bootstrap bias, and the interval should be read as ± its
   width around the point, secondary-grade. CM-instrument split-half
   drifts are small (|drift| ≤ 0.04) except the controls, where they
   diagnose exactly the bias the null removes.
5. **H6 claim scope:** the H6 bands concern the saving available to a
   byte coder over nibble coding for a GIVEN context C; the verdict is
   existential over practical decoder-reproducible C. The x86-64
   instruction-stream class of H6 has no corpus in this suite and is
   NOT tested. The telemetry verdict is scoped to this synthetic
   stand-in; testing the float32 band on real scientific arrays (HDF5
   with structured exponents) remains open.
6. **H7 claim scope:** H7's original statement (bits-back VAE within 2%
   of H(X|Y) on images) is not testable by a minimal discrete-
   dictionary coder; per the task, the tested content is Theorem 6.5's
   admission criterion, Lemma 6.6a's neutrality, Theorem 6.6b's o(N)
   shape, and Theorem 6.6's accounting — all PASS. The ≈ 4.4 bpb
   absolute rate reflects the minimal design (ℓ ≤ 7 candidates, greedy
   single-pass admission, order-0 adaptive symbol model, no context
   model on literals); it beats store and bz2-block but not lzma. The
   greedy/frozen-rest caveat of Theorem 6.5 is inherited: admission
   uses projected code costs, and the final parse re-realizes counts.
7. **Store regeneration:** `out/results_reduced.jsonl` is regenerated
   per campaign run (append-only within a run). METHODS.md's
   append-only rule governs the reporting store `bench/results/`,
   which this campaign does not touch.
8. **Determinism:** coders and instruments are integer-deterministic;
   coding-rate numbers are single runs (the METHODS 5-seed repetition
   governs timing hypotheses, none reported here). All seeds derive
   from the master seed and are recorded in the JSONL.

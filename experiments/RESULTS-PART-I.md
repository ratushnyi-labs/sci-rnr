# RNR Part I — Consolidated Empirical Results

<!-- sources: methods, exp-design -->

Synthesis report, 2026-07-10. This document consolidates the five Part-I
findings files into one hypothesis-by-hypothesis reconciliation against
(a) the pre-registered targets and falsification thresholds of the
experimental-design companion (`tex/rnr_experimental_design.tex`, §2) and
(b) the concrete predictions of Part I §11 (`tex/papers/core/rnr_core.tex`).

**No new measurements were taken for this report.** Every number below is
copied verbatim (or truncated, never re-rounded) from a cited findings
file; the gate `experiments/run_checks.py` verifies mechanically that each
quoted number appears in its cited source and that the mandatory sections
of this report are present.

**Statistical status (per `bench/METHODS.md` v1.0, frozen 2026-07-10).**
Primary-grade numbers in this report are: (i) the H2/H3/H4 paired
coding-gap BCa intervals (frozen procedure, B = 10000, master seed
20260710, one designated primary comparison per hypothesis), and (ii) the
conformance exact counts (REPRO family: binary pass/fail, excluded from
p-value correction). Everything else — the precondition screen, the
H5/H6/H7 campaign, the H1 context-mixing instrument check, all TC/CMI
uncertainty intervals, and the Clopper–Pearson divergence bounds — is
**secondary/labeled** per METHODS §1 and is marked as such where quoted.

Findings files consolidated (the only measurement sources):

| key | file | scope |
|---|---|---|
| preconditions | `experiments/preconditions/findings.md` | corpus qualification screen; P-ADV; sync-overhead datum |
| h2_h4 | `experiments/h2_h4_typeii/findings.md` | H2/H3/H4 Type-II separation, TC identity, permutation control |
| conformance | `experiments/conformance/findings.md` | H8/H9/H14 determinism and cross-ISA conformance |
| impl_cm | `impl_cm/findings.md` | H1 first measurement, context-mixing instrument |
| h5_h7 | `experiments/h5_h7_modes/findings.md` | H5/H6/H7 mode selection, nibble CMI, dictionary coding |

Verdict vocabulary: **CONFIRMED** (pre-registered statement met at the
tested scope) / **SUPPORTED** (conformance-style pass on the tested
subset) / **INSTRUMENT-SCOPED** (mechanism verified, literal claim awaits
the stated instrument) / **BELOW-BAND-NOT-ADVERSE** (below the predicted
band, above the falsification floor) / **ADVERSE** (a pre-declared
threshold crossed against the hypothesis, scope stated) /
**UNTESTED-EXTERNAL** (requires equipment/scale outside this program).

## 1. Hypothesis-by-hypothesis verdict table

<!-- sources: preconditions, h2_h4, conformance, impl_cm, h5_h7, exp-design, core-s11, methods -->

| Hypothesis (family) | Exp-design target → falsification | Key measurements (95% CI where applicable) | Statistical grade | Verdict | Source |
|---|---|---|---|---|---|
| Precondition screen P-STAT / P-MIX / P-PRED | pre-declared: spread ≤ 0.05 / ≤ 0.2; far-\|acf\| ≤ 0.05 / ≤ 0.15; relative drop ≥ 0.15 / ≥ 0.05 | enwik8 and scripts_src qualify (spread 0.97% / 1.18%; drop 35.9% / 36.9%); calgary / canterbury fail P-STAT (46.66% / 30.10% — qualifying for H5); telemetry_grid fails P-PRED (drop 1.3%) → H4 descoped by manifest patch | secondary/labeled screen (BCa, B = 10000) | CONFIRMED (as qualification instrument; per-corpus scoping delivered) | preconditions |
| Precondition P-ADV (Def. 2.1 advantage) | qualifies iff paired-diff 95% BCa CI < 0 AND full-archive ratio_block < 1 | every natural corpus fails: paired diff vs zstd-block from +0.2609 [+0.1842, +0.3738] (calgary) to +0.8423 [+0.8126, +0.8802] (scripts_src); enwik8 +0.6208 [+0.5768, +0.6704]; ratio_block 1.0827–1.5743; controls corrected to marginal | secondary/labeled; Holm within screen family only | ADVERSE at counting-instrument grade (scoped: the paper's claim is neural-predictor-scoped; see §3) | preconditions |
| H1 Type-I rate (RATE, quantifying) | 1.5–2.0 bpb with the 8–16 MB transformer; LZMA-max 2.1–2.5, NNCP 1.3–1.6 → falsified if > 2.5 bpb on enwik9 or > NNCP + 0.5 after E12-C | CM instrument, 1 MiB dev slices: 2.907 / 2.762 / 2.888 / 2.589 bpb (enwik8_dev / calgary_text / canterbury_text / scripts_src); beats the plain order-W predictor by 0.73–0.94 bpb everywhere; coder redundancy ≤ 7.5e-5 bpb; container+sync ~0.007 bpb; brotli-11 leads by 0.35–0.36 bpb | instrument check scale; single-seed deterministic rates, labeled | INSTRUMENT-SCOPED (band untested; full claim UNTESTED-EXTERNAL pending the transformer predictor) | impl_cm |
| H2 Type-II separation (RATE, quantifying) | Type-II coder 10–20% below LZMA-max on DB pages / records / instruction streams → falsified if rate ≥ LZMA on DB pages | joint beats factorized at matched W: sqlite_hdr +0.142 [0.116, 0.171] bpb (corrected +0.796); record_log +1.371 [1.364, 1.378]; columnar_int16 +2.414 [2.405, 2.422]; Cohen d 2.5 / 90.3 / 137.4; lzma control confirms joint layout structure (0.777 vs 1.072 planar on columnar) | paired BCa B = 10000, frozen (primary comparison: sqlite_hdr structured gap) | mechanism CONFIRMED; literal rate-vs-LZMA claim INSTRUMENT-SCOPED (needs the semi-non-factorized Type-II coder) | h2_h4 |
| H3 TC = Θ(N) (INFO, confirming) | TC per byte ≥ 0.1 bits over 1–16 KiB blocks → falsified per class if TC = o(N) | sqlite_hdr ≥ 0.760 (chain-MI lower bound); record_log 0.947–1.045 (analytic 0.9444); columnar 2.379 (analytic 2.3789); calib_dup 1.000 (analytic 1.0000); prefix-TC linear: columnar 5.02 → 16.03 → 27.05 → 38.06 bits at r = 4/8/12/16; telemetry_expo TC 0.069 < 0.1 | gap CIs frozen BCa; TC intervals secondary (percentile-multinomial-grassberger) | CONFIRMED for the structured-binary classes; ADVERSE for the telemetry_expo (b3,b2) decomposition (threshold not met); 1–16 KiB block scaling still open (small-N proxy) | h2_h4 |
| H4 manifold restriction (INFO, confirming) | Levina–Bickel intrinsic dimension d ≪ N → falsified if d within factor 2 of N | operational permutation control: destroying joint structure drives chain-MI to ≤ 0.0006 bpb on synthetics (vs 0.95–2.19 structured) and collapses the coding gap below the structured gap with separated CIs on every corpus and seed; permuted gap ≤ 0 everywhere | paired BCa frozen (primary: sqlite_hdr control separation) | CONFIRMED in the task-operational (permutation-control) form; the pre-registered Levina–Bickel estimator form is untested (exp-design §9) | h2_h4 |
| H5 Type-III-A scheduling (RATE, quantifying) | within 5% of per-tile oracle; ε̄ < 0.1 bits/position → falsified if regret > 20% | prefix selector regret +0.38% [0.09%, 0.77%] vs oracle 6.6323 bpb, beating every fixed mode (best fixed lzma-9e +1.20%); causal selector +2.45% [0.99%, 4.35%]; causal ε̄ = 0.076 bpb, p95 = 0.233 (meets the 0.1 target); prefix ε̄ = 1.136 (cheap-instrument inflation, documented); Theorem 6.1's 2ε bound holds on every tile | secondary/labeled campaign, reduced scale (256 KiB composite, one seed) | CONFIRMED at tested scale | h5_h7 |
| H6 I(H;L\|C) bands (INFO, quantifying) | ASCII 0.3–0.7 bpb; x86-64 0.4–0.9; float32 1.0–2.0 → falsified if < 0.1 on any class | ASCII in/above band: enwik8 0.940 [0.937, 0.950] (prev1) / 0.714 [0.711, 0.746] (cm); canterbury 0.611 / 0.646; scripts_src 0.875 / 0.692; float32 stand-in below band: telemetry 0.259 (prev1) → 0.283 (cm), best plane 0.440 — above the 0.1 floor; controls at the ~0 anchor | secondary/labeled; debiased Miller–Madow + conditional-permutation null | ASCII CONFIRMED (above-band favorable); float32 BELOW-BAND-NOT-ADVERSE (synthetic stand-in); x86-64 class untested (no corpus) | h5_h7 |
| H7 Type-III-C bits-back (RATE, confirming) | rate within 2% of H(X\|Y) with bits-back VAE on images → falsified if > +5% | bits-back identity (6.4) exact to 3.6e-15 bits on the enumerable micro model; greedy parse KL term 0.126 bpb (≈ 2.9%) measured exactly; dictionary overhead o(N): 0.0400 → 0.0063 over 8× N (K sublinear 215 → 281); index exactly rate-neutral (0.0e+00 ideal); minimal coder absolute rate ≈ 4.3–4.4 bpb | secondary/labeled campaign | mechanism CONFIRMED (Thm 6.5 / 6.6 / 6.6a / 6.6b accounting); the original VAE-on-images 2% claim INSTRUMENT-SCOPED / untested | h5_h7 |
| H8 byte-identical archives (REPRO, confirming) | 6 platforms × 100-block corpus, all archive hashes equal → any byte divergence fails | 107/107 blocks byte-identical across arm64-macOS vs x86-64-Linux legs, independent processes, three hash-seed regimes; cross-decode matrix 107/107 both directions; Clopper–Pearson 95% upper bound 2.8% per block (secondary, descriptive) | exact counts (REPRO family, primary grade) | SUPPORTED on the 2 of 6 platforms tested | conformance |
| H9 online-adaptation determinism (REPRO, confirming) | byte-identical archives after 1000 adaptation steps, same 6 platforms | adaptive-state hash ladder identical across all legs, 403 rungs; encoder-vs-decoder state agreement on all 228 coded sub-blocks on both ISAs; up to 67,584 adaptation steps per sub-block (67× the H9 threshold) | exact counts (REPRO family) | SUPPORTED for the count-based compatible-online path; integer SGD/signSGD/Adam variants untested | conformance |
| H14 inference-path portability (REPRO, confirming) | 100% byte-identity on a 1000-input set across targets (a)–(f) → any divergence falsifies | complete integer inference stream bit-identical arm64 vs x86-64 on 827,070 positions, 107/107 blocks (827× the pre-registered scale on the CPU targets); per-position CP bound ≈ 3.6×10⁻⁶ (secondary) | exact counts (REPRO family) | SUPPORTED for CPU scalar/SIMD targets; GPU / ANE / TPU legs UNTESTED-EXTERNAL | conformance |

Notes to the table.

1. The H10(c) sync-overhead datum (6.904–10.629% at cold-reset instrument
   grade vs the 2% target) is not a Part-I-scoped hypothesis row; it is
   carried in full in the ledger, §3, item 2.
2. H15 (Chinchilla-class enwik9) and the H1 full claim are
   UNTESTED-EXTERNAL; see §3 items 5–6. H16 is touched only through the
   precondition controls (see §2, item 11.6).

## 2. Section-11 prediction reconciliation

<!-- sources: core-s11, exp-design, preconditions, h2_h4, impl_cm, h5_h7 -->

Each item quotes the concrete Part I §11 prediction our measurements
touch, then states agree / partial / open with the numbers.

**§11 preamble (H10–H11 cross-reference).** Quote: "the central
value-proposition claim that the side-information+repair recasting beats
existing dictionary-based seekable formats on compression ratio remains an
empirical conjecture pending the measurement specified there." —
**OPEN, with a first adverse-at-instrument datum.** The P-ADV screen is
the first instrumented measurement touching this claim: at
counting-predictor grade, the paired per-block difference vs zstd-22 at
matched 64 KiB cold-start granularity is positive (RNR worse) on every
natural corpus, +0.2609 to +0.8423 bpb, ratio_block 1.0827–1.5743
(preconditions §1). Scoped: the claim is conditioned on decoder-reproducible
neural side information with amortized L(M), which this instrument does
not carry (preconditions caveat 2); the conjecture stands, now with a
measured lower bar for the neural instrument to clear.

**§11.1 (H1).** Quote: "predicted to achieve 1.5 to 2.0 bits per byte
before model amortization, competitive with LZMA at maximum preset (2.1 to
2.5 bits per byte) and behind NNCP (1.3 to 1.6 bits per byte on similar
corpora)." — **OPEN (instrument-scoped progress).** The context-mixing
instrument reaches 2.589–2.907 bpb on 1 MiB text/code dev slices —
outside the predicted band, but the gap decomposition shows the coding
chain adds ≤ 7.5e-5 bpb over the model's own cross-entropy (2.899 / 2.755
/ 2.879 / 2.579 bpb), so the remaining distance to the band is entirely
predictor quality, the exact term the 8–16 MB transformer upgrade targets
(impl_cm). Local baselines on the same slices: lzma-9e 2.930 / 2.686 /
2.955 / 2.384, brotli-11 2.547–2.680. No falsification condition is
touched (the rule is stated on enwik9 with the named predictor).

**§11.2 (H2/H3/H4).** Quote: "predicted to achieve coding rates 10 to 20%
below LZMA, because LZMA exploits string repetition that is weak in these
classes while Type-II exploits positional regularity." — **PARTIAL:
mechanism agrees, literal comparison open.** The positional-regularity
mechanism is confirmed at matched context: joint-over-factorized gaps
+0.142 [0.116, 0.171] bpb on real DB page headers, +1.371 on fixed-record
logs, +2.414 on columnar int16, with the permutation control collapsing
each gap (h2_h4). Positional periodicity is present in the wild at the
predicted strides (lags 4096 = SQLite page, 8192 = telemetry row, 77 =
base64 line; preconditions §3). The literal 10–20%-below-LZMA statement
needs the semi-non-factorized Type-II coder, which the Type-I reference
instrument is not (h2_h4 verdict). Supporting conditions: TC = Θ(N)
confirmed on the structured-binary classes (H3), and the prediction's own
caveat class — data whose structure a sequential model cannot express —
is instantiated honestly by Theorem 5.4's boundary (h2_h4 note 7).

**§11.3 (H5).** Quote: "predicted to outperform any single fixed mode on
archives with heterogeneous content ... provided the per-mode prediction
error and the per-tile mode-flag overhead together remain smaller than the
heterogeneity advantage." — **AGREE at tested scale.** On the 256 KiB
heterogeneous composite the prefix selector (flag bits charged) codes at
6.6576 bpb vs oracle 6.6323, regret +0.38% [0.09%, 0.77%], and beats every
fixed mode outright — best fixed mode lzma-9e sits at +1.20% (h5_h7). The
proviso is measured, not assumed: the decoder-reproducible causal
instrument has ε̄ = 0.076 bpb (< the 0.1 calibration target), and the
Theorem 6.1 2ε bound holds tile-by-tile. Reduced scale, one composition
seed; the composition-can-lose branch of the prediction was not exercised.

**§11.4 (H6).** Quote: "predicted to recover 0.3 to 0.7 bits per byte over
independent nibble coding on ASCII text ... On structured numeric arrays
(e.g., float32 columns), Type-III-B may recover larger gains." —
**AGREE on ASCII (favorably above band); float32 part NOT REALIZED on the
synthetic stand-in.** Measured I(H;L|C) on ASCII text: canterbury 0.611
(in band), enwik8 0.940 and scripts_src 0.875 (above band — favorable
direction); under the stronger CM context 0.646–0.714 (h5_h7). On the
float32 stand-in the coupling is 0.259 (prev1) / 0.283 (cm), diluted ~4×
by near-uniform low-mantissa planes (best plane, sign+exponent under CM:
0.440) — below the exp-design 1.0–2.0 band, 2.8× above the 0.1
falsification floor. The stand-in is synthetic smoothed noise; the band on
real scientific float32 arrays remains open (h5_h7 caveat 5).

**§11.5 (H7).** Quote: "predicted to achieve coding rates competitive with
the best dictionary-based methods (zstd at maximum level, LZMA at maximum
preset), with the advantage growing as the corpus size grows and the
dictionary amortizes." — **PARTIAL: amortization shape agrees;
competitive-rate open.** The dictionary-overhead fraction strictly
decreases along both ladders (enwik8 0.0400 → 0.0063 over 8× in N;
scripts_src 0.0285 → 0.0158) with K sublinear (215 → 281) — the predicted
amortization realized (h5_h7, gate G6). The absolute rate of the
deliberately minimal coder (≈ 4.3–4.4 bpb) beats store and bz2-block but
not lzma, so the competitive-rate half of the prediction is untested at
this instrument grade (h5_h7 caveat 6). The underlying accounting —
admission rule, rate-neutral index, bits-back identity to 3.6e-15 bits —
is exact.

**§11.6 (negative controls, H16 touch).** Quote: "On random or encrypted
data, all RNR types are predicted to revert to baseline plus overhead ...
On already-compressed data (gzip, xz, zstd output), RNR is predicted to be
slightly worse than passthrough." — **AGREE at instrument grade
(secondary datum).** On ctrl_urandom and ctrl_zstd the reference coder
emits 8.0043 bpb total with store fraction 1.00 — raw plus container
overhead, slightly worse than passthrough exactly as predicted; ladder
drop 0.000 confirms incompressibility (preconditions §4–5). One scoped
observation: on ctrl_base64 (a calibration control, not an H16 class) the
store escape also triggered (8.0043 bpb against H0 = 6.0222), leaving
base64 structure unexploited at this instrument grade — fail-safe behavior
held, compression did not (paired diff +1.9591 vs zstd-block). The H16
primary (SAFETY family) remains with the frozen conformance program.

**§11.7 (H15).** Quote: "encoder output is ≈ 79 MB for the 1 GB enwik9
source, with practical CLT-based concentration ± 0.04 MB at δ = 10⁻⁶ ...
a 12.0× bitstream-rate compression." — **UNTESTED-EXTERNAL.** Requires a
Chinchilla-70B-class predictor (H_M' ≈ 0.664 bpb); no local instrument
approaches this regime. enwik9 itself was deliberately not fetched
(preconditions caveat 6), so the H15 scope is untouched by this program.

## 3. Adverse / open ledger

<!-- sources: preconditions, h2_h4, h5_h7, impl_cm, conformance, exp-design, core-s11 -->

Honest list of every adverse, below-band, or open item, with scope.

1. **P-ADV counting-instrument outcome — ADVERSE (scoped).** Every
   natural corpus fails the Def. 2.1 advantage screen: paired per-block
   diff vs zstd-22 positive with 95% BCa CIs above 0 on all of
   enwik8/calgary/canterbury/scripts_src/sqlite_synth/telemetry_grid
   (+0.2609 to +0.8423 bpb), full-archive ratio_block 1.0827–1.5743, and
   the post-correction rule moves even the store-mode controls to
   marginal (preconditions §1 + correction note). Scope: the instrument
   is a deliberately weak order-W counting predictor without the
   amortized neural L(M) the paper's conditional claim assumes
   (preconditions caveat 2); a P-ADV failure does not falsify the scoped
   claim, but the advantage is *not demonstrated* at this grade, and the
   burden now sits with the neural instrument.
2. **Sync-point overhead 6.904–10.629% vs the 2% H10(c) target —
   ADVERSE at cold-reset instrument grade (scoped).** At K = 64 KiB the
   coded corpora show total-archive overhead 6.904% (enwik8), 6.934%
   (scripts_src), 10.604% (canterbury), 10.629% (calgary) — the last two
   above H10's 10% falsification threshold — while sqlite_synth at
   1.442% already meets the target even under this instrument
   (preconditions §2). Scope, stated precisely: H10(c) is a claim about
   archives built with a **W-bounded transformer predictor with fixed
   weights** (W = 4 KiB), where a sync point costs only the excess code
   length of the ≤ W bytes coded with truncated context plus an index
   entry — the weights stay warm. The counting instrument instead zeroes
   its *entire learned state* at each sync (its count tables are its
   weights), so each sync pays a full re-learning transient: the
   measured 6.904–10.629% is an upper bound from the worst-case
   cold-context regime (preconditions §2 note). **What remains to be
   shown:** for a fixed-weight W-bounded predictor at (W = 4 KiB,
   K = 64 KiB), the truncated-context excess over the first W bytes
   after each sync, plus index bytes, stays below 2% of the no-sync
   archive size — a context-truncation measurement, not a state
   re-learning one.
3. **Telemetry below-band datum (H6 float32) — BELOW-BAND-NOT-ADVERSE.**
   0.259 → 0.283 bpb (prev1 → cm), best byte-plane 0.440, against the
   1.0–2.0 band; 2.8× above the 0.1 falsification floor. Scoped to a
   synthetic smoothed-noise stand-in whose mantissa bytes are near
   incompressible; not yet evidence against the band on real scientific
   float32 data (h5_h7).
4. **H3 telemetry_expo decomposition — ADVERSE for that class.**
   TC = 0.069 bpb < the 0.1 threshold for the (b3,b2) sign/exponent pair;
   Theorem 5.5's separation regime is weak there and the measurement
   shows it (structured gap negative; only the corrected gap +0.054 ≈ TC
   0.069 is positive). The precondition screen independently descoped
   telemetry from H4 (manifest patch) (h2_h4; preconditions §6).
5. **Literal H2 — UNTESTED.** The 10–20%-below-LZMA statement requires a
   semi-non-factorized Type-II coder with a positional model; the Type-I
   reference coder is not one. Mechanism evidence only (h2_h4).
6. **Neural-instrument externals — UNTESTED-EXTERNAL.** (a) H1's full
   claim (1.5–2.0 bpb with the shared 8–16 MB transformer; NNCP
   comparison); (b) H15 (Chinchilla-class enwik9, 79 ± 0.04 MB, ≤ 50 ms
   random-access on a B100-class GPU). Neither is falsified nor supported
   by anything measured here; the impl_cm gap decomposition shows
   precisely which term (predictor cross-entropy) the upgrade must move.
7. **H4 pre-registered form — UNTESTED.** The Levina–Bickel
   intrinsic-dimension estimation (exp-design §9, with Two-NN /
   Grassberger–Procaccia cross-checks) was not run; only the operational
   permutation control was (h2_h4).
8. **H6 x86-64 instruction-stream class — UNTESTED.** No corpus in the
   suite exercises the 0.4–0.9 band (h5_h7 caveat 5).
9. **H7 original claim — UNTESTED.** Bits-back VAE within 2% of H(X|Y)
   on image data is not testable by the minimal discrete-dictionary
   instrument (h5_h7 caveat 6).
10. **Conformance coverage — OPEN remainder.** 2 of 6 exp-design §8.2
    platforms tested (AVX-512, RISC-V/QEMU, Windows/MSVC, Asahi legs
    open; both tested legs little-endian); H9's integer SGD/signSGD/Adam
    variants unimplemented; H14's accelerator targets (GPU int8 tensor
    cores, ANE, TPU) out of reach here (conformance caveats 1–3).
11. **Scale caveats.** All coding-rate work at dev/reduced scale
    (preconditions ≤ 1024 KiB coding slices; h5_h7 256 KiB composite,
    one seed; impl_cm 1 MiB slices; h2_h4 records 2–16 bytes vs the
    pre-registered 1–16 KiB blocks). Full-corpus campaigns are wired
    (`--full` modes) but not run.

## 4. Instrument inventory

<!-- sources: preconditions, h2_h4, conformance, impl_cm, h5_h7, exp-design -->

| Instrument | What it is | Used for | What it can / cannot decide |
|---|---|---|---|
| Counting n-gram predictor (`impl/rnr1.py`) | order-W adaptive integer count model behind the reference Type-I container, K-KiB sync, store escape | preconditions screen (P-ADV, sync overhead, indicators); both arms of H2/H3/H4; an H5 mode; the conformance corpus | can verify container accounting, separation mechanism, determinism; cannot demonstrate the Def. 2.1 advantage (no amortized neural L(M)) nor warm-predictor sync overhead (cold state reset per sync) |
| Context-mixing predictor (`impl_cm/cm_predictor.py`) | integer-deterministic logistic mixing of order-0..W byte models, fixed-point, spec R-3.1–R-3.8 conformant, same Predictor interface | H1 first measurement; an H5 mode; the strong H6 context instrument | closes 0.73–0.94 bpb of the n-gram gap and isolates the remaining gap as predictor quality; cannot reach the H1 band (2.589–2.907 bpb measured) |
| Minimal Type-III-C dictionary coder (`experiments/h5_h7_modes/h7_dict.py`) | greedy MDL admission (Thm 6.5 frozen-rest rule), single adaptive stream, byte-exact decoder | H7 mechanism: admission, o(N) shape, rate-neutral index, bits-back accounting | verifies identities and shapes exactly; absolute rate (≈ 4.3–4.4 bpb) deliberately non-competitive |
| Estimator stack (`bench/metrics.py`, experiment-local estimators) | frozen BCa/paired-bootstrap machinery; Grassberger + chain-MI TC estimators with closed-form anchors; Miller–Madow + conditional-permutation-null CMI; held-out entropy ladder | all uncertainty intervals; H3/H6 information measurements | anchor-validated on closed forms (e.g. calib_dup TC 1.0001 vs analytic 1.0000); TC/CMI intervals remain secondary-grade |
| Conformance harness (`experiments/conformance/`) | dual-ISA legs (arm64 native, x86-64 containerized), full integer inference-stream and state-ladder hashing | H8/H9/H14 | exact-count byte-identity on the tested platform pair; cannot speak to unimplemented optimizer variants or accelerator targets |

**What the transformer upgrade would decide** (the single most valuable
next instrument — the shared 8–16 MB W-bounded byte-level transformer of
exp-design §4.3): (i) H1's 1.5–2.0 bpb band, and with it the §11.1
prediction; (ii) the P-ADV / Def. 2.1 advantage screen in its intended
amortized-L(M) regime — the value-proposition datum; (iii) the H10(c)
warm-predictor sync overhead (ledger item 2's remaining measurement);
(iv) the neural-conformance leg of H8/H14 (R-3.x for a neural inference
path, conformance caveat 7); and (v) whether the H6 CM-vs-neural trend
(stronger context absorbing nibble coupling, 0.940 → 0.714 on enwik8)
continues, which calibrates the Type-III-B advantage under deployed
predictors.

## 5. Follow-ups

<!-- sources: preconditions, h2_h4, conformance, impl_cm, h5_h7, exp-design -->

One line per discrepancy; each is a self-contained actionable task.

1. Build or adapt the 8–16 MB W-bounded transformer predictor behind the
   `Predictor` interface and rerun the P-ADV screen + H1 band on the full
   corpora (decides ledger items 1 and 6a).
2. Measure sync overhead at (W = 4 KiB, K = 64 KiB) with a fixed-weight
   W-bounded predictor — truncated-context excess only — against the 2%
   H10(c) target (ledger item 2).
3. Implement the semi-non-factorized Type-II coder (positional model) and
   run the literal H2 comparison vs lzma-9e on sqlite/parquet-class
   corpora (ledger item 5).
4. Source a real scientific float32 corpus (HDF5 with structured exponent
   dynamics) to replace the telemetry stand-in; re-run the H6 float32
   band and the H3/H4 telemetry scoping (ledger items 3–4).
5. Run the exp-design §9 intrinsic-dimension estimation (Levina–Bickel +
   Two-NN cross-check) on the qualified corpora for the pre-registered H4
   form (ledger item 7).
6. Add an x86-64 instruction-stream corpus and test the H6 0.4–0.9 band
   (ledger item 8).
7. Scale the h2_h4 record framework to 1–16 KiB blocks for the
   pre-registered Θ(N) test (ledger item 11 / h2_h4 note 3).
8. Run the h5_h7 `--full` campaign (2 × 1 MiB composites, more
   composition seeds, B = 1000, 1 MiB H7 ladder, 8 MiB H6 slices).
9. Extend conformance to the four remaining §8.2 platforms and implement
   one integer-optimizer H9 variant (ledger item 10).
10. Build the bits-back VAE instrument for H7's original
    2%-of-H(X|Y)-on-images claim (CIFAR-10 per exp-design §3.3)
    (ledger item 9).
11. Investigate the ctrl_base64 store-mode escape (8.0043 bpb emitted
    against H0 = 6.0222): determine why the coded path lost to store at
    this instrument grade and whether the mode-selection margin needs a
    finer threshold.
12. External: the H15 enwik9 Chinchilla-class campaign remains outside
    this program's equipment envelope; keep it flagged external, do not
    fold into local milestones.

## 6. Proposed section-11 caveat lines (for the author; tex/ untouched)

<!-- sources: core-s11, exp-design, preconditions, h2_h4, impl_cm, h5_h7 -->

This report does not edit `tex/`. The following are the exact caveat
sentences Part I §11 should gain, one per touched prediction, written in
the paper's register and consistent with the findings files.

**§11 preamble** (append after the sentence naming H10–H12): "A first
instrumented screen (order-W counting predictor as side information,
dev-slice scale) did not demonstrate the per-block advantage of
Definition 2.1 — the paired per-block difference against zstd at matched
64 KiB cold-start granularity was positive on every natural corpus tested
(+0.2609 to +0.8423 bits per byte) — consistent with, but not yet
resolving, the neural-predictor scoping of the claim; the same screen
measured sync-point overhead of 6.904–10.629% at K = 64 KiB under
worst-case cold-context resets (full predictor-state reset per sync
point), above Hypothesis H10(c)'s 2% target and, on two corpora, its 10%
falsification threshold, so the H10(c) claim now rests specifically on
the fixed-weight W-bounded predictor regime, where a sync point costs
only truncated context over at most W bytes rather than a re-learning
transient."

**§11.1** (append): "A first instrument-grade measurement (an integer
context-mixing predictor inside the reference container, 1 MiB dev
slices) achieves 2.589–2.907 bits per byte with coder redundancy below
1e-4 bits per byte over the model's own cross-entropy, isolating the
remaining distance to the predicted band entirely in predictor quality;
the 1.5–2.0 band itself remains untested pending the 8–16 MB transformer
predictor."

**§11.2** (append): "The separation mechanism has been confirmed at
instrument grade — joint-over-factorized coding gaps of +0.142 to +2.414
bits per byte at matched context on structured-binary corpora, with 95%
confidence intervals excluding zero and a marginal-preserving permutation
control collapsing the gap — while the literal 10–20%-below-LZMA
comparison awaits a Type-II semi-non-factorized implementation; the
supporting H3 condition held on the structured-binary classes but failed
for a float32 sign/exponent decomposition (total correlation 0.069 bits
per byte, below the 0.1 threshold)."

**§11.3** (append): "At reduced scale (a 256 KiB heterogeneous
composite), a prefix-probe Type-III-A selector realized regret +0.38%
against the per-tile oracle and outperformed every fixed single mode,
with the decoder-reproducible variant's prediction error ε̄ = 0.076 bits
per byte meeting the calibration target — consistent with this
prediction, at one composition seed."

**§11.4** (append): "Measured I(H;L|C) on ASCII text is 0.611–0.940 bits
per byte (in or above the predicted band); on a synthetic float32
stand-in it is 0.259–0.283 (below the 1.0–2.0 band of the companion's
H6 but above its 0.1 falsification floor, the shortfall traced to
near-incompressible mantissa byte planes), leaving the numeric-array part
of this prediction open on real scientific data."

**§11.5** (append): "A minimal dictionary coder has confirmed the
admission rule, the o(N) dictionary-overhead shape (overhead fraction
0.0400 → 0.0063 over an 8× corpus-size ladder), and exact bits-back
accounting; the competitive-rate half of this prediction remains
untested at that instrument grade."

**§11.6** (append): "Confirmed at instrument grade: on random and
precompressed controls the reference coder emits 8.0043 bits per byte
with store fraction 1.00 — raw plus container overhead, slightly worse
than passthrough exactly as predicted."

**§11.7** (append): "The prediction remains external to the local
program: enwik9 was deliberately not fetched and no local instrument
approaches the Chinchilla-class regime, so nothing measured here bears on
it in either direction."

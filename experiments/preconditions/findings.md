# Precondition statistics — qualification of the corpus suite

Campaign scale: **full** (coding slice 1024 KiB, W in [2, 3, 4], K = 64 KiB; indicator slice 8 MiB; front-of-file slices). Generated 2026-07-10T15:36:40+00:00 by `experiments/preconditions/measure_preconditions.py`.

Resolves, at dev-slice scale, the measurements named external by `docs/backlog/blocked/BUG-001-D` (conditional-advantage ratio of the Abstract's scoped claim, Definition 2.1 / inequality (10.1)), `BUG-003-D` (sync-point ratio overhead at typical (K, W), exp-design H10(c)), and `BUG-010-D` (which corpora satisfy the section-10.2.1 precondition). Statistics follow the frozen methodology (`bench/METHODS.md` v1.0): BCa bootstrap, B = 10000, derived seeds; this screen is **secondary/labeled** relative to the frozen H1–H18 primary families, so no amendment is required.

## 1. Conditional-advantage ratio (BUG-001-D instrument)

Instrument: reference Type-I coder (`impl/rnr1.py`), order-W adaptive counting predictor as side information, sync K = 64 KiB. Primary comparison: paired per-block repair bpb vs per-block zstd-22 on the same 64 KiB blocks (both coders restart cold at each block — the in-process proxy for seekable formats at matched granularity). `ratio_block` charges the FULL RNR archive (header+index+hashes) against bare baseline blobs (conservative); `ratio_stream` compares against the best whole-slice streaming baseline (harder). Def. 2.1 advantage needs ratio < 1.

| corpus | best W | RNR total bpb | zstd-block bpb | lzma-block bpb | paired diff bpb [95% BCa] | p (boot, secondary) | Holm reject | ratio_block | stream best | ratio_stream |
|---|---|---|---|---|---|---|---|---|---|---|
| enwik8 | 2 | 3.4241 | 2.7990 | 2.7142 | +0.6208 [+0.5768, +0.6704] | 0.0002 | yes | 1.2233 | brotli 2.2358 | 1.5315 |
| calgary | 2 | 3.3437 | 3.0785 | 3.0306 | +0.2609 [+0.1842, +0.3738] | 0.0002 | yes | 1.0861 | bz2 2.3768 | 1.4068 |
| canterbury | 2 | 2.2269 | 1.4145 | 1.1961 | +0.8081 [+0.5871, +0.9471] | 0.0002 | yes | 1.5743 | lzma 1.0714 | 2.0784 |
| scripts_src | 2 | 3.1279 | 2.2813 | 2.2396 | +0.8423 [+0.8126, +0.8802] | 0.0002 | yes | 1.3711 | brotli 1.5791 | 1.9808 |
| sqlite_synth | 2 | 4.4663 | 3.9584 | 3.2369 | +0.5037 [+0.4708, +0.6096] | 0.0002 | yes | 1.1283 | lzma 2.8244 | 1.5813 |
| telemetry_grid | 2 | 8.0043 | 7.3927 | 6.5476 | +0.6073 [+0.5932, +0.6206] | 0.0002 | yes | 1.0827 | lzma 6.6242 | 1.2083 |
| ctrl_urandom | 2 | 8.0043 | 8.0012 | 8.0078 | -0.0012 [-0.0012, -0.0012] | 0.0002 | yes | 1.0004 | brotli 8.0000 | 1.0005 |
| ctrl_zstd | 2 | 8.0043 | 8.0012 | 8.0078 | -0.0012 [-0.0012, -0.0012] | 0.0002 | yes | 1.0004 | brotli 8.0000 | 1.0005 |
| ctrl_base64 | 2 | 8.0043 | 6.0409 | 6.1886 | +1.9591 [+1.9583, +1.9598] | 0.0002 | yes | 1.3250 | brotli 6.0280 | 1.3278 |

### 1a. Inequality (10.1) term breakdown at best W (bits/byte)

| corpus | (a) H(X|Y) proxy | (b) model | (c) rho | (d) nu | (e) mu | total | coded frac |
|---|---|---|---|---|---|---|---|
| enwik8 | 3.4198 | 0.000061 | 0.00007 | 0.00104 | 0.00317 | 3.4241 | 1.00 |
| calgary | 3.3394 | 0.000061 | 0.00006 | 0.00104 | 0.00317 | 3.3437 | 1.00 |
| canterbury | 2.2225 | 0.000061 | 0.00006 | 0.00104 | 0.00317 | 2.2269 | 1.00 |
| scripts_src | 3.1235 | 0.000061 | 0.00007 | 0.00104 | 0.00317 | 3.1279 | 1.00 |
| sqlite_synth | 4.4620 | 0.000061 | 0.00006 | 0.00104 | 0.00317 | 4.4663 | 1.00 |
| telemetry_grid | n/a | 0.000061 | n/a | 0.00104 | 0.00317 | 8.0043 | 0.00 |
| ctrl_urandom | n/a | 0.000061 | n/a | 0.00104 | 0.00317 | 8.0043 | 0.00 |
| ctrl_zstd | n/a | 0.000061 | n/a | 0.00104 | 0.00317 | 8.0043 | 0.00 |
| ctrl_base64 | n/a | 0.000061 | n/a | 0.00104 | 0.00317 | 8.0043 | 0.00 |

Term (b) is only the 8-byte model-description hash: the counting predictor carries no trained parameters, so the neural L(M)/mN term of (10.1) is NOT exercised by this instrument (it would only make the inequality harder). (a)+(c) sum to the repair rate over coded blocks; store-mode (raw) blocks pay 8 bpb.

## 2. Sync-point overhead at K = 64 KiB, W = 3 (BUG-003-D datum)

| corpus | window | overhead (total archive) | overhead (repair only) | note |
|---|---|---|---|---|
| enwik8 | 128 KiB | 6.904% | 6.854% | coded; cold-context + index cost per sync |
| calgary | 128 KiB | 10.629% | 10.580% | coded; cold-context + index cost per sync |
| canterbury | 128 KiB | 10.604% | 10.553% | coded; cold-context + index cost per sync |
| scripts_src | 128 KiB | 6.934% | 6.881% | coded; cold-context + index cost per sync |
| sqlite_synth | 128 KiB | 1.442% | 1.403% | coded; cold-context + index cost per sync |
| telemetry_grid | 128 KiB | 0.024% | 0.000% | store-mode data; overhead = index entries only (no cold-context term) |
| ctrl_urandom | 128 KiB | 0.024% | 0.000% | store-mode data; overhead = index entries only (no cold-context term) |
| ctrl_zstd | 128 KiB | 0.024% | 0.000% | store-mode data; overhead = index entries only (no cold-context term) |
| ctrl_base64 | 128 KiB | 0.024% | 0.000% | store-mode data; overhead = index entries only (no cold-context term) |

H10(c) target: < 2% at K = 64 KiB; falsification threshold: > 10%. Measured on a 2-block window (the per-sync cold-context cost is deterministic and scale-free per sync point); the adaptive counting predictor restarts from empty counts at each sync, which is the worst-case cold-context regime for this instrument.

## 3. Stationarity / mixing indicators

| corpus | per-third spread (max of H0,H2) | verdict | max abs acf (lags >= 256) | periodic lags | verdict |
|---|---|---|---|---|---|
| enwik8 | 0.97% | qualifies | 0.0174 | — | qualifies |
| calgary | 46.66% | fails | 0.3104 | [77, 4096, 8192] | fails |
| canterbury | 30.10% | fails | 0.3753 | [77, 4096, 8192] | fails |
| scripts_src | 1.18% | qualifies | 0.1058 | [77] | marginal |
| sqlite_synth | 12.18% | marginal | 0.1494 | [77, 4096, 8192] | marginal |
| telemetry_grid | 0.23% | qualifies | 0.3961 | [8192] | fails |
| ctrl_urandom | 0.00% | qualifies | 0.0008 | — | qualifies |
| ctrl_zstd | 0.00% | qualifies | 0.0005 | — | qualifies |
| ctrl_base64 | 0.01% | qualifies | 0.0019 | [77] | qualifies |

A large per-third spread is DISqualifying for stationary-source hypotheses (H1/H3) but QUALIFYING for the heterogeneous-scheduler hypothesis (H5). Likewise, strong far-lag autocorrelation (positional periodicity: lag 4096 = SQLite page, 8192 = telemetry row stride, 77 = base64 line) is qualifying evidence for Type-II positional structure (H2). Full ACF table: `out/indicators_full.csv`.

## 4. Entropy-rate ladder (held-out order-0..4) and I(H;L|C)

| corpus | H0 | min CE @ order | drop bpb | drop % | coverage@min | I(H;L|C) |
|---|---|---|---|---|---|---|
| enwik8 | 5.1130 | 3.2750 @ 3 | 1.8381 | 35.9% | 0.989 | 0.9397 |
| calgary | 6.6697 | 4.2174 @ 2 | 2.4522 | 36.8% | 0.878 | 0.8312 |
| canterbury | 5.1254 | 3.5906 @ 3 | 1.5348 | 29.9% | 0.925 | 0.5902 |
| scripts_src | 5.0932 | 3.2154 @ 2 | 1.8778 | 36.9% | 0.998 | 0.8661 |
| sqlite_synth | 7.0803 | 4.8412 @ 2 | 2.2391 | 31.6% | 0.996 | 0.7182 |
| telemetry_grid | 7.3620 | 7.2678 @ 1 | 0.0942 | 1.3% | 1.000 | 0.2608 |
| ctrl_urandom | 8.0000 | 8.0000 @ 4 | 0.0000 | 0.0% | 0.001 | 0.0000 |
| ctrl_zstd | 7.9993 | 7.9993 @ 0 | 0.0000 | 0.0% | 1.000 | 0.0002 |
| ctrl_base64 | 6.0222 | 6.0222 @ 0 | 0.0000 | 0.0% | 1.000 | 0.4002 |

Held-out cross-entropy (counts from the first half, Laplace-smoothed evaluation on the second half) is a real code length, hence an honest UPPER bound on the order-k entropy rate — no plug-in sparsity under-bias. drop = H0 − min_k CE_k lower-bounds the total-correlation rate: TC(D_N)/N >= H0 − h_k − o(1), the H3 screening quantity. Per-order values: `out/ladder_full.csv`.

## 5. QUALIFICATION TABLE (corpus x precondition)

Thresholds (pre-declared): P-ADV qualifies iff the paired-diff 95% BCa CI lies below 0 (fails iff above 0); P-STAT spread <= 0.05 / <= 0.2; P-MIX far-|acf| <= 0.05 / <= 0.15; P-PRED relative drop >= 0.15 / >= 0.05.

| corpus | class | P-ADV | P-STAT | P-MIX | P-PRED | per-hypothesis verdicts |
|---|---|---|---|---|---|---|
| enwik8 | 3.1-text-source | fails | qualifies | qualifies | qualifies | H1: **qualifies**, H3: **qualifies**, H4: **screen-pass**, H6: **qualifies**, H8: **qualifies** |
| calgary | 3.4-heterogeneous-negative | fails | fails | fails | qualifies | H1: **marginal**, H5: **qualifies**, H8: **qualifies** |
| canterbury | 3.1-text-source | fails | fails | fails | qualifies | H1: **marginal**, H3: **qualifies**, H4: **screen-pass**, H6: **qualifies**, H8: **qualifies** |
| scripts_src | 3.1-text-source | fails | qualifies | marginal | qualifies | H1: **qualifies**, H3: **qualifies**, H4: **screen-pass**, H6: **qualifies** |
| sqlite_synth | 3.2-structured-binary | fails | marginal | marginal | qualifies | H2: **qualifies**, H3: **qualifies**, H4: **screen-pass**, H8: **qualifies** |
| telemetry_grid | 3.3-image-numeric | fails | qualifies | fails | fails | H4: **fails**, H6: **qualifies** |
| ctrl_urandom | 3.4-heterogeneous-negative | qualifies | qualifies | qualifies | fails | H16: **qualifies** |
| ctrl_zstd | 3.4-heterogeneous-negative | qualifies | qualifies | qualifies | fails | H16: **qualifies** |
| ctrl_base64 | 3.4-heterogeneous-negative | fails | qualifies | qualifies | fails | H1: **qualifies** |

## 6. Per-hypothesis scoping conclusions

**enwik8** — prior `qualifies_for`: ['H1', 'H3', 'H4', 'H6', 'H8']
- H1: qualifies — held-out drop 36% and per-third spread 1.0%
- H3: qualifies — TC-rate lower bound H0-minCE = 1.84 bpb (TC(D_N) >= 1.84*N bits; H3 needs Theta(N))
- H4: screen-pass — intrinsic-dimension estimator not run here; predictive-structure screen qualifies
- H6: qualifies — I(H;L|C)=0.940 bpb, predicted band 0.3-0.7 (ASCII), falsification floor 0.1
- H8: qualifies — byte-identity conformance corpus role is data-agnostic
- MANIFEST patch: no change (all prior scopes survive the screen)

**calgary** — prior `qualifies_for`: ['H1', 'H5', 'H8']
- H1: marginal — predictive structure qualifies, stationarity fails
- H5: qualifies — heterogeneity confirmed: per-third spread 46.7% (scheduler target)
- H8: qualifies — byte-identity conformance corpus role is data-agnostic
- MANIFEST patch: no change (all prior scopes survive the screen)

**canterbury** — prior `qualifies_for`: ['H1', 'H3', 'H4', 'H8']
- H1: marginal — predictive structure qualifies, stationarity fails
- H3: qualifies — TC-rate lower bound H0-minCE = 1.53 bpb (TC(D_N) >= 1.53*N bits; H3 needs Theta(N))
- H4: screen-pass — intrinsic-dimension estimator not run here; predictive-structure screen qualifies
- H6: qualifies — I(H;L|C)=0.590 bpb, predicted band 0.3-0.7 (ASCII), falsification floor 0.1
- H8: qualifies — byte-identity conformance corpus role is data-agnostic
- MANIFEST patch: remove nothing, add ['H6'] -> proposed ['H1', 'H3', 'H4', 'H8', 'H6']

**scripts_src** — prior `qualifies_for`: ['H1', 'H3', 'H4']
- H1: qualifies — held-out drop 37% and per-third spread 1.2%
- H3: qualifies — TC-rate lower bound H0-minCE = 1.88 bpb (TC(D_N) >= 1.88*N bits; H3 needs Theta(N))
- H4: screen-pass — intrinsic-dimension estimator not run here; predictive-structure screen qualifies
- H6: qualifies — I(H;L|C)=0.866 bpb, predicted band 0.3-0.7 (ASCII), falsification floor 0.1
- MANIFEST patch: remove nothing, add ['H6'] -> proposed ['H1', 'H3', 'H4', 'H6']

**sqlite_synth** — prior `qualifies_for`: ['H2', 'H3', 'H4', 'H8']
- H2: qualifies — positional periodicity at lags [77, 4096, 8192], held-out drop 2.24 bpb
- H3: qualifies — TC-rate lower bound H0-minCE = 2.24 bpb (TC(D_N) >= 2.24*N bits; H3 needs Theta(N))
- H4: screen-pass — intrinsic-dimension estimator not run here; predictive-structure screen qualifies
- H8: qualifies — byte-identity conformance corpus role is data-agnostic
- MANIFEST patch: no change (all prior scopes survive the screen)

**telemetry_grid** — prior `qualifies_for`: ['H4', 'H6']
- H4: fails — no low-dimensional predictive structure signal (screen); LB estimator not run here
- H6: qualifies — I(H;L|C)=0.261 bpb, predicted band 1.0-2.0 (float32), falsification floor 0.1
- MANIFEST patch: remove ['H4'], add nothing -> proposed ['H6']

**ctrl_urandom** — prior `qualifies_for`: ['H16']
- H16: qualifies — incompressible confirmed: best baseline 8.000 bpb, ladder drop 0.000; RNR store-mode fraction 1.00
- MANIFEST patch: no change (all prior scopes survive the screen)

**ctrl_zstd** — prior `qualifies_for`: ['H16']
- H16: qualifies — incompressible confirmed: best baseline 8.000 bpb, ladder drop 0.000; RNR store-mode fraction 1.00
- MANIFEST patch: no change (all prior scopes survive the screen)

**ctrl_base64** — prior `qualifies_for`: ['H1']
- H1: qualifies — calibration control: measured H0=6.022 vs closed-form ~5.99 bpb (inside [5.9,6.1])
- MANIFEST patch: no change (all prior scopes survive the screen)

Patch file: `out/manifest_qualifies_patch.json` (data/ is NOT edited directly).

## 7. Honest caveats

1. **Scale**: dev slices (coding <= 1024 KiB, indicators <= 8 MiB, front-of-file). BUG-001-D/BUG-010-D name FULL-corpus campaigns with a neural predictor in the amortized regime; this campaign is the local, instrumented version, and the external full-scale items remain open. Front-of-file slices bias heterogeneous tars (calgary: the coding slice covers only the first member(s); indicator slices cover the whole file for corpora <= the indicator budget).
2. **Instrument strength**: the order-W counting predictor is a deliberately weak side-information instrument. A P-ADV failure does NOT falsify the paper's conditional claim (the claim is scoped to decoder-reproducible NEURAL side information with amortized L(M)); it means the advantage is not demonstrated by this instrument at this scale. A P-ADV pass, conversely, is strong evidence the precondition is satisfiable.
3. **Determinism**: coder and indicators are deterministic; coding-rate numbers are single runs (the METHODS 5-seed repetition governs TIMING hypotheses, none reported here). All bootstrap seeds derive from the master seed per METHODS section 4 and are recorded in the JSONL.
4. **Secondary status**: the precondition screen is not one of the frozen H1–H18 primary tests; all p-values are labeled secondary, with Holm-Bonferroni applied within this screening family (9 corpora) for transparency only.
5. **Plug-in bias**: per-third ladder values are plug-in estimates (biased low at high orders/sparse contexts); qualification uses held-out cross-entropy where bias direction is safe, and per-third plug-in values only through the RELATIVE spread.
6. **enwik9**: deferred payload, not fetched, not measured; the H15 scope is untouched.
7. **Store-mode corpora**: for incompressible controls the sync-overhead datum reduces to index bytes (no cold-context term); reported as such, not blended with coded corpora.


## Correction (same run, before any downstream use)

The P-ADV verdict rule was tightened to the Def-2.1 full-cost form — `qualifies` requires BOTH the paired per-block 95% CI below 0 AND full-archive `ratio_block` < 1 — because the store-mode negative controls otherwise 'qualify' by undercutting the baseline's per-block framing overhead with their raw fallback (8.0000 bpb/block vs zstd-block ~8.001).

Corrected P-ADV verdicts (measured numbers unchanged): **ctrl_urandom**: qualifies -> marginal, **ctrl_zstd**: qualifies -> marginal. `qualification.csv` / `out/qualification_full.csv` are rewritten with the corrected labels; for the corpora listed here the P-ADV cells of the section-5 table above are superseded by this note.

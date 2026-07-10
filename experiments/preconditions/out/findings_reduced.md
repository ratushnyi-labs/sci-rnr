# Precondition statistics — qualification of the corpus suite

Campaign scale: **reduced** (coding slice 256 KiB, W in [2, 3], K = 64 KiB; indicator slice 2 MiB; front-of-file slices). Generated 2026-07-10T15:48:35+00:00 by `experiments/preconditions/measure_preconditions.py`.

Resolves, at dev-slice scale, the measurements named external by `docs/backlog/blocked/BUG-001-D` (conditional-advantage ratio of the Abstract's scoped claim, Definition 2.1 / inequality (10.1)), `BUG-003-D` (sync-point ratio overhead at typical (K, W), exp-design H10(c)), and `BUG-010-D` (which corpora satisfy the section-10.2.1 precondition). Statistics follow the frozen methodology (`bench/METHODS.md` v1.0): BCa bootstrap, B = 10000, derived seeds; this screen is **secondary/labeled** relative to the frozen H1–H18 primary families, so no amendment is required.

## 1. Conditional-advantage ratio (BUG-001-D instrument)

Instrument: reference Type-I coder (`impl/rnr1.py`), order-W adaptive counting predictor as side information, sync K = 64 KiB. Primary comparison: paired per-block repair bpb vs per-block zstd-22 on the same 64 KiB blocks (both coders restart cold at each block — the in-process proxy for seekable formats at matched granularity). `ratio_block` charges the FULL RNR archive (header+index+hashes) against bare baseline blobs (conservative); `ratio_stream` compares against the best whole-slice streaming baseline (harder). Def. 2.1 advantage needs ratio < 1.

| corpus | best W | RNR total bpb | zstd-block bpb | lzma-block bpb | paired diff bpb [95% BCa] | p (boot, secondary) | Holm reject | ratio_block | stream best | ratio_stream |
|---|---|---|---|---|---|---|---|---|---|---|
| enwik8 | 2 | 3.3544 | 2.6896 | 2.6091 | +0.6594 [+0.5377, +0.7892] | 0.0002 | yes | 1.2472 | brotli 2.2371 | 1.4994 |
| calgary | 2 | 3.3010 | 2.9413 | 2.8732 | +0.3544 [+0.1558, +0.5529] | 0.0002 | yes | 1.1223 | bz2 2.4320 | 1.3573 |
| canterbury | 2 | 3.1255 | 2.9074 | 2.8690 | +0.2127 [+0.1329, +0.2604] | 0.0002 | yes | 1.0750 | bz2 2.3859 | 1.3100 |
| scripts_src | 2 | 3.2232 | 2.3650 | 2.3232 | +0.8528 [+0.8421, +0.8597] | 0.0002 | yes | 1.3629 | brotli 2.0512 | 1.5714 |
| sqlite_synth | 2 | 4.6026 | 4.0444 | 3.3501 | +0.5528 [+0.4321, +0.7926] | 0.0002 | yes | 1.1380 | lzma 3.1002 | 1.4846 |
| telemetry_grid | 2 | 8.0054 | 7.3701 | 6.4677 | +0.6299 [+0.6213, +0.6384] | 0.0002 | yes | 1.0862 | lzma 6.3862 | 1.2535 |
| ctrl_urandom | 2 | 8.0054 | 8.0012 | 8.0078 | -0.0012 [-0.0012, -0.0012] | 0.0002 | yes | 1.0005 | brotli 8.0002 | 1.0007 |
| ctrl_zstd | 2 | 8.0054 | 8.0012 | 8.0078 | -0.0012 [-0.0012, -0.0012] | 0.0002 | yes | 1.0005 | brotli 8.0002 | 1.0007 |
| ctrl_base64 | 2 | 8.0054 | 6.0407 | 6.1882 | +1.9593 [+1.9570, +1.9611] | 0.0002 | yes | 1.3252 | brotli 6.0287 | 1.3279 |

### 1a. Inequality (10.1) term breakdown at best W (bits/byte)

| corpus | (a) H(X|Y) proxy | (b) model | (c) rho | (d) nu | (e) mu | total | coded frac |
|---|---|---|---|---|---|---|---|
| enwik8 | 3.3489 | 0.000244 | 0.00006 | 0.00122 | 0.00391 | 3.3544 | 1.00 |
| calgary | 3.2956 | 0.000244 | 0.00007 | 0.00122 | 0.00391 | 3.3010 | 1.00 |
| canterbury | 3.1201 | 0.000244 | 0.00005 | 0.00122 | 0.00391 | 3.1255 | 1.00 |
| scripts_src | 3.2177 | 0.000244 | 0.00010 | 0.00122 | 0.00391 | 3.2232 | 1.00 |
| sqlite_synth | 4.5972 | 0.000244 | 0.00003 | 0.00122 | 0.00391 | 4.6026 | 1.00 |
| telemetry_grid | n/a | 0.000244 | n/a | 0.00122 | 0.00391 | 8.0054 | 0.00 |
| ctrl_urandom | n/a | 0.000244 | n/a | 0.00122 | 0.00391 | 8.0054 | 0.00 |
| ctrl_zstd | n/a | 0.000244 | n/a | 0.00122 | 0.00391 | 8.0054 | 0.00 |
| ctrl_base64 | n/a | 0.000244 | n/a | 0.00122 | 0.00391 | 8.0054 | 0.00 |

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
| enwik8 | 3.65% | qualifies | 0.0226 | — | qualifies |
| calgary | 26.43% | fails | 0.1099 | [4096, 8192] | marginal |
| canterbury | 50.52% | fails | 0.4234 | [77, 4096, 8192] | fails |
| scripts_src | 2.50% | qualifies | 0.0884 | [77] | marginal |
| sqlite_synth | 0.66% | qualifies | 0.1092 | [77, 4096, 8192] | marginal |
| telemetry_grid | 0.53% | qualifies | 0.3977 | [8192] | fails |
| ctrl_urandom | 0.02% | qualifies | 0.0009 | — | qualifies |
| ctrl_zstd | 0.05% | qualifies | 0.0004 | — | qualifies |
| ctrl_base64 | 0.01% | qualifies | 0.0021 | [77] | qualifies |

A large per-third spread is DISqualifying for stationary-source hypotheses (H1/H3) but QUALIFYING for the heterogeneous-scheduler hypothesis (H5). Likewise, strong far-lag autocorrelation (positional periodicity: lag 4096 = SQLite page, 8192 = telemetry row stride, 77 = base64 line) is qualifying evidence for Type-II positional structure (H2). Full ACF table: `out/indicators_reduced.csv`.

## 4. Entropy-rate ladder (held-out order-0..4) and I(H;L|C)

| corpus | H0 | min CE @ order | drop bpb | drop % | coverage@min | I(H;L|C) |
|---|---|---|---|---|---|---|
| enwik8 | 5.0888 | 3.5350 @ 2 | 1.5537 | 30.5% | 0.997 | 0.9424 |
| calgary | 6.9681 | 4.9374 @ 2 | 2.0308 | 29.1% | 0.796 | 0.9059 |
| canterbury | 5.7036 | 3.7436 @ 2 | 1.9599 | 34.4% | 0.954 | 0.6113 |
| scripts_src | 5.0586 | 3.4235 @ 2 | 1.6351 | 32.3% | 0.994 | 0.8754 |
| sqlite_synth | 7.1094 | 4.5438 @ 2 | 2.5656 | 36.1% | 0.979 | 0.7752 |
| telemetry_grid | 7.3721 | 7.3139 @ 1 | 0.0582 | 0.8% | 1.000 | 0.2593 |
| ctrl_urandom | 8.0002 | 8.0000 @ 4 | 0.0002 | 0.0% | 0.000 | 0.0002 |
| ctrl_zstd | 7.9990 | 7.9990 @ 0 | 0.0000 | 0.0% | 1.000 | 0.0005 |
| ctrl_base64 | 6.0224 | 6.0224 @ 0 | 0.0000 | 0.0% | 1.000 | 0.3997 |

Held-out cross-entropy (counts from the first half, Laplace-smoothed evaluation on the second half) is a real code length, hence an honest UPPER bound on the order-k entropy rate — no plug-in sparsity under-bias. drop = H0 − min_k CE_k lower-bounds the total-correlation rate: TC(D_N)/N >= H0 − h_k − o(1), the H3 screening quantity. Per-order values: `out/ladder_reduced.csv`.

## 5. QUALIFICATION TABLE (corpus x precondition)

Thresholds (pre-declared): P-ADV qualifies iff the paired-diff 95% BCa CI lies below 0 (fails iff above 0); P-STAT spread <= 0.05 / <= 0.2; P-MIX far-|acf| <= 0.05 / <= 0.15; P-PRED relative drop >= 0.15 / >= 0.05.

| corpus | class | P-ADV | P-STAT | P-MIX | P-PRED | per-hypothesis verdicts |
|---|---|---|---|---|---|---|
| enwik8 | 3.1-text-source | fails | qualifies | qualifies | qualifies | H1: **qualifies**, H3: **qualifies**, H4: **screen-pass**, H6: **qualifies**, H8: **qualifies** |
| calgary | 3.4-heterogeneous-negative | fails | fails | marginal | qualifies | H1: **marginal**, H5: **qualifies**, H8: **qualifies** |
| canterbury | 3.1-text-source | fails | fails | fails | qualifies | H1: **marginal**, H3: **qualifies**, H4: **screen-pass**, H6: **qualifies**, H8: **qualifies** |
| scripts_src | 3.1-text-source | fails | qualifies | marginal | qualifies | H1: **qualifies**, H3: **qualifies**, H4: **screen-pass**, H6: **qualifies** |
| sqlite_synth | 3.2-structured-binary | fails | qualifies | marginal | qualifies | H2: **qualifies**, H3: **qualifies**, H4: **screen-pass**, H8: **qualifies** |
| telemetry_grid | 3.3-image-numeric | fails | qualifies | fails | fails | H4: **fails**, H6: **qualifies** |
| ctrl_urandom | 3.4-heterogeneous-negative | marginal | qualifies | qualifies | fails | H16: **qualifies** |
| ctrl_zstd | 3.4-heterogeneous-negative | marginal | qualifies | qualifies | fails | H16: **qualifies** |
| ctrl_base64 | 3.4-heterogeneous-negative | fails | qualifies | qualifies | fails | H1: **qualifies** |

## 6. Per-hypothesis scoping conclusions

**enwik8** — prior `qualifies_for`: ['H1', 'H3', 'H4', 'H6', 'H8']
- H1: qualifies — held-out drop 31% and per-third spread 3.6%
- H3: qualifies — TC-rate lower bound H0-minCE = 1.55 bpb (TC(D_N) >= 1.55*N bits; H3 needs Theta(N))
- H4: screen-pass — intrinsic-dimension estimator not run here; predictive-structure screen qualifies
- H6: qualifies — I(H;L|C)=0.942 bpb, predicted band 0.3-0.7 (ASCII), falsification floor 0.1
- H8: qualifies — byte-identity conformance corpus role is data-agnostic
- MANIFEST patch: no change (all prior scopes survive the screen)

**calgary** — prior `qualifies_for`: ['H1', 'H5', 'H8']
- H1: marginal — predictive structure qualifies, stationarity fails
- H5: qualifies — heterogeneity confirmed: per-third spread 26.4% (scheduler target)
- H8: qualifies — byte-identity conformance corpus role is data-agnostic
- MANIFEST patch: no change (all prior scopes survive the screen)

**canterbury** — prior `qualifies_for`: ['H1', 'H3', 'H4', 'H8']
- H1: marginal — predictive structure qualifies, stationarity fails
- H3: qualifies — TC-rate lower bound H0-minCE = 1.96 bpb (TC(D_N) >= 1.96*N bits; H3 needs Theta(N))
- H4: screen-pass — intrinsic-dimension estimator not run here; predictive-structure screen qualifies
- H6: qualifies — I(H;L|C)=0.611 bpb, predicted band 0.3-0.7 (ASCII), falsification floor 0.1
- H8: qualifies — byte-identity conformance corpus role is data-agnostic
- MANIFEST patch: remove nothing, add ['H6'] -> proposed ['H1', 'H3', 'H4', 'H8', 'H6']

**scripts_src** — prior `qualifies_for`: ['H1', 'H3', 'H4']
- H1: qualifies — held-out drop 32% and per-third spread 2.5%
- H3: qualifies — TC-rate lower bound H0-minCE = 1.64 bpb (TC(D_N) >= 1.64*N bits; H3 needs Theta(N))
- H4: screen-pass — intrinsic-dimension estimator not run here; predictive-structure screen qualifies
- H6: qualifies — I(H;L|C)=0.875 bpb, predicted band 0.3-0.7 (ASCII), falsification floor 0.1
- MANIFEST patch: remove nothing, add ['H6'] -> proposed ['H1', 'H3', 'H4', 'H6']

**sqlite_synth** — prior `qualifies_for`: ['H2', 'H3', 'H4', 'H8']
- H2: qualifies — positional periodicity at lags [77, 4096, 8192], held-out drop 2.57 bpb
- H3: qualifies — TC-rate lower bound H0-minCE = 2.57 bpb (TC(D_N) >= 2.57*N bits; H3 needs Theta(N))
- H4: screen-pass — intrinsic-dimension estimator not run here; predictive-structure screen qualifies
- H8: qualifies — byte-identity conformance corpus role is data-agnostic
- MANIFEST patch: no change (all prior scopes survive the screen)

**telemetry_grid** — prior `qualifies_for`: ['H4', 'H6']
- H4: fails — no low-dimensional predictive structure signal (screen); LB estimator not run here
- H6: qualifies — I(H;L|C)=0.259 bpb, predicted band 1.0-2.0 (float32), falsification floor 0.1
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

1. **Scale**: dev slices (coding <= 256 KiB, indicators <= 2 MiB, front-of-file). BUG-001-D/BUG-010-D name FULL-corpus campaigns with a neural predictor in the amortized regime; this campaign is the local, instrumented version, and the external full-scale items remain open. Front-of-file slices bias heterogeneous tars (calgary: the coding slice covers only the first member(s); indicator slices cover the whole file for corpora <= the indicator budget).
2. **Instrument strength**: the order-W counting predictor is a deliberately weak side-information instrument. A P-ADV failure does NOT falsify the paper's conditional claim (the claim is scoped to decoder-reproducible NEURAL side information with amortized L(M)); it means the advantage is not demonstrated by this instrument at this scale. A P-ADV pass, conversely, is strong evidence the precondition is satisfiable.
3. **Determinism**: coder and indicators are deterministic; coding-rate numbers are single runs (the METHODS 5-seed repetition governs TIMING hypotheses, none reported here). All bootstrap seeds derive from the master seed per METHODS section 4 and are recorded in the JSONL.
4. **Secondary status**: the precondition screen is not one of the frozen H1–H18 primary tests; all p-values are labeled secondary, with Holm-Bonferroni applied within this screening family (9 corpora) for transparency only.
5. **Plug-in bias**: per-third ladder values are plug-in estimates (biased low at high orders/sparse contexts); qualification uses held-out cross-entropy where bias direction is safe, and per-third plug-in values only through the RELATIVE spread.
6. **enwik9**: deferred payload, not fetched, not measured; the H15 scope is untouched.
7. **Store-mode corpora**: for incompressible controls the sync-overhead datum reduces to index bytes (no cold-context term); reported as such, not blended with coded corpora.


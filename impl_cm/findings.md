# H1 first measurement — context-mixing instrument (check scale)

Instrument: `cm_predictor.py`, an integer-deterministic context-mixing predictor
(logistic mixing of order-0..W byte context models in fixed-point arithmetic,
spec R-3.1–R-3.8 conformant: integer-only inference, table-based stretch/squash,
fixed reduction order, overflow guards) behind the same `Predictor` interface as
`impl/rnr1.py`, coded through the unmodified rnr1 container (K = 64 KiB, W = 3).
This is the local hybrid instrument of exp-design §2.1; the 8–16 MB transformer
predictor remains the external upgrade — the H1 verdict below is scoped
accordingly.

## Measured rates (bits/byte, 1 MiB dev slices, seeds/shas in results/h1_check.jsonl)

| corpus | rnr1-cm | rnr1-ngram | gzip | zstd-22 | bz2 | lzma-9e | brotli-11 |
|---|---|---|---|---|---|---|---|
| enwik8_dev | **2.907** | 3.775 | 3.166 | 3.021 | 2.866 | 2.930 | 2.547 |
| calgary_text | **2.762** | 3.600 | 2.871 | 2.748 | 2.701 | 2.686 | 2.426 |
| canterbury_text | **2.888** | 3.830 | 3.144 | 3.012 | 2.912 | 2.955 | 2.680 |
| scripts_src | **2.589** | 3.319 | 2.514 | 2.419 | 2.442 | 2.384 | 2.238 |

## Gap decomposition (Theorem 4.2 chain, per corpus)

- **Model cross-entropy** dominates: CM CE = 2.899 / 2.755 / 2.879 / 2.579 bpb —
  the achieved rate tracks the model's own CE to within ~0.01 bpb everywhere.
- **Coder redundancy** (measured − ideal assignment cost): ≤ 7.5e-5 bpb —
  the arithmetic-coding chain is essentially lossless against its model, as
  Theorem 4.2(A) requires.
- **Container + sync overhead**: ~0.007 bpb at K = 64 KiB on these slices.

## Verdict (H1, quantifying — instrument-scoped)

- The CM instrument improves on the plain order-W predictor by **0.73–0.94 bpb**
  on every text/code corpus (gate G4, strict).
- Against baselines: rnr1-cm beats gzip and zstd on the natural-text corpora and
  bz2/lzma on two of four; brotli-11 leads everywhere (0.35–0.36 bpb ahead).
  Since the coding chain adds < 1e-4 bpb over model CE, the remaining gap to
  the strongest baselines is **entirely predictor quality** — exactly the term
  the exp-design's neural instrument targets. H1's full claim (RNR with a
  matched neural predictor at or below context-mixing baselines) is therefore
  neither confirmed nor falsified at this instrument grade; the measurement
  isolates where the remaining bits are.

## Gate

`run_checks.py`: G1 table sanity, G2 round-trip bit-exactness (+ random access,
+ store mode), G3 two-process determinism, G4 CM strictly beats n-gram CE on
every corpus, G5 coder redundancy < 0.01 bpb — **5/5 PASS** (172.7 s).

Caveats: check-scale slices (1 MiB, front-of-file); single-seed rate rows
(deterministic coders; timing rows are single-run and marked as deviations in
the JSONL per METHODS); pure-Python throughput (~12 s/MiB) makes wall-clock
comparisons against C baselines meaningless — rate columns only.

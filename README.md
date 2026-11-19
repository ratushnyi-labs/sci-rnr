# Ratushnyi Neural Repair Coding

Research manuscripts and reference implementations for Ratushnyi Neural Repair Coding (RNR). The series studies source coding with shared or reproducible predictors, the cost of transmitting repair information, and the trade-offs introduced by random access and finite blocks.

## Papers

| Part | Subject | LaTeX source |
| --- | --- | --- |
| I | Framework, constructions and complexity | [rnr_core.tex](tex/papers/core/rnr_core.tex) |
| II | Random access, finite-block surprisal and the deviation hierarchy | [rnr_random_access.tex](tex/papers/random_access/rnr_random_access.tex) |
| III | Operational rate–distortion dispersion | [rnr_dispersion.tex](tex/papers/dispersion/rnr_dispersion.tex) |

The [combined manuscript](tex/rnr_coding.tex) is the source of truth for the split. See the [series guide](tex/papers/README.md) for the numbering conventions. Build all three PDFs from the repository root with Make and Tectonic.

```sh
make papers
```

## Reference coder

The Python Type-I coder defines the [RNR1 archive format](impl/FORMAT.md). The [C implementation](impl_fast/README.md) and [Rust CLI](rnr-cli/README.md) use the same format. The [two-pass implementation](impl_2pass/README.md) and [hashed context mixer](impl_hashed/) are experimental.

The reference coder requires Python 3 and NumPy. This example packs the format specification and checks that unpacking reproduces it byte for byte.

```sh
mkdir -p build
python3 impl/rnr1.py pack impl/FORMAT.md build/demo.rnr --W 3 --K 65536
python3 impl/rnr1.py info build/demo.rnr
python3 impl/rnr1.py unpack build/demo.rnr build/restored-format.md
cmp impl/FORMAT.md build/restored-format.md
```

## Reproduction

- `python3 scripts/verify/crossref_integrity.py` checks numbered references in the manuscripts and explicitly marked links to companion papers.
- `python3 impl/run_checks.py` exercises reference-coder round-trip, determinism and random access.
- The [benchmark methods](bench/METHODS.md) define the timing and statistical protocol. The [corpus manifest guide](data/README.md) gives the dataset sources, digests and fetch instructions.

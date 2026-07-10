# Three-paper split of `rnr_coding.tex`

`../rnr_coding.tex` (the monolith) remains the source of truth. These are the three
**self-contained, independently compiling** papers of the unified series, split per
`docs/backlog/SPLIT-PLAN.md`. Each carries the shared preamble, a fresh title/abstract
scoped to its own content, a prerequisites section restating (statement-only, with
companion provenance) the results it depends on, companion-citation repair for
cross-paper references, and its own reference list. Original theorem/section numbering
is preserved across the series (a note in each paper explains this); companion results
are cited `[RNR-I]` / `[RNR-II]` / `[RNR-III]`.

| Folder | File | Part | Body | Pages |
|---|---|---|---|---|
| `core/` | `rnr_core.tex` | I — Framework, Constructions, Complexity | §1–6 taxonomy + five-type matrix + constructions + hardness, §7 lead, §8–13 | 201 |
| `random_access/` | `rnr_random_access.tex` | II — Random Access & the Deviation Hierarchy | §7.2–7.33, §7.35–7.36, §10.7–10.8 | 132 |
| `dispersion/` | `rnr_dispersion.tex` | III — Operational RD Dispersion | the §7.34 family incl. the Route B replica spectral theorem | 78 |

## Correctness verified

- **Compiles clean**: each builds with `tectonic` to a valid PDF — exit 0, zero
  undefined references, zero hard LaTeX errors (Overfull/Underfull hbox warnings are
  inherited verbatim from the monolith's dense displays).
- **Dependency policy**: Part I is fundamental except for the §§10.7–10.8
  random-access blocks, which it restates statement-only in a dedicated
  "P. Prerequisites (restated from Part II)" section; Part II depends only on
  Part I; Part III on Parts I/II the same way. Every cross-Part mention carries
  an [RNR-x] tag within the gate's window: `scripts/verify/crossref_integrity.py`
  gates each paper independently — Part I and Part III currently PASS; Part II
  has 13 residual bare references (audit findings 11–13, its prerequisite
  restatements of Theorem 10.1/10.2 and the Lemma 5.6b–f cluster are the
  remaining repair, tracked in `docs/backlog/AUDIT-SPLIT-2026-07-07.md`).
- **Label containment**: the corpus' only `\label` (`eq:G1chain`) lives in §7.34, so
  it is confined to Part III — no cross-file reference resolution needed.
- **Stealth-clean**: no tool/model/assistant names in any file.

Build products (`*.pdf`, `*.aux`, `*.log`, …) are git-ignored (`.gitignore`).
Rebuild any paper with `cd <folder> && tectonic -X compile <file>.tex`.

## Polish applied

- Reference lists pruned to drop entries with no author or title footprint in the body,
  then renumbered contiguously: Part I 135→134, Part II 135→132, Part III 160→137.
  The prune is deliberately conservative (a false drop of a cited work is worse than an
  untidy list): it removes only clearly off-topic entries (grammar-compression,
  bits-back neural coders, PCP-hardness, crypto, graphical-model works) that a
  dispersion / random-access paper never cites, and keeps every entry whose author is
  named in the text — including classics cited without a year (Shannon, Blahut, Marton).
- All three recompile clean after the prune.

## Remaining polish (optional editorial judgment)

- A *lean* bibliography for Part III (down to the ~40 works it actually cites) would
  require per-entry review of ~100 borderline master-list entries whose common author
  surnames (Chen, Gray, White, …) appear coincidentally — separating those from genuine
  classic-work citations needs reading, not automation, so it is left as a judgment call.
- Part III's ~25 added references render a few venues best-known without invented page
  numbers — verify against the primary sources.
- `scripts/verify/` partitions ~1:1 by filename; retarget the document-invariant probes
  (`bug_001`, `bug_010`, `crossref_integrity`) to the three files (`crossref_integrity`
  retargeted 2026-07-07; `bug_001`/`bug_010` remain).
- The companion docs (`rnr_summary`, `rnr_engineering_spec`, `rnr_experimental_design`)
  reference "the main paper" — sweep their pointers to the correct Part.

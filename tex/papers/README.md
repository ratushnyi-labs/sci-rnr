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
- **Self-contained**: cross-paper theorem references are repaired to companion
  citations; each dependent result is restated (statement-only) in the citing paper's
  prerequisites section.
- **Label containment**: the corpus' only `\label` (`eq:G1chain`) lives in §7.34, so
  it is confined to Part III — no cross-file reference resolution needed.
- **Stealth-clean**: no tool/model/assistant names in any file.

Build products (`*.pdf`, `*.aux`, `*.log`, …) are git-ignored (`.gitignore`).
Rebuild any paper with `cd <folder> && tectonic -X compile <file>.tex`.

## Remaining polish (optional, non-blocking)

- Part III's reference list adds ~25 information-theory / probability works cited by
  name in §7.34 that were absent from the master list; a few venues are rendered
  best-known without invented page numbers — verify against the primary sources.
- The full master reference list is carried in each paper (over-inclusion is harmless);
  prune to each paper's cited subset if desired.
- `scripts/verify/` partitions ~1:1 by filename; retarget the document-invariant probes
  (`bug_001`, `bug_010`, `bug_003`, `crossref_integrity`) to the three files.
- The companion docs (`rnr_summary`, `rnr_engineering_spec`, `rnr_experimental_design`)
  reference "the main paper" — sweep their pointers to the correct Part.

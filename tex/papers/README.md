# Three-paper split of `rnr_coding.tex`

`../rnr_coding.tex` (the monolith) remains the source of truth. These three
papers compile independently and follow the split in
`docs/backlog/SPLIT-PLAN.md`. Each has its own title, abstract, prerequisites
and reference list. Prerequisites include statement-only restatements of
selected results. The papers preserve the original theorem and section
numbering. Companion results are cited as `[RNR-I]`, `[RNR-II]` or `[RNR-III]`.

| Folder | File | Part | Body |
|---|---|---|---|
| `core/` | `rnr_core.tex` | I — Framework, Constructions, Complexity | §1–6 taxonomy + five-type matrix + constructions + hardness, §7 lead, §8–13 |
| `random_access/` | `rnr_random_access.tex` | II — Random Access & the Deviation Hierarchy | §7.2–7.33, §7.35–7.36, §10.7–10.8 |
| `dispersion/` | `rnr_dispersion.tex` | III — Operational RD Dispersion | the §7.34 family incl. the Route B replica spectral theorem |

## Correctness verified

- **Build.** Each paper builds with `tectonic`. The exported layout
  still produces overfull and underfull boxes, especially in dense displays.
- **Cross-paper references.** Part I restates selected random-access
  results from Part II. Part II points to Part I for the theory it uses,
  including Lemmas 5.6b–f and Theorems 10.1 and 10.2. Part III
  draws on both companion papers. Run `python3 scripts/verify/crossref_integrity.py`
  from the root to check numbered references and nearby companion tags.
  The [split audit](../../docs/backlog/AUDIT-SPLIT-2026-07-07.md) records the
  broader editorial and proof-dependency review.
- **Label containment**: the corpus' only `\label` (`eq:G1chain`) lives in §7.34, so
  it is confined to Part III — no cross-file reference resolution needed.

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
- Retarget the remaining monolith-specific document probes (`bug_001` and
  `bug_010`) to the split papers.
- The companion docs (`rnr_summary`, `rnr_engineering_spec`, `rnr_experimental_design`)
  reference "the main paper" — sweep their pointers to the correct Part.

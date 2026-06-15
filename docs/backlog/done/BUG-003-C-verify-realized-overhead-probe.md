# BUG-003-C - verify-realized-overhead-probe

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** scripts/verify/ (new probe); tex/rnr_coding.tex (§10.7, read-only anchor)
**Business Rule Conflict Check:** NOT_CHECKED
**Conflict Resolution Reference:** N/A
**Atomic Check:** ATOMIC
**Budget Forecast:** TBD
**Priority:** MEDIUM
**Lawbook Version (intake):** 0.44.0
**Applicable Lawbook Version:** TBD
**Created:** 2026-06-15  **Updated:** 2026-06-15
**Recovered:** N/A

## Raw Request
> mark each gap as bug and follow the global rules
>
> (Filed gap 003 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

(Child sub-task created at section 12.7 decomposition of BUG-003 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: §10.7 (the "approximately one percent" heuristic conjecture, lines ~15966--15970).

The "~1% overhead at typical parameters" figure rests on "empirical observations of cold-vs-warm
rates for small distilled transformers" with no executable artifact in the repo backing the
ratio claim. The project verification rule ("Практична верифікація") requires every
lemma/conjecture-support to have an executable check in scripts/verify/. There IS an existing
scripts/verify/lemma_5_1c_tight_cold_context_markov.py, but it verifies the EXACT identity
`δ_W' = H(X^W') − W'·h(X)` -- it does NOT compute the *ratio* of total sync overhead to total
warm-context length at typical (K, W) parameters, which is what the "~1%" figure actually is.
So the conjecture's order-of-magnitude (≈1%) has no numerical support artifact distinct from
the tight-identity check.

**FIX TARGET (this leaf):** add a verification probe (e.g.
scripts/verify/syncpoint_overhead_ratio_probe.py) that, on small synthetic sources
(stationary order-W' Markov chains plus at least one non-Markov / longer-range source as a
stress case), computes for a sweep of (K, W, η) the realized total sync overhead
`Σ δ_realized` (or the floor-bounded `m·S·W·log₂(1/η)`) as a FRACTION of the total
warm-context coded length `K·h(X)·m`, and reports the ratio in the ~1% regime where Markov
order-W' is matched (K ≫ W'). The probe must PRINT PASS/FAIL: PASS = (i) realized ratio
sits below the loose `S·W·log₂(1/η)`-per-sync ceiling for every swept point, and (ii) at
typical K = 64 KiB-analogue (K/W large) the realized ratio falls in the conjectured
low-single-digit-percent band for the Markov sources, FAIL otherwise. This is HEURISTIC
SUPPORT for the conjecture on synthetic data -- it is explicitly NOT a measurement on a named
real corpus (that is BUG-003-D, closeable_in_repo=NO).

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class repo; the probe is the project-mandated executable check for
  a conjecture-support claim. Closeable in-repo NOW (pure synthetic-data Python; no external
  data, no reference implementation).
- Run locally with /Users/para/.venvs/rnr/bin/python (per project memory; CI runs the
  scripts/verify job). gh cannot observe this repo's CI.
- The probe supports the order-of-magnitude conjecture on synthetic sources; it does NOT
  discharge the empirical-on-named-corpus residual (BUG-003-D). Keep the two strictly
  separate so the conjecture is not overclaimed as "measured."
- Distinct from lemma_5_1c probe: that one checks the exact tight identity; this one checks
  the OVERHEAD-AS-FRACTION (the ~1%) across a (K, W, η) sweep including a non-Markov stress case.

## Refined Description
**Scope:** one new scripts/verify/ probe computing realized-overhead-as-fraction over a
(K, W, η) sweep on synthetic sources, with PASS/FAIL on (loose-ceiling-respected) +
(typical-K ratio in conjectured band).  **Non-goals:** any real-corpus measurement, any
reference codec, any change to existing lemma_5_1c probe.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-003; conceptually anchored by BUG-003-B's accounting ladder
(the probe verifies the realized number sits inside that ladder) -- can proceed in parallel
but the PASS criteria should reference B's three ceilings.
**Risks / Open Questions:** must not overstate -- synthetic Markov sources are favorable to
the conjecture; include a non-Markov stress case so the probe can FAIL honestly and the
result is reported as "supports ~1% on Markov synthetic; does not establish real-corpus value."
**Expected Code Change:** new Python file in scripts/verify/ (+ CI job inclusion).
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | executable heuristic support for ~1% on synthetic; not a real-corpus measurement |
| ESTIMATION | DONE | ABOVE_EASY (synthetic Markov entropy bookkeeping + sweep + PASS/FAIL) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | source set (Markov orders + 1 non-Markov), (K, W, η) sweep, PASS thresholds |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | implement probe; print PASS/FAIL |
| TESTING | PENDING | run via /Users/para/.venvs/rnr/bin/python; add to CI scripts/verify job |
| POST_AUDIT | PENDING | adversarial re-review: probe cannot trivially pass; non-Markov stress case can fail |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh as a decomposition child of BUG-003; no probe
  is written before this file.
- Class: GOVERNANCE/DOCUMENTATION; runtime product stages N/A, but scripts/verify executable
  checks ARE in-scope per the project verification rule.
- Closeable in-repo: YES (synthetic-data probe). The named-corpus measurement is BUG-003-D.


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-003 A/B/C (commit pending). A: Abstract '~1%' scoped as an explicit heuristic CONJECTURE (rests on the warm-rate assumption of Cor 10.2a, not the bare floor) + crosslink to §10.7 ladder + to H10(c) measurement. B: Corollary 10.2a 'assumption-graded per-sync overhead ladder' + 3-row table consolidating the three PROVED ceilings -- rung1 loose S*W*log2(1/eta) [floor only; Lemma 5.1a], rung2 warm-rate-refined [floor + uniform warm LOWER bound], rung3 exact-Markov delta_W' [stationary order-W' Markov; Lemma 5.1c eqn 5.9]; states the fully-general tight bound REMAINS OPEN. C: scripts/verify/bug_003_syncpoint_overhead_fraction.py (realized overhead FRACTION = cold-restart excess / no-sync warm baseline, matching H10(c); rung-1 ceiling respected at 360 points; natural sources in ~1% band; ~1/K scaling; counterweights legitimately exceed, probe can fail honestly) + CI job. Hostile-referee POST_AUDIT clean (3 rungs verified against sources, no hypothesis flipped; probe non-vacuous via injected broken source).

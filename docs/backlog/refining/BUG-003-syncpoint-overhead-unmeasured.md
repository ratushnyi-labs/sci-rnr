# BUG-003 - syncpoint-overhead-unmeasured

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DECOMPOSED
**Complexity:** HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§10.7; Abstract)
**Business Rule Conflict Check:** NOT_CHECKED
**Conflict Resolution Reference:** N/A
**Atomic Check:** ATOMIC
**Budget Forecast:** TBD
**Priority:** MEDIUM
**Lawbook Version (intake):** 0.44.0
**Applicable Lawbook Version:** TBD
**Created:** 2026-06-14  **Updated:** 2026-06-15
**Recovered:** YES

## Raw Request
> mark each gap as bug and follow the global rules
>
> (Filed gap 003 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

## Gap (defect statement; document-grounded)
Source: §10.7; Abstract.

The random-access sync-point overhead is conjectured '~1% for natural data at typical parameters' but 'neither a tight upper bound in general nor a measured value on a specific corpus is provided.' The only proven bound is S*W*log2(1/eta) under a probability-floor assumption. FIX TARGET: a tight general upper bound, OR a measured value on a named corpus, OR scope the 1% as a conjecture.

## Clarification / Assumptions
- This is a GOVERNANCE/DOCUMENTATION-class defect in the theory paper
  `tex/rnr_coding.tex` (no runtime product). "Fix" means: close the proof/empirical
  gap, OR scope the paper's claim to match what is actually established.
- Severity reflects how load-bearing the gap is for the document's MAIN IDEA
  (side-information+repair: neural ratio AND polylog random access).

## Refined Description
**Scope:** the gap above and its FIX TARGET.  **Non-goals:** unrelated §-content.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** none
**Risks / Open Questions:** see Gap; some fixes need work outside this repo (experiments).
**Expected Code Change:** UNKNOWN
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | scope the claim vs close the gap |
| ESTIMATION | PENDING | set Complexity/Budget |
| DECOMPOSITION | DONE | split into BUG-003-A, BUG-003-B, BUG-003-C, BUG-003-D |
| DESIGN | PENDING | proof strategy or experiment design |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex edit / proof / probe (or N/A if empirical-only) |
| TESTING | PENDING | scripts/verify probe or measured datum |
| POST_AUDIT | PENDING | adversarial re-review per project review protocol |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): the gap analysis was performed in the immediately
  preceding turn (user request "list the available gaps"), BEFORE this backlog file
  existed; the CodexOfLaws §12 format was read this turn to file compliantly. Flagged
  `Recovered: YES` per §12.1.3; no code/proof change was made to the paper before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.

## Decomposition (section 12.7)
Split on 2026-06-15 into four blocking child leaves along the doc-scoping / theory-bound /
verification-script / external-measurement boundary. The three in-repo leaves (A, B, C) form
the shippable core; the genuinely-external residual (D) is isolated and clearly flagged so the
parent need not wait on a reference implementation. Each child is ATOMIC.

- **BUG-003-A** (scope-abstract-conjecture) — Make the Abstract's "~1%" conjecture sentence
  self-consistent and traceable: explicit conjecture label + cross-link to the §10.7 warm-rate
  assumption it rests on + cross-link to the §11 / experimental-design H10 measurement.
  Complexity EASY. Closeable in-repo: YES.
- **BUG-003-B** (tighten-bound-accounting) — Consolidate the three already-proved per-sync
  ceilings (loose floor-only `S·W·log₂(1/η)`; warm-rate-refined; exact-Markov `δ_W'`) into one
  assumption-graded ladder (corollary/remark + 3-row table) in §10.7; the fully-general tight
  bound stays open. Complexity ABOVE_EASY. Closeable in-repo: YES.
- **BUG-003-C** (verify-realized-overhead-probe) — New scripts/verify/ probe computing realized
  sync overhead AS A FRACTION over a (K, W, η) sweep on synthetic Markov + a non-Markov stress
  source, PASS/FAIL on loose-ceiling-respected and typical-K ratio in the conjectured band;
  heuristic support on synthetic data, NOT a real-corpus measurement. Complexity ABOVE_EASY.
  Closeable in-repo: YES.
- **BUG-003-D** (measure-named-corpus) — The genuinely-external residual: a real measured
  sync-point ratio-overhead datum on a NAMED corpus at typical (K, W) per experimental-design
  H10(c), requiring a reference RNR implementation and a compute run outside this repo.
  Complexity HARD. Closeable in-repo: NO (external work; in-repo follow-up is only a one-line
  claim update once a datum exists).

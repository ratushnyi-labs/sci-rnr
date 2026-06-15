# BUG-007 - bitexact-reproduction-spec-only

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DECOMPOSED
**Complexity:** HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Thm 10.1; §10.1)
**Business Rule Conflict Check:** NOT_CHECKED
**Conflict Resolution Reference:** N/A
**Atomic Check:** ATOMIC
**Budget Forecast:** TBD
**Priority:** HIGH
**Lawbook Version (intake):** 0.44.0
**Applicable Lawbook Version:** TBD
**Created:** 2026-06-14  **Updated:** 2026-06-15
**Recovered:** YES

## Raw Request
> mark each gap as bug and follow the global rules
>
> (Filed gap 007 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

## Gap (defect statement; document-grounded)
Source: Thm 10.1; §10.1.

The whole side-information recasting REQUIRES the decoder to reproduce the model output bit-for-bit across platforms; this is reduced to a precision/determinism SPECIFICATION (Thm 10.1) but not closed -- 'an implementer must verify against the full specification on a chosen target.' FIX TARGET: a reference deterministic-inference implementation + a cross-platform bit-exactness test demonstrating the spec is satisfiable end-to-end.

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
**Expected Code Change:** YES
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | scope the claim vs close the gap |
| ESTIMATION | PENDING | set Complexity/Budget |
| DECOMPOSITION | DONE | split into BUG-007-A, BUG-007-B, BUG-007-C, BUG-007-D |
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

## Decomposition (§12.7)
Split on 2026-06-15 into four BLOCKING child leaves, separating the part that can be
honestly closed IN-REPO now (claim scoping, math verification probe, cross-doc accounting)
from the genuinely EXTERNAL residual (reference implementation + real cross-platform test).
The parent's HARD complexity drops to a shippable in-repo core (A + B + C, all EASY /
ABOVE_EASY) plus one clearly-flagged external leaf (D, closeable_in_repo = NO).

| Child ID | Objective (one line) | Complexity | Closeable in-repo |
|---|---|---|---|
| BUG-007-A | Audit/tighten §10.1 / §10.5 / Thm 10.1 prose so it claims exactly the proven SUFFICIENT CONDITION and explicitly defers the reference impl (§13.3); no overclaim of end-to-end demonstration. | EASY | YES |
| BUG-007-B | Add `scripts/verify/thm_10_1_integer_reorder_bitexact.py` demonstrating Thm 10.1's math core: per-layer bit-width no-overflow bound (with undersized negative control) + integer reduction-order invariance ⇒ bit-identity. | ABOVE_EASY | YES |
| BUG-007-C | Make the §13.3 deferral precise: cross-link Thm 10.1 / §10.5 / Thm 10.6 to the engineering-spec conformance suite (§3 R-3.x, §12) and to any bit-exactness item in the experimental-design doc, so "satisfiable end-to-end" is scoped as conditional-on-certified-conformance with a named artifact. | EASY | YES |
| BUG-007-D | EXTERNAL: build the reference deterministic-integer inference codec (engineering-spec §3) and run the §12 conformance suite across real hardware targets to demonstrate cross-platform byte-identity end-to-end. | ABOVE_EASY | NO |

Recommended first leaf: BUG-007-A (scope the claim honestly; quickest, unblocks the
cross-link wording in C). BUG-007-B is independent and can proceed in parallel.

# BUG-006 - type3c-builder-universality-open

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Rem 6.6b; Thm 6.6b/6.6c)
**Business Rule Conflict Check:** NOT_CHECKED
**Conflict Resolution Reference:** N/A
**Atomic Check:** ATOMIC
**Budget Forecast:** TBD
**Priority:** MEDIUM
**Lawbook Version (intake):** 0.44.0
**Applicable Lawbook Version:** TBD
**Created:** 2026-06-14  **Updated:** 2026-06-14
**Recovered:** YES

## Raw Request
> mark each gap as bug and follow the global rules
>
> (Filed gap 006 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

## Gap (defect statement; document-grounded)
Source: Rem 6.6b; Thm 6.6b/6.6c.

Type-III-C is 'asymptotically rate-optimal' (Nh+o(N), converse-matched) only UNDER a universal dictionary builder; whether the concrete greedy marginal-MDL builder (Thm 6.5) is itself a universal transform is left OPEN (Rem 6.6b). So the headline III-C optimality rests on an unproven property of the actual builder. FIX TARGET: prove the greedy-MDL builder is universal, or replace it with a proven-universal one (LZ78/irreducible-grammar) as the default.

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
| CLARIFICATION | DONE | partial-close (augmented universality) + honest residual localization, not bare-rule resolution |
| ESTIMATION | DONE | HARD; single proof-lemma + probe |
| DECOMPOSITION | DONE | one lemma (6.6d) + one scoping remark (6.6e); no empirical split needed |
| DESIGN | DONE | guard-against-universal-baseline: features (A) entropy-coded literal fallback + (B) min{learned,empty} baseline guard |
| FRONTEND | N/A | no UI |
| BACKEND | DONE | Lemma 6.6d + Remark 6.6e in tex/rnr_coding.tex (~5720-5876); Rem 6.6b + §13.2 updated; bib [134] Ryabko 1984, [135] Willems et al. 1995 |
| TESTING | DONE | scripts/verify/lemma_6_6d_mdl_entropy_fallback_universal.py (D1-D4 PASS, cold-run); CI job added to build.yml |
| POST_AUDIT | DONE | hostile-referee re-review PASSED ("sound and honestly scoped"): proof upper/converse/a.s. branches non-circular; no overclaim (abstract left saying open; §13.2 bare-rule (i)-(iii) open); fixed 2 dangling citations |

## Resolution
Partial-close per FIX TARGET's "OR scope the claim" branch, upgraded to a positive
universality result for the augmented builder. **Lemma 6.6d:** the greedy marginal-MDL
builder of Thm 6.5 is universal (rate -> h(X), a.s. for ergodic) for EVERY stationary
source once its entropy stage satisfies (A) entropy-coded literal fallback under a
universal predictor Q and (B) a one-bit baseline guard min{learned dict, empty dict} ---
because the empty dict already entropy-codes at Nh+o(N), the guard forbids worse, the
converse (Thm 6.6c) forbids better, so the admission decisions are rate-irrelevant.
**Remark 6.6e** localizes the genuinely-open residual to the *bare* flat-8*ell literal
rule of Thm 6.5 (where the empty-dict baseline is 8N >> Nh, so the guard is useless and
the burden falls entirely on the admissions = the smallest-grammar-adjacent question);
§13.2 sub-questions (i)-(iii) for the bare rule remain open. The MAIN-IDEA load is
discharged: any implementation wanting a provable rate guarantee adds (A)+(B) and may use
the greedy dictionary freely for its structural (random-access) benefits with no rate penalty.

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** drafted 2026-06-14 (prior session, uncommitted); recovered, POST_AUDIT'd, and closed 2026-06-15 (model: claude-opus-4-8; hostile-referee subagent for POST_AUDIT)

## Compliance Notes
- §12.1.1 (Backlog-File-First): the gap analysis was performed in the immediately
  preceding turn (user request "list the available gaps"), BEFORE this backlog file
  existed; the CodexOfLaws §12 format was read this turn to file compliantly. Flagged
  `Recovered: YES` per §12.1.3; no code/proof change was made to the paper before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.

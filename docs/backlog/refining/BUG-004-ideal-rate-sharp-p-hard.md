# BUG-004 - ideal-rate-sharp-p-hard

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE (fully closed; no external leaf)
**Complexity:** VERY_HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Rem 7.15g; §7.15 open middle)
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
> (Filed gap 004 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

## Gap (defect statement; document-grounded)
Source: Rem 7.15g; §7.15 open middle.

The exact ideal RNR rate L*_RNR = H(X^N) is #P-hard / FP^#P-complete, and the exact finite-block entropy oracle is #P-hard (only eps-approximable under filter stability, 7.15c). So the optimum the framework targets is not efficiently computable in general. FIX TARGET (mostly characterization): make explicit, wherever the 'ideal rate' is invoked, that it is an intractable target and the operative object is the eps-approximation / achievable scheme.

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
**Expected Code Change:** NO
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | scope the claim vs close the gap |
| ESTIMATION | PENDING | set Complexity/Budget |
| DECOMPOSITION | DONE | split into BUG-004-A, BUG-004-B, BUG-004-C |
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
Split on 2026-06-15 into three blocking child leaves. The underlying hardness theorem (Rem
7.15g: exact ideal rate `L*_RNR = H(X^N)` is #P-hard / FP^#P-complete) is ALREADY proven and
committed; the §7.15g Corollary (~line 9070) already states the operational reading locally.
The parent's residual VERY_HARD framing was an over-estimate: the remaining work is a
documentation-scoping audit + surgical cross-link edits + a small verification probe, all
honestly closeable in-repo. No genuinely-hard research residual remains.

| Child ID | Objective (one line) | Complexity | Closeable in-repo |
|---|---|---|---|
| BUG-004-A-audit-ideal-rate-invocations | Inventory every `L*_RNR = H(X^N)` / entropy-floor invocation across tex/rnr_coding.tex and verify which carry / imply / lack the intractability caveat (audit only, no edits). | EASY | YES |
| BUG-004-B-scope-ideal-rate-as-intractable | At each flagged site, add a minimal truthful cross-link: exact ideal is #P-hard (Rem 7.15g), operative object is the eps-approx / achievable scheme (Rem 7.15b/c); do not weaken the converse. | ABOVE_EASY | YES |
| BUG-004-C-verify-ideal-vs-achievable-gap-probe | Small scripts/verify/ probe on a tiny finite-state source: poly one-pass realized rate vs exact ideal `H(X^N)` (the #P-hard target) + the converse inequality; verification anchor for B. | EASY | YES |

Sequencing: A -> B (B consumes A's inventory); C can run in parallel and pairs with B as its
green-CI anchor. Recommended first leaf: BUG-004-A.


## Resolution (2026-06-15)
FULLY CLOSED in-repo (A/B/C -> done/; no D-leaf). The exact ideal RNR rate H(X^N) is now
scoped as an intractable (#P-hard/FP^#P-complete, Rem 7.15g) BENCHMARK -- not an efficiently
computable encoder target -- at every site that presented it as the operative target (Thm
6.6c.1, §7.1b, §7.27'), with the operative object being the eps-approximation tractable under
filter stability (Rem 7.15c) / the achievable scheme; converse inequalities preserved.
Gap probe bug_004_ideal_vs_achievable_gap.py (intractable exact vs poly eps-approx vs
achievable-realizes-rate; OVERALL -> PASS) + CI job. POST_AUDIT (hostile referee) clean;
both probes (bug_004 + the cited remark_7_15g) pass cold. Nothing parked -- bug complete.
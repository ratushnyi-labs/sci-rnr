# BUG-001 - empirical-ratio-unmeasured

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE (in-repo); D external-blocked
**Complexity:** UNESTIMATABLE
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Abstract; §1.1; §11)
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
> (Filed gap 001 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

## Gap (defect statement; document-grounded)
Source: Abstract; §1.1; §11.

The central value proposition -- that RNR's side-information+repair recasting beats existing seekable formats on ratio -- is EMPIRICAL and unmeasured. The paper states it makes "no quantitative claim about absolute compression ratios on any specific corpus" and that the ratio comparison "is an empirical conjecture (§11), not a theorem." No corpus result confirms the main idea. FIX TARGET: run the §11 / experimental-design protocol (or scope the claim to 'conditional, pending measurement').

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
| DECOMPOSITION | DONE | split into BUG-001-A, BUG-001-B, BUG-001-C, BUG-001-D |
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
Split into four BLOCKING child leaves. The §12.7 split isolates the part that can be
HONESTLY CLOSED NOW IN-REPO (scope the claim, cross-link the pre-registered protocol,
guard with a verify probe) from the IRREDUCIBLE external measurement. After A/B/C land,
the parent's in-repo posture is a precisely-scoped conditional claim traceable to its
pre-registered test; D is the clearly-flagged deferred residual needing out-of-repo work.

- **BUG-001-A — scope-claim-conditional**: make the absolute-ratio / beats-seekable-format
  claim uniformly "conditional, pending measurement" across Abstract, §1.1, §11 intro,
  §§11.1--11.6 (one canonical phrasing; excludes §11.7). Complexity EASY; closeable-in-repo YES.
- **BUG-001-B — crosslink-preregistered-protocol**: verify + add the §11-prediction <->
  companion-hypothesis (H1--H7, H15) mapping and cross-references so the claim is traceable
  to the exact pre-registered experiment, including the seekable-format baseline coverage.
  Complexity ABOVE_EASY; closeable-in-repo YES.
- **BUG-001-C — conditional-status-verify-probe**: add a `scripts/verify/` PASS/FAIL probe
  (+ CI job) asserting the document-invariants from A/B (no unconditional ratio claim; no
  orphan in the §11<->companion mapping). Complexity ABOVE_EASY; closeable-in-repo YES.
- **BUG-001-D — run-empirical-ratio-measurement**: execute the pre-registered protocol
  (reference implementation + baselines incl. bgzip/zstd-seekable + corpora + compute) and
  report measured bits-per-byte, then reconcile §11. Complexity VERY_HARD;
  closeable-in-repo NO (external work; isolated research residual).


## Resolution (2026-06-15)
A/B/C closed in-repo (-> done/): §11.1-11.6 uniformly scoped conditional/pre-registered
+ §11-intro canonical paragraph (Abstract/§1.1 were already scoped); §11.x<->H1-H7/H15-H16
mapping with seekable-baseline coverage note; document-invariant probe
bug_001_ratio_claim_conditional.py (inv1 conditional-scoping + inv2 no-orphan, self-tested
non-vacuous) + CI job. POST_AUDIT (hostile referee) clean; probe hardened (per-paragraph
locality, dominate-needs-baseline). HONEST RESIDUAL feeding D: §11.5 repetitive-ratio and
the beats-seekable value-prop have only partial pre-registered coverage (no dedicated
bgzip/zstd-seekable RATIO hypothesis); flagged in-text, not invented. D (run the measurement)
-> blocked/ (external).
# BUG-001-A - scope-claim-conditional

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Abstract; §1.1; §11 intro + §§11.1--11.6)
**Business Rule Conflict Check:** NOT_CHECKED
**Conflict Resolution Reference:** N/A
**Atomic Check:** ATOMIC
**Budget Forecast:** TBD
**Priority:** HIGH
**Lawbook Version (intake):** 0.44.0
**Applicable Lawbook Version:** TBD
**Created:** 2026-06-15  **Updated:** 2026-06-15
**Recovered:** N/A

## Raw Request
> mark each gap as bug and follow the global rules
>
> (Filed gap 001 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at §12.7 decomposition of BUG-001 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: Abstract; §1.1; §11 (intro and §§11.1--11.6).

The closeable-in-repo half of BUG-001: make the ratio claim's CONDITIONAL status
uniform and unambiguous across every place it appears. The paper already disclaims
in some spots ("we make no quantitative claim about absolute compression ratios on
any specific corpus"; "is an empirical conjecture (§11), not a theorem"; predictions
"may be off by 10-30%"), but the disclaimers are scattered and the §§11.1--11.5
prediction statements are phrased assertively ("is predicted to achieve ...",
"is predicted to dominate", "outperform any single fixed mode") without a per-claim
"conditional, pending measurement" tag. FIX TARGET (this leaf): edit the Abstract,
§1.1, and the §11 intro plus §§11.1--11.6 so that every absolute-ratio /
beats-existing-seekable-format claim is explicitly scoped as a falsifiable empirical
conjecture pending measurement, with one canonical phrasing reused consistently and a
single forward pointer to the pre-registered protocol (see sibling BUG-001-B). No
new measurement, no proof, no §11.7 change (the §11.7 number is theorem-backed and
is BUG-002's scope, not this leaf's).

## Clarification / Assumptions
- This is a GOVERNANCE/DOCUMENTATION-class defect in the theory paper
  `tex/rnr_coding.tex` (no runtime product). "Fix" means here: scope the paper's
  claim to match what is actually established (the SCOPE-the-claim branch of
  BUG-001's FIX TARGET, not the run-the-protocol branch).
- This leaf is fully closeable in-repo: it is a wording/scoping edit only.
- Non-goal: §11.7 / Theorem 7.15's 0.664 bpb number (theorem-backed; BUG-002).
- Non-goal: model-inclusive accounting (BUG-002), sync-point overhead (BUG-003).

## Refined Description
**Scope:** uniform conditional-scoping wording for the absolute-ratio claim in
Abstract, §1.1, §11 intro, §§11.1--11.6.  **Non-goals:** unrelated §-content; §11.7;
the actual measurement (BUG-001-D).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-001
**Risks / Open Questions:** wording must not silently weaken a claim that is in fact
theorem-backed elsewhere; cross-check against §11.7 before retagging anything.
**Expected Code Change:** YES (tex wording)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | confirm canonical conditional phrasing |
| ESTIMATION | DONE | EASY (localized wording edit, 3 sites) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | choose one reusable "conditional, pending measurement" sentence |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex wording edits in Abstract/§1.1/§11 |
| TESTING | PENDING | manual re-read; covered by BUG-001-C probe |
| POST_AUDIT | PENDING | adversarial re-review per project review protocol |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at §12.7 decomposition of BUG-001;
  no code/proof change made to the paper before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- §12.7: this is one BLOCKING child leaf of BUG-001 (the closeable-now scoping layer).


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-001 A/B/C (commit pending). A: §11.1-11.6 per-prediction conditional scoping ('falsifiable empirical conjecture, conditional and pending measurement, not a theorem') + canonical §11-intro paragraph; Abstract/§1.1 already scoped (left intact). B: §11.x<->companion-hypothesis mapping added (H1-H7, H15-H16) incl. seekable-baseline coverage note (bgzip/zstd-seekable = H10/H11 random-access TIME baselines; squashfs+zstd = H12 ratio baseline). C: scripts/verify/bug_001_ratio_claim_conditional.py (inv1 conditional-scoping + inv2 no-orphan mapping, non-vacuity self-test; per-paragraph locality + dominate-needs-baseline matcher hardened by supervisor; OVERALL -> PASS) + CI job. Hostile-referee POST_AUDIT clean. HONEST RESIDUAL (feeds D): §11.5 repetitive-ratio-vs-zstd/LZMA target + the beats-seekable value-prop have only PARTIAL pre-registered coverage (H7 = image bits-back; no dedicated bgzip/zstd-seekable RATIO hypothesis) -- flagged in-text, no H-number invented.

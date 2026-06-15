# BUG-001-C - conditional-status-verify-probe

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** scripts/verify/ (new probe); tex/rnr_coding.tex (Abstract; §1.1; §11) as read-only fixture
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
> (Filed gap 001 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at §12.7 decomposition of BUG-001 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: Abstract; §1.1; §11 (as the artifact under test) + project rule "Практична
верифікація — кожен крок як робоча програма" (every claim gets an executable check).

The closeable-now scoping work (BUG-001-A) and cross-link work (BUG-001-B) need a
GUARD that prevents future regressions of the conditional status: a `scripts/verify/`
probe that mechanically asserts the paper does NOT make an unconditional absolute-ratio
claim and that the §11<->companion mapping is complete. Per the project's verification
discipline, a wording/governance fix without an executable check is not fully closed.

FIX TARGET (this leaf): add a Python probe under `scripts/verify/`
(e.g. `bug_001_conditional_ratio_claim.py`) that, against the committed tex sources,
asserts PASS/FAIL on document-grounded invariants -- for example: (a) every §11.1--11.6
prediction paragraph co-occurs with a conditional/falsifiable marker or a companion
cross-reference; (b) the Abstract/§1.1 absolute-ratio sentences retain their
"no quantitative claim ... empirical conjecture (§11)" disclaimer; (c) the
§11<->companion hypothesis mapping (from BUG-001-B) has no orphan on either side.
The probe prints PASS/FAIL to stdout and is wired as a CI job per the project's
verify-job convention. It is a string/structure invariant check, NOT a measurement
of compression ratio.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class; the "executable check" is a documentation-invariant
  probe over the tex sources, not a runtime product test.
- Fully closeable in-repo (the probe and its fixtures are committed files).
- Verify locally with the project interpreter (/Users/para/.venvs/rnr/bin/python),
  per the memory note that CI cannot be observed from this account.
- Non-goal: asserting any numeric bpb result (that requires BUG-001-D's external run).

## Refined Description
**Scope:** a PASS/FAIL `scripts/verify/` probe guarding the conditional-status and
mapping invariants established by BUG-001-A/B, plus its CI wiring.
**Non-goals:** measuring ratio; editing the claim text (that is BUG-001-A/B).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-001; sequentially after BUG-001-A and BUG-001-B (the
probe encodes the invariants those leaves establish)
**Risks / Open Questions:** probe must be robust to benign rewording (assert on
stable anchor phrases / labels, not brittle full-sentence matches) to avoid false CI
failures.
**Expected Code Change:** YES (new verify script + CI job entry)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | fix the exact invariants to assert |
| ESTIMATION | DONE | ABOVE_EASY (single probe over tex fixtures + CI wiring) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | choose stable anchor phrases / labels to match |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | scripts/verify/bug_001_conditional_ratio_claim.py |
| TESTING | PENDING | probe self-verifies (PASS/FAIL); run locally then CI |
| POST_AUDIT | PENDING | adversarial re-review per project review protocol |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at §12.7 decomposition of BUG-001;
  no code/proof change made before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
  The probe is a documentation-invariant check, consistent with the project's
  "every claim gets an executable verification" rule for a no-runtime-product repo.
- §12.7: BLOCKING child leaf of BUG-001 (the verification-script layer).


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-001 A/B/C (commit pending). A: §11.1-11.6 per-prediction conditional scoping ('falsifiable empirical conjecture, conditional and pending measurement, not a theorem') + canonical §11-intro paragraph; Abstract/§1.1 already scoped (left intact). B: §11.x<->companion-hypothesis mapping added (H1-H7, H15-H16) incl. seekable-baseline coverage note (bgzip/zstd-seekable = H10/H11 random-access TIME baselines; squashfs+zstd = H12 ratio baseline). C: scripts/verify/bug_001_ratio_claim_conditional.py (inv1 conditional-scoping + inv2 no-orphan mapping, non-vacuity self-test; per-paragraph locality + dominate-needs-baseline matcher hardened by supervisor; OVERALL -> PASS) + CI job. Hostile-referee POST_AUDIT clean. HONEST RESIDUAL (feeds D): §11.5 repetitive-ratio-vs-zstd/LZMA target + the beats-seekable value-prop have only PARTIAL pre-registered coverage (H7 = image bits-back; no dedicated bgzip/zstd-seekable RATIO hypothesis) -- flagged in-text, no H-number invented.

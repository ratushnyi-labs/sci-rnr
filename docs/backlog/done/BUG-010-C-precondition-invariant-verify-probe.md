# BUG-010-C - precondition-invariant-verify-probe

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** scripts/verify/ (new probe); tex/rnr_coding.tex (Abstract; §1.1; Def 2.1; Theorem 8.1; §10.2; §10.3) as read-only fixture
**Business Rule Conflict Check:** NOT_CHECKED
**Conflict Resolution Reference:** N/A
**Atomic Check:** ATOMIC
**Budget Forecast:** TBD
**Priority:** LOW
**Lawbook Version (intake):** 0.44.0
**Applicable Lawbook Version:** TBD
**Created:** 2026-06-15  **Updated:** 2026-06-15
**Recovered:** N/A

## Raw Request
> mark each gap as bug and follow the global rules
>
> (Filed gap 010 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at section 12.7 decomposition of BUG-010 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: Abstract; §1.1; Def 2.1; Theorem 8.1; §10.2; §10.3 (as the artifact under test) +
project rule "Практична верифікація — кожен крок як робоча програма" (every claim gets an
executable check).

The closeable-now scoping (BUG-010-A) and the precondition-characterization inequality
(BUG-010-B) need a GUARD that prevents future regression of the precondition's conditional
status and term-completeness: a `scripts/verify/` probe that mechanically asserts, against
the committed tex sources, that (i) the Abstract/§1.1 precondition sentences remain
conditional and carry the cross-references to Def 2.1 / Thm 8.1 / §10.2 that BUG-010-A
establishes; (ii) the paper makes NO unconditional universal-dominance claim (the Thm 8.1
invariant -- no sentence asserts RNR shortens every input or dominates every baseline);
(iii) the BUG-010-B precondition inequality lists every component cost term
(H(X|Y), L(M)/m, repair residual, verification overhead, baseline rate) and each term has
a resolvable §-anchor in the document. Per the project's verification discipline, a
governance/characterization fix without an executable check is not fully closed.

FIX TARGET (this leaf): add a Python probe under `scripts/verify/`
(e.g. `bug_010_precondition_characterization.py`) that prints PASS/FAIL to stdout on the
document-grounded invariants above, asserting on stable anchor phrases / labels (Def 2.1,
Thm 8.1, §10.2 labels; the precondition-inequality term names) rather than brittle
full-sentence matches, and is wired as a CI job per the project's verify-job convention.
It is a string/structure invariant check over the tex sources, NOT a measurement of
compression ratio and NOT a numeric evaluation of the precondition inequality on any
corpus (that is BUG-010-D).

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class; the "executable check" is a documentation-invariant
  probe over the tex sources, not a runtime product test and not a corpus measurement.
- Fully closeable in-repo (the probe and its fixtures are committed files).
- Verify locally with the project interpreter (/Users/para/.venvs/rnr/bin/python), per the
  memory note that CI cannot be observed from this account.
- Non-goal: asserting any numeric value of H(X|Y), L(M)/m, or ratio (those need BUG-010-D's
  external run); editing the claim text (BUG-010-A) or writing the inequality (BUG-010-B).

## Refined Description
**Scope:** a PASS/FAIL `scripts/verify/` probe guarding the conditional-precondition,
no-universal-dominance, and term-completeness invariants established by BUG-010-A/B, plus
its CI wiring.  **Non-goals:** measuring or numerically evaluating the precondition
(BUG-010-D); editing claim text (BUG-010-A); writing the inequality (BUG-010-B).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-010; sequentially after BUG-010-A and BUG-010-B (the probe
encodes the invariants those leaves establish).
**Risks / Open Questions:** probe must be robust to benign rewording (match stable anchor
labels / term names, not brittle full sentences) to avoid false CI failures; coordinate
with the BUG-001-C / BUG-002-C probes so the verify jobs do not duplicate or conflict.
**Expected Code Change:** YES (new verify script + CI job entry)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | fix the exact invariants and anchor labels to assert |
| ESTIMATION | DONE | ABOVE_EASY (single probe over tex fixtures + CI wiring) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | choose stable anchor phrases / labels (Def 2.1, Thm 8.1, §10.2, term names) |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | scripts/verify/bug_010_precondition_characterization.py |
| TESTING | PENDING | probe self-verifies (PASS/FAIL); run locally then CI |
| POST_AUDIT | PENDING | adversarial re-review per project review protocol |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at section 12.7 decomposition of BUG-010;
  no code/proof change made before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
  The probe is a documentation-invariant check, consistent with the project's
  "every claim gets an executable verification" rule for a no-runtime-product repo.
- §12.7: BLOCKING child leaf of BUG-010 (the verification-script layer).

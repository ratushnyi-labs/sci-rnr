# BUG-007-A - scope-thm-10-1-claim-vs-deferral

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Thm 10.1; §10.1; §10.5; cross-refs to §13.3)
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
> (Filed gap 007 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at §12.7 decomposition of BUG-007 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: Thm 10.1; §10.1; §10.5.

The whole side-information recasting requires the decoder to reproduce the model
output bit-for-bit. The paper REDUCES this to a precision/determinism specification
(Theorem 10.1, §10.5) but the surrounding prose must not overclaim. This leaf audits
ONLY the editorial honesty of the claim boundary: does §10.1/§10.5/Thm 10.1 state
exactly what is PROVEN (a sufficient condition: integer arithmetic + bit-width bound
+ deterministic non-linearity tables + pinned reduction order ⇒ bit-identity across
conforming implementations) and explicitly defer what is NOT proven (that a concrete
implementation exists and that named hardware actually conforms — currently routed to
§13.3)? Any sentence asserting end-to-end satisfiability has been DEMONSTRATED, rather
than SPECIFIED-and-deferred, is an overclaim to be tightened.

## Refined Description
**Scope (this leaf only):** Read §10.1, §10.5 (Thm 10.1 + its proof + trailing
discussion at lines ~15705-15711), and the §13.3 deferral bullet on "Reference
implementation conformance". Confirm the language is honest: Thm 10.1 is a SUFFICIENT
CONDITION (proven), not a claim that a conforming implementation has been built or that
end-to-end bit-exactness has been measured. If any prose implies the latter, refine the
wording to "specifies a sufficient condition / deferred to companion artifact" and ensure
the deferral pointer (§13.3) is present and accurate at the point of claim.

**Fix target (this leaf):** a precise scoping edit (or a confirmation that no edit is
needed) so the paper claims exactly Thm 10.1's proven content and explicitly flags the
reference-implementation + cross-platform-test as deferred engineering.

**Non-goals:** writing/running any verification script (BUG-007-B); strengthening the
spec/experimental-design cross-links (BUG-007-C); the external implementation itself
(BUG-007-D). No change to the mathematical statement or proof of Thm 10.1.

**Closeable in-repo NOW:** YES — pure documentation scoping inside `tex/rnr_coding.tex`.

**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-007.
**Risks / Open Questions:** the wording may already be honest (lines 15705-15711 and
§13.3 read as deferral, not demonstration); if so, this leaf closes by confirming no
overclaim and recording the citation, no edit needed.
**Expected Code Change:** YES (tex prose; possibly zero-diff confirm-only)
**Documentation Recovery Required:** NO

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class defect in the theory paper `tex/rnr_coding.tex`
  (no runtime product). "Fix" = scope the paper's claim to match what is established.
- Severity inherits HIGH from parent: the bit-exact reproduction premise is load-bearing
  for the whole side-information recasting (main idea).

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = audit claim-vs-deferral honesty only |
| ESTIMATION | DONE | EASY (prose scoping audit) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | identify exact sentences to confirm/tighten |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex prose edit (or confirm zero-diff) |
| TESTING | N/A | no executable artifact in this leaf |
| POST_AUDIT | PENDING | adversarial re-read: any residual overclaim? |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7: child leaf of BUG-007; created fresh this turn (Recovered: N/A).
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Honest-step principle: this leaf carries the part that CAN be closed in-repo now
  (precise claim scoping), keeping the external implementation as a separate flagged leaf.

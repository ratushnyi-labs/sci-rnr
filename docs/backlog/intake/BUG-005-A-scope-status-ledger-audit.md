# BUG-005-A - scope-status-ledger-audit

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§13.2; §13.4.1/13.4.2/13.4.4; OP2(a),OP2(b),OP4)
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
> (Filed gap 005 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at section 12.7 decomposition of BUG-005 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: §13.2; §13.4.1/13.4.2/13.4.4; OP2(a),OP2(b),OP4.

The parent gap states that "'optimal' RNR is not known to be efficiently
approximable" because the unrestricted bit-cost SGP (OP2a discrete/ceiling),
the free-predictor joint Type-III-C (OP2b-APX), and general support-selection
approximability (OP4, the E_12 constant gap) are OPEN. The paper has, in fact,
already CLOSED several sub-cases of this gap (continuous-bit-cost OP2(a) via
Lemma 6.8e; continuous OP2(b) in P via Remark 6.9a; fixed-field Type-III-C
APX-hard 5/4 via Lemma 6.9; additive rule-cost SGP via Lemmas 6.8k / 6.8k-approx;
OP4' continuous via Lemma 6.8e). What is missing is a SINGLE precise ledger that
states, per sub-question, exactly which variant is CLOSED, which is OPEN, the
named structural blocker for each open variant, and where in the paper / which
verify script witnesses each entry. Absent this ledger the headline "optimal RNR
is not efficiently approximable" reads as a blanket open problem when it is in
fact a small, precisely-pinned residual set.

## Clarification / Assumptions
- This is the IN-REPO, CLOSEABLE-NOW scoping leaf of BUG-005: it neither closes
  nor attempts any open hardness reduction. It scopes the paper's claim to match
  what is actually established, so the headline complexity of the parent drops to
  a shippable core (this audit) plus a clearly-flagged research residual
  (sibling BUG-005-D).
- GOVERNANCE/DOCUMENTATION-class defect in the theory paper `tex/rnr_coding.tex`
  (no runtime product). "Fix" here = make the claim precise / build the ledger,
  NOT prove an open theorem.
- Severity inherits parent MEDIUM (load-bearing for the "optimal" framing of the
  side-information+repair main idea, but the headline rate/RA results do not
  depend on optimality of the encoder).

## Refined Description
**Scope:** Produce a single per-sub-question status ledger covering OP2(a)
{continuous = CLOSED via 6.8e; discrete/ceiling = OPEN, blocker = ceiling
discontinuity at power-of-two |Σ|+K boundaries / needs sharper-than-Charikar
8569/8568 SGP constant}, OP2(b) {continuous = in P via 6.9a; single-source
discrete = OPEN, blocker = predictor-absorption Lemma 13.4.2a; amortized family
= OPEN, PRG route gives decision-only Theorem 13.4.2b, not APX}, OP4 {E_12
constant gap = OPEN, blocker = partial-cube no-robust-NO promise per Remark 4.3d;
fixed-byte = OPEN, byte-encoding dilution; Type-III-C fixed-field = CLOSED
APX-hard 5/4 via Lemma 6.9; count/entropy objectives all dilute, gap lives in
Type-III-C log-det}. For each row: status, paper location, named blocker (if
open), witnessing verify script (if any). Then confirm the §13.2 / §13.4 /
§1.3 / abstract prose states the residual at exactly this granularity (no
blanket overclaim, no understatement of what is closed).
**Non-goals:** any new reduction, any attempt at OP2a-discrete / OP2b-APX / OP4
(those are sibling BUG-005-D); verification-script consolidation (sibling
BUG-005-B); dead-end ledger cross-linking of memory journals (sibling
BUG-005-C).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-005.
**Risks / Open Questions:** the ledger must not claim a closure the paper does
not actually establish (the op_progress_2026_06 journal records that OP2(b)-APX
was briefly mis-recorded as closed; the correct status is OPEN — "both have a
poly witness" ≠ "optima equal"). Cold re-read against live HEAD required before
asserting any row CLOSED.
**Expected Code Change:** YES (tex prose/table edit only; no math)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | confirm ledger granularity = per sub-question variant |
| ESTIMATION | PENDING | EASY: read §13.2/§13.4, tabulate existing closures+blockers |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | ledger as a §13.2 status table (mirror existing status-update prose) |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | add/refresh the per-sub-question status table in §13.2 |
| TESTING | PENDING | none needed (doc-only); rows cross-ref verify scripts owned by sibling B |
| POST_AUDIT | PENDING | cold re-read + independent adversarial re-review per project review protocol (each CLOSED row must survive "Attack failed; holds") |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 DECOMPOSITION child of BUG-005. This is the honestly-closeable-now,
  in-repo scoping leaf: it lowers the parent's effective complexity by separating
  the precise-claim work from the genuinely-open research residual.
- §12.1.1 (Backlog-File-First): created fresh this turn as a decomposition child;
  no paper change precedes this file. `Recovered: N/A`.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Per project rules (Truth > peremoga): do NOT record any sub-question as CLOSED
  unless the cited lemma + verify script actually establish it on live HEAD.

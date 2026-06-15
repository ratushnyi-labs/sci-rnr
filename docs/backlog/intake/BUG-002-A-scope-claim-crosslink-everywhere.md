# BUG-002-A - scope-claim-crosslink-everywhere

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Abstract ll.202-204; §10.2 Model amortization regimes ll.15573-15582; Theorem 7.15 caveat ll.8573-8606; §11.7 ll.16814-16835)
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
> (Filed gap 002 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at §12.7 decomposition of BUG-002 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: parent BUG-002 FIX TARGET, clause "...state the break-even amortization
(blocks/shared-model) explicitly **wherever the number appears**."

The model-overhead caveat and the break-even amortization volume ($V^* \approx
152\,\mathrm{GB}$) already exist, but only in two of the four places the headline
number is stated:
- **Theorem 7.15 "Side-information framing and model overhead (caveat)"** (ll.8573-8606)
  carries the full caveat + negative-savings ($\approx140\times$ expansion for a
  single archive) + the break-even formula $V^*\cdot(\log_2|\Sigma|-H_{M'})=140\,\mathrm{GB}\cdot\log_2|\Sigma|\Rightarrow V^*\approx152\,\mathrm{GB}$.
- **§11.7** (ll.16826-16831) cross-links the caveat and restates $V^*\gtrsim152\,\mathrm{GB}$.

But the **Abstract** (ll.202-204) carries only a one-clause caveat ("excludes the
$\approx$140 GB model") with **no break-even pointer**, and **§10.2 "Model
amortization regimes"** (ll.15573-15582) discusses single-block / multi-block /
streaming regimes and $L(M)/m$ amortization but **never states the break-even
volume nor cross-links the T7.15 caveat**, so the "wherever the number appears"
audit is incomplete and the document is internally under-cross-linked.

## Refined Description
**One-leaf scope:** a pure documentation cross-link / consistency pass. Make the
break-even amortization explicit (or an explicit pointer to the T7.15 caveat that
states it) at EVERY site where the headline 79 MB / $12.0\times$ / 0.664 bpb number
appears: Abstract, §10.2, and any other occurrence surfaced by an audit grep. Do
NOT re-derive the number; only ensure each occurrence is accompanied by (i) the
model-inclusive (self-contained) framing and (ii) a break-even pointer.

**Fix target:** add a short break-even cross-reference clause to the Abstract
caveat and to §10.2; verify T7.15 and §11.7 remain consistent. No new theorem,
no numeric re-derivation (that is BUG-002-B / BUG-002-C scope).

**Closeable in-repo now:** YES (paper-scoping / doc edit only).

**Non-goals:** the consolidated accounting table (BUG-002-B); the verification
probe (BUG-002-C); the empirical enwik9 measurement (BUG-002-D).
**Impacted UCs:** N/A (theory paper)  **Impacted BR/WF:** N/A
**Dependencies:** parent BUG-002. Soft-ordered after BUG-002-B (the table gives a
canonical wording/anchor to point to), but can proceed independently with a
pointer to the existing T7.15 caveat.
**Risks / Open Questions:** the supervisor edits tex/; this leaf only specifies and
stages the precise wording/sites. NOTE: actual tex/rnr_coding.tex edits are made by
the supervisor on reconciliation, not in the intake directory.
**Expected Code Change:** YES (tex prose, supervisor-applied)
**Documentation Recovery Required:** NO

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class defect in the theory paper; "fix" = scope the
  claim precisely / restore internal cross-link consistency.
- The break-even number itself ($V^*\approx152\,\mathrm{GB}$) is taken as given
  from the T7.15 caveat; its correctness is audited in BUG-002-C, not here.

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = cross-link consistency, no re-derivation |
| ESTIMATION | DONE | EASY (prose cross-link pass) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | list exact insertion sites + clause wording |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex prose clauses (supervisor-applied) |
| TESTING | PENDING | grep-audit that every headline-number site now carries break-even pointer |
| POST_AUDIT | PENDING | adversarial re-read: did any occurrence get missed? |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7: created fresh as a BLOCKING child leaf of BUG-002 on 2026-06-15; not recovered.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Stays in `docs/backlog/intake/`; no tex/ edits performed in this turn.

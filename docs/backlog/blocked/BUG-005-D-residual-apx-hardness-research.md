# BUG-005-D - residual-apx-hardness-research

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** BLOCKED
**Complexity:** HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Conjecture 6.8d; Conjecture 4.3e; Theorem 13.4.2b; OP2(a)-discrete, OP2(b)-APX, OP4 E_12)
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
Source: §13.2; §13.4; Conjecture 6.8d, Conjecture 4.3e, Theorem 13.4.2b.

This is the IRREDUCIBLE research residual of BUG-005, after the scoping (A),
verification (B), and dead-end-cross-link (C) leaves carve off everything
honestly closeable in-repo. Three named open questions remain, all sharing the
SAME structural obstacle pattern (a hardness gap that DILUTES under the
encoding / objective transform, leaving no robust YES/NO promise):
1. **OP2(a) discrete/ceiling bit-cost SGP** — Conjecture 6.8d (weak form, "some
   ε>0"). Residual blocker: the ceiling discontinuity at power-of-two |Σ|+K
   boundaries; closing the strong γ_SGP form needs a sharper-than-Charikar
   8569/8568 SGP inapproximability constant (itself a ~20-yr open SGP-literature
   problem).
2. **OP2(b)-APX free-predictor / amortized-family joint Type-III-C** — OPEN
   (continuous case is in P, Remark 6.9a). Residual blocker: predictor-absorption
   (Lemma 13.4.2a is an UPPER bound, not a closure); the PRG construction
   (Theorem 13.4.2b) yields average-case DECISION hardness only, never
   multiplicative APX.
3. **OP4 E_12 constant-gap inapproximability** — Conjecture 4.3e. Residual
   blocker: partial-cube no-robust-NO promise (Remark 4.3d); a constant gap needs
   a NEW max-dilation-1 / min-cubical-deletion gap theorem (none in the
   literature). Count and entropy objectives all dilute; the gap, where it does
   exist, lives in Type-III-C log-det (Lemma 6.9), already closed.
FIX TARGET (parent's words): close ONE open direction OR pin the inapproximability
constant.

## Clarification / Assumptions
- This leaf is the genuinely-hard, NOT-closeable-in-repo-now research residual.
  It is kept SINGLE (not split into three hard leaves) because the three open
  questions share one obstacle pattern and the project's waterfall rule forbids
  parallelising multiple open theorems; whichever direction is attacked, only ONE
  is taken at a time, and "weakened theorem holds" / "Attack failed; holds" /
  honest non-closure are all acceptable terminal states.
- closeable_in_repo = NO: per the OP2/OP4 attempt journals and op_progress
  2026-05/06, every standard reduction route has been tried and pinned to a named
  blocker; genuine closure needs either new mathematical machinery or a
  sharper external (SGP-literature / PCP) constant — explicitly flagged as
  beyond an in-session closure by the journals.
- The closeable-in-repo FALLBACK for this leaf (if no closure is achieved) is to
  do NOTHING new in the paper beyond what siblings A/C already scope — i.e. the
  honest terminal state is "residual remains open, precisely pinned." Do not
  manufacture a borderline-tautological closure (op2_attempts journal: Observation
  6.8f excluded the obstacle by fiat — not acceptable).

## Refined Description
**Scope:** Attempt to close exactly ONE of {OP2(a)-discrete, OP2(b)-APX, OP4
E_12} OR to pin a concrete inapproximability constant, by a route that ADDRESSES
the recorded blockers (a sharper SGP constant > Charikar's for OP2a; a
gap-preserving multiplicative reduction surviving predictor-absorption for OP2b;
a max-dilation-1/min-cubical-deletion robust-gap theorem for OP4). Any attempt
MUST follow the project 3-step post-fix protocol (cold re-review +
independent adversarial re-review + adversarial disprove "PROVE wrong / Attack failed")
before any commit, and the "small honest steps" rule (a weakened but proven
sub-claim beats a strong hand-wave).
**Non-goals:** the in-repo scoping table (sibling A), verify consolidation
(sibling B), and failed-path cross-linking (sibling C) — those carry no open
math and must be completed first/independently. Do NOT re-attempt any journalled
dead end (density-shifted BK, Dinur-Safra FGLSS, K-escape budget, Label-Cover
K-uniform, multiway-cut/heavy-edge, bit-cut/spectral, cube-host-QAP layout
transfer, Liu-Pass-as-APX, padding-bit multiplicative) without first defeating
its named blocker.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-005; sequentially AFTER siblings BUG-005-A and
BUG-005-C (the precise status + faithful dead-end ledger must exist before
attacking the residual, so the attempt starts from a clean obstacle map).
**Risks / Open Questions:** HIGH research risk — the journals document this as
genuinely research-level open, plausibly needing months of focused proof
engineering or new machinery. The most likely terminal state is an honest
non-closure (residual remains, more sharply pinned), which is an acceptable
outcome under the project rules and does NOT block shipping siblings A/B/C.
**Expected Code Change:** UNCERTAIN (NO if non-closure; a new lemma + verify
script if a sub-case genuinely closes)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | pick ONE direction (waterfall: one open theorem at a time) |
| ESTIMATION | PENDING | HARD/research; outcome may be honest non-closure |
| DECOMPOSITION | N/A | atomic residual; cannot be honestly split further (shared obstacle) |
| DESIGN | PENDING | route must defeat the named blocker for the chosen direction |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | new lemma + proof IF closed; else no paper change beyond A/C scoping |
| TESTING | PENDING | constructive/reduction verify script IF a sub-case closes |
| POST_AUDIT | PENDING | mandatory 3-step protocol (cold + independent adversarial re-review + adversarial disprove) before any commit |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 DECOMPOSITION child of BUG-005. This is the SINGLE genuinely-hard
  research residual; siblings A/B/C carry all the in-repo-closeable work so the
  parent ships a precise core plus this clearly-flagged open residual.
- closeable_in_repo = NO: needs external/new mathematics (sharper SGP constant,
  PCP-style robust-gap source, or new partial-cube embedding gap theorem) per the
  recorded blockers; an honest non-closure that more sharply pins the residual is
  an acceptable terminal state and must not be replaced by a tautological closure.
- §12.1.1 (Backlog-File-First): created fresh this turn; no paper change precedes
  this file. `Recovered: N/A`.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Mandatory: project 3-step post-fix protocol + adversarial disprove before any
  commit; "Truth > peremoga"; verify LOCALLY (gh cannot observe this repo's CI).

## Status (2026-06-15; HEAD ec42c03)
BLOCKED / OPEN -- intentionally NOT worked. This is the research residual; per the
parent decomposition it stays open. Siblings A (status ledger), B (index probe +
CI), and C (failed-path blocker crosslink + Remark-4.3d framing correction) are
DONE and shipped, which precisely pins this residual to the set
{ OP2(a)-discrete (Conj 6.8d, ceiling discontinuity), OP2(b)-discrete-APX
(predictor-absorption upper-bound-only / PRG reduces to OP2(a)-APX), OP4 Type-II
approximability, OP4 E_12 (Conj 4.3e, missing dilation-1/robust-NO source),
OP4 fixed-byte }. The op2/op4 attempt-journals document DEAD-END routes; do NOT
re-attempt them. Moved to blocked/.

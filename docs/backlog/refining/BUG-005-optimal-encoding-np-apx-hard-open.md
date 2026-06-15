# BUG-005 - optimal-encoding-np-apx-hard-open

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DECOMPOSED
**Complexity:** VERY_HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§13.2; OP2(a),OP2(b),OP4)
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
> (Filed gap 005 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

## Gap (defect statement; document-grounded)
Source: §13.2; OP2(a),OP2(b),OP4.

Computing the OPTIMAL RNR encoding (support/grammar selection) is NP-/APX-hard in restricted forms; the unrestricted bit-cost SGP (OP2a integer/ceiling variant), the free-predictor joint (OP2b), and general support-selection approximability (OP4) are OPEN. So 'optimal' RNR is not known to be efficiently approximable. NOTE: op2/op4 attempt-journals document dead-end routes -- do not re-attempt without addressing the recorded blockers. FIX TARGET: close one open direction or pin the inapproximability constant.

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
| DECOMPOSITION | DONE | split into BUG-005-A, BUG-005-B, BUG-005-C, BUG-005-D (see Decomposition section) |
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
Split on 2026-06-15 into four BLOCKING children. The closeable-now-in-repo work
(precise scoping, executable-verification coverage, faithful dead-end ledger) is
carved off into EASY/ABOVE_EASY leaves so the parent's complexity drops to a
shippable core plus ONE clearly-flagged genuinely-open research residual.

- **BUG-005-A — scope-status-ledger-audit** (EASY; closeable in-repo: YES).
  Build a single per-sub-question status ledger (OP2a continuous CLOSED / discrete
  OPEN; OP2b continuous in-P / single+amortized OPEN; OP4 E_12 OPEN / Type-III-C
  fixed-field CLOSED), each row with status + paper location + named blocker +
  witnessing script; make §13.2 / §13.4 / §1.3 / abstract prose match this exact
  granularity (no blanket overclaim, no understatement of what is closed).
- **BUG-005-B — verify-probe-consolidation** (ABOVE_EASY; closeable in-repo: YES).
  Cross-check every §13.4 closed attack vector against a present, PASSing
  `scripts/verify/` script; add one consolidated index/status probe asserting
  presence+PASS; confirm the CI `verify` job covers the OP2/OP4 closures (verify
  LOCALLY).
- **BUG-005-C — dead-end-blocker-crosslink** (ABOVE_EASY; closeable in-repo: YES).
  Transcribe the validated OP2/OP4 attempt-journal blockers into the §13.4
  failed-path roadmap at re-attempt-prevention granularity; correct any stale
  obstacle framing (e.g. "row-stochastic F is THE obstacle" → partial-cube
  no-robust-NO promise per Remark 4.3d); doc-only, no new math.
- **BUG-005-D — residual-apx-hardness-research** (HARD; closeable in-repo: NO).
  The irreducible research residual: close ONE of {OP2a-discrete (Conj 6.8d),
  OP2b-APX, OP4 E_12 (Conj 4.3e)} OR pin an inapproximability constant, by a route
  that defeats the recorded blocker. Kept SINGLE (shared obstacle pattern;
  waterfall = one open theorem at a time). Honest non-closure is an acceptable
  terminal state and does not block shipping A/B/C; mandatory 3-step adversarial
  protocol before any commit.

Recommended first leaf: **BUG-005-A** (precise status ledger; unblocks C and D by
giving a clean obstacle map, and is independently shippable).

## BUILD outcome (2026-06-15; HEAD ec42c03)
In-repo leaves A/B/C DONE (moved to done/); research residual D BLOCKED/OPEN
(moved to blocked/) -- intentionally not worked, per decomposition.
- A: 9-row per-sub-question status ledger table added at head of §13.2
  (tex/rnr_coding.tex). Card granularity MATCHES the paper; no card-vs-paper
  status disagreement. Abstract + §1.3 prose re-read, already precise/consistent
  (no blanket "optimal RNR is NP-hard" overclaim, no understatement of closures).
- B: consolidated index probe scripts/verify/bug_005_optimal_encoding_status_index.py
  (asserts CLOSED/in-P witnessing scripts present + exit 0; OPEN rows carried with
  NO fabricated PASS; OVERALL -> PASS). Wired into CI verify job (build.yml).
  Full local CI verify set (16 probes) all OVERALL -> PASS.
- C: §13.4.1/13.4.2 failed-path lists already faithful; corrected the stale
  "row-stochastic F is THE obstacle" framing in the §13.2 E_12 paragraph to the
  partial-cube / no-robust-NO-promise obstruction (Remark 4.3d, verified to say
  this), and cross-linked it into the §13.4.4 OP4 roadmap.
Compiles clean (tectonic). Remaining open residual (D): { OP2(a)-discrete
(Conj 6.8d), OP2(b)-discrete-APX, OP4 Type-II approx, OP4 E_12 (Conj 4.3e),
OP4 fixed-byte }. Parent stays in refining/ with this clearly-flagged residual,
mirroring the BUG-004 pattern.

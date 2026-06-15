# BUG-002 - headline-bpb-excludes-model

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE (in-repo); D external-blocked
**Complexity:** MEDIUM
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Abstract (Thm 7.15 synthesis); §10.2)
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
> (Filed gap 002 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

## Gap (defect statement; document-grounded)
Source: Abstract (Thm 7.15 synthesis); §10.2.

The headline end-to-end number (~0.664 bpb / ~79 MB on enwik9) carries an explicit caveat that it EXCLUDES the ~140 GB model. For a single archive L(M) dominates, so the stated 'compression' does not account for model cost. FIX TARGET: present a model-inclusive accounting alongside the amortized one, and state the break-even amortization (blocks/shared-model) explicitly wherever the number appears.

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
**Expected Code Change:** YES
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | collate+crosslink+verify in-repo (A/B/C); empirical measurement (D) external |
| ESTIMATION | PENDING | set Complexity/Budget |
| DECOMPOSITION | DONE | split into BUG-002-A, BUG-002-B, BUG-002-C, BUG-002-D |
| DESIGN | DONE | consolidated table in §10.2; break-even V* identity; arithmetic probe |
| FRONTEND | N/A | no UI |
| BACKEND | DONE | §10.2 consolidated accounting table + Abstract & §10.2 break-even crosslinks (A,B) |
| TESTING | DONE | scripts/verify/bug_002_model_inclusive_accounting.py (C), OVERALL -> PASS + CI job |
| POST_AUDIT | DONE | hostile-referee PASSED: no-new-claims, arithmetic (V*=152.67/12.048x/140.08x/T7.21), crosslinks, compile clean |

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
Split on 2026-06-15 into four BLOCKING child leaves separating the part that can be
HONESTLY CLOSED NOW IN-REPO (precise scoping, an accounting table, an executable
arithmetic probe) from the IRREDUCIBLE EXTERNAL measurement residual. Audit of the
paper found the model-overhead caveat and break-even $V^*\approx152\,\mathrm{GB}$
already exist at T7.15 (ll.8573-8606) and §11.7, with the model-encoded total-cost
formula at T7.21 (ll.9750-9782) and the regime taxonomy at §10.2 (ll.15573-15582);
the gaps are (i) the break-even is not stated/cross-linked at the Abstract and §10.2,
(ii) no single consolidated model-inclusive-vs-amortized table, (iii) no executable
verification of the accounting arithmetic, and (iv) no actual measured end-to-end run.

- **BUG-002-A** (`scope-claim-crosslink-everywhere`) — EASY — closeable in-repo: YES.
  Cross-link consistency pass: state the break-even amortization (or a pointer to the
  T7.15 caveat) at EVERY site where the 79 MB / 12.0× / 0.664 bpb number appears
  (Abstract ll.202-204 and §10.2 currently lack it).
- **BUG-002-B** (`model-inclusive-accounting-table`) — ABOVE_EASY — closeable in-repo: YES.
  Synthesize the scattered T7.15 / T7.21 / §10.2 numbers into ONE consolidated
  model-inclusive-vs-amortized accounting table (regime, bitstream cost, model
  included, net savings/expansion, break-even); collation only, no new claim.
- **BUG-002-C** (`breakeven-accounting-verify-probe`) — EASY — closeable in-repo: YES.
  Add/extend a `scripts/verify/` probe that recomputes the self-contained-vs-amortized
  arithmetic, the $V^*\approx152\,\mathrm{GB}$ break-even, and T7.21's $L^*(N)$ optimum,
  emitting PASS/FAIL (internal-consistency check only; does NOT measure the empirical datum).
- **BUG-002-D** (`empirical-enwik9-measurement-external`) — HARD — closeable in-repo: NO.
  The irreducible external residual: actually build an RNR archive with a frontier
  predictor on enwik9 and measure bitstream + self-contained sizes (H15, experimental-design
  §2.11). Needs external compute + reference implementation (§13.3); the in-repo
  scoping sliver is delegated to BUG-002-A so this stays a clean external residual.

Recommended first leaf: **BUG-002-B** (canonical accounting table the others anchor to).


## Resolution (2026-06-15)
A/B/C closed in-repo (moved to done/): consolidated model-inclusive-vs-amortized
accounting table in §10.2 (collation of T7.15/T7.21/§10.2, every cell cited, no new
claims), break-even V*~152 GB crosslinks at the Abstract and §10.2, and the arithmetic
verify probe + CI job. POST_AUDIT (hostile referee) passed clean. D (empirical enwik9
measurement) is the irreducible external residual -> blocked/ (needs GPU + frontier LM +
reference impl, §13.3); NOT claimed closed -- the headline number stays scoped as a
prediction/conjecture, not a measured datum.

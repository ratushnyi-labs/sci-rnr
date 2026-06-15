# BUG-010 - conditional-scope-no-universal-dominance

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DECOMPOSED
**Complexity:** MEDIUM
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Abstract; §1.1)
**Business Rule Conflict Check:** NOT_CHECKED
**Conflict Resolution Reference:** N/A
**Atomic Check:** ATOMIC
**Budget Forecast:** TBD
**Priority:** LOW
**Lawbook Version (intake):** 0.44.0
**Applicable Lawbook Version:** TBD
**Created:** 2026-06-14  **Updated:** 2026-06-15
**Recovered:** YES

## Raw Request
> mark each gap as bug and follow the global rules
>
> (Filed gap 010 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

## Gap (defect statement; document-grounded)
Source: Abstract; §1.1.

By construction the framework disclaims universal dominance (no lossless code shortens every input) and is 'evaluated conditionally on data classes for which side information reduces residual uncertainty by more than the model and repair overhead.' Which real corpora actually qualify -- the data-class precondition -- is itself the unmeasured empirical question, and only the model-amortization regime (§10.2) can win. FIX TARGET: a precise, checkable characterization of the qualifying data-class precondition (tie to BUG-001/002).

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
**Dependencies:** BUG-001 / BUG-002
**Risks / Open Questions:** see Gap; some fixes need work outside this repo (experiments).
**Expected Code Change:** NO
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope-the-claim branch (A/B/C in-repo) vs close-the-gap (D, external); the conditional claim now traces to its formal predicate |
| ESTIMATION | DONE | A EASY, B/C ABOVE_EASY (in-repo, DONE); D VERY_HARD (external, BLOCKED) |
| DECOMPOSITION | DONE | split into BUG-010-A, BUG-010-B, BUG-010-C, BUG-010-D |
| DESIGN | DONE (A/B/C) | A: Abstract+§1.1 cross-refs to Def 2.1/Thm 8.1/§10.2; B: §10.2.1 precondition inequality (10.1); C: invariant probe. D: external measurement, blocked |
| FRONTEND | N/A | no UI |
| BACKEND | DONE (A/B/C) | tex: Abstract+§1.1 conditional cross-links, new §10.2.1 inequality (10.1), and a new \textbf{Definition 2.1} header (the cross-ref target did not exist — caught in cold review); scripts/verify/bug_010_precondition_invariant.py (inv1-inv4 + 6 self-tests); CI step. D: BLOCKED |
| TESTING | DONE (A/B/C) | bug_010_precondition_invariant.py inv1-inv4 OVERALL -> PASS; lualatex 371pp EXIT=0, no undefined refs; sibling BUG-001/002 probes still PASS |
| POST_AUDIT | DONE (A/B/C) | hostile-referee audit PASSED; PLUS my cold re-read caught a real defect BOTH agents missed — "Definition 2.1" was cross-referenced 6x with no such header (same dangling-prose-ref class as the 7.34c fix). Fixed by promoting the §2.1 predicate to a labelled Definition 2.1; hardened the probe with inv4 (anchor-resolution) to prevent recurrence |

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
Decomposed on 2026-06-15 into four BLOCKING child leaves. The split separates the part
that can be honestly closed now in-repo (scope the conditional claim against its formal
predicate; write the precondition as an auditable inequality; guard it with a verify
probe) from the genuinely-external empirical residual (which real corpora actually
satisfy the precondition). Atomic Check remains ATOMIC: this is a decomposition into
children, not an umbrella peer-split.

- **BUG-010-A** (scope-precondition-claim-crosslink) — edit the Abstract + §1.1
  conditional-claim sentences to point explicitly at Def 2.1 (conditionally-advantageous
  predicate), Theorem 8.1 (no universal dominance), and §10.2 (model-amortization regime),
  uniformly phrased as conditional pending the empirical question. Complexity: EASY.
  Closeable in-repo: YES.
- **BUG-010-B** (precondition-accounting-characterization) — state the qualifying data-class
  precondition as one precise inequality with each component cost term (H(X|Y), L(M)/m,
  repair residual, §10.3 verification overhead, baseline rate) named and cross-referenced,
  naming candidate structural data-class families as conjecture, reconciled with
  BUG-001/002 accounting. Complexity: ABOVE_EASY. Closeable in-repo: YES.
- **BUG-010-C** (precondition-invariant-verify-probe) — add a `scripts/verify/` PASS/FAIL
  probe asserting the document invariants: precondition stays conditional with its
  cross-references, no unconditional universal-dominance sentence exists, and the
  BUG-010-B inequality's terms are complete and §-anchored; wire as a CI job.
  Complexity: ABOVE_EASY. Closeable in-repo: YES.
- **BUG-010-D** (which-corpora-qualify-measurement-external) — the irreducible external
  residual: MEASURE which real corpora satisfy the BUG-010-B precondition inequality in
  the §10.2 amortized regime, reconcile with the §1.1 plausibly-qualifying-families
  conjecture; rides on the BUG-001-D / BUG-002-D campaign. Complexity: VERY_HARD.
  Closeable in-repo: NO (needs reference implementation + corpora + compute).

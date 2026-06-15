# BUG-010-A - scope-precondition-claim-crosslink

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Abstract; §1.1) with read-only anchors §2.1 (Def 2.1 conditionally-advantageous), Theorem 8.1, §10.2
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
Source: Abstract ("the framework is evaluated conditionally on data classes for which
side information reduces residual uncertainty by more than the model and repair overhead
costs"); §1.1 ("The central claim is conditional: for data classes in which neural side
information reduces residual uncertainty by more than the cost of describing the model,
the repair, the verification, and the metadata, RNR can in principle outperform ...").

The closeable-in-repo wording half of BUG-010. The conditional disclaimer is stated in
the Abstract and §1.1, but in prose only: the "data-class precondition" is named without
being tied to the formal predicate that already exists in the paper -- Definition 2.1
("conditionally advantageous on D" iff E[L(Enc(X))] < E[L(B(X))]), the no-universal-
dominance Theorem 8.1 ("the claim is conditional and distributional: for data classes
whose conditional entropy under deterministic side information is below the baseline rate
by more than the model and overhead cost"), and the §10.2 model-amortization regime where
L(M)/m governs whether the cost term is payable. A reader cannot today trace the
Abstract/§1.1 word "conditionally" to its operational definition.

FIX TARGET (this leaf): edit the Abstract and §1.1 conditional-claim sentences so the
word "conditionally" / "data-class precondition" carries an explicit forward pointer to
(a) Definition 2.1 (the conditionally-advantageous predicate), (b) Theorem 8.1 (no
universal dominance, the conditional-and-distributional framing), and (c) §10.2 (the
model-amortization regime that makes the L(M) cost term payable -- the "only the model-
amortization regime can win" observation), and so the claim is uniformly phrased as
conditional pending the empirical data-class question (cross-link sibling BUG-010-D and,
where relevant, BUG-001/002). No new measurement, no proof, no change to Def 2.1 / Thm 8.1
/ §10.2 themselves (read-only anchors).

## Clarification / Assumptions
- This is a GOVERNANCE/DOCUMENTATION-class defect in the theory paper
  `tex/rnr_coding.tex` (no runtime product). "Fix" means here: scope the paper's claim to
  match what is actually established (the SCOPE-the-claim branch of BUG-010's FIX TARGET).
- This leaf is fully closeable in-repo: it is a wording/cross-reference edit only.
- Non-goal: writing the precondition accounting inequality (BUG-010-B), the verify probe
  (BUG-010-C), and the empirical which-corpora-qualify question (BUG-010-D).
- Aligns with BUG-001-A (uniform conditional scoping) without overlapping it: this leaf is
  specifically the no-universal-dominance / data-class-precondition sentences, not the
  §11 absolute-ratio predictions.

## Refined Description
**Scope:** Abstract + §1.1 conditional-claim sentences gain explicit pointers to Def 2.1,
Thm 8.1, §10.2 and to the empirical residual (BUG-010-D).  **Non-goals:** the accounting
table (BUG-010-B); the probe (BUG-010-C); the measurement (BUG-010-D); §11 predictions
(BUG-001).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-010
**Risks / Open Questions:** the cross-reference must not silently strengthen the claim
(it remains conditional + in-principle); verify the cited labels (Def 2.1, Thm 8.1, §10.2)
resolve before pointing at them.
**Expected Code Change:** YES (tex wording / cross-references)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | confirm canonical "conditional, see Def 2.1 / Thm 8.1 / §10.2" phrasing |
| ESTIMATION | DONE | EASY (localized wording + cross-ref edit, 2 sites) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | choose one reusable conditional-precondition sentence with the three anchors |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex wording edits in Abstract + §1.1 |
| TESTING | PENDING | manual re-read; mechanically guarded by BUG-010-C probe |
| POST_AUDIT | PENDING | adversarial re-review per project review protocol |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at section 12.7 decomposition of BUG-010;
  no code/proof change made to the paper before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- §12.7: this is one BLOCKING child leaf of BUG-010 (the closeable-now scoping layer).

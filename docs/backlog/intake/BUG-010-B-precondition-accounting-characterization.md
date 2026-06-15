# BUG-010-B - precondition-accounting-characterization

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Abstract; §1.1; §2.1 Def 2.1; Theorem 8.1; §10.2; §10.3 overhead accounting) -- a new precondition-characterization passage anchored here
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
Source: Abstract / §1.1 ("data classes for which side information reduces residual
uncertainty by more than the model and repair overhead costs"); Def 2.1 (conditionally-
advantageous predicate); Theorem 8.1 ("conditional entropy ... below the baseline rate by
more than the model and overhead cost"); §10.2 (L(M)/m amortization); §10.3 (verification
overhead accounting, the per-subblock hash + offset-table bits).

BUG-010's FIX TARGET is "a precise, checkable CHARACTERIZATION of the qualifying data-class
precondition." Today the precondition is described in words across the Abstract, §1.1, Thm
8.1, and §10.2, but it is never written down as a single explicit inequality with every
term named and matched to where that term is defined/bounded in the paper. The closeable-
in-repo half here is to STATE the precondition precisely: assemble the terms the paper
already defines into one characterization -- in symbols, the qualifying-precondition
inequality is (schematically)

  H(X | Y) + L(M)/m + E[L(repair residual)] + verification_overhead(§10.3) < baseline_rate(B),

i.e. Def 2.1's E[L(Enc(X))] < E[L(B(X))] expanded into its component costs -- with each
term cross-referenced to its definition/bound (H(X|Y) per §2.2 conditional source coding;
L(M)/m per §10.2; verification_overhead per §10.3; baseline B per Def 2.1) and with the
§10.2 observation made explicit that the L(M)/m term is what forces the
"only the model-amortization regime can win" conclusion (single-block: L(M) dominates;
streaming m large: L(M)/m vanishes). This is an ACCOUNTING / characterization artifact,
not a new theorem and not a measurement -- it makes "the precondition" auditable.

FIX TARGET (this leaf): add a precise precondition-characterization passage (a labelled
inequality + a term-by-term accounting, mirroring §10.3's overhead bookkeeping) that
operationalizes Def 2.1 for the data-class precondition, names which structural data-class
FAMILIES the paper conjectures plausibly satisfy it (positional / geometric / multi-scale
/ non-adjacent-coupled, per §1.1) WITHOUT claiming any specific corpus does (that is the
empirical BUG-010-D), and cross-links BUG-001 (ratio claim) / BUG-002 (model-inclusive
accounting) so the cost terms reconcile with their accounting. The characterization must
stay strictly within terms already defined in the paper; if any term lacks a paper-side
definition or bound, flag it as a residual rather than inventing one.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class; the artifact is a paper-content characterization
  (inequality + accounting table), not a runtime test or a numeric measurement.
- Fully closeable in-repo: every term is sourced from existing paper definitions; no
  external data needed. (The QUESTION of which real corpora satisfy the inequality is the
  separate external leaf BUG-010-D.)
- ABOVE_EASY rather than EASY because it requires faithfully reconciling cost terms across
  §2.2 / §10.2 / §10.3 / Def 2.1 and aligning with BUG-002's model-inclusive accounting,
  not merely rewording one sentence.
- Non-goal: the wording/cross-link-only edit (BUG-010-A); the probe (BUG-010-C); the
  measurement (BUG-010-D).

## Refined Description
**Scope:** one precise precondition inequality + term-by-term accounting (each term
cross-referenced), naming candidate structural data-class families, reconciled with
BUG-001/002 accounting.  **Non-goals:** measuring which corpora qualify (BUG-010-D);
the conditional-scoping wording (BUG-010-A); the verify probe (BUG-010-C).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-010; coheres with BUG-010-A (the wording points here);
cross-checks BUG-002-B (model-inclusive accounting table) for term consistency.
**Risks / Open Questions:** the accounting must not over-claim -- listing a data-class
family as "plausibly qualifying" is a conjecture, must be tagged as such; the repair-
residual term may not have a closed paper-side bound (flag as residual if so).
**Expected Code Change:** YES (tex: new characterization passage / inequality + accounting)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | fix the exact term set and which §-anchor bounds each |
| ESTIMATION | DONE | ABOVE_EASY (multi-section reconciliation, one passage) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | choose the inequality form + accounting-table layout (mirror §10.3) |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex passage near §1.1/§8/§10.2 stating the precondition inequality + accounting |
| TESTING | PENDING | internal-consistency re-read; mechanically guarded by BUG-010-C probe |
| POST_AUDIT | PENDING | adversarial re-review per project review protocol (no over-claim of qualifying corpora) |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at section 12.7 decomposition of BUG-010;
  no code/proof change made to the paper before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- §12.7: this is one BLOCKING child leaf of BUG-010 (the characterization/accounting layer
  that turns the prose precondition into an auditable inequality).

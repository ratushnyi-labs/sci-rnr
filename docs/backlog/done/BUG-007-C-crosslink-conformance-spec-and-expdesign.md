# BUG-007-C - crosslink-conformance-spec-and-expdesign

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Thm 10.1; §10.5; Thm 10.6; §13.3) cross-referenced against tex/rnr_engineering_spec.tex (§3, §12 conformance suite) and tex/rnr_experimental_design.tex
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
> (Filed gap 007 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at §12.7 decomposition of BUG-007 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: Thm 10.1; §10.5; Thm 10.6; §13.3.

The FIX TARGET names "a reference deterministic-inference implementation + a
cross-platform bit-exactness test demonstrating the spec is satisfiable end-to-end."
That artifact's SPECIFICATION already partly exists: `rnr_engineering_spec.tex` defines
integer-arithmetic / determinism requirements (R-3.x), per-hardware-target conformance
(§4), and a conformance test suite with acceptance criteria (§12, 10,000-input test
vectors, per-target certification). The main paper defers to "a companion artifact
(§13.3)" but the pointer is generic. This leaf makes the deferral PRECISE: cross-link
Thm 10.1 / §10.5 / §13.3 to the exact conformance-suite section of the engineering spec,
and to any bit-exactness validation hypothesis in the experimental-design doc, so the
"satisfiable end-to-end" claim is scoped as CONDITIONAL-ON-CONFORMANCE with a named,
existing companion specification rather than an unspecified future deliverable.

## Refined Description
**Scope (this leaf only):** verify and, where missing, add precise cross-references so the
reader can trace the bit-exactness deferral to its specification:
1. confirm `rnr_engineering_spec.tex` §3 (integer/determinism R-3.x) and §12 (conformance
   test suite) constitute the "full specification on a chosen target" referenced by §10.1
   / §10.5 (the prose currently says "an implementer must verify against the full
   specification on a chosen target"); make the §13.3 deferral cite that spec section
   explicitly (not just "companion artifact");
2. check `rnr_experimental_design.tex` for a cross-platform bit-exactness / reproducibility
   validation item; if one exists, cross-link it from §13.3 (empirical validation of
   satisfiability end-to-end); if none exists, NOTE the absence in this leaf's findings so
   it can be routed to the experimental-design doc owner (not fixed here);
3. ensure the accounting is consistent: Thm 10.1 (proven sufficient condition) → spec §3
   requirements → spec §12 conformance test (per-target certification) → end-to-end
   satisfiability (empirical, external = BUG-007-D).

**Fix target (this leaf):** precise, bidirectional cross-links (paper ↔ engineering spec
↔ experimental design) so the deferred FIX TARGET has a named home, and the paper's claim
is scoped as "specified and conditional on certified conformance."

**Non-goals:** authoring NEW conformance requirements or NEW experimental-design
hypotheses; writing verify scripts (BUG-007-B); the external implementation/test itself
(BUG-007-D); main-paper claim-honesty wording, which is BUG-007-A (this leaf is the
SPEC/cross-doc accounting layer, A is the in-paragraph claim layer — kept distinct).

**Closeable in-repo NOW:** YES — documentation cross-referencing across the three existing
tex documents; no external work. (Edits land only in tex/rnr_coding.tex per the
constraint that BUG-007-* touch the paper's claims; spec/exp-design edits, if any are
needed beyond confirming pointers, are recorded as findings for their doc owners.)

**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-007; coordinates with BUG-007-A (claim wording) — A scopes
the in-paragraph claim, C wires the deferral to the named spec; do A first if editing the
same paragraph to avoid churn.
**Risks / Open Questions:** the engineering spec already covers most of the "full
specification"; the residual is whether the experimental-design doc has a bit-exactness
validation hypothesis at all — if not, that is a finding, not an in-leaf fix.
**Expected Code Change:** YES (tex cross-refs; mostly main paper §13.3 / §10.5)
**Documentation Recovery Required:** NO

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class defect; "fix" = make the deferral precise and traceable
  so the claim boundary is auditable.
- Priority MEDIUM (refined down from parent HIGH): this is an accounting/traceability
  improvement; the load-bearing honesty edit is BUG-007-A and the math probe is BUG-007-B.

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = cross-doc traceability of the deferral |
| ESTIMATION | DONE | EASY (cross-reference audit + pointer edits) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | map paper claim → spec §3/§12 → exp-design item |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | add precise §13.3 / §10.5 cross-refs in main paper |
| TESTING | N/A | no executable artifact |
| POST_AUDIT | PENDING | adversarial: is the end-to-end claim now fully traceable? |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7: child leaf of BUG-007; created fresh this turn (Recovered: N/A).
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Honest-step principle: this leaf closes the cross-doc accounting NOW; the actual
  end-to-end demonstration remains external (BUG-007-D).


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-007 A/B/C (commit pending). A: §10.5/§10.8/§13.3 prose tightened to the proven sufficient-condition + explicit external deferral (no end-to-end-demonstration overclaim). B: scripts/verify/thm_10_1_integer_reorder_bitexact.py (no-overflow bit-width bound + integer reduction-order invariance + saturation/wrap negative control; OVERALL -> PASS) + CI job. C: §13.3 deferral crosslinked to engineering-spec conformance suite (R-3.x/§5/§12 R-12.2.3/§14) + experimental-design H8/H9/H14 (all anchors verified real). Hostile-referee POST_AUDIT clean.

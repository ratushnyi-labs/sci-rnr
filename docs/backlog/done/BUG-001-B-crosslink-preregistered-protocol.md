# BUG-001-B - crosslink-preregistered-protocol

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§1.1; §11) <-> tex/rnr_experimental_design.tex (H1--H7, H15)
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
> (Filed gap 001 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at §12.7 decomposition of BUG-001 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: §1.1 ("the experimental-design companion document specifies the pre-registered
measurement protocol"); §11 intro and §§11.1--11.5; cross-document target
`tex/rnr_experimental_design.tex` (H1--H7 quantifying/confirming hypotheses; H15
the enwik9 end-to-end test; §§3/5/6 corpora, metrics, analysis, pre-registration).

The empirical ratio claim is the load-bearing main-idea promise, and the companion
document already operationalizes it as a pre-registered protocol (falsifiable
bit-per-byte targets, named baselines bgzip/zstd-seekable/LZMA/NNCP, corpora,
multiple-comparison correction). What is MISSING in-repo is an explicit, audited
MAPPING that ties each §11 prediction to the specific hypothesis (and falsification
criterion) in the companion that would confirm-or-refute it, so a reader can trace
"central value proposition -> the exact experiment that settles it." Today §1.1 names
the companion only generically and §11 lists predictions without per-prediction
companion-hypothesis cross-references.

FIX TARGET (this leaf): (1) verify the claim<->hypothesis correspondence is COMPLETE
(every §11.1--11.6 prediction has a matching companion hypothesis with a stated
falsification rule; flag any prediction with no companion counterpart, and any
companion hypothesis with no §11 prediction); (2) add precise cross-references in §11
(and the §1.1 sentence) pointing each prediction at its companion hypothesis number;
(3) make the "beats existing dictionary-based SEEKABLE formats" claim specifically
checkable -- confirm the companion's baseline set and seekable-format comparison
actually cover bgzip / zstd-seekable, the formats named in the Abstract/§1.1, and add
the cross-reference (or flag the residual if the seekable comparison is absent).

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class; "fix" = scope/cross-link the claim so the in-repo
  artifacts honestly reflect what is established vs. deferred.
- Closeable in-repo: the cross-reference audit and the §11/§1.1 reference edits live
  entirely in the two committed tex documents. The leaf produces a verified mapping
  table and the cross-references; it does NOT run any experiment.
- Edits to `tex/rnr_experimental_design.tex` are limited to cross-reference / pointer
  additions if the companion lacks a back-pointer; the parent task's constraint
  forbidding edits to `tex/rnr_coding.tex` belongs to THIS intake turn only -- the
  eventual BUILD of this leaf does edit `tex/rnr_coding.tex` cross-references.
- Non-goal: changing any predicted number; non-goal: §11.7 (theorem-backed).

## Refined Description
**Scope:** complete + verified §11-prediction <-> companion-hypothesis mapping, and
the cross-references that make the claim traceable to its pre-registered test.
**Non-goals:** running the protocol (BUG-001-D); model accounting (BUG-002).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-001; soft-ordered after BUG-001-A (scoping wording
should land first so cross-references attach to the conditionalized claims)
**Risks / Open Questions:** a prediction with no companion hypothesis, or a "seekable
format" comparison absent from the companion baselines, would be a residual to flag
(possibly feeding BUG-001-D), not silently paper over.
**Expected Code Change:** YES (tex cross-references in §11/§1.1; possibly companion back-pointer)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | confirm which companion hypotheses are in scope (H1--H7, H15) |
| ESTIMATION | DONE | ABOVE_EASY (mapping audit + cross-ref edits across two docs) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | mapping-table format; cross-reference convention |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | §11/§1.1 cross-references; optional companion back-pointer |
| TESTING | PENDING | mapping completeness check (also asserted by BUG-001-C probe) |
| POST_AUDIT | PENDING | adversarial re-review per project review protocol |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at §12.7 decomposition of BUG-001;
  no code/proof change made before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- §12.7: BLOCKING child leaf of BUG-001 (the cross-link/traceability layer).


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-001 A/B/C (commit pending). A: §11.1-11.6 per-prediction conditional scoping ('falsifiable empirical conjecture, conditional and pending measurement, not a theorem') + canonical §11-intro paragraph; Abstract/§1.1 already scoped (left intact). B: §11.x<->companion-hypothesis mapping added (H1-H7, H15-H16) incl. seekable-baseline coverage note (bgzip/zstd-seekable = H10/H11 random-access TIME baselines; squashfs+zstd = H12 ratio baseline). C: scripts/verify/bug_001_ratio_claim_conditional.py (inv1 conditional-scoping + inv2 no-orphan mapping, non-vacuity self-test; per-paragraph locality + dominate-needs-baseline matcher hardened by supervisor; OVERALL -> PASS) + CI job. Hostile-referee POST_AUDIT clean. HONEST RESIDUAL (feeds D): §11.5 repetitive-ratio-vs-zstd/LZMA target + the beats-seekable value-prop have only PARTIAL pre-registered coverage (H7 = image bits-back; no dedicated bgzip/zstd-seekable RATIO hypothesis) -- flagged in-text, no H-number invented.

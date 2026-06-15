# BUG-010-D - which-corpora-qualify-measurement-external

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** BLOCKED
**Complexity:** VERY_HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** EXTERNAL (reference implementation + corpus runs); tex/rnr_coding.tex (Abstract; §1.1; Def 2.1; Theorem 8.1; §10.2) and tex/rnr_experimental_design.tex (per-Type ratio + amortization hypotheses) as the protocol anchor
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
Source: Abstract / §1.1 (the data-class precondition); the parent BUG-010 observation that
"Which real corpora actually qualify -- the data-class precondition -- is itself the
unmeasured empirical question, and only the model-amortization regime (§10.2) can win";
companion `tex/rnr_experimental_design.tex` (per-Type ratio targets and amortization-
regime hypotheses).

This is the IRREDUCIBLE external half of BUG-010. The characterization leaf (BUG-010-B)
writes the precondition inequality precisely; this leaf MEASURES which real corpora
actually satisfy it. Confirming "this corpus qualifies" requires evaluating the component
terms of the precondition inequality on real data -- H(X|Y) under a decoder-reproducible
neural side-information model, the realized repair-residual cost, L(M)/m at the operating
block count m (the §10.2 amortization regime), and the verification overhead -- and
comparing against the baseline rate. That needs a working reference RNR coder, the named
corpora, and compute, and must be run in the model-amortized regime (large m) where §10.2
says a win is even possible. It CANNOT be closed inside this docs/theory repo.

FIX TARGET (this leaf): execute the precondition-evaluation on the pre-registered corpora
(per the companion's per-Type ratio + amortization hypotheses), report which data-class
families measurably satisfy the BUG-010-B inequality (and at what m / amortization regime),
and feed the outcome back to retire/condition the §1.1 claim that structural data classes
plausibly qualify (confirm specific families, refute -> downgrade the conjecture, or
quantify the margin). Until then the in-repo posture is the scoped-conditional claim from
BUG-010-A plus the auditable precondition inequality from BUG-010-B. Coordinate with
BUG-001-D (ratio measurement) / BUG-002-D (enwik9 model-inclusive measurement) -- this
leaf reads the precondition-term breakdown those runs would produce; do not duplicate the
campaign.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class file, but the WORK is external: this leaf is flagged
  closeable_in_repo = NO precisely because measuring per-corpus precondition terms needs a
  reference implementation + corpora + compute that do not exist in this repo.
- This leaf is the genuinely-hard residual the section 12.7 split isolates so the rest of
  BUG-010 (A/B/C) can be honestly closed now. It is intentionally NOT split further:
  building/obtaining a reference coder and running per-corpus precondition evaluation in
  the amortized regime is a single large external program, not a stack of EASY in-repo
  steps. It overlaps the BUG-001-D / BUG-002-D campaign and should ride on it rather than
  spawn a parallel one.
- Stays open (or moves to a tracking/external-work state) after A/B/C land; it is the
  documented residual, clearly flagged as needing out-of-repo work.
- Non-goal: any in-repo wording / inequality / probe (those are A/B/C).

## Refined Description
**Scope:** measure per-corpus satisfaction of the BUG-010-B precondition inequality in the
§10.2 amortized regime; reconcile §1.1's "plausibly-qualifying structural families"
conjecture with measurement.  **Non-goals:** in-repo scoping (BUG-010-A), the inequality
itself (BUG-010-B), the verify probe (BUG-010-C).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-010; reads the precondition-term breakdown produced by
BUG-001-D / BUG-002-D's runs; consumes the BUG-010-B inequality as the evaluation target;
does not block A/B/C.
**Risks / Open Questions:** requires a reference implementation, corpus access, compute,
and runs in the large-m amortized regime; results may show NO common corpus family
satisfies the precondition outside narrow regimes -- that is an acceptable, honest outcome
forcing a §1.1 conjecture downgrade, not a failure of this leaf.
**Expected Code Change:** NO (in THIS repo); external implementation + data are required
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | freeze the corpus set + the amortization regime (m) to evaluate at |
| ESTIMATION | DONE | VERY_HARD; closeable_in_repo = NO (external work) |
| DECOMPOSITION | N/A | atomic residual leaf (single large external program; rides on BUG-001-D/002-D, not further splittable into EASY in-repo steps) |
| DESIGN | PENDING | precondition-term measurement plan layered on the companion protocol |
| FRONTEND | N/A | no UI |
| BACKEND | N/A (this repo) | measurement is external; no tex/proof change here |
| TESTING | PENDING (external) | measured per-term breakdown per (corpus, m); inequality satisfied yes/no |
| POST_AUDIT | PENDING | reconcile measured qualifying-families with §1.1 conjecture; adversarial re-review |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD (external; large)
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at section 12.7 decomposition of BUG-010;
  no code/proof change made before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- §12.7: BLOCKING child leaf of BUG-010 -- the isolated external/research residual.
  Flagged closeable_in_repo = NO so the parent's complexity collapses to the shippable
  A/B/C core plus this clearly-marked deferred measurement.

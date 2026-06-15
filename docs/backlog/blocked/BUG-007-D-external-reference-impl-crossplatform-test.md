# BUG-007-D - external-reference-impl-crossplatform-test

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** BLOCKED (external)
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** EXTERNAL artifact (reference codec + cross-platform CI); specified by tex/rnr_engineering_spec.tex §3, §12; deferred from tex/rnr_coding.tex §13.3
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
Source: Thm 10.1; §13.3; engineering spec §3, §12.

This leaf carries the part of the FIX TARGET that genuinely CANNOT be closed in this
theory repo: build "a reference deterministic-inference implementation + a cross-platform
bit-exactness test demonstrating the spec is satisfiable end-to-end." This is the
engineering deliverable the paper explicitly defers (§13.3: "software engineering work,
not mathematics ... a separate code release accompanied by a conformance specification").
It requires real hardware targets (e.g., scalar CPU, CPU SIMD, GPU integer kernels, an
NPU such as Apple Neural Engine), an actual integer-arithmetic inference engine
conforming to engineering-spec R-3.x, and execution of the §12 conformance test suite
(≥10,000-input test vectors, per-target certification) to demonstrate byte-identity of
`M(C)` across platforms.

## Refined Description
**Scope (this leaf only):** the EXTERNAL artifact: a reference codec implementing
deterministic integer inference per engineering-spec §3 (R-3.x), plus a multi-target
conformance run per engineering-spec §12 producing a cross-platform bit-exactness result
(byte-identical predictor output / archive across at least two distinct conforming
targets). Deliverable is a separate code release + CI, NOT a change to the paper.

**Fix target (this leaf):** demonstrated end-to-end satisfiability of Thm 10.1 / Thm 10.6
conformance conditions on real hardware — i.e., the empirical confirmation that the
specification is realizable.

**Non-goals:** any change inside `tex/rnr_coding.tex` (the paper-side scoping is BUG-007-A
and BUG-007-C; the in-repo math probe is BUG-007-B). This leaf does not modify the paper.

**Closeable in-repo NOW:** NO — this requires EXTERNAL work (a reference implementation
and real cross-platform hardware testing). It cannot be discharged inside this theory
repo. Flagged explicitly so the parent's shippable core (A + B + C) is separated from this
external residual. The repo-side honest statement (scope the claim, cross-link the spec)
is already carried by BUG-007-A / BUG-007-C; this leaf is the measurement/implementation
itself and stays OPEN-EXTERNAL.

**Impacted UCs:** N/A (theory paper; deliverable is a separate code release)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-007; specification provided by engineering spec §3/§12
(traceability ensured by BUG-007-C). Logically downstream of A/B/C but executed outside
this repo.
**Risks / Open Questions:** whether any specific named hardware class ACTUALLY conforms is
an engineering claim about its instruction set, kernel library, and serialization
(per Thm 10.6's conditional caveat) — discovering a non-conforming target is a valid,
informative outcome, not a failure of the theory.
**Expected Code Change:** NO (in this repo) — external code release.
**Documentation Recovery Required:** NO

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class theory repo has no runtime product; this deliverable
  lives in a SEPARATE engineering artifact. The theory paper already (correctly) defers it
  to §13.3; this backlog leaf records the external work item so the deferral is tracked,
  not silently dropped.
- Priority MEDIUM: necessary for a DEPLOYABLE system but not for the paper's validity
  (the paper claims a sufficient condition, proven; deployment is out of paper scope).

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = external reference impl + cross-platform test |
| ESTIMATION | DONE | ABOVE_EASY as execution (spec exists); EXTERNAL, not in-repo |
| DECOMPOSITION | N/A | atomic leaf (further engineering breakdown belongs to the code release) |
| DESIGN | DEFERRED | governed by engineering spec §3/§12 (external) |
| FRONTEND | N/A | no UI |
| BACKEND | N/A (external) | reference codec lives outside this repo |
| TESTING | DEFERRED (external) | §12 conformance suite on real targets |
| POST_AUDIT | DEFERRED | conformance certification per target |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD (external)
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7: child leaf of BUG-007; created fresh this turn (Recovered: N/A).
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Honest-step principle: this is the explicitly-flagged EXTERNAL residual
  (closeable_in_repo = NO). Separating it lets the in-repo core (BUG-007-A scoping,
  BUG-007-B math probe, BUG-007-C cross-links) close now while this stays openly deferred.


## Status (2026-06-15)
IRREDUCIBLE EXTERNAL residual: reference deterministic-integer inference codec (engineering-spec §3) + §12 conformance suite across real hardware for end-to-end cross-platform byte-identity. Not closeable in this theory repo; the in-repo prose now correctly defers it (BUG-007-A). Parked in blocked/.

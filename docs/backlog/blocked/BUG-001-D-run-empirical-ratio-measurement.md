# BUG-001-D - run-empirical-ratio-measurement

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** BLOCKED (external)
**Complexity:** VERY_HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** EXTERNAL (reference implementation + corpus runs); tex/rnr_coding.tex §11 / tex/rnr_experimental_design.tex H1--H7,H15 as the protocol
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
Source: Abstract; §1.1; §11 ("Empirical validation against NNCP, PixelCNN++, LZMA at
maximum preset, and zstd at maximum level on the corpora named below is the natural
next step"); companion protocol H1--H7 (per-Type ratio targets) and H15 (enwik9
end-to-end).

This is the IRREDUCIBLE external half of BUG-001: actually MEASURING the central value
proposition. Closing the empirical gap (as opposed to scoping it, BUG-001-A/B/C)
requires building or obtaining a reference RNR coder, running it and the named
baselines (NNCP, PixelCNN++, LZMA -9, zstd -22, and the SEEKABLE formats bgzip /
zstd-seekable) on the pre-registered corpora, and reporting measured bits-per-byte
with the companion's statistical methodology (5 seeds, effect sizes, multiple-comparison
correction). Only a measured datum can confirm "RNR's side-information+repair recasting
beats existing seekable formats on ratio." This CANNOT be closed inside this docs/theory
repo: it needs a working implementation, corpus data, and compute.

FIX TARGET (this leaf): execute the pre-registered protocol (H1--H7, H15) and produce
measured ratio results, then feed the outcome back to retire/condition the §11
predictions (confirm, refute -> downgrade, or quantify the off-by-X). Until then the
in-repo posture is the scoped-conditional claim from BUG-001-A and the traceable
protocol pointer from BUG-001-B.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class file, but the WORK is external: this leaf is flagged
  closeable_in_repo = NO precisely because the measurement needs a reference
  implementation + corpora + compute that do not exist in this repo.
- This leaf is the genuinely-hard residual the §12.7 split isolates so the rest of
  BUG-001 (A/B/C) can be honestly closed now. It is intentionally NOT split further:
  building a reference coder and running a 14+ hypothesis pre-registered campaign is a
  single large external program of work, not a stack of EASY in-repo steps.
- Stays open (or moves to a tracking/external-work state) after A/B/C land; it is the
  documented residual, clearly flagged as needing out-of-repo work.
- Non-goal: any in-repo wording/probe (those are A/B/C).

## Refined Description
**Scope:** run the pre-registered ratio protocol and report measured results;
reconcile §11 predictions with measurement.  **Non-goals:** in-repo scoping/cross-link/
probe (BUG-001-A/B/C); §11.7 theorem verification (already done in scripts/verify).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-001; protocol traceability from BUG-001-B (so the run
targets the exact pre-registered hypotheses); does not block A/B/C.
**Risks / Open Questions:** requires a reference implementation (engineering-spec
realization), corpus access (enwik9 etc.), compute, and seekable-format baselines;
results may refute predictions and force a §11 downgrade -- that is an acceptable,
honest outcome, not a failure of this leaf.
**Expected Code Change:** NO (in THIS repo); external implementation + data are required
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | freeze protocol scope (H1--H7, H15) and baseline set |
| ESTIMATION | DONE | VERY_HARD; closeable_in_repo = NO (external work) |
| DECOMPOSITION | N/A | atomic residual leaf (single large external program; not further splittable into EASY in-repo steps) |
| DESIGN | PENDING | reference-implementation + measurement plan (engineering spec + companion) |
| FRONTEND | N/A | no UI |
| BACKEND | N/A (this repo) | measurement is external; no tex/proof change here |
| TESTING | PENDING (external) | measured bpb per (coder, corpus) with 5-seed methodology |
| POST_AUDIT | PENDING | reconcile measured results with §11 predictions; adversarial re-review |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD (external; large)
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at §12.7 decomposition of BUG-001;
  no code/proof change made before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- §12.7: BLOCKING child leaf of BUG-001 -- the isolated external/research residual.
  Flagged closeable_in_repo = NO so the parent's complexity collapses to the
  shippable A/B/C core plus this clearly-marked deferred measurement.


## Status (2026-06-15)
IRREDUCIBLE EXTERNAL residual: run the pre-registered §11/experimental-design protocol (reference RNR impl + baselines incl. bgzip/zstd-seekable + squashfs+zstd + corpora + compute) and report measured bpb; reconcile §11. Not closeable in this theory repo. ALSO ABSORBS the in-repo-flagged coverage gap: add a dedicated quantifying hypothesis for the §11.5 repetitive-structured ratio-vs-dictionary / beats-bgzip-zstd-seekable RATIO target when the measurement is designed. Parked in blocked/.

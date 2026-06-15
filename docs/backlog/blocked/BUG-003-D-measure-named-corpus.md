# BUG-003-D - measure-named-corpus

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** BLOCKED (external)
**Complexity:** HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Abstract; §10.7; §11) -- claim consumes the result;
tex/rnr_experimental_design.tex (H10, read-only) -- the measurement protocol
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
> (Filed gap 003 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

(Child sub-task created at section 12.7 decomposition of BUG-003 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: Abstract; §10.7; §11.

This is the residual external leaf of BUG-003: the "measured value on a specific corpus" fix
target. The Abstract conjectures ~1% overhead for natural data "at typical parameters" but
"a measured value on a specific corpus is [not] provided here." A real measurement requires a
working RNR encoder with sync points (a reference implementation -- see rnr_engineering_spec.tex)
run on a NAMED corpus (e.g. enwik9, or a Linux distribution ISO as in experimental-design
H12), with a `W`-bounded predictor at `K = 64 KiB`, comparing archive size with vs. without
sync points (the H10(c) "ratio overhead ... below 2 percent" hypothesis). This CANNOT be
closed inside this docs/theory repo: it needs an implementation, a corpus download, and a
compute run -- all external work. The in-repo half of this fix (precisely scoping the claim
as a conjecture and cross-linking the H10 measurement protocol) is handled by BUG-003-A; the
in-repo executable heuristic support on synthetic data is BUG-003-C. What remains HERE is only
the genuinely-external measurement itself.

**FIX TARGET (this leaf):** obtain a real measured sync-point ratio-overhead datum on a named
corpus at typical (K, W) parameters, per experimental-design H10(c), and report it back so the
Abstract can either (i) cite the measured value, downgrading the conjecture to an empirical
fact with a corpus name, or (ii) record a falsification if overhead exceeds the H10
falsification threshold (>10% at K = 64 KiB). Until then the conjecture stands scoped (via
BUG-003-A) and synthetic-supported (via BUG-003-C).

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class repo, but THIS leaf's work product (a measurement) lives
  OUTSIDE the repo: it needs the reference implementation (rnr_engineering_spec.tex) and a
  named corpus. closeable_in_repo = NO.
- The measurement protocol already exists as experimental-design H10 (tex/rnr_experimental_design.tex
  lines ~218--248); this leaf executes that protocol externally, it does not invent it.
- Once measured, the only in-repo follow-up is a one-line Abstract/§10.7/§11 update citing the
  datum -- that follow-up edit can be a trivial addendum tracked under BUG-003-A's umbrella or
  filed separately; it is NOT the hard part. The hard part is the external run.
- This leaf is intentionally left at HARD: it is the one genuinely-external residual that
  cannot be split further without an implementation in hand.

## Refined Description
**Scope:** a real measured sync-point ratio-overhead value on a NAMED corpus at typical
(K, W) per H10(c), executed against a reference RNR implementation.  **Non-goals:** the
in-repo scoping (BUG-003-A), the assumption-graded bound ladder (BUG-003-B), the synthetic
probe (BUG-003-C).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-003; experimental-design H10 protocol; a reference implementation
(rnr_engineering_spec.tex) -- EXTERNAL, not present in this repo.
**Risks / Open Questions:** requires building/obtaining an RNR encoder with sync points; result
may falsify the ~1% conjecture (then BUG-003-A's scoped claim must be revised to the measured
band). Blocked on external implementation + compute; not schedulable inside this repo's turn.
**Expected Code Change:** NONE in this repo until a datum exists; then a one-line claim update.
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | external measurement leaf; closeable_in_repo = NO |
| ESTIMATION | DONE | HARD (needs reference impl + named corpus + compute run; cannot split further now) |
| DECOMPOSITION | N/A | atomic leaf (the only honest external residual) |
| DESIGN | N/A | protocol already specified as experimental-design H10 |
| FRONTEND | N/A | no UI |
| BACKEND | BLOCKED | external: requires RNR reference implementation (rnr_engineering_spec.tex) |
| TESTING | BLOCKED | external: run on named corpus (enwik9 / Linux ISO), measure size-with vs size-without sync |
| POST_AUDIT | PENDING | when a datum exists: feed back to BUG-003-A claim wording (cite or falsify) |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD (dominated by external implementation, not estimable in-repo)
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh as a decomposition child of BUG-003.
- Class: GOVERNANCE/DOCUMENTATION; the measurement work product is external to this repo by
  nature (no runtime product lives here). closeable_in_repo = NO is the honest status.
- This leaf isolates the genuinely-external residual so the parent's in-repo core (A, B, C)
  is shippable without waiting on an implementation.


## Status (2026-06-15)
IRREDUCIBLE EXTERNAL residual: a real measured sync-point ratio-overhead datum on a NAMED corpus at typical (K,W) per H10(c), needing a reference RNR impl + compute run. Not closeable in this theory repo; in-repo follow-up is a one-line claim update once a datum exists. Parked in blocked/. (Also: the fully-general TIGHT per-sync bound is a separate OPEN theory residual, noted in Corollary 10.2a.)

# BUG-008 - random-access-parameter-conditional

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE (in-repo); D external-blocked
**Complexity:** MEDIUM
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Thm 10.3; Rem 10.3a)
**Business Rule Conflict Check:** NOT_CHECKED
**Conflict Resolution Reference:** N/A
**Atomic Check:** ATOMIC
**Budget Forecast:** TBD
**Priority:** MEDIUM
**Lawbook Version (intake):** 0.44.0
**Applicable Lawbook Version:** TBD
**Created:** 2026-06-14  **Updated:** 2026-06-15
**Recovered:** YES

## Raw Request
> mark each gap as bug and follow the global rules
>
> (Filed gap 008 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

## Gap (defect statement; document-grounded)
Source: Thm 10.3; Rem 10.3a.

The polylog random-access claim holds only once sync-spacing K, read window k, predictor cost(M), and log L(R) are treated as fixed/polylog -- but cost(M) (a real neural eval) is non-universal and possibly large, and whole-block Type-II/III-B partition-rank modes get only sub-block-granular access (Rem 10.3a). FIX TARGET: state the cost(M) dependence prominently and quantify the access granularity per mode; ideally a measured seek-latency.

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
| CLARIFICATION | PENDING | scope the claim vs close the gap |
| ESTIMATION | PENDING | set Complexity/Budget |
| DECOMPOSITION | DONE | split into BUG-008-A, BUG-008-B, BUG-008-C, BUG-008-D |
| DESIGN | PENDING | proof strategy or experiment design |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex edit / proof / probe (or N/A if empirical-only) |
| TESTING | PENDING | scripts/verify probe or measured datum |
| POST_AUDIT | PENDING | adversarial re-review per project review protocol |

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
Split into four BLOCKING child leaves. The §12.7 split isolates the part that can be
HONESTLY CLOSED NOW IN-REPO (surface the `cost(M)` dependence, tabulate per-mode access
granularity, guard the cost accounting with a verify probe) from the IRREDUCIBLE external
measurement (a real wall-clock seek-latency). After A/B/C land, the parent's in-repo posture
is a precisely-scoped random-access claim — `cost(M)` shown as a non-universal dominant
neural-eval factor, per-mode granularity consolidated, and the cost formula guarded by an
executable probe; D is the clearly-flagged deferred residual needing out-of-repo work.

- **BUG-008-A — scope-cost-m-dependence**: surface the non-universal, dominant `cost(M)`
  dependence of the $O(\log L(R)+(K+k)\cdot\mathrm{cost}(M))$ bound at the Abstract / §1.1
  headline sites and the Theorem 10.3 scope sentence, in one canonical phrasing aligned with
  Lemma 5.1b. Complexity EASY; closeable-in-repo YES.
- **BUG-008-B — per-mode-granularity-table**: add one consolidated table near Remark 10.3a
  giving, per coding mode, the access granularity (byte vs sub-block), the per-query cost
  formula, the wasted-decode overhead, and source cross-refs — faithful to Thm 10.3 / Rem
  10.3a / §6.3 / §10.7. Complexity ABOVE_EASY; closeable-in-repo YES.
- **BUG-008-C — randomaccess-cost-verify-probe**: add a `scripts/verify/` PASS/FAIL probe
  (+ CI job) modelling the seek + (K+k) decode-window accounting and asserting the Thm 10.3 /
  Rem 10.3a cost relations, regime boundaries ($K=\Theta(|X|)$ degeneracy; simultaneous-polylog
  impossibility), and sub-block-granular cost. Complexity ABOVE_EASY; closeable-in-repo YES.
- **BUG-008-D — measure-seek-latency-external**: execute the pre-registered seek-latency
  protocol (companion §2.8 H10/H11, optionally §2.9 H12/H13) with a reference coder across
  archive sizes and hardware paths, then confirm/condition the Rem 10.3a 650 ms and §11
  20--200 ms / sub-ms estimates. Complexity HARD; closeable-in-repo NO (external work;
  isolated measurement residual).


## Resolution (2026-06-15)
A/B/C closed in-repo (-> done/): cost(M) non-universal/dominant scoping at the headline
random-access sites + Thm 10.3 (polylog in archive params, cost(M) fixed, not wall-clock);
6-row per-mode access-granularity table near Rem 10.3a (byte vs sub-block, per-query cost,
wasted decode, sources); cost-accounting probe bug_008_random_access_cost.py (decode-window,
K=Theta(|X|) degeneracy, simultaneous-polylog impossibility, sub-block cost) + CI job.
POST_AUDIT (hostile referee) clean. D (measured seek latency) -> blocked/ (external).
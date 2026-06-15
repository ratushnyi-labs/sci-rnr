# BUG-008-B - per-mode-granularity-table

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Theorem 10.3 scope clause ~ll.15972--15996; Remark 10.3a ~ll.16070--16121; §6.3 Type-III-B strategies (i)/(ii)/(iii); §10.7 "Type-II semi-non-factorized" trade-off bullet)
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
> (Filed gap 008 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at section 12.7 decomposition of BUG-008 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: Thm 10.3 (per-position sequential modes, byte-granular); Rem 10.3a (Type-II
partition modes, sub-block-granular); §6.3 Type-III-B strategies; §10.7 trade-off bullet.

The closeable-in-repo half of BUG-008's second FIX TARGET: "quantify the access
granularity per mode." The information IS present in the paper, but it is dispersed:
Theorem 10.3 lists which modes are byte-granular (Type-I, per-byte sequential Type-III-B
strategy (iii), tile-aligned Type-III-A) and which are excluded (Type-II enumerative rank,
Type-III-B strategies (i)/(ii)); Remark 10.3a derives the sub-block-granular cost
$O(K\cdot\mathrm{unrank}(K))$ per requested sub-block for the partition-rank modes. A reader
must reconstruct the per-mode access-granularity map from prose scattered across Thm 10.3,
Rem 10.3a, §6.3, and §10.7.

FIX TARGET (this leaf): add ONE compact accounting table (or a tightly-structured
enumerated paragraph) in or adjacent to Remark 10.3a that, per coding mode, states:
(i) the access granularity (byte-granular vs sub-block-granular), (ii) the per-query
random-access cost formula in that mode ($O(\log L(R)+(K+k)\cdot\mathrm{cost}(M))$ for
per-position sequential; $O(\log L(R)+K\cdot\mathrm{unrank}(K))$ per affected sub-block for
partition-rank), (iii) the wasted-decode overhead (sub-block modes decode and discard up to
a full sub-block per request; byte-granular decode at most $K+k$ positions of which $\le K$
are warmup), and (iv) the cross-reference to where each mode is defined (§6.1--§6.3) and
analyzed. Every row must be a faithful restatement of an existing result — no new claim. The
table is the single canonical "access-granularity per mode" object the Abstract / §1.1 can
point to instead of re-deriving the split inline.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class defect in `tex/rnr_coding.tex` (no runtime product).
  "Fix" here = make the already-established per-mode granularity explicit and consolidated,
  matching what the theorems/remarks already prove.
- Fully closeable in-repo: a consolidation/restatement edit (a table plus cross-refs). No
  new theorem, no proof, no measurement.
- ABOVE_EASY (not EASY): requires faithfully reconciling four scattered sources (Thm 10.3
  inclusion/exclusion list, Rem 10.3a cost, §6.3 strategy taxonomy, §10.7 trade-off bullet)
  and getting the per-mode cost formulas and $S$ values exactly right, without introducing a
  granularity claim the proofs do not support.
- Non-goal: the `cost(M)` headline-scoping wording (BUG-008-A); the verify probe (BUG-008-C);
  the measured seek-latency (BUG-008-D).

## Refined Description
**Scope:** one consolidated per-mode access-granularity + per-query-cost table near
Remark 10.3a, faithful to Thm 10.3 / Rem 10.3a / §6.3 / §10.7.  **Non-goals:** cost(M)
headline wording (A); probe (C); measurement (D); any new granularity result.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-008; coheres with BUG-008-A's canonical phrasing (the table
should use the same cost-formula notation A standardizes).
**Risks / Open Questions:** must not over-claim — e.g. Type-III-A is byte-granular ONLY
under the tile/sync alignment condition Thm 10.3 states; the table must carry that caveat.
The $S$ (symbols-per-byte) value differs per mode ($S=1$ Type-I/direct byte, $S=2$ per-byte
Type-III-B) and must match Definition 10.3 / §10.7.
**Expected Code Change:** YES (tex table + cross-refs)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | enumerate the exact mode set and confirm each row against its source |
| ESTIMATION | DONE | ABOVE_EASY (consolidate 4 scattered sources into one faithful table) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | choose table columns: mode / granularity / per-query cost / wasted-decode / source-ref |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex table near Rem 10.3a + cross-references from Abstract/§1.1 |
| TESTING | PENDING | manual re-read; row-consistency invariant covered by BUG-008-C probe |
| POST_AUDIT | PENDING | adversarial re-review (no over-claimed granularity row) per project review protocol |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at §12.7 decomposition of BUG-008;
  no code/proof change made to the paper before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- §12.7: BLOCKING child leaf of BUG-008 (the closeable-now per-mode granularity-accounting layer).

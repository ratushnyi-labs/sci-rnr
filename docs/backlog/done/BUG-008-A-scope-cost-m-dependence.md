# BUG-008-A - scope-cost-m-dependence

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Abstract random-access sentence ~ll.209--222; §1.1 contributions bullet ~ll.610--667; Theorem 10.3 statement ~ll.15972--16026; Remark 10.3a tail ~ll.16101--16121)
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
Source: Thm 10.3; Rem 10.3a; the Abstract and §1.1 restatements of the polylog
random-access claim.

The closeable-in-repo half of BUG-008's first FIX TARGET: make the
`cost(M)`-dependence of the random-access bound PROMINENT and unambiguous wherever the
"polylogarithmic in archive size" claim appears. The cost formula
$O(\log L(R) + (K+k)\cdot\mathrm{cost}(M))$ is correct, and Theorem 10.3 already carries
the conditional clause "polylogarithmic in archive size $|X|$ only when $K$, $k$, and
$\mathrm{cost}(M)$ are held fixed or grow at most polylogarithmically in $|X|$". But
$\mathrm{cost}(M)$ is a REAL neural-predictor evaluation: it is non-universal (depends on
the architecture and $W$), it is not asymptotically negligible (a windowed transformer pays
$\Theta(W)=\Theta(\log N)$ per token; Lemma 5.1b already records this), and it is the term
that actually dominates the per-query wall-clock cost (the $(K+k)\cdot\mathrm{cost}(M)$
decode-window term swamps $\log L(R)$ in any realistic regime). In the Abstract (ll.211--221)
and the §1.1 contributions bullet (ll.610--622) `cost(M)` appears inside the big-O but the
"polylog" headline is stated first and the conditionality on `cost(M)` reads as a fixed
constant rather than as a non-universal, possibly-large neural-eval factor.

FIX TARGET (this leaf): edit the Abstract random-access sentence, the §1.1 contributions
bullet, and (if needed for emphasis) one clarifying clause in the Theorem 10.3 statement or
its trailing scope sentence, so that the `cost(M)` dependence is surfaced explicitly: state
once, in a canonical reusable phrasing, that (a) `cost(M)` is a per-position neural-predictor
evaluation that is architecture-dependent and not universal; (b) it is the dominant term of
the per-query cost in any realistic regime; (c) "polylog in $|X|$" is a statement about the
DEPENDENCE ON ARCHIVE SIZE holding $\mathrm{cost}(M)$ (and $K$, $k$, $\log L(R)$) fixed or
polylog, NOT a claim that the per-query cost is small in absolute terms. No new theorem, no
measurement, no change to the cost formula itself (which is already correct).

## Clarification / Assumptions
- This is a GOVERNANCE/DOCUMENTATION-class defect in the theory paper
  `tex/rnr_coding.tex` (no runtime product). "Fix" means here: scope the paper's claim to
  match what is actually established (the SURFACE-the-dependence branch of BUG-008's FIX
  TARGET, not the measure-it branch).
- Fully closeable in-repo: a wording/scoping edit only. The cost formula and its
  conditional clause already exist; this leaf makes the `cost(M)` dependence visible at the
  headline sites and unifies the phrasing.
- Non-goal: the per-mode access-granularity accounting table (BUG-008-B).
- Non-goal: any verify probe (BUG-008-C) or measured seek-latency (BUG-008-D).

## Refined Description
**Scope:** surface and uniformly phrase the non-universal `cost(M)` dependence of the
random-access bound at the Abstract / §1.1 headline sites and the Theorem 10.3 scope
sentence.  **Non-goals:** granularity table (B); probe (C); measurement (D); unrelated
§-content.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-008
**Risks / Open Questions:** wording must not contradict Lemma 5.1b's existing
per-evaluation-cost discussion (transformer $\Theta(W)$ vs lookup-table $\Theta(\sigma^W)$
model storage); reuse that lemma's characterization rather than restating it differently.
**Expected Code Change:** YES (tex wording)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | confirm canonical "cost(M) is a non-universal dominant neural-eval" phrasing |
| ESTIMATION | DONE | EASY (localized wording edit, ~3 sites) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | choose one reusable sentence; align with Lemma 5.1b |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex wording edits in Abstract / §1.1 / Thm 10.3 scope sentence |
| TESTING | PENDING | manual re-read; invariant covered by BUG-008-C probe |
| POST_AUDIT | PENDING | adversarial re-review per project review protocol |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh at §12.7 decomposition of BUG-008;
  no code/proof change made to the paper before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- §12.7: this is one BLOCKING child leaf of BUG-008 (the closeable-now `cost(M)`-scoping layer).


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-008 A/B/C (commit pending). A: cost(M) surfaced as a non-universal, dominant per-position neural-eval term at Abstract/§1.1/§10.7/Thm 10.3 (polylog = in archive params with cost(M) fixed, not wall-clock; aligned w/ Lemma 5.1b). B: 6-row per-mode access-granularity table near Rem 10.3a (byte vs sub-block; per-query cost; wasted decode; sources Thm 10.3/Rem 10.3a/§6.3/§10.7; every cell cited, no new claim). C: scripts/verify/bug_008_random_access_cost.py (decode-window count p+k'-P<=K+k; K=Theta(|X|) degeneracy; simultaneous-polylog impossibility via interior K*; sub-block-granular cost; OVERALL -> PASS) + CI job. Hostile-referee POST_AUDIT clean (table byte-vs-sub-block split verified both directions). Supervisor: shortened the Type-III-C cost cell to a §10.7 pointer; table's 39.81pt overfull is house-style (188 such in the paper, many larger).

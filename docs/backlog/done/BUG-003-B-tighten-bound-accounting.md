# BUG-003-B - tighten-bound-accounting

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§10.7; Lemma 5.1a; Lemma 5.1c)
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
Source: §10.7; Lemma 5.1a; Lemma 5.1c.

Parent BUG-003's "tight general upper bound" fix target is genuinely hard in full generality
(the abstract says "neither a tight upper bound in general ... is provided"). But the paper
already contains THREE strictly different per-sync ceilings at different assumption strengths,
and they are scattered: (1) the loose floor-only bound `Δ_sync ≤ S·W·log₂(1/η)` (§10.7,
Lemma 5.1a); (2) the warm-rate-refined bound `Δ_sync ≤ W·(S·log₂(1/η) − r_warm_min)` under a
uniform per-byte warm-cost LOWER bound `r_warm_min` (§10.7 lines ~15945--15952); (3) the
EXACT Markov per-sync excess `δ_W' = H(X^W') − W'·h(X)`, total `(m−1)·δ_W'`, for stationary
order-W' sources with W ≥ W' and K ≥ W' (Lemma 5.1c). There is no single accounting
statement that (a) orders these three by assumption strength, (b) states precisely which
assumption buys which sharpness, and (c) shows that the realized natural-data figure the
Abstract conjectures lives in the gap between (2)/(3) and (1).

**FIX TARGET (this leaf):** add a consolidated **per-sync overhead accounting** (a short
corollary or remark + a 3-row table) in §10.7 that lays out the loose / warm-rate / Markov-tight
ceilings side by side, with each row labelled by its exact hypothesis and the resulting
per-sync and total-`m` form. This converts the diffuse "no tight bound in general" status into
an explicit honest ladder: tight under named conditions (Markov, warm-rate), loose otherwise.
It does NOT claim a tight bound in full generality (that residual stays open) -- it scopes
precisely how tight we can be and under what hypotheses.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class defect; this is a consolidation/sharpening of already-proved
  material (Lemma 5.1a worst-case, §10.7 warm-rate refinement, Lemma 5.1c exact-Markov), not
  a new open-research result. Closeable in-repo NOW.
- The genuinely-general tight bound remains open and is explicitly out of scope here; this
  leaf only assembles the assumption-graded ladder from existing proved pieces.
- Sequentially blocked by nothing, but should land before BUG-003-A's cross-link wording is
  finalized so the Abstract can point at the consolidated table.

## Refined Description
**Scope:** a §10.7 accounting corollary/remark + 3-row assumption-vs-ceiling table over the
loose, warm-rate, and Markov-tight bounds; verify each row matches its source statement.
**Non-goals:** a tight bound in full generality (stays open); any empirical number.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-003.
**Risks / Open Questions:** must verify the warm-rate refinement's direction (LOWER bound on
warm cost yields an UPPER bound on overhead -- §10.7 already flags this) and that Lemma 5.1c's
hypotheses (stationary order-W', W ≥ W', K ≥ W', N = mK) are reproduced exactly; do not
silently relax them.
**Expected Code Change:** new short remark/corollary + table in §10.7 of tex/rnr_coding.tex.
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | consolidate three proved ceilings into one assumption-graded ladder |
| ESTIMATION | DONE | ABOVE_EASY (assembles existing proofs; no new proof, careful hypothesis bookkeeping) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | choose corollary form + table columns (hypothesis, per-sync, total-m) |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex edit: new remark/corollary + table in §10.7 |
| TESTING | PENDING | cross-check each row vs Lemma 5.1a / §10.7 warm-rate / Lemma 5.1c; build green |
| POST_AUDIT | PENDING | adversarial re-review: no hypothesis silently dropped; ladder honest |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh as a decomposition child of BUG-003; no paper
  edit precedes this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Closeable in-repo: YES (consolidates already-proved ceilings; general-case tightness stays open).


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-003 A/B/C (commit pending). A: Abstract '~1%' scoped as an explicit heuristic CONJECTURE (rests on the warm-rate assumption of Cor 10.2a, not the bare floor) + crosslink to §10.7 ladder + to H10(c) measurement. B: Corollary 10.2a 'assumption-graded per-sync overhead ladder' + 3-row table consolidating the three PROVED ceilings -- rung1 loose S*W*log2(1/eta) [floor only; Lemma 5.1a], rung2 warm-rate-refined [floor + uniform warm LOWER bound], rung3 exact-Markov delta_W' [stationary order-W' Markov; Lemma 5.1c eqn 5.9]; states the fully-general tight bound REMAINS OPEN. C: scripts/verify/bug_003_syncpoint_overhead_fraction.py (realized overhead FRACTION = cold-restart excess / no-sync warm baseline, matching H10(c); rung-1 ceiling respected at 360 points; natural sources in ~1% band; ~1/K scaling; counterweights legitimately exceed, probe can fail honestly) + CI job. Hostile-referee POST_AUDIT clean (3 rungs verified against sources, no hypothesis flipped; probe non-vacuous via injected broken source).

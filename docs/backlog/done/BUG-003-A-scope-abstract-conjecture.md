# BUG-003-A - scope-abstract-conjecture

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Abstract; §10.7 cross-link); tex/rnr_experimental_design.tex (H10 cross-link, read-only)
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
Source: §10.7; Abstract.

Parent BUG-003 offers three fix targets (tight general upper bound, OR measured value, OR
scope the 1% as a conjecture). This leaf is the **scope-as-conjecture** target and is the
closeable-now core. The Abstract (tex/rnr_coding.tex lines ~229--234) already states
"We conjecture the actual overhead is on the order of one percent ... but neither a tight
upper bound in general nor a measured value ... is provided here," and §10.7 (lines
~15966--15970) already labels the "approximately one percent" figure a heuristic conjecture.
The remaining defect is purely epistemic-hygiene: the Abstract conjecture sentence does not
cross-link (i) the §10.7 sharper-bound *assumption* that the 1% depends on (the uniform
warm-rate lower bound `r_warm_min`), nor (ii) the experimental-design hypothesis that would
measure it (H10(c) in tex/rnr_experimental_design.tex, "ratio overhead from sync points ...
below 2 percent for K = 64 KiB"). A reader cannot, from the Abstract alone, trace the
conjecture to its assumption-set or its falsification plan.

**FIX TARGET (this leaf):** make the Abstract's 1% conjecture sentence self-consistent and
traceable -- explicit conjecture label (verify it is already unambiguous), a parenthetical
pointer to the §10.7 warm-rate assumption the figure rests on, and a pointer to the §11 /
experimental-design H10 measurement that would confirm/falsify it. No new claim is asserted;
the epistemic status is only made precise and navigable.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class defect; "fix" = scope/cross-link the paper claim to match
  what is actually established, no new theorem.
- Closeable in-repo NOW (text-only edit to tex/rnr_coding.tex Abstract + at most a one-line
  forward pointer; tex/rnr_experimental_design.tex is read-only here, used only to cite H10).
- This leaf does NOT tighten the bound (that is BUG-003-B) and does NOT measure anything
  (that is BUG-003-D). It only ensures the conjecture is labelled and traceable.

## Refined Description
**Scope:** the Abstract 1% conjecture sentence and its cross-links to §10.7 assumption and
H10 measurement.  **Non-goals:** the bound itself, any numeric estimate, any measurement.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-003.
**Risks / Open Questions:** must not weaken or strengthen the existing claim; pure scoping.
**Expected Code Change:** small text edit to tex/rnr_coding.tex Abstract (no math change).
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope-as-conjecture target; closeable in-repo |
| ESTIMATION | DONE | EASY (text-only Abstract edit + cross-links) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | choose exact cross-link wording (§10.7 assumption + H10) |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex edit to Abstract conjecture sentence |
| TESTING | PENDING | re-read Abstract+§10.7 for self-consistency; xelatex build green |
| POST_AUDIT | PENDING | adversarial re-review: claim neither weakened nor overclaimed |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): created fresh as a decomposition child of BUG-003; no paper
  edit precedes this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Closeable in-repo: YES (text-only scoping edit).


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-003 A/B/C (commit pending). A: Abstract '~1%' scoped as an explicit heuristic CONJECTURE (rests on the warm-rate assumption of Cor 10.2a, not the bare floor) + crosslink to §10.7 ladder + to H10(c) measurement. B: Corollary 10.2a 'assumption-graded per-sync overhead ladder' + 3-row table consolidating the three PROVED ceilings -- rung1 loose S*W*log2(1/eta) [floor only; Lemma 5.1a], rung2 warm-rate-refined [floor + uniform warm LOWER bound], rung3 exact-Markov delta_W' [stationary order-W' Markov; Lemma 5.1c eqn 5.9]; states the fully-general tight bound REMAINS OPEN. C: scripts/verify/bug_003_syncpoint_overhead_fraction.py (realized overhead FRACTION = cold-restart excess / no-sync warm baseline, matching H10(c); rung-1 ceiling respected at 360 points; natural sources in ~1% band; ~1/K scaling; counterweights legitimately exceed, probe can fail honestly) + CI job. Hostile-referee POST_AUDIT clean (3 rungs verified against sources, no hypothesis flipped; probe non-vacuous via injected broken source).

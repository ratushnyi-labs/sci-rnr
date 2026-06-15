# BUG-004-B - scope-ideal-rate-as-intractable

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (ideal-rate invocation sites identified by BUG-004-A)
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
> (Filed gap 004 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at section 12.7 decomposition of BUG-004 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: parent BUG-004 (Rem 7.15g; §7.15 open middle).

This is the EDIT half of the parent FIX TARGET: at each ideal-rate invocation site that
BUG-004-A flags as "caveat missing" or "caveat only locally implied", add a precise,
minimal cross-reference making explicit that the exact ideal rate `L*_RNR = H(X^N)` is an
**intractable target** (exact value / coefficient query is #P-hard = FP^#P-complete, Rem
7.15g) and that the **operative object** is the eps-approximation (poly-eps under filter
stability, Rem 7.15c) or the **achievable scheme** (realized causal codelength, one forward
pass per symbol, Rem 7.15b). The §7.15g Corollary (~line 9070) already states the operational
reading locally ("which is precisely why RNR reports achievable rates and bounds rather than
the exact ideal"); this leaf propagates that one-line cross-link to the operationally central
invocation sites (notably the §7.1 rate-optimality theorem and §6.6c converse) so a reader
encountering `H(X^N)` as a target is not misled into thinking it is efficiently computable.

The edits must be MINIMAL and TRUTHFUL: a parenthetical or a single sentence pointing to
Rem 7.15g (hardness) + 7.15b/7.15c (operative surrogate). Do NOT weaken any correct converse
(the floor `E[L] >= H(X^N)` is a valid information-theoretic bound regardless of its
computational cost; the caveat is about COMPUTING the optimum, not about the bound's validity).

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class defect; closeable IN-REPO NOW (paper-text edits + recompile).
- This leaf does NOT modify tex/rnr_coding.tex during INTAKE; it is the scoped work item whose
  BUILD stage will edit the .tex. (The intake author is forbidden from touching the .tex this turn.)
- The distinction to preserve: the entropy floor is a sound CONVERSE; its EXACT computation is
  #P-hard; the achievable scheme prices the optimum approximately. The edit clarifies the third
  point without disturbing the first.
- Inherits MEDIUM priority from parent; ABOVE_EASY because each site needs a context-faithful
  one-liner (not a blind find-replace) and the result must survive an adversarial re-read.

## Refined Description
**Scope:** add the intractability + operative-object cross-link at the flagged sites.
**Non-goals:** the audit itself (BUG-004-A), the verification probe (BUG-004-C), any new proof,
any change to the §7.15g hardness statement (already correct).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** BUG-004-A (must consume its inventory table) — SEQUENTIALLY BLOCKED on A.
**Risks / Open Questions:** over-caveating could clutter the prose or wrongly suggest the
converse bound is invalid — keep edits surgical. If BUG-004-A surfaces a site that genuinely
NEEDS the exact ideal (not a surrogate), that is a substantive finding to escalate, not a
one-liner.
**Expected Code Change:** YES (LaTeX prose edits at BUILD; none at INTAKE)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = surgical cross-link edits, converse untouched |
| ESTIMATION | DONE | ABOVE_EASY (context-faithful one-liners at ~4-6 sites) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | one-line template: "(exact value #P-hard, Rem 7.15g; operative object is the achievable rate, Rem 7.15b/c)" |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | apply the cross-link at each flagged site; recompile PDF |
| TESTING | PENDING | PDF builds clean; BUG-004-C probe referenced where the operational reading is asserted |
| POST_AUDIT | PENDING | cold re-read + adversarial pass: caveat is correct, converse not weakened |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 child of BUG-004 (decomposition on 2026-06-15). Created fresh this turn (Recovered: N/A).
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Sequentially blocked on BUG-004-A (needs the inventory). This is the in-repo-closeable-now
  edit core; no genuinely-hard residual (the hard math, Rem 7.15g, is already proven).


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-004 (FULLY closed -- no external leaf). A: audited all ideal-rate/H(X^N)-target sites; B: scoped the 3 unscoped computable-target sites (Thm 6.6c.1 floor, §7.1b converse, §7.27' ideal rate) with one canonical clause -- H(X^N) is an info-theoretic BENCHMARK, exact value #P-hard/FP^#P-complete (Rem 7.15g), poly-eps-approximable only under filter stability (Rem 7.15c); operative object = achievable scheme; converse inequalities preserved verbatim. Sites with a non-finite-block sense (§7.36 drain rate, §7.20 LD rate function, §7.32 empirical Hhat) correctly left alone. C: scripts/verify/bug_004_ideal_vs_achievable_gap.py (G1 intractable exact target recovers #SAT; G2 poly eps-approx under FS geometric A_d->H; G3 achievable one-pass realizes the rate, exact optimum never computed) + CI job. Hostile-referee POST_AUDIT clean.

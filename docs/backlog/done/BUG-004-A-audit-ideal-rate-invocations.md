# BUG-004-A - audit-ideal-rate-invocations

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (all invocations of the ideal rate L*_RNR = H(X^N))
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

The parent FIX TARGET says: "make explicit, **wherever** the 'ideal rate' is invoked, that
it is an intractable target and the operative object is the eps-approximation / achievable
scheme." This leaf is the AUDIT half: enumerate every place in `tex/rnr_coding.tex` where the
exact ideal rate `L*_RNR = H(X^N)` (or the entropy floor `H(X^N)`, the chain-rule
`sum_t H(X_t | X_{<t})`) is invoked as the *target* the framework optimises toward, and
record for each whether the intractability caveat (exact value is #P-hard / FP^#P-complete;
Rem 7.15g) is present, locally implied, or missing.

Known invocation sites to inspect (from intake grep; verify line numbers live before editing):
- §6.6c converse Theorem 6.6c.1: `E[L_IIIC] >= H(X^N) >= N h(X)` (~line 5666) — entropy floor as target.
- §7.1 RNR rate-optimality: `E[L_RNR] >= H(X^N)`, `(1/N) E[L_RNR] -> h(X)` (~lines 6777, 6795).
- §7.15g Corollary "the exact ideal RNR rate is #P-hard" (~line 9070) — the caveat ALREADY lives here.
- §7.27'/neural bisection Remark: "ideal rate is `H(X^N) = N h + O(1)`" (~line 11764).
- §7.5 sub-block "draining at the matched ideal rate h" (~line 15362).
- §11.7 enwik9 numerical prediction (~line 8463): uses predictor cross-entropy `H_{M'}` (an
  achievable/realized object) — confirm it does NOT silently equate to the exact ideal.

This leaf produces only the AUDIT TABLE (a finding artifact recorded in the backlog item's
Execution Tracking / a follow-up note), NOT the .tex edits. The edits are BUG-004-B.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class defect in the theory paper `tex/rnr_coding.tex`; no runtime
  product. This leaf is closeable IN-REPO NOW (a read-only inventory pass).
- The underlying hardness math (Rem 7.15g) is already proven and committed; this leaf does NOT
  re-derive or re-verify it — it only locates where the result must be cross-referenced.
- "Fix" for this leaf = produce a complete, line-anchored inventory of ideal-rate invocations
  with a present/implied/missing caveat verdict for each. No claim is closed or downgraded here.

## Refined Description
**Scope:** the inventory of ideal-rate invocation sites and their caveat status. **Non-goals:**
the actual scoping edits (BUG-004-B), the verification probe (BUG-004-C), any new proof.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-004.
**Risks / Open Questions:** line numbers drift; the audit must grep live HEAD (parallel §7.34
development per project memory). Risk that an invocation is found that genuinely needs the exact
ideal (not the achievable surrogate) — if so, flag it for BUG-004-B as a substantive, not
cosmetic, scoping fix.
**Expected Code Change:** NO (audit only; the edit is BUG-004-B)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = inventory only, no edits |
| ESTIMATION | DONE | EASY (read-only grep + classify) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | grep patterns: `L_{\mathrm{RNR}}`, `H(X^N)`, `ideal rate`, `chain rule`, `optimum` |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | produce inventory table (no .tex change) |
| TESTING | PENDING | cross-check each site's surrounding paragraph for the caveat link |
| POST_AUDIT | PENDING | hand the table to BUG-004-B; confirm completeness |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 child of BUG-004 (decomposition on 2026-06-15). Created fresh this turn (Recovered: N/A).
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- This leaf is the in-repo-closeable-now audit core; it carries no genuinely-hard residual.


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-004 (FULLY closed -- no external leaf). A: audited all ideal-rate/H(X^N)-target sites; B: scoped the 3 unscoped computable-target sites (Thm 6.6c.1 floor, §7.1b converse, §7.27' ideal rate) with one canonical clause -- H(X^N) is an info-theoretic BENCHMARK, exact value #P-hard/FP^#P-complete (Rem 7.15g), poly-eps-approximable only under filter stability (Rem 7.15c); operative object = achievable scheme; converse inequalities preserved verbatim. Sites with a non-finite-block sense (§7.36 drain rate, §7.20 LD rate function, §7.32 empirical Hhat) correctly left alone. C: scripts/verify/bug_004_ideal_vs_achievable_gap.py (G1 intractable exact target recovers #SAT; G2 poly eps-approx under FS geometric A_d->H; G3 achievable one-pass realizes the rate, exact optimum never computed) + CI job. Hostile-referee POST_AUDIT clean.

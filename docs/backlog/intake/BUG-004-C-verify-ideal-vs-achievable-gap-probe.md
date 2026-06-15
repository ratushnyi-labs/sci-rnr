# BUG-004-C - verify-ideal-vs-achievable-gap-probe

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** scripts/verify/ (probe separating exact-ideal H(X^N) from achievable realized rate)
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

Per the project rule "each theorem/remark has an executable check in scripts/verify/", the
scoping claim that "the exact ideal rate is the intractable target while the operative object
is the achievable scheme" deserves a small executable PROBE that DEMONSTRATES the gap on a
tiny finite-state source. The existing `scripts/verify/remark_7_15g_hmm_entropy_sharp_p_hard.py`
already verifies the #SAT -> coefficient-query hardness reduction (V1-V7). This leaf adds the
OPERATIONAL companion probe asserting, on a small explicit HMM / finite-state neural-surrogate:

1. the **achievable / realized** causal codelength `sum_t -log2 Q(x_t | x_{<t})` is computed in
   ONE forward pass (poly), matching what the deployed encoder ships (Rem 7.15b);
2. the **exact ideal** `L*_RNR = H(X^N) = sum_t H(X_t | X_{<t})` equals the brute-force block
   entropy on the small instance (identity check), confirming the object that 7.15g proves
   #P-hard is precisely the optimisation target;
3. the realized rate >= H(X^N) (the converse) and the two coincide only in expectation under
   the true source — i.e. the achievable surrogate is the operative object, the exact ideal is
   the intractable target.

This is a small-N, exact-fraction probe (no external data, no experiment); fully closeable
IN-REPO. It does NOT require any new mathematics — it instruments already-proven facts so the
scoping edits of BUG-004-B carry a green-CI verification anchor.

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class; the probe is the project's "executable check" form for a
  characterization/scoping claim (identity + poly-vs-target separation on a small instance).
- Closeable IN-REPO NOW: a self-contained Python script under `scripts/verify/`, PASS/FAIL to
  stdout, includable as a CI job (per the project verification rule). Verify locally with the project venv
  `/Users/para/.venvs/rnr/bin/python` (gh cannot observe this repo's CI per memory).
- Does NOT re-prove 7.15g; reuses/extends its small HMM. May share a fixture with the existing
  7.15g script but stays a distinct file (one file per claim).
- This leaf does NOT create the script during INTAKE; it is the scoped work item whose BUILD
  stage writes the probe. (The .tex is untouched; only scripts/verify/ is in scope.)

## Refined Description
**Scope:** a small executable probe in scripts/verify/ exhibiting (poly achievable rate) vs
(exact ideal = H(X^N), the #P-hard target) on a tiny finite-state source, with the converse
inequality. **Non-goals:** the hardness reduction (already in the 7.15g script), the .tex edits
(BUG-004-B), the audit (BUG-004-A), any large-scale / external experiment.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-004; logically pairs with BUG-004-B (the probe is the verification
anchor for B's operational-reading edits) but can be authored in parallel with A/B.
**Risks / Open Questions:** keep N tiny so brute-force H(X^N) is exact; use exact fractions
(Fraction) to avoid float drift, matching the 7.15g script's convention.
**Expected Code Change:** YES (new scripts/verify/ probe at BUILD; none at INTAKE)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = small operational probe, reuses 7.15g fixture |
| ESTIMATION | DONE | EASY (tiny exact-fraction script, no new math) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | checks: P1 one-pass realized length; P2 H(X^N) identity vs brute force; P3 converse |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | write scripts/verify/ probe; PASS/FAIL stdout |
| TESTING | PENDING | run with /Users/para/.venvs/rnr/bin/python; all checks PASS; wire into CI job |
| POST_AUDIT | PENDING | confirm the probe verifies the POLY claim itself, not a 2^N proxy (per 7.15f lesson 5) |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 child of BUG-004 (decomposition on 2026-06-15). Created fresh this turn (Recovered: N/A).
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- In-repo-closeable-now verification core; no external work, no genuinely-hard residual.

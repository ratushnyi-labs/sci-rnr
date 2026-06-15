# BUG-007-B - verify-thm-10-1-integer-reorder-probe

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Thm 10.1) → scripts/verify/ (new probe)
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
> (Filed gap 007 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at §12.7 decomposition of BUG-007 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: Thm 10.1.

Thm 10.1 has no executable verification in `scripts/verify/`, unlike the project's
per-theorem convention (project rule: each theorem/lemma has a runnable PASS/FAIL probe).
The theorem's MATHEMATICAL core is two independently checkable claims: (1) under the
per-layer signed bit-width bound `b_ℓ ≥ ⌈log2(d_ℓ·W_ℓ·I_{ℓ-1}+1)⌉+1`, no accumulator
overflows for any conforming integer inference graph; (2) integer add/multiply is
associative/commutative, so two implementations using different (but each pinned)
reduction orders produce a bit-identical result. These can be demonstrated in-repo with
a small construction WITHOUT any real hardware — closing the verification-script gap for
Thm 10.1 even though the full cross-platform implementation (BUG-007-D) is external.

## Refined Description
**Scope (this leaf only):** add `scripts/verify/thm_10_1_integer_reorder_bitexact.py`
that, for randomized small integer inference graphs (a few layers, integer weights and
inputs within declared magnitude bounds, integer-to-integer lookup-table non-linearities):
1. computes the per-layer accumulator bit-width bound of Thm 10.1 and asserts the
   worst-case accumulator magnitude is within the signed range (no overflow) — and as a
   negative control, demonstrates that a deliberately undersized bit-width DOES overflow,
   matching the int16-fan-in numerical sanity check in the §10.5 discussion (`b=32`
   overflows at fan-in ≥ 2 under int16 bounds; `b≥44` needed);
2. evaluates each matrix-vector reduction under several distinct summation orders
   (sequential left-to-right, reversed, randomized pairwise tree) and asserts the exact
   integer outputs are bit-identical across orders;
3. prints `PASS`/`FAIL` to stdout in the project's convention; PASS confirms the
   no-overflow + reorder-invariance claims of Thm 10.1 as stated.

**Fix target (this leaf):** one executable verify script demonstrating Thm 10.1's
arithmetic core, suitable for inclusion in the CI verify job.

**Non-goals:** real GPU/NPU/CPU-SIMD hardware (BUG-007-D); paper prose scoping
(BUG-007-A); engineering-spec / experimental-design cross-links (BUG-007-C). This probe
validates the THEOREM AS STATED (integer associativity + bit-width bound), not end-to-end
cross-platform conformance of a real model.

**Closeable in-repo NOW:** YES — pure-Python integer arithmetic, no external deps beyond
the project venv (`/Users/para/.venvs/rnr/bin/python`); verify locally per project rule.

**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-007. (Logically reads Thm 10.1; not blocked by BUG-007-A.)
**Risks / Open Questions:** the probe demonstrates the theorem's claim within Python's
arbitrary-precision integers; the bit-width bound is checked symbolically (worst-case
magnitude vs declared `b`), which is exactly what Thm 10.1 asserts — the script must make
clear it tests the MATH, not a hardware target.
**Expected Code Change:** YES (new script under scripts/verify/)
**Documentation Recovery Required:** NO

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class defect; the "executable verification" is the project's
  standard practical-check artifact (project rules: "кожен крок як робоча програма" —
  every step as a runnable program).
- Follows existing naming/format convention of `scripts/verify/` (one file per theorem,
  PASS/FAIL to stdout, CI-includable).

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = arithmetic-core probe only |
| ESTIMATION | DONE | ABOVE_EASY (small construction + negative control) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | choose graph shape, magnitude bounds, reduction orders |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | write thm_10_1_integer_reorder_bitexact.py |
| TESTING | PENDING | run locally; PASS expected; add to CI verify job |
| POST_AUDIT | PENDING | adversarial: does PASS actually exercise overflow + reorder? |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7: child leaf of BUG-007; created fresh this turn (Recovered: N/A).
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Honest-step principle: closes the IN-REPO verification of the theorem's math now;
  the external cross-platform demonstration stays a separate flagged leaf (BUG-007-D).
- Verify LOCALLY (`/Users/para/.venvs/rnr/bin/python`); CI is account-mismatched and
  cannot be observed from this repo.

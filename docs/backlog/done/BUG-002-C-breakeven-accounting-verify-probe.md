# BUG-002-C - breakeven-accounting-verify-probe

**Type:** TASK
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Theorem 7.15 caveat ll.8593-8598 break-even; T7.21 ll.9765-9782); scripts/verify/
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
> (Filed gap 002 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at §12.7 decomposition of BUG-002 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: parent BUG-002 FIX TARGET (executable-verification discipline per the
project local-rules file: "кожна теорема... повинна мати виконувану перевірку у `scripts/verify/`").

The model-inclusive accounting arithmetic — self-contained vs amortized size, the
break-even volume, and the per-regime net savings — is currently stated in prose
(T7.15 caveat, T7.21, §10.2) with NO standalone executable probe that recomputes
it and emits PASS/FAIL. The break-even identity
$V^*\cdot(\log_2|\Sigma|-H_{M'})=140\,\mathrm{GB}\cdot\log_2|\Sigma|$ giving
$V^*\approx140\cdot8/(8-0.664)\approx152\,\mathrm{GB}$ (ll.8595-8597), the
$\approx140\times$ single-archive expansion factor, and T7.21's
$\mathrm{TotalCost}(L,N)=L+N\,h(X)+\Theta(N^{1/(1+\alpha)})$ with
$L^*=(NC\alpha)^{1/(1+\alpha)}$ are unverified by code. Existing
`scripts/verify/thm_7_15_neural_rnr_end_to_end.py` and
`thm_7_21_predictor_archive_tradeoff.py` exist but should be checked for whether
they already cover the model-inclusive/break-even arithmetic; if not, a small probe
is needed.

## Refined Description
**One-leaf scope:** add (or extend) a small `scripts/verify/` probe that
recomputes the model-inclusive accounting and break-even arithmetic and asserts it
matches the paper's stated figures:
1. self-contained single-archive size vs raw source $\Rightarrow$ $\approx140\times$
   expansion (negative savings) at $N=1\,\mathrm{GB}$, $140\,\mathrm{GB}$ model;
2. break-even $V^*\approx152\,\mathrm{GB}$ from the stated identity;
3. T7.21 $L^*(N)$ optimum and the sub-linear $\Theta(N^{1/(1+\alpha)})$ total
   overhead (calculus optimum reproduced numerically for the cited $\alpha$ values);
4. monotone "net savings becomes positive once amortized volume $V>V^*$" check.
Emit PASS/FAIL to stdout per project convention; suitable for the CI verify job.

**Fix target:** one self-contained Python file under `scripts/verify/`
(e.g. `bug_002_model_inclusive_accounting.py` or an extension of the two existing
scripts), runnable with the project venv `/Users/para/.venvs/rnr/bin/python`.

**Closeable in-repo now:** YES (executable verification of the paper's own
arithmetic — no external data needed; this is a pure numeric-identity / calculus
check on stated constants).

**Non-goals:** the prose table (BUG-002-B); the cross-link pass (BUG-002-A); the
REAL empirical enwik9 measurement (BUG-002-D) — this probe verifies the paper's
arithmetic, NOT the empirical $0.664$ bpb datum.
**Impacted UCs:** N/A (theory paper)  **Impacted BR/WF:** N/A
**Dependencies:** parent BUG-002; numbers should agree with BUG-002-B's table.
**Risks / Open Questions:** must NOT overclaim — the probe verifies internal
arithmetic consistency only; the empirical $H_{M'}=0.664$ bpb is an input constant
(Delétang et al. 2024), not something this probe measures.
NOTE: this leaf writes to `scripts/verify/` (allowed: it is the project's
verification directory); the supervisor reconciles/commits.
**Expected Code Change:** YES (verify script)
**Documentation Recovery Required:** NO

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class project; the executable-probe discipline applies
  to numeric claims per the project local-rules "Практична верифікація".
- Inputs are the paper's stated constants ($H_{M'}=0.664$, $|\Sigma|=256$,
  model $\approx140\,\mathrm{GB}$, $\alpha\in\{0.34,0.5,1\}$); the probe checks the
  derived quantities, treating these as given.

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = executable arithmetic/break-even probe, internal-consistency only |
| ESTIMATION | DONE | EASY (closed-form arithmetic + calculus optimum) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | check existing two scripts; decide extend-vs-new |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | write/extend scripts/verify probe |
| TESTING | PENDING | run under /Users/para/.venvs/rnr/bin/python; PASS/FAIL; CI verify job |
| POST_AUDIT | PENDING | adversarial re-read: does the probe avoid overclaiming the empirical datum? |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7: created fresh as a BLOCKING child leaf of BUG-002 on 2026-06-15; not recovered.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A,
  but the project mandates an executable verify probe for numeric claims.
- Verify LOCALLY (`/Users/para/.venvs/rnr/bin/python`); CI runs the verify job.


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-002 A/B/C (commit pending). Consolidated model-inclusive-vs-amortized table in §10.2 + break-even V*~152 GB crosslinks at Abstract & §10.2 + arithmetic probe scripts/verify/bug_002_model_inclusive_accounting.py (OVERALL -> PASS) + CI job. Collation only, no new claims; hostile-referee POST_AUDIT clean.

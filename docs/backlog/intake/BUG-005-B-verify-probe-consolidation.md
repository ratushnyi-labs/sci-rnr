# BUG-005-B - verify-probe-consolidation

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** scripts/verify/ (op2a_*, op4_*, lemma_6_8*, lemma_6_9*, remark_4_3d*) against tex/rnr_coding.tex §13.4
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
> (Filed gap 005 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at section 12.7 decomposition of BUG-005 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: §13.4.1/13.4.2/13.4.4 attack vectors; their cited verify scripts.

The §13.4 attack-vector roadmap closes numerous attack vectors NEGATIVELY (e.g.
density-shifted Berman-Karpinski, Dinur-Safra FGLSS, K-escape budget, Label-Cover
K-uniform, DkS affine dilution, non-count entropy objectives) and CLOSES sub-cases
POSITIVELY (continuous-bit-cost 6.8e/6.8b, fixed-field Type-III-C 6.9, additive
rule-cost 6.8k). Each closure cites a verify script in `scripts/verify/`
(op2a_dinur_safra_attack_check.py, op2a_k_escape_budget_check.py,
op2a_label_cover_kuniform_check.py, op2a_additive_rule_cost_check.py,
op2a_sgp_constant_survey_check.py, op4_dksh_affine_dilution.py,
op4_noncount_entropy_objective.py, lemma_6_8b/6_8k/6_9*, remark_4_3d*). What is
missing: (1) a check that EVERY closed attack vector named in §13.4 has a
present, runnable, PASSing verify script (no orphan claims); and (2) a single
consolidated status-probe that asserts the "closed-vector ledger" the paper
states matches what the scripts actually establish (so the CI `verify` job
covers the OP2/OP4 closures, not only the §7.34 RD-dispersion probes).

## Clarification / Assumptions
- This is the IN-REPO VERIFICATION leaf of BUG-005: it does NOT close any open
  hardness question. It guarantees the EXISTING closures are executable and that
  the paper's closed-vector claims are not orphaned.
- GOVERNANCE/DOCUMENTATION-class; the "executable verification" here is per the
  project rule "each theorem/lemma/proposition has a runnable PASS/FAIL probe and
  is wired into CI (green CI = all formal+executable checks pass)".
- closeable_in_repo = YES: the relevant scripts already exist; this leaf audits
  coverage and adds at most one consolidated index/status probe + (if any gap is
  found) wires the missing OP2/OP4 probes into the CI `verify` job.

## Refined Description
**Scope:** (a) Inventory the §13.4 closed attack vectors (positive and negative)
and cross-check each against a present `scripts/verify/` script; flag any orphan
(claim with no script) or stale (script absent / not PASSing) entry. (b) Run the
OP2/OP4 verify scripts locally (project venv /Users/para/.venvs/rnr/bin/python)
and confirm PASS. (c) Add a single consolidated index/status probe (e.g.
`scripts/verify/op2_op4_closed_vector_ledger_check.py`) that enumerates the
closed-vector ledger and asserts each named script exists and prints PASS,
emitting one OVERALL -> PASS/FAIL line in the CI style. (d) If the CI `verify`
job (.github/workflows/build.yml) does not already run the OP2/OP4 probes, note
the wiring gap for the build (the build.yml edit itself is the BUILD step, not
intake).
**Non-goals:** any new hardness reduction (sibling BUG-005-D); the §13.2 status
ledger prose/table (sibling BUG-005-A); memory-journal dead-end cross-linking
(sibling BUG-005-C). No modification to tex/rnr_coding.tex.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-005; sibling BUG-005-A (the ledger this probe
encodes is the one A makes precise — B may proceed in parallel but should align
its row set with A's table before final commit).
**Risks / Open Questions:** the consolidated probe must assert script PRESENCE +
PASS, not re-derive the math (avoid duplicating each script's own assertions).
gh cannot observe this repo's CI (account mismatch) — verify LOCALLY only.
**Expected Code Change:** YES (new verify script; possibly CI job wiring noted)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | confirm ledger row set with sibling A |
| ESTIMATION | PENDING | ABOVE_EASY: inventory + run scripts + write index probe |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | consolidated index probe asserts presence+PASS, no math re-derivation |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | new scripts/verify/op2_op4_closed_vector_ledger_check.py |
| TESTING | PENDING | run all OP2/OP4 probes locally; confirm OVERALL -> PASS |
| POST_AUDIT | PENDING | confirm CI `verify` job covers OP2/OP4 (or flag wiring gap) |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 DECOMPOSITION child of BUG-005. This is the executable-verification leaf:
  it lowers parent complexity by ensuring all already-closed sub-cases are
  CI-witnessed, isolating the genuinely-open residual to sibling BUG-005-D.
- §12.1.1 (Backlog-File-First): created fresh this turn; no script/paper change
  precedes this file. `Recovered: N/A`.
- Class: GOVERNANCE/DOCUMENTATION; runtime test stages are class-level N/A; the
  scripts/verify probes are the class-appropriate executable checks.
- Verify LOCALLY (gh cannot observe this repo's CI).

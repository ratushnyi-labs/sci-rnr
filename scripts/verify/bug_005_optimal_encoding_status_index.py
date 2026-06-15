#!/usr/bin/env python3
"""
BUG-005  --  consolidated status-index probe for the "optimal RNR encoding"
             complexity ledger (OP2(a) / OP2(b) / OP4 support/grammar selection).

Scope (governance/documentation; NO new mathematics). This probe instruments the
§13.2 status ledger so the BUG-005-A/B/C scoping edits carry a green-CI anchor. It
does NOT re-derive any hardness or fabricate a PASS for an OPEN sub-case. It
asserts:

  (1) for each CLOSED / in-P ledger row, the WITNESSING verify script is present
      in scripts/verify/ and RUNS to a green exit (return code 0). These scripts
      carry their own internal assertions; this index re-runs them and confirms
      each one exits 0 (its own PASS/verdict line is its proof, not re-derived here);
  (2) for each OPEN ledger row, NO closure script is claimed -- the row is printed
      as OPEN and the index intentionally provides no PASS for it (a fabricated PASS
      on an open sub-case would be a Truth>peremoga violation);
  (3) the ledger printed here matches the §13.2 table granularity exactly
      (per-sub-question status / location / blocker-or-proof / script).

Final line: 'OVERALL -> PASS' iff every CLOSED/in-P witnessing script is present and
exits 0, and the OPEN rows are correctly carried without a fabricated PASS.

The witnessing scripts (sources of truth, each verified standalone):
  OP2(a) continuous bit-cost        CLOSED  -> lemma_6_8b_bit_cost_asymptotic.py
                                               op2a_k_escape_budget_check.py (alphabet floor)
  OP2(a) additive rule-cost         CLOSED+ -> lemma_6_8k_approx_additive_rule_cost.py
                                               op2a_additive_rule_cost_check.py
  OP2(a) discrete ceiling residual  OPEN    -> op2a_sgp_constant_survey_check.py (negative route)
  OP2(b) free-field continuous      in P    -> remark_6_9a_continuous_op2b_easy.py
  OP4 Type-II general               OPEN    -> op4_dksh_affine_dilution.py (negative route)
                                               op4_noncount_entropy_objective.py (negative route)
  OP4 Type-III-C fixed-field        CLOSED  -> lemma_6_9_typeiiic_mesp_support_selection.py
  OP4 E_12 constant-gap             OPEN    -> remark_4_3d_e12_constant_gap_dilution.py (negative route)
  OP4 fixed-byte                    OPEN    -> (none; no gap-preserving construction)

Note: the §13.4.1 negative-control attack-vector scripts (op2a_dinur_safra_attack_check.py,
op2a_label_cover_kuniform_check.py) verify CLOSED-NEGATIVELY routes; they support OPEN-sub-case
re-attempt prevention, NOT a closure, so they are checked for presence+exit-0 but the
sub-cases they touch stay OPEN.
"""

import os
import sys
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))


def run_script(name):
    """Run a sibling verify script; return (exit_code, present_bool)."""
    path = os.path.join(HERE, name)
    if not os.path.isfile(path):
        return (None, False)
    proc = subprocess.run(
        [sys.executable, path],
        capture_output=True, text=True, cwd=os.path.dirname(HERE) or ".",
    )
    return (proc.returncode, True)


# Ledger rows: (sub_question, status, location, blocker_or_proof, [witness scripts], must_pass)
# must_pass = True  -> CLOSED / in-P: script(s) must be present AND exit 0 (their own proof).
# must_pass = False -> OPEN: script(s) (if any) are negative-route controls; presence+exit-0
#                      checked, but the SUB-CASE remains OPEN and gets NO fabricated PASS.
LEDGER = [
    ("OP2(a) continuous bit-cost SGP", "CLOSED",
     "Lemma 6.8e/6.8b; §13.4.1",
     "inherits Charikar gamma_SGP=8569/8568 with O(1/(|V|log|V|)) correction",
     ["lemma_6_8b_bit_cost_asymptotic.py", "op2a_k_escape_budget_check.py"], True),
    ("OP2(a) additive rule-cost L_alpha", "CLOSED+",
     "Lemma 6.8k / 6.8k-approx; §13.4.1",
     "APX-hard for every alpha>=0, ceiling-free gap monotone in alpha",
     ["lemma_6_8k_approx_additive_rule_cost.py", "op2a_additive_rule_cost_check.py"], True),
    ("OP2(a) discrete integer/ceiling SGP", "OPEN",
     "Conj 6.8d; §13.4.1",
     "ceiling discontinuity at 2^k boundaries; 1/8568 < 1/log2 N for all N<2^8568",
     ["op2a_sgp_constant_survey_check.py"], False),
    ("OP2(b) free-field continuous (Gaussian)", "in P",
     "Remark 6.9a; §13.4.4",
     "sample-covariance ML plug-in; iid Gaussian cannot hide poly-estimable cov",
     ["remark_6_9a_continuous_op2b_easy.py"], True),
    ("OP2(b) single/amortised discrete (free M)", "OPEN",
     "Lemma 13.4.2a; Thm 13.4.2b; §13.4.2",
     "predictor-absorption is UPPER BOUND only; PRG reduces to OP2(a)-APX (Conj 6.8d)",
     [], False),
    ("OP4 Type-II general support selection", "OPEN",
     "Lemma 5.7c/5.7d; §13.4.4",
     "affine n_0(S) dilution 1+Theta(D_YES/m); natural A_s saturates at 1+O(1/log s)",
     ["op4_dksh_affine_dilution.py", "op4_noncount_entropy_objective.py"], False),
    ("OP4 Type-III-C fixed-field (log-det)", "CLOSED (Partial)",
     "Lemma 6.9; §13.4.4",
     "Schur-complement log-det = MESP; APX-hard 5/4 (Ohsaka) under P!=NP, high-SNR survives",
     ["lemma_6_9_typeiiic_mesp_support_selection.py"], True),
    ("OP4 E_12 constant-gap labeling", "OPEN",
     "Conj 4.3e; Remark 4.3d",
     "missing dilation-1/robust-NO promise; NOT row-stochasticity, NOT affine dilution",
     ["remark_4_3d_e12_constant_gap_dilution.py"], False),
    ("OP4 fixed-byte restricted Type-III-C", "OPEN",
     "§13.2",
     "byte-encoding adds Theta(log|V|) multiplicative overhead diluting the gap",
     [], False),
]

# §13.4.1 negative-control attack-vector scripts (re-attempt prevention for OPEN OP2(a)).
NEG_CONTROL = [
    "op2a_dinur_safra_attack_check.py",
    "op2a_label_cover_kuniform_check.py",
]

print("=" * 78)
print("BUG-005 optimal-encoding status index (OP2/OP4 support/grammar selection)")
print("=" * 78)

all_ok = True
closed_count = 0
open_count = 0

for (sq, status, loc, blocker, scripts, must_pass) in LEDGER:
    tag = "CLOSED/in-P" if must_pass else "OPEN"
    if must_pass:
        closed_count += 1
    else:
        open_count += 1
    print(f"\n[{status}] {sq}")
    print(f"    location: {loc}")
    label = "proof" if must_pass else "blocker"
    print(f"    {label}: {blocker}")

    if must_pass:
        # CLOSED / in-P: every witnessing script must be present and exit 0.
        for s in scripts:
            rc, present = run_script(s)
            if not present:
                print(f"    MISSING witness script: {s}  -> FAIL")
                all_ok = False
            elif rc != 0:
                print(f"    witness {s}: exit {rc}  -> FAIL")
                all_ok = False
            else:
                print(f"    witness {s}: present, exit 0  -> PASS")
    else:
        # OPEN: NO fabricated PASS. Negative-route controls (if any) are checked for
        # presence+exit-0 only; the sub-case status stays OPEN.
        if scripts:
            for s in scripts:
                rc, present = run_script(s)
                if not present:
                    print(f"    (negative-route control absent: {s})")
                elif rc != 0:
                    print(f"    (negative-route control {s}: exit {rc} -- non-zero)")
                    all_ok = False
                else:
                    print(f"    negative-route control {s}: present, exit 0 (route closed negatively; sub-case OPEN)")
        else:
            print(f"    (no closure script -- correctly carried as OPEN, no PASS fabricated)")

# §13.4.1 re-attempt-prevention negative controls (presence + exit-0).
print("\n" + "-" * 78)
print("§13.4.1 negative-control attack-vector scripts (re-attempt prevention):")
for s in NEG_CONTROL:
    rc, present = run_script(s)
    if not present:
        print(f"    MISSING: {s}  -> FAIL")
        all_ok = False
    elif rc != 0:
        print(f"    {s}: exit {rc}  -> FAIL")
        all_ok = False
    else:
        print(f"    {s}: present, exit 0  -> PASS (route closed negatively)")

print("\n" + "=" * 78)
print(f"SUMMARY: {closed_count} CLOSED/in-P rows (all witnessing scripts present & exit 0); "
      f"{open_count} OPEN rows (no fabricated PASS).")
print("  Open residual = { OP2(a)-discrete (Conj 6.8d), OP2(b)-discrete-APX,")
print("                    OP4 Type-II approx, OP4 E_12 (Conj 4.3e), OP4 fixed-byte }.")
print(f"OVERALL -> {'PASS' if all_ok else 'FAIL'}")
sys.exit(0 if all_ok else 1)

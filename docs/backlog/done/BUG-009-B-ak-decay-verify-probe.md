# BUG-009-B - ak-decay-verify-probe

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§7.34: Remark 7.34m$''$ E8 Chebyshev-tail closure; Remark 7.34o; Remark 7.34q); scripts/verify/
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
> (Filed gap 009 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at section 12.7 decomposition of BUG-009 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: §7.34 (Remark 7.34m$''$ / Remark 7.34o / Remark 7.34q).

The achievability residual reduces (Remark 7.34m$''$, "All-$D$ via a Chebyshev
tail") to the numerically-observed coefficient decay $a_k(D)=O(D^{k+1})$ for
$k\ge2$ (alternating signs, geometric ratio $\approx0.9D$, $a_2=2\kappa_A^2 D^3$),
which makes the $u$-independent lower bound
$\mathrm{LB}(D)=4a_2-\sum_{k\ge3}|a_k|\tfrac13 k^2(k^2-1)$ positive on the memory
bulk. The paper currently cites the probe
`probe_7_34m_convexity_endpoint.py` (E1–E8) and the beyond-Gray probe
`probe_7_34o_beyond_gray.py` (O1–O3) and the distortion-type probe
`probe_7_34q_lee_general_distortion.py` (Q1–Q4). What this leaf supplies is the
EXECUTABLE-VERIFICATION coverage for the precise quantity the hard sibling
(BUG-009-D) must bound rigorously: a dedicated probe that measures the
Chebyshev-coefficient decay ORDER (the $D^{k+1}$ scaling and the alternating-sign
cancellation that the crude Bauer–Fike / Bernstein bound misses by one $D$-power),
and a consolidated check that every E/O/Q assertion the achievability claim depends
on is present, PASSes, and is covered by the CI `verify` job (run LOCALLY — `gh`
cannot observe this repo's CI per the account-mismatch note).

## Clarification / Assumptions
- This is an IN-REPO, CLOSEABLE-NOW verification leaf. It does NOT prove the
  $a_k$ bound (that is the analytic capstone, sibling BUG-009-D); it makes the
  numerical SIGN of the residual reproducible and CI-covered, and isolates the
  exact $D^{k+1}$ cancellation quantity so the hard leaf has a clean numerical
  target / falsifier.
- GOVERNANCE/DOCUMENTATION-class defect in `tex/rnr_coding.tex`. Per project rules
  every theorem/remark must have an executable PASS/FAIL check in `scripts/verify/`
  included in CI.
- Severity inherits parent MEDIUM.

## Refined Description
**Scope:**
  1. Add (or extend the existing 7.34m$''$ probe with) an explicit `a_k`-decay
     diagnostic: fit the scaling exponent of $|a_k(D)|$ vs $D$ for $k=2,\dots,6$
     (expect order $k+1$: 2.98, 3.89, 4.74, ... from the journal), confirm the
     alternating signs ($a_2>0,a_3<0,a_4>0,\dots$) and the geometric ratio
     $|a_{k+1}/a_k|\approx0.9D$, and EXPLICITLY contrast against the crude
     Bauer–Fike/Bernstein majorant $|a_k|\le 2M e^{-k\tau^*}$ with
     $\tau^*\approx\ln(1/D)$ to demonstrate the one-$D$-power looseness
     (the documented negative finding) — so the probe witnesses BOTH the true
     cancellation and why the crude bound fails. Use high precision (mpmath,
     45–60 dps) since the margin is $\sim D^3$ and float-diff cannot resolve it.
  2. Confirm `probe_7_34m_convexity_endpoint.py` E1–E8, `probe_7_34o_beyond_gray.py`
     O1–O3, `probe_7_34q_lee_general_distortion.py` Q1–Q4 all exist and PASS on
     live HEAD; add a one-line consolidated index/status assertion if absent.
  3. Confirm the CI `verify` job (.github/workflows) covers these probes; verify
     LOCALLY via `/Users/para/.venvs/rnr/bin/python`.
**Non-goals:** the rigorous analytic $a_k$ bound or beyond-Gray exact
$V_{\mathrm{op}}$ (sibling BUG-009-D); the prose status ledger (sibling BUG-009-A);
the beyond-Gray finite-$n$ map sharpening (sibling BUG-009-C).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-009; benefits from (not blocked by) BUG-009-A's ledger
to know which rows each probe witnesses.
**Risks / Open Questions:** the $a_k$ extraction is numerically delicate — the
journal records that a naive 4-point float $g_4$ ($s^4$) extraction gave a SPURIOUS
$g_4<v/12$ at strongly-asymmetric corners (sign-flipping non-monotone = noise tell);
use `mpmath` `mp.diff`, and use the 5×5 quotient operator (not the full $A^3$, on
which `mp.eig` QR failed at 70 dps). Watch DCT aliasing in any Chebyshev transform
(use enough nodes). The probe must MEASURE the residual's sign, not assert closure.
**Expected Code Change:** YES (new/extended `scripts/verify/` probe + possibly a
one-line tex cross-reference to the new probe name)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | confirm probe target = a_k decay ORDER + crude-bound contrast, not a closure |
| ESTIMATION | PENDING | ABOVE_EASY: reuse 7.34m'' machinery; mpmath high-precision a_k fit |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | a_k-decay probe (scaling fit + sign + ratio + Bauer–Fike contrast) |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | add/extend scripts/verify probe; tex cross-ref if new file |
| TESTING | PENDING | PASS/FAIL in stdout; included in CI verify job (checked LOCALLY) |
| POST_AUDIT | PENDING | cold re-run + reproduce independently; confirm no float-precision artifact |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 DECOMPOSITION child of BUG-009. Closeable-now executable-verification leaf;
  isolates the exact $D^{k+1}$ cancellation quantity so the research residual
  (BUG-009-D) has a reproducible numerical target and falsifier.
- §12.1.1 (Backlog-File-First): created fresh this turn; no paper change precedes
  this file. `Recovered: N/A`.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Verify LOCALLY (`/Users/para/.venvs/rnr/bin/python`); `gh` cannot observe this
  repo's CI (account mismatch) — confirm the green CI gate by reading the workflow,
  not by querying remote runs.

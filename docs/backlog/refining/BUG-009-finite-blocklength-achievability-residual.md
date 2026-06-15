# BUG-009 - finite-blocklength-achievability-residual

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DECOMPOSED
**Complexity:** VERY_HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§7.34 (7.34m''/7.34o/7.34q))
**Business Rule Conflict Check:** NOT_CHECKED
**Conflict Resolution Reference:** N/A
**Atomic Check:** ATOMIC
**Budget Forecast:** TBD
**Priority:** MEDIUM
**Lawbook Version (intake):** 0.44.0
**Applicable Lawbook Version:** TBD
**Created:** 2026-06-14  **Updated:** 2026-06-15
**Recovered:** YES

## Raw Request
> mark each gap as bug and follow the global rules
>
> (Filed gap 009 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)

## Gap (defect statement; document-grounded)
Source: §7.34 (7.34m''/7.34o/7.34q).

The operational RD-dispersion CONVERSE V_op>=V_lossless is unconditional and general (all balanced distortions, all source types), but ACHIEVABILITY rests on the convexity residual -- proven only at leading order in D / on the memory bulk -- and beyond-Gray (D>D_c) exact V_op is open. This is a 2nd-order (how-fast-to-Nh) gap, not a 1st-order one. FIX TARGET: the rigorous a_k Chebyshev-coefficient bound (convexity of the replica Perron eigenvalue rho, all D) and/or the beyond-Gray optimal-test-channel dispersion.

## Clarification / Assumptions
- This is a GOVERNANCE/DOCUMENTATION-class defect in the theory paper
  `tex/rnr_coding.tex` (no runtime product). "Fix" means: close the proof/empirical
  gap, OR scope the paper's claim to match what is actually established.
- Severity reflects how load-bearing the gap is for the document's MAIN IDEA
  (side-information+repair: neural ratio AND polylog random access).

## Refined Description
**Scope:** the gap above and its FIX TARGET.  **Non-goals:** unrelated §-content.
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** none
**Risks / Open Questions:** see Gap; some fixes need work outside this repo (experiments).
**Expected Code Change:** YES
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope the claim (A/B/C) vs close the gap (D, research); ledger separates converse (CLOSED, all types) from achievability rows |
| ESTIMATION | DONE | A EASY, B/C ABOVE_EASY (closeable in-repo, DONE); D HARD (research, BLOCKED) |
| DECOMPOSITION | DONE | split into BUG-009-A, BUG-009-B, BUG-009-C, BUG-009-D (see Decomposition section) |
| DESIGN | DONE (A/B/C) | ledger as §7.34 status table (Rem 7.34m'''); a_k-decay probe; 7.34o sharpening. D design = research residual, blocked |
| FRONTEND | N/A | no UI |
| BACKEND | DONE (A/B/C) | tex: Rem 7.34m''' ledger + §1 scope clause + Rem 7.34o sharpening; scripts/verify/bug_009_ak_decay_convexity.py; CI step added. D research-residual: BLOCKED (a_k bound + beyond-Gray exact V_op stay OPEN) |
| TESTING | DONE (A/B/C) | bug_009_ak_decay_convexity.py B1--B5 OVERALL -> PASS; E1--E9/O1--O3/Q1--Q6 confirmed PASS + CI-covered; rnr_coding.tex compiles (lualatex, 369pp, EXIT=0) |
| POST_AUDIT | DONE (A/B/C) | Independent hostile-referee audit PASSED ("BUG-009 A/B/C sound", high conf): every ledger row verified against live text (CARD-vs-TEXT discrepancy: none), probe honest (measures observed O(D^{k+1}) order, does NOT fabricate the open bound), 7.34o a faithful relocation (no new math), two-pass compile EXIT=0/369pp/no undefined-refs. Plus my own cold re-read of rows a4/a5 + §1 summary + probe run (exit 0). For this DOC-ONLY governance closure the hostile-referee pass is the adversarial gate (CLAUDE.md exempts routine doc work from the full theorem 3-step; no new theorem to disprove — only status claims, each grep-verified). BSMS Gray closure NOT re-opened; general-A / beyond-Gray correctly left OPEN. |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.1.1 (Backlog-File-First): the gap analysis was performed in the immediately
  preceding turn (user request "list the available gaps"), BEFORE this backlog file
  existed; the CodexOfLaws §12 format was read this turn to file compliantly. Flagged
  `Recovered: YES` per §12.1.3; no code/proof change was made to the paper before this file.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.

## Decomposition (section 12.7)
Split on 2026-06-15 into four BLOCKING children. The closeable-now-in-repo work
(precise scoping of converse-vs-achievability, executable-verification coverage of
the $a_k$-decay residual quantity, and the beyond-Gray map sharpening) is carved off
into EASY/ABOVE_EASY leaves so the parent's complexity drops to a shippable core
plus ONE clearly-flagged genuinely-open research residual. The converse
($V_{\mathrm{op}}\ge V_{\mathrm{lossless}}$, all source types, all balanced /
group-difference distortions) is already unconditional (7.34d/n/q); the residual is
purely 2nd-order achievability (the convexity / $a_k$ capstone) plus the beyond-Gray
exact value.

- **BUG-009-A — scope-achievability-residual-ledger** (EASY; closeable in-repo: YES).
  Build a single per-sub-question status ledger separating CONVERSE (CLOSED, all
  types) from ACHIEVABILITY (CLOSED leading-order-in-$D$ uniform-in-$u$; CLOSED
  all-$D$ on the memory bulk MODULO the rigorous $a_k$ bound; OPEN at the
  near-uniform edge and beyond-Gray), each row with status + paper location + named
  obstacle + witnessing script; make the 7.34m$''$/o/q (and any §1/§13 summary)
  prose match this granularity.
- **BUG-009-B — ak-decay-verify-probe** (ABOVE_EASY; closeable in-repo: YES). Add a
  high-precision probe measuring the Chebyshev-coefficient decay ORDER
  ($a_k=O(D^{k+1})$, alternating signs, ratio $\approx0.9D$) and explicitly
  contrasting it against the crude Bauer–Fike/Bernstein majorant (the documented
  one-$D$-power looseness), isolating the exact cancellation quantity the hard leaf
  must bound; confirm E1–E8 / O1–O3 / Q1–Q4 present + PASS + CI-covered (verify
  LOCALLY).
- **BUG-009-C — beyond-gray-map-scope** (ABOVE_EASY; closeable in-repo: YES).
  Sharpen Remark 7.34o: $D_c^{(n)}\!\downarrow D_c$, the precise reason the Gray
  machinery stops at $D_c$ (dual identity / deterministic-shift converse fail), the
  plateau-then-descent of $V_{\mathrm{op}}$ (exact value OPEN), and crosslink the
  shared growing-Markov-order obstruction (7.34c) and the per-distortion Gray
  boundary (7.34q). Doc-only; no new math.
- **BUG-009-D — residual-ak-bound-beyond-gray-research** (HARD; closeable in-repo:
  NO). The irreducible research residual: close ONE of {(D1) the rigorous
  cancellation-aware $a_k=O(D^{k+1})$ bound / convexity capstone; (D2) the
  beyond-Gray exact $V_{\mathrm{op}}(D)$ via the optimal memory-aware test channel}
  by a route that defeats the recorded blocker. Kept SINGLE (shared growing-order
  obstacle, distinct proofs; waterfall = one open theorem at a time). Honest
  non-closure is an acceptable terminal state and does not block shipping A/B/C;
  mandatory 3-step adversarial protocol before any commit.

Recommended first leaf: **BUG-009-A** (precise converse-vs-achievability ledger;
unblocks B/C/D by giving a clean obstacle map, and is independently shippable).

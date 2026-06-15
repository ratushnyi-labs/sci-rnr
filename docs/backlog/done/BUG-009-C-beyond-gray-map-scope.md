# BUG-009-C - beyond-gray-map-scope

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§7.34: Remark 7.34o beyond-Gray $D>D_c$; cross-ref Remark 7.34c growing-Markov-order; Remark 7.34q Gray-region boundary)
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
Source: §7.34 (Remark 7.34o — the beyond-Gray regime $D>D_c$).

The exact $V_{\mathrm{op}}(D)$ for $D>D_c$ is open (it needs the optimal non-SLB,
memory-aware test channel), but the QUALITATIVE map is established and is closeable
in-repo at a sharper grade than the headline "open": the per-$n$ deconvolution
threshold $D_c^{(n)}\!\downarrow D_c=\inf_n D_c^{(n)}$ (plateau ends EXACTLY at
$D_c$ in the operational limit); for $D>D_c$ the dual identity $M=\gamma^n P_X$ and
the deterministic-shift $j_n=i_n-nc$ FAIL, so the converse argument dies; with the
Blahut–Arimoto optimal output $V_{\mathrm{op}}(D)$ is flat $=V_{\mathrm{lossless}}$
on $[0,D_c]$ and DECREASES below it for $D>D_c$, continuously at $D_c$. The memory
journal (`blackwell_ac_lead_dead_skewproduct_pivot.md`, §7.34o + the off-Gray
$\rho(D)$ characterization) further records the UNIFYING insight that the SAME
growing-Markov-order obstruction (optimal output needs growing Markov order;
Eswaran–Gastpar 2022) walls BOTH in-Gray achievability AND the off-Gray
converse-ratio's non-closed-form $\rho(D)=\lim\mathrm{Var}(j_n)/\mathrm{Var}(i_n)$,
and that NO elementary closed form for $\rho(D)$ exists (don't re-hunt one). This
leaf's job is to ensure Remark 7.34o states the map, the precise reason the Gray
machinery stops at $D_c$, and the exact-value residual at this granularity — and to
cross-link the growing-order obstruction shared with 7.34c — WITHOUT attempting the
exact $V_{\mathrm{op}}(D)$.

## Clarification / Assumptions
- This is an IN-REPO, CLOSEABLE-NOW scoping/sharpening leaf for the beyond-Gray
  sub-question. The finite-$n$ map and the plateau-then-descent are already
  validated (probe `probe_7_34o_beyond_gray.py` O1–O3); this leaf scopes the
  CLAIM precisely and crosslinks the shared obstacle. The exact $V_{\mathrm{op}}(D)$
  for $D>D_c$ remains the genuinely-open residual (sibling BUG-009-D).
- GOVERNANCE/DOCUMENTATION-class defect in `tex/rnr_coding.tex`. "Fix" = precise
  claim + obstacle crosslink, NOT closing the open value.
- Severity inherits parent MEDIUM. Beyond-Gray is a regime where the headline
  Gray-region results explicitly do not apply; scoping it correctly protects the
  in-Gray results from over-extension.

## Refined Description
**Scope:**
  1. Confirm Remark 7.34o states, at the right granularity: (a) $D_c^{(n)}\downarrow
     D_c$ monotone, plateau ends exactly at $D_c$ in the $n\to\infty$ limit;
     (b) the precise mechanism by which the Gray machinery stops ($P_{Y^*}$ not a
     valid distribution, $K_\nu K^{-1}=\gamma I$ / dual identity fails,
     deterministic-shift converse no longer applies); (c) $V_{\mathrm{op}}$ flat on
     $[0,D_c]$, decreasing and continuous beyond — exact value OPEN.
  2. Cross-link the UNIFYING growing-Markov-order obstruction (optimal output needs
     growing Markov order, Eswaran–Gastpar 2022) shared with Remark 7.34c, noting it
     walls both the in-Gray achievability capstone AND the off-Gray exact value; and
     record the established negative that $\rho(D)=\lim\mathrm{Var}(j_n)/
     \mathrm{Var}(i_n)$ has NO elementary closed form (it is the growing-order
     Green–Kubo / $\chi''(0)$ rate of the $\lambda^*$-tilted joint
     (source, optimal-output) transfer operator) — so the residual is the exact
     test channel, not a missing closed-form.
  3. Make sure the beyond-Gray boundary statement is consistent with Remark 7.34q
     (each balanced / group-difference distortion has its OWN Gray region with its
     own $D_c$; the beyond-Gray map is per-distortion).
**Non-goals:** computing the exact $V_{\mathrm{op}}(D)$ for $D>D_c$, or any attempt
at the optimal memory-aware test channel (sibling BUG-009-D); the in-Gray $a_k$
convexity bound (BUG-009-D / BUG-009-B); the converse-side prose ledger (BUG-009-A).
Do NOT re-hunt an elementary closed form for $\rho(D)$ (journal: provably none).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-009; coordinate with BUG-009-A so the beyond-Gray rows
of the converse/achievability ledger and the 7.34o prose agree.
**Risks / Open Questions:** the descent numbers in 7.34o are finite-$n$ (e.g.
$p=0.2,n=7$); state them as illustrative of the plateau-then-descent SHAPE, not as
the operational $V_{\mathrm{op}}(D)$ (which is the $n\to\infty$ limit). Do not let
the sharpening imply the exact curve is known.
**Expected Code Change:** YES (tex prose / cross-reference edit only; no math; the
O1–O3 probe already exists)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | confirm scope = map + obstacle crosslink, exact value stays open |
| ESTIMATION | PENDING | ABOVE_EASY: prose precision + crosslink; probe already present |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | tighten 7.34o; crosslink 7.34c growing-order + 7.34q per-distortion boundary |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex prose edit in §7.34o (+ cross-refs); no new math |
| TESTING | PENDING | none new (O1–O3 already PASS); confirm rows match BUG-009-A ledger |
| POST_AUDIT | PENDING | cold re-read + adversarial check that the map is not overstated as the exact value |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 DECOMPOSITION child of BUG-009. Closeable-now scoping/crosslink leaf for the
  beyond-Gray sub-question; separates the precise map (in-repo) from the open exact
  value (research residual, BUG-009-D).
- §12.1.1 (Backlog-File-First): created fresh this turn; no paper change precedes
  this file. `Recovered: N/A`.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Per project rules (Truth > peremoga): the finite-$n$ descent is illustrative, not
  the operational $V_{\mathrm{op}}(D)$; do not record the beyond-Gray exact value as
  anything but OPEN.

# BUG-009-A - scope-achievability-residual-ledger

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (§7.34: Remark 7.34m$''$, Remark 7.34o, Remark 7.34q; Theorems 7.34m/7.34n)
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

The parent gap states that the RD-dispersion CONVERSE $V_{\mathrm{op}}\ge
V_{\mathrm{lossless}}$ is unconditional and general (all balanced / group-difference
distortions, all source types), but ACHIEVABILITY rests on the convexity residual
of Remark 7.34m$''$ (the binding $s=\pi$ point) which is proven only (i) at leading
order in $D$ uniformly in $u$, and (ii) on the memory bulk ($p$ in the memory
regime) via the Chebyshev-tail lower bound $\mathrm{LB}(D)>0$ — and that all-$D$
bulk closure is itself stated "modulo making the $O(D^{k+1})$ coefficient bounds
rigorous". Separately, beyond-Gray ($D>D_c$) the exact $V_{\mathrm{op}}(D)$ is OPEN
(Remark 7.34o). What is MISSING is a single precise ledger that states, per
sub-question, exactly what is CLOSED (converse all types; achievability
leading-order-in-$D$ uniform-in-$u$; achievability all-$D$ on the memory bulk
modulo $a_k$), what is OPEN (the rigorous cancellation-aware $a_k$ bound; the
near-uniform / memoryless edge $p\to(A-1)/A$; the beyond-Gray exact $V_{\mathrm{op}}$),
the named structural obstacle for each open item, and where in the paper / which
verify script witnesses each entry. Absent this ledger the headline could read as a
flat "achievability open" when in fact it is a small, precisely-pinned residual set
sitting on top of an unconditionally-proven converse and a memory-bulk achievability
closure.

## Clarification / Assumptions
- This is the IN-REPO, CLOSEABLE-NOW scoping leaf of BUG-009. It neither proves
  the rigorous $a_k$ bound nor closes beyond-Gray; it scopes the paper's
  achievability claim to match exactly what is established, so the headline
  complexity of the parent drops to a shippable core (this audit) plus a
  clearly-flagged research residual (sibling BUG-009-D).
- GOVERNANCE/DOCUMENTATION-class defect in the theory paper `tex/rnr_coding.tex`
  (no runtime product). "Fix" here = make the claim precise / build the ledger,
  NOT prove an open lemma.
- Severity inherits parent MEDIUM. The 2nd-order achievability (how fast to $nR(D)$)
  is a refinement; the 1st-order rate $R(D)$ and the converse dispersion bound do
  NOT depend on this residual, so the side-information+repair main idea is not
  load-bearing on it.

## Refined Description
**Scope:** Produce a single per-sub-question status ledger for §7.34 RD-dispersion,
distinguishing CONVERSE from ACHIEVABILITY and tabulating, per row, {status, paper
location, named obstacle if open, witnessing verify script}:
  - CONVERSE $V_{\mathrm{op}}\ge V_{\mathrm{lossless}}$: CLOSED unconditionally, all
    source types (7.34d/i/m/n), all balanced / group-difference distortions on each
    distortion's own Gray region (7.34q); source-agnostic dual identity
    $K_\nu K^{-1}=\gamma I$.
  - ACHIEVABILITY, in-Gray, leading order in $D$: CLOSED, uniform in $u$, constant
    $8\kappa_A^2 D^3$ (7.34m$''$); the prefactor $|C_s/C_0|^2$ is affine in $u$
    (carries no curvature) so the residual is purely spectral ($\rho(W_A)$).
  - ACHIEVABILITY, in-Gray, all $D\le D_c$ on the MEMORY BULK: CLOSED MODULO the
    rigorous $a_k=O(D^{k+1})$ Chebyshev-coefficient bound (the cancellation-aware
    analyticity capstone); the crude Bauer–Fike/Bernstein bound is one $D$-power too
    loose (obstacle = the harmonic cancellation, $a_2\sim D^3$ not $D^2$).
  - ACHIEVABILITY, near-uniform / memoryless edge $p\to(A-1)/A$: OPEN (crude
    $T_k''(1)$ over-counts; true $G''>0$ via alternating-$a_k$ cancellation) — but
    this is the $V_{\mathrm{lossless}}\to0$ limit, outside the memory regime of
    interest.
  - ACHIEVABILITY & exact $V_{\mathrm{op}}$, beyond-Gray $D>D_c$: OPEN (needs the
    optimal non-SLB memory-aware test channel); finite-$n$ map (plateau-then-descent,
    continuous at $D_c$) is established (7.34o).
Then confirm the §7.34 prose (the "Status" paragraph of 7.34m$''$; 7.34o; 7.34q;
and any §1/§13 summary that references the dispersion result) states the residual at
exactly this granularity — no blanket "achievability open" overclaim, no
understatement of the memory-bulk closure.
**Non-goals:** the rigorous $a_k$ bound or beyond-Gray exact $V_{\mathrm{op}}$
(sibling BUG-009-D); verification-script consolidation / the new $a_k$-decay probe
(sibling BUG-009-B); the beyond-Gray finite-$n$ map sharpening (sibling BUG-009-C).
**Impacted UCs:** N/A (theory paper)
**Impacted BR/WF:** N/A
**Dependencies:** parent BUG-009.
**Risks / Open Questions:** the ledger must not record the all-$D$ memory-bulk
achievability as a clean THEOREM — it is explicitly "modulo making the
$O(D^{k+1})$ coefficient bounds rigorous" (memory journal
`blackwell_ac_lead_dead_skewproduct_pivot.md`: the crude Bauer–Fike bound FAILS;
the cancellation $a_k\sim D^{k+1}$ is essential). Cold re-read against live HEAD
before asserting any row CLOSED. (Note: the SYMMETRIC operational dispersion
$V_{\mathrm{op}}=V_{\mathrm{lossless}}$ on the closed Gray region was separately
protocol-closed via the quenched-tail + KP + GAP-1 + R2 chain — verify which of the
A-ary / asymmetric achievability rows that closure already covers vs which inherit
the convexity residual, and record both faithfully.)
**Expected Code Change:** YES (tex prose / status table edit only; no math)
**Documentation Recovery Required:** NO

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | PENDING | confirm ledger granularity = converse-vs-achievability × {leading-order / all-D-bulk / edge / beyond-Gray} |
| ESTIMATION | PENDING | EASY: read 7.34m''/o/q, tabulate existing closures + obstacles |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | ledger as a §7.34 status table (mirror existing status-paragraph prose) |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | add/refresh the per-sub-question status table in §7.34 |
| TESTING | PENDING | none needed (doc-only); rows cross-ref verify scripts owned by sibling B |
| POST_AUDIT | PENDING | cold re-read + adversarial re-review per project review protocol (each CLOSED row must survive "Attack failed; holds") |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7 DECOMPOSITION child of BUG-009. This is the honestly-closeable-now, in-repo
  scoping leaf: it lowers the parent's effective complexity by separating the
  precise-claim work from the genuinely-open research residual.
- §12.1.1 (Backlog-File-First): created fresh this turn as a decomposition child;
  no paper change precedes this file. `Recovered: N/A`.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Per project rules (Truth > peremoga): do NOT record any sub-question as CLOSED
  unless the cited remark/theorem + verify script actually establish it on live HEAD;
  the all-$D$ memory-bulk row is "closed MODULO $a_k$", not closed.

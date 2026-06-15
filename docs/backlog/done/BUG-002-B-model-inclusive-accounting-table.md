# BUG-002-B - model-inclusive-accounting-table

**Type:** BUG
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** DONE
**Complexity:** ABOVE_EASY
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Theorem 7.15 caveat ll.8573-8606; Theorem 7.21 total-cost formula ll.9750-9814; §10.2 Model amortization regimes ll.15573-15582)
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
Source: parent BUG-002 FIX TARGET, clause "present a **model-inclusive accounting
alongside the amortized one**."

The model-inclusive ingredients exist but are SCATTERED across three locations and
never assembled into one side-by-side accounting:
- **Amortized / bitstream-only:** $N\cdot H_{M'}\approx79\,\mathrm{MB}$, $12.0\times$
  (T7.15 boxed bound + worked example, ll.8441-8515).
- **Self-contained / model-inclusive (single archive):** $\approx140\,\mathrm{GB}$
  model + $1\,\mathrm{GB}$ source $\Rightarrow$ negative savings, $\approx140\times$
  expansion (T7.15 caveat, ll.8581-8586).
- **Amortized total cost with model encoded:** $\mathrm{TotalCost}(L,N)=L+N\cdot h(X)+O(\sqrt N)$,
  optimal $L^*(N)=(NC\alpha)^{1/(1+\alpha)}$ (T7.21, ll.9765-9782).
- **Regime boundaries:** single-block / multi-block / streaming, $L(M)/m$ (§10.2, ll.15575-15582).

No single table puts the amortized bitstream figure, the single-archive
self-contained figure, the break-even volume $V^*\approx152\,\mathrm{GB}$, and the
per-regime amortized cost side by side, so a reader cannot see the full accounting
in one place.

## Refined Description
**One-leaf scope:** synthesize the existing scattered numbers into ONE consolidated
"model-inclusive vs amortized accounting" table (placed at the T7.15 caveat or
§10.2), with columns for: regime (single-archive / multi-archive volume $V$ /
streaming $V\to\infty$); bitstream cost ($N\cdot H_{M'}$); model cost included
(yes/no, and the $\approx140\,\mathrm{GB}$ or distilled $L^*$ figure); net
savings / expansion factor; and the governing break-even condition. The table
collates and cross-references EXISTING established quantities (T7.15, T7.21,
§10.2); it introduces NO new claim and re-derives NO number from scratch.

**Fix target:** add the consolidated accounting table + a one-paragraph reading
guide; cross-link T7.21's $L^*(N)$ as the "distilled-model" row of the table.

**Closeable in-repo now:** YES (paper-scoping / doc synthesis only). It is
ABOVE_EASY rather than EASY because it requires carefully reconciling units
(bytes vs bits, GB vs GiB) across three sources and laying out a correct table.

**Non-goals:** the cross-link pass at every occurrence (BUG-002-A); the executable
arithmetic check of the table's numbers (BUG-002-C); empirical measurement (BUG-002-D).
**Impacted UCs:** N/A (theory paper)  **Impacted BR/WF:** N/A
**Dependencies:** parent BUG-002. Naturally precedes BUG-002-A (gives the canonical
table to point to) and precedes BUG-002-C (the script checks this table's cells).
**Risks / Open Questions:** unit-reconciliation errors across the three sources;
the break-even $V^*$ figure must agree with what BUG-002-C recomputes.
NOTE: actual tex/rnr_coding.tex edits are made by the supervisor on reconciliation.
**Expected Code Change:** YES (tex table + prose, supervisor-applied)
**Documentation Recovery Required:** NO

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class defect; "fix" = present the accounting precisely
  by collating already-established quantities, not by establishing new ones.
- $H_{M'}\approx0.664$ bpb, $\approx140\,\mathrm{GB}$ model, $V^*\approx152\,\mathrm{GB}$
  are inherited from T7.15/T7.21 as the source of truth.

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = one consolidated accounting table, collation only |
| ESTIMATION | DONE | ABOVE_EASY (unit reconciliation across 3 sources) |
| DECOMPOSITION | N/A | atomic leaf |
| DESIGN | PENDING | table columns/rows + placement (T7.15 caveat vs §10.2) |
| FRONTEND | N/A | no UI |
| BACKEND | PENDING | tex table + reading-guide paragraph (supervisor-applied) |
| TESTING | PENDING | numbers cross-checked against BUG-002-C probe output |
| POST_AUDIT | PENDING | adversarial re-read: are all four regimes + break-even row present and unit-consistent? |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7: created fresh as a BLOCKING child leaf of BUG-002 on 2026-06-15; not recovered.
- Class: GOVERNANCE/DOCUMENTATION; §6 runtime stages / §7 coverage are class-level N/A.
- Stays in `docs/backlog/intake/`; no tex/ edits performed in this turn.


## Resolution (2026-06-15, opus-4-8; POST_AUDIT passed)
Closed in-repo as part of BUG-002 A/B/C (commit pending). Consolidated model-inclusive-vs-amortized table in §10.2 + break-even V*~152 GB crosslinks at Abstract & §10.2 + arithmetic probe scripts/verify/bug_002_model_inclusive_accounting.py (OVERALL -> PASS) + CI job. Collation only, no new claims; hostile-referee POST_AUDIT clean.

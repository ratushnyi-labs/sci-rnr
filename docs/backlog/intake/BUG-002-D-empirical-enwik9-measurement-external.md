# BUG-002-D - empirical-enwik9-measurement-external

**Type:** TASK
**Entry Mode:** CODE_FIRST
**User Mode:** EXPERT
**Guided Flow Stage:** N/A
**Status:** INTAKE
**Complexity:** HARD
**Audit Type:** DOCUMENTATION_AUDIT
**Scan Depth:** N/A
**Audit Scope:** tex/rnr_coding.tex (Theorem 7.15 numerical prediction ll.8463-8515; §11.7 ll.16814-16835); tex/rnr_experimental_design.tex (H15 §2.11 ll.292-315; §5.1 model-inclusive metric l.519)
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
> (Filed gap 002 of the document gap-analysis (preceding review) from the immediately
> preceding turn; the gap text below is that analysis, not the user's words.)
>
> (Child sub-task created at §12.7 decomposition of BUG-002 on 2026-06-15; scope below is this leaf only.)

## Gap (defect statement; document-grounded)
Source: parent BUG-002 — the empirical core of the FIX TARGET that CANNOT be closed
in this repo. The headline $\approx79\,\mathrm{MB}$ / $12.0\times$ figure is a
*prediction* derived from a *cited third-party* per-position cross-entropy
($H_{M'}\approx0.664$ bpb, Delétang et al. 2024 Table 1), NOT a measured RNR
end-to-end archive. The genuinely-external work is to actually build an RNR archive
with a Chinchilla-class predictor on enwik9 and MEASURE both the bitstream size and
the model-inclusive (self-contained) size, confirming or falsifying H15 (experimental
-design §2.11, ll.294-315) including falsification clause (b): the
$\approx140\,\mathrm{GB}$ self-contained expansion at single-document scale.

This requires a frontier-class LM, bit-exact integer-arithmetic inference (Theorem
10.1/10.6 conformance), a reference RNR encoder/decoder, and GPU compute — none of
which exists in this theory repo (no runtime product; companion artifact deferred
to §13.3).

## Refined Description
**One-leaf scope (the in-repo half is DELIBERATELY carved out to BUG-002-A/B/C):**
this leaf is the IRREDUCIBLE external-measurement residual. Its only in-repo action
is to ensure the claim is precisely scoped as a *prediction* and cross-linked to
the pre-registered experimental-design hypothesis H15 (which already exists at §2.11
and already lists the model-inclusive metric at §5.1 l.519). The MEASUREMENT itself
— a real RNR end-to-end run on enwik9 with a frontier predictor, reporting bitstream
and self-contained sizes with confidence bands — needs external compute and a
reference implementation and CANNOT be performed in this repo.

**Fix target:** (in-repo, trivial) confirm Abstract/§11.7/T7.15 each label the
number a *prediction* and point to H15; (external, the real work) execute H15 and
report measured bitstream + model-inclusive figures.

**Closeable in-repo now:** NO — flagged as needing EXTERNAL work (frontier LM,
reference RNR implementation, GPU compute). The in-repo scoping/cross-link sliver is
already covered by BUG-002-A; this leaf exists to make the external residual
explicit and separately trackable so the parent's shippable core (A/B/C) is not
blocked on it.

**Non-goals:** any of the in-repo doc/table/probe work (BUG-002-A/B/C). This leaf
must NOT be marked done by a paper edit; only a real measurement (or an explicit
project decision to leave H15 pre-registered-but-unrun) closes it.
**Impacted UCs:** N/A (theory paper)  **Impacted BR/WF:** N/A
**Dependencies:** parent BUG-002; the §13.3 companion reference-implementation
artifact; experimental-design H15 (§2.11).
**Risks / Open Questions:** genuinely-hard / out-of-repo; depends on a reference
implementation that does not yet exist. Honest status: this stays OPEN as the
external residual after A/B/C close the in-repo core.
**Expected Code Change:** NO (in this repo) — measurement output lives outside.
**Documentation Recovery Required:** NO

## Clarification / Assumptions
- GOVERNANCE/DOCUMENTATION-class project with no runtime product; runtime/experiment
  execution is class-level out-of-repo.
- Per parent's Clarification: "some fixes need work outside this repo (experiments)."
  This leaf IS that part, isolated so the closeable-now leaves can ship independently.
- The honest in-repo posture is "claim scoped as prediction + cross-linked to H15";
  the empirical confirmation is external and may remain open indefinitely.

## Lifecycle Coverage Map
| Stage | Status | Notes |
|---|---|---|
| CLARIFICATION | DONE | scope = external measurement residual; in-repo sliver delegated to BUG-002-A |
| ESTIMATION | DONE | HARD (frontier LM + reference impl + GPU; out-of-repo) |
| DECOMPOSITION | N/A | atomic leaf; cannot be split further without an actual implementation |
| DESIGN | N/A | already designed as H15 in experimental-design §2.11 |
| FRONTEND | N/A | no UI |
| BACKEND | BLOCKED | needs §13.3 reference implementation (does not exist) |
| TESTING | BLOCKED | the measurement IS the test; needs external compute |
| POST_AUDIT | PENDING | on any future run: check measured self-contained size matches H15(b) |

## Execution Tracking (§12.11)
**Estimate (hours, before BUILD):** TBD (external; unbounded without reference impl)
**Start Timestamp:** TBD
**End Timestamp:** TBD
**Token Stats (per work session):** input: UNKNOWN (not yet worked) ; output: UNKNOWN ; model: UNKNOWN ; level: UNKNOWN

## Compliance Notes
- §12.7: created fresh as a BLOCKING child leaf of BUG-002 on 2026-06-15; not recovered.
- Class: GOVERNANCE/DOCUMENTATION; the measurement is class-level out-of-repo (no
  runtime product). Flagged closeable_in_repo=NO per the parent's experimental caveat.
- The in-repo scoping/cross-link sliver is intentionally NOT duplicated here — it is
  BUG-002-A — so this leaf is a clean external residual.

#!/usr/bin/env python3
"""One-shot correction pass: re-apply the Def-2.1 full-cost guard to the
P-ADV verdicts of an already-written campaign run.

The guard (added to measure_preconditions.measure_corpus during the
first full campaign, before any result was consumed): a P-ADV
"qualifies" from the paired per-block CI is downgraded to "marginal"
when the full-archive ratio_block >= 1, because Definition 2.1 compares
the COMPLETE encoder output with the baseline. Without the guard the
store-mode negative controls (raw fallback = exactly 8.0000 bpb per
block) would "qualify" merely by undercutting the baseline's per-block
framing overhead.

Measured numbers are untouched: only the verdict LABEL in
summary_full.json and the qualification CSVs is corrected, and a
labeled correction note is appended to findings.md (append-only, in
the spirit of the results-store correction convention).
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"


def main() -> int:
    spath = OUT / "summary_full.json"
    summary = json.loads(spath.read_text(encoding="utf-8"))
    corrected = []
    for r in summary["results"]:
        a = r["adv"]
        if a["verdict"] == "qualifies" and a["ratio_block"] >= 1.0:
            a["verdict"] = "marginal"
            a["note"] = ("per-block repair rate beats the baseline, but the "
                         "full-archive Def-2.1 ratio is >= 1 (container "
                         "overhead); downgraded")
            corrected.append(r["corpus"])
    if not corrected:
        print("no verdicts needed correction")
        return 0
    spath.write_text(json.dumps(summary, indent=1, sort_keys=True),
                     encoding="utf-8")

    by = {r["corpus"]: r["adv"]["verdict"] for r in summary["results"]}
    for csv_path in (OUT / "qualification_full.csv", HERE / "qualification.csv"):
        if not csv_path.exists():
            continue
        import csv as _csv

        rows = list(_csv.DictReader(csv_path.open(encoding="utf-8")))
        for row in rows:
            row["P_ADV"] = by[row["corpus"]]
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            w = _csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    findings = HERE / "findings.md"
    note = [
        "",
        "## Correction (same run, before any downstream use)",
        "",
        "The P-ADV verdict rule was tightened to the Def-2.1 full-cost "
        "form — `qualifies` requires BOTH the paired per-block 95% CI "
        "below 0 AND full-archive `ratio_block` < 1 — because the "
        "store-mode negative controls otherwise 'qualify' by "
        "undercutting the baseline's per-block framing overhead with "
        "their raw fallback (8.0000 bpb/block vs zstd-block ~8.001).",
        "",
        "Corrected P-ADV verdicts (measured numbers unchanged): "
        + ", ".join(f"**{c}**: qualifies -> marginal" for c in corrected)
        + ". `qualification.csv` / `out/qualification_full.csv` are "
        "rewritten with the corrected labels; for the corpora listed "
        "here the P-ADV cells of the section-5 table above are "
        "superseded by this note.",
        "",
    ]
    with findings.open("a", encoding="utf-8") as f:
        f.write("\n".join(note))
    print("corrected:", ", ".join(corrected))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

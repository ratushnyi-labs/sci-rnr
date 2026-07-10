#!/usr/bin/env python
"""Consistency gate for experiments/RESULTS-PART-I.md.

Verifies, in under a minute and with no new measurements:
  G1  the report exists and every mandatory section heading is present;
  G2  the verdict table carries a row for every Part-I-scoped hypothesis
      (precondition program, H1-H9, H14) and each row states a verdict
      from the declared vocabulary;
  G3  every source key declared in the report resolves to an existing file;
  G4  every quoted number in the report appears verbatim (substring) in at
      least one of the sources declared for its section or table row --
      the traceability rule: numbers only from the findings files;
  G5  the adverse/open ledger names the mandatory items.

Run: /Users/para/.venvs/rnr/bin/python experiments/run_checks.py
"""

import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "experiments" / "RESULTS-PART-I.md"

SOURCES = {
    "preconditions": ROOT / "experiments" / "preconditions" / "findings.md",
    "h2_h4": ROOT / "experiments" / "h2_h4_typeii" / "findings.md",
    "conformance": ROOT / "experiments" / "conformance" / "findings.md",
    "impl_cm": ROOT / "impl_cm" / "findings.md",
    "h5_h7": ROOT / "experiments" / "h5_h7_modes" / "findings.md",
    "exp-design": ROOT / "tex" / "rnr_experimental_design.tex",
    "core-s11": ROOT / "tex" / "papers" / "core" / "rnr_core.tex",
    "methods": ROOT / "bench" / "METHODS.md",
}
# Targets/thresholds and paper predictions may be cited anywhere.
ALWAYS_OK = {"exp-design", "core-s11", "methods"}

MANDATORY_HEADINGS = [
    "## 1. Hypothesis-by-hypothesis verdict table",
    "## 2. Section-11 prediction reconciliation",
    "## 3. Adverse / open ledger",
    "## 4. Instrument inventory",
    "## 5. Follow-ups",
    "## 6. Proposed section-11 caveat lines",
]

REQUIRED_ROWS = [
    r"\|\s*Precondition screen P-STAT",
    r"\|\s*Precondition P-ADV",
    r"\|\s*H1 ",
    r"\|\s*H2 ",
    r"\|\s*H3 ",
    r"\|\s*H4 ",
    r"\|\s*H5 ",
    r"\|\s*H6 ",
    r"\|\s*H7 ",
    r"\|\s*H8 ",
    r"\|\s*H9 ",
    r"\|\s*H14 ",
]

VERDICT_WORDS = [
    "CONFIRMED",
    "SUPPORTED",
    "INSTRUMENT-SCOPED",
    "BELOW-BAND-NOT-ADVERSE",
    "ADVERSE",
    "UNTESTED-EXTERNAL",
]

LEDGER_MUSTS = [
    "P-ADV",
    "6.904",
    "10.629",
    "H10(c)",
    "2%",
    "telemetry",
    "literal H2",
    "H15",
    "Levina",
]

# Numeric tokens worth tracing: comma-grouped ints, scientific notation,
# decimals, and integers of >= 3 digits. 1-2 digit bare integers (section
# numbers, KiB sizes, percents like 2%) are structural, not measurements.
NUM_RE = re.compile(
    r"\d{1,3}(?:,\d{3})+(?:\.\d+)?"  # 67,584  827,070
    r"|\d+\.\d+e[+-]?\d+"            # 3.6e-15
    r"|\d+e[+-]?\d+"                 # 1e-4
    r"|\d+\.\d+"                     # 0.142  6.904
    r"|\d{3,}"                       # 403  10000  20260710
)

SOURCES_DECL_RE = re.compile(r"<!--\s*sources:\s*([^>]*?)\s*-->")


def parse_keys(text):
    return {k for k in SOURCES if re.search(r"\b" + re.escape(k) + r"\b", text)}


def main():
    t0 = time.time()
    results = []

    def gate(name, ok, detail=""):
        results.append((name, ok, detail))
        print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))

    # ---- G1: report + headings -------------------------------------------
    if not REPORT.exists():
        gate("G1 report-exists", False, f"missing {REPORT}")
        print("FAIL (1 gate)")
        return 1
    text = REPORT.read_text(encoding="utf-8")
    missing = [h for h in MANDATORY_HEADINGS if h not in text]
    gate("G1 mandatory-sections", not missing,
         f"missing: {missing}" if missing else f"{len(MANDATORY_HEADINGS)} headings present")

    # Split into sections keyed by their '## ' heading for scoped checks.
    sections = {}
    current = "_header"
    sections[current] = []
    for line in text.splitlines():
        if line.startswith("## "):
            current = line.strip()
            sections[current] = []
        sections[current].append(line)

    def section_text(prefix):
        for k, v in sections.items():
            if k.startswith(prefix):
                return "\n".join(v)
        return ""

    # ---- G2: verdict table rows + vocabulary -----------------------------
    tbl = section_text("## 1.")
    missing_rows = [p for p in REQUIRED_ROWS if not re.search(p, tbl)]
    bad_verdict = []
    for p in REQUIRED_ROWS:
        m = re.search(p + r".*", tbl)
        if m and not any(w in m.group(0) for w in VERDICT_WORDS):
            bad_verdict.append(p)
    ok2 = not missing_rows and not bad_verdict
    gate("G2 verdict-table-rows", ok2,
         (f"missing rows: {missing_rows} " if missing_rows else "") +
         (f"rows without a vocabulary verdict: {bad_verdict}" if bad_verdict else
          f"{len(REQUIRED_ROWS)} rows, all with vocabulary verdicts"))

    # ---- G3: source declarations resolve ---------------------------------
    declared = set()
    for m in SOURCES_DECL_RE.finditer(text):
        for piece in re.split(r"[,\s]+", m.group(1)):
            if piece:
                declared.add(piece)
    unknown = declared - set(SOURCES)
    missing_files = [k for k in SOURCES if not SOURCES[k].exists()]
    gate("G3 sources-resolve", not unknown and not missing_files,
         f"unknown keys: {unknown}; missing files: {missing_files}"
         if unknown or missing_files else
         f"{len(declared)} keys declared, all resolve; {len(SOURCES)} source files exist")

    # ---- G4: numeric traceability ----------------------------------------
    src_text = {k: p.read_text(encoding="utf-8") for k, p in SOURCES.items() if p.exists()}

    def token_in(tok, keys):
        for k in keys:
            body = src_text.get(k, "")
            if tok in body:
                return True
            if "," in tok and tok.replace(",", "") in body.replace(",", ""):
                return True
        return False

    failures = []
    scope = set()
    checked = 0
    for lineno, line in enumerate(text.splitlines(), 1):
        decl = SOURCES_DECL_RE.search(line)
        if decl:
            scope = {p for p in re.split(r"[,\s]+", decl.group(1)) if p} & set(SOURCES)
            continue
        if "<!--nocheck-->" in line:
            continue
        allowed = set(scope) | ALWAYS_OK
        stripped = line.strip()
        if stripped.startswith("|"):
            # Table row: the last cell may name this row's findings source(s).
            cells = [c for c in stripped.replace("\\|", "\x00").split("|") if c.strip()]
            if cells:
                row_keys = parse_keys(cells[-1])
                if row_keys:
                    allowed = row_keys | ALWAYS_OK
        for tok in NUM_RE.findall(line):
            checked += 1
            if not token_in(tok, allowed):
                failures.append((lineno, tok, sorted(allowed)))
    for lineno, tok, allowed in failures[:20]:
        print(f"    line {lineno}: number {tok!r} not found in sources {allowed}")
    gate("G4 numbers-traceable", not failures,
         f"{len(failures)} untraceable numbers" if failures
         else f"{checked} numeric tokens all trace to their cited sources")

    # ---- G5: ledger mandatory items --------------------------------------
    ledger = section_text("## 3.")
    absent = [s for s in LEDGER_MUSTS if s.lower() not in ledger.lower()]
    gate("G5 ledger-mandatory-items", not absent,
         f"missing: {absent}" if absent else f"{len(LEDGER_MUSTS)} mandatory items present")

    n_fail = sum(1 for _, ok, _ in results if not ok)
    dt = time.time() - t0
    print(f"{'PASS' if n_fail == 0 else 'FAIL'} "
          f"({len(results) - n_fail}/{len(results)} gates, {dt:.1f}s)")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

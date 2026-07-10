#!/usr/bin/env python3
r"""
Document-invariant probe for BUG-001 (empirical-ratio-unmeasured).

GOVERNANCE/DOCUMENTATION guard, NOT a measurement. It parses the committed
LaTeX sources and asserts two structural invariants that BUG-001-A (uniform
conditional scoping) and BUG-001-B (the section 11 <-> companion-hypothesis
cross-link) establish, so that a future edit that silently removes a
qualifier or breaks the mapping turns the CI job red.

  inv1 (conditional scoping): in the Abstract, section 1.1, and the section 11
       region of the checked document, every occurrence of the absolute-ratio /
       beats-seekable-format claim co-occurs -- within a bounded character
       neighborhood -- with at least one conditional/falsifiable marker from
       {conditional, pending, conjecture, "Section 11"/"section 11"/"S11",
        "not a theorem", unmeasured, "empirical", "in principle",
        "pre-registered", "predicted"}. The paper must never assert the ratio
       claim unconditionally. (Section 11.7 is EXCLUDED: it is theorem-backed
       -- Theorem 7.15 -- and verified against a cited measurement, so it is
       deliberately phrased as a derived prediction, not a conjecture.)

  inv2 (no-orphan mapping): every H-number that the section 11 predictions
       (11.1-11.6) cite as a pre-registered hypothesis EXISTS as a defined
       hypothesis ("\\textbf{H<k> {[}...{]}:}") in
       tex/rnr_experimental_design.tex; and conversely every ratio/predictor-
       relevant companion hypothesis (the H1-H7, H15, H16 set named in the
       section 1.1 cross-link) is reachable from the main paper (cited in
       section 11.1-11.7 or in the section 1.1 cross-link enumeration). No
       dangling H-reference on either side.

DOCUMENTS CHECKED (post three-paper split; the invariant must hold in every
document it appears in):

  * tex/rnr_coding.tex (monolith, still the source of truth): inv1 + inv2 +
    self-test, unchanged semantics.
  * tex/papers/core/rnr_core.tex (Part I): carries the Abstract, section 1.1,
    and section 11 (11.1-11.7 all stayed in Part I), so the full inv1 + inv2 +
    self-test battery runs on it identically (the companion
    tex/rnr_experimental_design.tex is shared, not split).
  * tex/papers/random_access/rnr_random_access.tex (Part II): the monolith
    Abstract's sync-overhead / ratio-comparison passage ("The ratio comparison
    against ... bgzip and zstd seekable ... is empirical") moved into Part II's
    section 10.7 comparison block (SPLIT-PLAN section 5.4: "Abstract lines
    230-252 -> moves to P2"). The probe follows the moved TEXT: it checks
    Part II's own Abstract plus the enclosing paragraph of the moved
    ratio-comparison claim for inv1 conditional scoping. inv2 does not apply
    (the section 1.1 <-> section 11 <-> companion mapping lives in Part I).

A negative self-test (--selftest, also run automatically) MUTATES the parsed
text -- strips a qualifier; injects a non-existent H-number -- and asserts the
checks then FAIL, proving the guard is not vacuous.

Stdlib only. Prints PASS/FAIL per check and "OVERALL -> PASS" on success.
Run locally:  /Users/para/.venvs/rnr/bin/python scripts/verify/bug_001_ratio_claim_conditional.py
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
MAIN = os.path.join(REPO, "tex", "rnr_coding.tex")
CORE = os.path.join(REPO, "tex", "papers", "core", "rnr_core.tex")
RA = os.path.join(REPO, "tex", "papers", "random_access", "rnr_random_access.tex")
COMPANION = os.path.join(REPO, "tex", "rnr_experimental_design.tex")

# Neighborhood (characters) within which a ratio-claim keyword must see a
# conditional marker. Chosen wide enough to span a sentence/paragraph but
# narrow enough that a genuinely unconditional claim elsewhere would fail.
NEIGHBORHOOD = 600

# Keywords that signal the absolute-ratio / beats-seekable value-proposition
# claim. Lower-cased substring match. Two tiers:
#
#  STRONG: inherently a ratio / RNR-superiority claim wherever it appears, so
#  it must always carry a nearby conditional marker.
RATIO_CLAIM_KEYWORDS_STRONG = [
    "compression ratio",
    "ratio comparison",
    "ratio advantage",
    "outperform",
    "bits per byte",
    "bit per byte",
    "bit-per-byte",
    "below lzma",
    "competitive with",
]
# AMBIGUOUS comparatives ("dominate" also means "X is dominated BY redundancy",
# a descriptive use that is NOT the ratio claim): count only when a FORMAT /
# baseline object co-occurs in the same sentence (RNR dominates <a seekable
# format>), excluding "redundancy is dominated by local repetition".
RATIO_CLAIM_KEYWORDS_AMBIGUOUS = [
    "dominate",
    "dominates",
]
#
#  FORMAT: a baseline-format NAME (seekable/bgzip/zstd/dictionary-based). A
#  bare descriptive mention of how existing formats work is NOT the guarded
#  claim; it only becomes the value-proposition claim when a comparative-
#  superiority term co-occurs in the same neighborhood (RNR beats/loses-no
#  advantage-over/recasts vs. those formats). Only then is a conditional
#  marker required.
RATIO_CLAIM_KEYWORDS_FORMAT = [
    "seekable",
    "bgzip",
    "zstd",
    "dictionary-based",
    "dictionary based",
]
# Matched with word boundaries (a comparative term must be a whole word, so
# "ratio" does not spuriously match inside "operationally", etc.).
COMPARATIVE_TERMS = [
    "beat",
    "beats",
    "outperform",
    "outperforms",
    "advantage",
    "better",
    "below",
    "competitive",
    "dominate",
    "dominates",
    "recast",
    "recasting",
    "without losing",
    "ratio",
]
_COMPARATIVE_RE = re.compile(
    r"|".join(r"\b" + re.escape(t) + r"\b" for t in COMPARATIVE_TERMS))

# Conditional / falsifiable markers that legitimize a ratio-claim occurrence.
CONDITIONAL_MARKERS = [
    "conditional",
    "pending",
    "conjecture",
    "section 11",
    "§11",          # §11
    "§§11",    # §§11 (rare)
    "s11",
    "not a theorem",
    "unmeasured",
    "empirical",
    "in principle",
    "pre-registered",
    "pre-register",
    "predicted",
    "prediction",
    "falsifiable",
    "not validated",
    "not been empirically validated",
    "no quantitative claim",
]


def read(path):
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


def line(ok, label, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {detail}")
    return ok


# --------------------------------------------------------------------------
# Region extraction
# --------------------------------------------------------------------------
def find_marker(text, marker, start=0):
    idx = text.find(marker, start)
    if idx < 0:
        raise RuntimeError(f"region marker not found: {marker!r}")
    return idx


def extract_regions(text):
    """Return the in-scope regions of the main paper as {name: substring}.

    In scope for inv1: Abstract, section 1.1, and section 11 intro + 11.1-11.6.
    section 11.7 is explicitly EXCLUDED (theorem-backed, BUG-002 territory).
    """
    abs_start = find_marker(text, "\\textbf{Abstract}")
    abs_end = find_marker(text, "\\textbf{Keywords:}", abs_start)

    s11_start = find_marker(text, "\\textbf{1. Introduction}")
    # section 1.1 spans from "1.1 Motivation" up to "1.2 Position"
    s11_mot = find_marker(text, "\\textbf{1.1 Motivation}", s11_start)
    s11_end = find_marker(text, "\\textbf{1.2 Position", s11_mot)

    sec11_start = find_marker(text, "\\textbf{11. Experimental predictions}")
    # 11.1-11.6 stop where 11.7 begins (11.7 is out of scope).
    sec11_7 = find_marker(text, "\\textbf{11.7 ", sec11_start)

    return {
        "Abstract": text[abs_start:abs_end],
        "section 1.1": text[s11_mot:s11_end],
        "section 11 (intro + 11.1-11.6)": text[sec11_start:sec11_7],
    }


def extract_regions_ra(text):
    """In-scope regions of Part II (random access) as {name: substring}.

    The monolith Abstract's sync-overhead / ratio-comparison passage moved
    into Part II's section 10.7 comparison block (SPLIT-PLAN section 5.4), so
    the invariant follows the moved TEXT: (a) Part II's own Abstract must
    never assert the ratio claim unconditionally, and (b) the enclosing
    paragraph of the moved ratio-comparison claim must keep its conditional
    markers.
    """
    abs_start = find_marker(text, "\\textbf{Abstract}")
    abs_end = find_marker(text, "\\textbf{Keywords:}", abs_start)

    s107 = find_marker(text, "\\textbf{10.7 Random access}")
    anchor = find_marker(text, "The ratio comparison", s107)
    para_lo = text.rfind("\n\n", 0, anchor)
    para_lo = 0 if para_lo < 0 else para_lo + 2
    para_hi = text.find("\n\n", anchor)
    para_hi = len(text) if para_hi < 0 else para_hi

    return {
        "Abstract (Part II)": text[abs_start:abs_end],
        "section 10.7 moved sync-overhead/ratio-comparison paragraph":
            text[para_lo:para_hi],
    }


# --------------------------------------------------------------------------
# inv1: every ratio-claim occurrence carries a nearby conditional marker
# --------------------------------------------------------------------------
def _all_indices(low, kw):
    out, start = [], 0
    while True:
        idx = low.find(kw, start)
        if idx < 0:
            break
        out.append(idx)
        start = idx + 1
    return out


def claim_occurrences(region_text):
    """Occurrences of the guarded value-proposition claim.

    STRONG keywords always count. FORMAT (baseline-name) keywords count only
    when a comparative-superiority term co-occurs within the neighborhood --
    so a purely descriptive mention of how bgzip/zstd-seekable work is not
    flagged, but an RNR-beats-seekable comparison is.
    """
    low = region_text.lower()
    hits = []
    for kw in RATIO_CLAIM_KEYWORDS_STRONG:
        for idx in _all_indices(low, kw):
            hits.append((idx, kw))
    for kw in RATIO_CLAIM_KEYWORDS_FORMAT:
        for idx in _all_indices(low, kw):
            # Comparative term must appear in the SAME SENTENCE as the format
            # name -- otherwise a purely descriptive mention of how existing
            # seekable formats work (no RNR comparison) is not the guarded
            # value-proposition claim and is not flagged.
            sent = _sentence_around(low, idx)
            if _COMPARATIVE_RE.search(sent):
                hits.append((idx, kw))
    for kw in RATIO_CLAIM_KEYWORDS_AMBIGUOUS:
        for idx in _all_indices(low, kw):
            # Counts only when a baseline/format object is in the SAME sentence
            # ("dominate <seekable format>"), so descriptive "dominated by local
            # repetition" (no baseline object) is not flagged.
            sent = _sentence_around(low, idx)
            if any(fmt in sent for fmt in RATIO_CLAIM_KEYWORDS_FORMAT):
                hits.append((idx, kw))
    return hits


# Sentence boundary: ". " or "; " or paragraph break. LaTeX uses literal
# periods; we split on a period/semicolon followed by whitespace or end.
_SENT_SPLIT = re.compile(r"(?<=[.;])\s")


def _sentence_around(low, idx):
    """The sentence (period/semicolon-delimited) containing position idx."""
    # left boundary
    lo = 0
    for m in _SENT_SPLIT.finditer(low, 0, idx):
        lo = m.end()
    # right boundary
    m = _SENT_SPLIT.search(low, idx)
    hi = m.start() if m else len(low)
    return low[lo:hi]


def has_nearby_marker(region_text, idx):
    low = region_text.lower()
    # Locality: clamp the marker search to the ENCLOSING PARAGRAPH (blank-line
    # bounded), so a stripped qualifier in one subsection cannot be masked by an
    # adjacent subsection's markers -- a paragraph break never crosses a
    # \textbf{} run-in subsection header. Then intersect with the +-NEIGHBORHOOD
    # cap for long paragraphs.
    para_lo = low.rfind("\n\n", 0, idx)
    para_lo = 0 if para_lo < 0 else para_lo + 2
    para_hi = low.find("\n\n", idx)
    para_hi = len(low) if para_hi < 0 else para_hi
    lo = max(para_lo, idx - NEIGHBORHOOD)
    hi = min(para_hi, idx + NEIGHBORHOOD)
    window = low[lo:hi]
    return any(m in window for m in CONDITIONAL_MARKERS)


def check_inv1(regions, verbose=True):
    if verbose:
        print("\ninv1: conditional scoping of the absolute-ratio / beats-seekable claim")
    ok = True
    total = 0
    for name, region in regions.items():
        unconditional = []
        for idx, kw in claim_occurrences(region):
            total += 1
            if not has_nearby_marker(region, idx):
                snippet = re.sub(r"\s+", " ", region[max(0, idx - 40): idx + 60]).strip()
                unconditional.append((kw, snippet))
        region_ok = not unconditional
        if verbose:
            line(region_ok, name,
                 "all claim occurrences carry a nearby conditional marker"
                 if region_ok else f"{len(unconditional)} UNCONDITIONAL: {unconditional[:3]}")
        ok &= region_ok
    if verbose and ok:
        print(f"    ({total} ratio-claim keyword occurrences checked, all scoped)")
    return ok


# --------------------------------------------------------------------------
# inv2: no orphan in the section 11 <-> companion-hypothesis mapping
# --------------------------------------------------------------------------
# H-numbers the section 1.1 cross-link enumerates as ratio/predictor-relevant.
RATIO_RELEVANT_HS = {"H1", "H2", "H3", "H4", "H5", "H6", "H7", "H15", "H16"}


def defined_hypotheses(companion_text):
    """H-numbers DEFINED in the companion: \\textbf{H<k> {[}...{]}:}."""
    return set(re.findall(r"\\textbf\{(H\d+)\s*\{\[\}", companion_text))


def cited_hypotheses_in_main(text):
    """H-numbers the main paper CITES as pre-registered hypotheses.

    Restrict to the section 11 region + the section 1.1 cross-link, where the
    mapping is asserted. Pattern: 'Hypothes(is|es) H<k>' or 'H<k>--H<k>' ranges
    or bare 'H<k>' immediately after 'Hypothes...'.
    """
    sec1 = find_marker(text, "\\textbf{1.1 Motivation}")
    sec11 = find_marker(text, "\\textbf{11. Experimental predictions}")
    sec12 = find_marker(text, "\\textbf{12. Related work}")
    scope = text[sec1: find_marker(text, "\\textbf{1.2 Position", sec1)] + text[sec11:sec12]

    cited = set()
    # ranges like H10--H11, H1--H7, H15--H16
    for a, b in re.findall(r"H(\d+)\s*--\s*H(\d+)", scope):
        for k in range(int(a), int(b) + 1):
            cited.add(f"H{k}")
    # singletons explicitly tagged as Hypothesis H<k>
    for k in re.findall(r"Hypothes[ie]s\s+H(\d+)", scope):
        cited.add(f"H{k}")
    # also "Hypotheses H3 and H4" style
    for k in re.findall(r"\bH(\d+)\b", scope):
        # only count if a "hypothes" word appears within 80 chars before it
        pass
    return cited, scope


def check_inv2(main_text, companion_text, verbose=True, ratio_relevant=RATIO_RELEVANT_HS):
    if verbose:
        print("\ninv2: no orphan in the section 11 <-> companion-hypothesis mapping")
    defined = defined_hypotheses(companion_text)
    cited, scope = cited_hypotheses_in_main(main_text)

    ok = True
    # (a) every cited H exists in the companion.
    dangling = sorted(cited - defined, key=lambda h: int(h[1:]))
    ok &= line(not dangling, "every cited H exists in companion",
               "no dangling H-reference" if not dangling
               else f"DANGLING (cited but undefined): {dangling}")

    # (b) every ratio-relevant H is reachable (cited) from the main paper.
    unreached = sorted(ratio_relevant - cited, key=lambda h: int(h[1:]))
    ok &= line(not unreached, "every ratio-relevant H is cited from the main paper",
               f"all of {sorted(ratio_relevant, key=lambda h: int(h[1:]))} reachable"
               if not unreached else f"ORPHAN companion H (defined+relevant, never cited): {unreached}")

    # (c) the seekable-format baseline coverage is explicitly named.
    low = scope.lower()
    seekable_named = ("bgzip" in low and "zstd" in low and "seekable" in low)
    ok &= line(seekable_named, "seekable-format baseline named in the cross-link",
               "bgzip / zstd-seekable named" if seekable_named
               else "seekable-format baseline NOT named in section 1.1/11 cross-link")

    if verbose:
        print(f"    (companion defines {len(defined)} hypotheses; "
              f"main paper cites {sorted(cited, key=lambda h: int(h[1:]))})")
    return ok


# --------------------------------------------------------------------------
# Non-vacuity self-test: mutate the inputs, assert the checks FAIL.
# --------------------------------------------------------------------------
def _strip_markers(regions):
    """Regions with every conditional marker replaced by 'XXX' (mutation)."""
    mutated = {}
    for name, region in regions.items():
        m = region
        for marker in CONDITIONAL_MARKERS:
            m = re.sub(re.escape(marker), "XXX", m, flags=re.IGNORECASE)
        mutated[name] = m
    return mutated


def selftest_marker_strip(regions):
    """inv1 must FAIL on marker-stripped regions (non-vacuity of inv1)."""
    inv1_on_mutated = check_inv1(_strip_markers(regions), verbose=False)
    return line(not inv1_on_mutated, "inv1 detects a stripped qualifier",
                "mutated (markers removed) text correctly FAILS inv1"
                if not inv1_on_mutated else "MUTATION NOT DETECTED -- inv1 is vacuous!")


def selftest(main_text, companion_text):
    print("\nself-test (non-vacuity): mutated inputs MUST fail the checks")
    ok = True

    # (1) inv1 must fail if we strip conditional markers from a region.
    ok &= selftest_marker_strip(extract_regions(main_text))

    # (2) inv2 must fail if the main paper cites a non-existent hypothesis.
    injected = main_text.replace(
        "pre-registered as Hypothesis H1 (quantifying)",
        "pre-registered as Hypothesis H99 (quantifying)", 1)
    inv2_on_injected = check_inv2(injected, companion_text, verbose=False)
    ok &= line(not inv2_on_injected, "inv2 detects a dangling H-reference",
               "injected H99 (undefined) correctly FAILS inv2"
               if not inv2_on_injected else "INJECTION NOT DETECTED -- inv2 is vacuous!")

    # (3) inv2's orphan side must fire if a hypothesis that SHOULD be cross-
    #     referenced is defined in the companion but never cited from the main
    #     paper. Exercise it by extending the required set with H17 (a real,
    #     defined-but-uncited companion hypothesis) and asserting inv2 fails.
    extended = set(RATIO_RELEVANT_HS) | {"H17"}
    assert "H17" in defined_hypotheses(companion_text), "self-test premise: H17 is defined"
    inv2_orphan = check_inv2(main_text, companion_text, verbose=False, ratio_relevant=extended)
    ok &= line(not inv2_orphan, "inv2 detects an uncited (orphan) required H",
               "a defined-but-uncited required H (H17) correctly FAILS inv2"
               if not inv2_orphan else "ORPHAN NOT DETECTED -- inv2 orphan-side is vacuous!")

    return ok


def check_full_document(label, text, companion_text):
    """Full battery (inv1 + inv2 + self-test) for a document that carries the
    Abstract / section 1.1 / section 11 regions (monolith and Part I)."""
    print("\n" + "-" * 70)
    print(f"Document: {label}")
    regions = extract_regions(text)
    return {
        f"[{label}] inv1 conditional scoping (Abstract/1.1/11)": check_inv1(regions),
        f"[{label}] inv2 no-orphan mapping (11 <-> companion)":
            check_inv2(text, companion_text),
        f"[{label}] self-test non-vacuity": selftest(text, companion_text),
    }


def check_ra_document(label, text):
    """Part II battery: inv1 on the regions the split moved there (its own
    Abstract + the migrated section 10.7 ratio-comparison paragraph), plus the
    marker-strip non-vacuity self-test on those regions. inv2 does not apply:
    the section 1.1/11 <-> companion mapping lives in Part I."""
    print("\n" + "-" * 70)
    print(f"Document: {label}")
    regions = extract_regions_ra(text)
    inv1_ok = check_inv1(regions)
    print("\nself-test (non-vacuity): mutated inputs MUST fail the checks")
    strip_ok = selftest_marker_strip(regions)
    return {
        f"[{label}] inv1 conditional scoping (Abstract + moved 10.7 passage)": inv1_ok,
        f"[{label}] self-test non-vacuity": strip_ok,
    }


def main():
    print("Verification: BUG-001 conditional-ratio-claim + mapping invariants")
    print("=" * 70)
    print("DOCUMENT-INVARIANT GUARD -- asserts the paper never makes the")
    print("absolute-ratio/beats-seekable claim unconditionally, and that the")
    print("section 11 <-> companion-hypothesis mapping has no orphan.")
    print("Checked in every document the invariant appears in after the split:")
    print("monolith (source of truth), Part I (core), and Part II (the moved")
    print("Abstract sync-overhead/ratio-comparison passage, now in its 10.7).")
    print("Does NOT measure any compression ratio (that is BUG-001-D, external).")

    companion_text = read(COMPANION)

    results = {}
    results.update(check_full_document(
        "monolith tex/rnr_coding.tex", read(MAIN), companion_text))
    results.update(check_full_document(
        "Part I tex/papers/core/rnr_core.tex", read(CORE), companion_text))
    results.update(check_ra_document(
        "Part II tex/papers/random_access/rnr_random_access.tex", read(RA)))

    print("\n" + "=" * 70)
    print("Summary:")
    all_ok = True
    for name, ok in results.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        all_ok &= ok

    print()
    if all_ok:
        print("OVERALL -> PASS")
        return 0
    print("OVERALL -> FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())

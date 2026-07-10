#!/usr/bin/env python3
r"""
Document-invariant probe for BUG-010 (conditional-scope-no-universal-dominance).

GOVERNANCE/DOCUMENTATION guard, NOT a measurement. It parses the committed
LaTeX sources and asserts the document invariants that BUG-010-A (the
Abstract/section 1.1 conditional cross-link) and BUG-010-B (the section 10.2.1
precondition inequality, accounting form) establish, so that a future edit that
silently strengthens the claim, drops a cross-reference, or removes a cost term
turns the CI job red.

  inv1 (conditional scoping with cross-references): in the Abstract and section
       1.1 of the checked document, the central no-universal-dominance / data-class
       precondition claim is CONDITIONAL and carries an explicit forward pointer
       to all three formal anchors -- Definition 2.1 (the conditionally-
       advantageous predicate), Theorem 8.1 (no universal dominance), and section
       10.2 (the model-amortization regime that makes L(M)/m payable). Each of the
       two regions must also contain a conditional marker (conditional / in
       principle / open empirical question / pending).

  inv2 (NO AFFIRMATIVE universal-dominance sentence): across the whole paper, no
       sentence AFFIRMATIVELY asserts that RNR (or "the framework") beats / shortens /
       dominates the baseline universally / always / on every input. This is the
       FALSE-POSITIVE-HAZARD invariant: the paper DELIBERATELY disclaims universal
       dominance, so negation/disclaimer phrasings ("no lossless code shortens
       every input", "we do not claim universal dominance", "not universal", "RNR
       does not assert that any fixed coder dominates any baseline on every input")
       WILL appear and are CORRECT. The probe flags a sentence ONLY when it has, in
       the SAME sentence: (1) an affirmative RNR subject, (2) a universal quantifier
       applied to inputs/baselines, and (3) NO negation/disclaimer cue. An explicit
       allowlist of the known disclaimer sentences provides a second layer of
       safety.

  inv3 (precondition inequality term-completeness + section-anchoring): the
       section 10.2.1 precondition passage exists, names every component cost term
       of inequality (10.1) -- residual uncertainty H(X|Y), amortized model L(M)/m,
       repair residual, verification overhead, metadata, and the baseline rate --
       and each term carries a resolvable section-anchor (a section number cited in
       the same item).

  inv4 (named cross-references resolve): every named-object cross-reference the
       claim relies on -- Definition 2.1, Theorem 8.1, and the section 10.2.1
       precondition block -- actually resolves to a header that exists in the
       document. Catches a FABRICATED pointer that inv1's presence-only check
       cannot (a "Definition 2.1" citation with no such labelled header anywhere).

A negative self-test (run automatically) exercises the guard on tiny inline
fixtures: a FABRICATED affirmative-universal sentence the checker MUST REJECT,
and a real-style DISCLAIMER sentence it MUST ACCEPT; plus mutations of the parsed
text (strip a cross-reference; delete a cost term) that the checks MUST then fail.

DOCUMENTS CHECKED (post three-paper split; the invariant must hold in every
document it appears in):

  * tex/rnr_coding.tex (monolith, still the source of truth): inv1-inv4 +
    self-test, unchanged semantics.
  * tex/papers/core/rnr_core.tex (Part I): carries the Abstract, section 1.1,
    section 10.2.1, Definition 2.1, and Theorem 8.1 (per SPLIT-PLAN section 2
    all of sections 1, 8, 10.1-10.5 stayed in Part I), so the full inv1-inv4 +
    self-test battery runs on it identically.
  * tex/papers/random_access/rnr_random_access.tex (Part II) and
    tex/papers/dispersion/rnr_dispersion.tex (Part III): inv1/inv3/inv4's
    anchors (Definition 2.1, Theorem 8.1, section 10.2.1) live in Part I, but
    inv2 is a WHOLE-PAPER scan and large blocks of the monolith's body
    (sections 7.2-7.36 and 10.7-10.8 -> Part II; the section 7.34 family ->
    Part III) moved into these documents, so the no-affirmative-universal-
    dominance scan follows the moved text: inv2 + its injection self-test run
    on each companion part.

Stdlib only. Prints PASS/FAIL per check and "OVERALL -> PASS" on success.
Run locally:  /Users/para/.venvs/rnr/bin/python scripts/verify/bug_010_precondition_invariant.py
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
MAIN = os.path.join(REPO, "tex", "rnr_coding.tex")
CORE = os.path.join(REPO, "tex", "papers", "core", "rnr_core.tex")
RA = os.path.join(REPO, "tex", "papers", "random_access", "rnr_random_access.tex")
DISP = os.path.join(REPO, "tex", "papers", "dispersion", "rnr_dispersion.tex")


def read(path):
    with open(path, "r", encoding="utf-8") as fh:
        return fh.read()


def line(ok, label, detail, verbose=True):
    if verbose:
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {detail}")
    return ok


def find_marker(text, marker, start=0):
    idx = text.find(marker, start)
    if idx < 0:
        raise RuntimeError(f"region marker not found: {marker!r}")
    return idx


# --------------------------------------------------------------------------
# Region extraction
# --------------------------------------------------------------------------
def extract_regions(text):
    """Return {name: substring} for the Abstract, section 1.1, and section 10.2.1."""
    abs_start = find_marker(text, "\\textbf{Abstract}")
    abs_end = find_marker(text, "\\textbf{Keywords:}", abs_start)

    s11_mot = find_marker(text, "\\textbf{1.1 Motivation}")
    s11_end = find_marker(text, "\\textbf{1.2 Position", s11_mot)

    s1021_start = find_marker(text, "\\textbf{10.2.1 ")
    s1021_end = find_marker(text, "\\textbf{10.3 ", s1021_start)

    return {
        "Abstract": text[abs_start:abs_end],
        "section 1.1": text[s11_mot:s11_end],
        "section 10.2.1": text[s1021_start:s1021_end],
    }


# --------------------------------------------------------------------------
# inv1: conditional scoping with the three cross-references
# --------------------------------------------------------------------------
# Each anchor: a list of acceptable surface forms (lower-cased substring).
ANCHOR_FORMS = {
    "Definition 2.1": ["definition 2.1", "def. 2.1", "def 2.1"],
    "Theorem 8.1": ["theorem 8.1", "thm 8.1", "thm. 8.1"],
    "section 10.2": ["10.2"],   # matches §10.2 and §10.2.1 alike
}
CONDITIONAL_MARKERS = [
    "conditional", "in principle", "open empirical", "pending",
    "not universal", "no claim of universal", "do not claim universal",
]


def check_inv1(regions, verbose=True):
    if verbose:
        print("\ninv1: Abstract & section 1.1 conditional claim carries Def 2.1 + Thm 8.1 + section 10.2")
    ok = True
    for name in ("Abstract", "section 1.1"):
        low = regions[name].lower()
        missing = [anchor for anchor, forms in ANCHOR_FORMS.items()
                   if not any(f in low for f in forms)]
        has_marker = any(m in low for m in CONDITIONAL_MARKERS)
        region_ok = (not missing) and has_marker
        detail = ("all three anchors present + conditional marker"
                  if region_ok else
                  f"missing anchors={missing} ; conditional_marker={has_marker}")
        ok &= line(region_ok, name, detail, verbose=verbose)
    return ok


# --------------------------------------------------------------------------
# inv2: NO AFFIRMATIVE universal-dominance sentence (false-positive-safe)
# --------------------------------------------------------------------------
# Sentence splitter: period/semicolon/colon followed by whitespace.
_SENT_SPLIT = re.compile(r"(?<=[.;:])\s")

# Affirmative RNR subject cues (the subject of a dominance claim).
RNR_SUBJECT = ["rnr", "the framework", "our coder", "our method", "our encoder",
               "the encoder \\mathsf{enc}", "\\mathsf{enc}"]
# Comparative-superiority verbs (the claim verb).
SUPERIORITY_VERB = ["dominate", "dominates", "beat", "beats", "outperform",
                    "outperforms", "shorten", "shortens", "compress", "compresses",
                    "wins", "win"]
# Universal quantifiers applied to inputs / baselines.
UNIVERSAL_QUANT = ["every input", "every baseline", "all inputs", "all baselines",
                   "any input", "any baseline", "every source", "all sources",
                   "universally", "always", "in all cases", "every string",
                   "every codeword", "universal dominance", "universal compression"]
# Negation / disclaimer cues that make a universal sentence a (correct) disclaimer.
NEGATION_CUES = ["no ", "not ", "never", "cannot", "can't", "do not", "does not",
                 "doesn't", "don't", "without", "disclaim", "rules out", "impossible",
                 "no claim", "fails to", "rare", "exponent", "pigeonhole",
                 "conditional", "n't"]

# Known, hand-verified disclaimer sentences (allowlist). A sentence whose
# normalized text contains any of these fragments is ALWAYS accepted as a
# (correct) disclaimer regardless of the heuristic. Keeps the heuristic from
# ever flagging the paper's deliberate negations.
DISCLAIMER_ALLOWLIST = [
    "we do not claim universal dominance",
    "no lossless code can map every input",
    "rnr does not assert that any fixed coder",
    "dominates any baseline on every input",
    "the claim is conditional and distributional",
    "no injective map",
    "assigns every input string of length",
    "not only does no code shorten",
    "there is no claim of universal dominance",
]


def _normalize(s):
    return re.sub(r"\s+", " ", s).strip().lower()


def is_affirmative_universal(sentence):
    """True iff `sentence` AFFIRMATIVELY claims RNR beats the baseline universally.

    False-positive-safe: requires (subject) AND (superiority verb) AND
    (universal quantifier) AND (NO negation/disclaimer cue) in the same sentence,
    and is overridden by the disclaimer allowlist.
    """
    low = _normalize(sentence)
    # Allowlisted disclaimers are never affirmative-universal.
    if any(frag in low for frag in DISCLAIMER_ALLOWLIST):
        return False
    # Any negation/disclaimer cue in the sentence => it is a disclaimer, not an
    # affirmative claim.
    if any(cue in low for cue in NEGATION_CUES):
        return False
    has_subject = any(s in low for s in RNR_SUBJECT)
    has_verb = any(v in low for v in SUPERIORITY_VERB)
    has_universal = any(q in low for q in UNIVERSAL_QUANT)
    return has_subject and has_verb and has_universal


def sentences(text):
    # Strip LaTeX comment lines to avoid commented-out fixtures tripping the check.
    body = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("%"))
    return _SENT_SPLIT.split(body)


def check_inv2(text, verbose=True):
    if verbose:
        print("\ninv2: no AFFIRMATIVE universal-dominance sentence (disclaimers are OK)")
    offenders = [_normalize(s)[:120] for s in sentences(text) if is_affirmative_universal(s)]
    ok = not offenders
    line(ok, "whole paper",
         "no affirmative universal-dominance sentence found"
         if ok else f"{len(offenders)} AFFIRMATIVE-UNIVERSAL: {offenders[:3]}",
         verbose=verbose)
    return ok


# --------------------------------------------------------------------------
# inv3: precondition inequality term-completeness + section-anchoring
# --------------------------------------------------------------------------
# Each cost term: (label, [surface forms], [acceptable section anchors]).
# The section anchor must appear in the section 10.2.1 region (the inequality
# item that defines the term cites where it is bounded).
COST_TERMS = [
    ("residual uncertainty H(X|Y)", ["h(x\\mid y)", "h(x \\mid y)", "residual uncertainty"], ["2.2"]),
    ("amortized model L(M)/m",      ["l(m)/m", "\\tfrac{l(m)}{m}"],                          ["10.2"]),
    ("repair residual",             ["repair residual", "\\rho(d)"],                         ["2.3", "5"]),
    ("verification overhead",       ["verification overhead"],                               ["10.3"]),
    ("metadata",                    ["metadata"],                                            ["10.4", "10.7"]),
    ("baseline rate",               ["baseline rate", "r_b(d)"],                             ["2.1", "def"]),
]


def check_inv3(regions, verbose=True):
    if verbose:
        print("\ninv3: section 10.2.1 precondition inequality -- all cost terms present + section-anchored")
    region = regions["section 10.2.1"]
    low = region.lower()
    ok = True
    # The labelled inequality itself must exist.
    has_eq = "\\tag{10.1}" in region or "(10.1)" in region
    ok &= line(has_eq, "labelled inequality (10.1)",
               "found \\tag{10.1}" if has_eq else "MISSING: no (10.1) labelled inequality in section 10.2.1",
               verbose=verbose)
    for label, forms, anchors in COST_TERMS:
        present = any(f in low for f in forms)
        anchored = any(a in low for a in anchors)
        term_ok = present and anchored
        ok &= line(term_ok, f"term: {label}",
                   "present + section-anchored"
                   if term_ok else f"present={present} anchored(any of {anchors})={anchored}",
                   verbose=verbose)
    return ok


# --------------------------------------------------------------------------
# inv4: cross-referenced named anchors RESOLVE to a header in the document
# --------------------------------------------------------------------------
# Guards against a FABRICATED cross-reference -- a "Definition 2.1" / "Theorem
# 8.1" pointer with no matching header anywhere. A presence-only check (inv1)
# cannot catch this: it confirms the Abstract/section 1.1 *cite* the anchor, not
# that the anchor *exists*. (This is the exact defect a cold re-read caught: the
# conditionally-advantageous predicate of section 2.1 was cross-referenced as
# "Definition 2.1" before any such labelled header existed.)
RESOLVABLE_ANCHORS = [
    ("Definition 2.1", ["\\textbf{Definition 2.1"]),
    ("Theorem 8.1",    ["\\textbf{Theorem 8.1"]),
    ("section 10.2.1", ["\\textbf{10.2.1 "]),
]


def check_inv4(text, verbose=True):
    if verbose:
        print("\ninv4: every named cross-reference resolves to a header in the document")
    ok = True
    for name, headers in RESOLVABLE_ANCHORS:
        resolved = any(h in text for h in headers)
        ok &= line(resolved, f"anchor: {name}",
                   "header present"
                   if resolved else f"DANGLING: cited but no header {headers} exists",
                   verbose=verbose)
    return ok


# --------------------------------------------------------------------------
# Non-vacuity self-test: inline fixtures + mutations
# --------------------------------------------------------------------------
def selftest_inv2_injection(text):
    """inv2 must FAIL when an affirmative-universal sentence is injected."""
    injected = text + "\n\nRNR outperforms every baseline universally and always wins.\n"
    inv2_inj = check_inv2(injected, verbose=False)
    return line(not inv2_inj, "inv2 detects an injected affirmative-universal sentence",
                "injected affirmative-universal correctly FAILS inv2"
                if not inv2_inj else "INJECTION NOT DETECTED -- inv2 is vacuous!")


def selftest(text):
    print("\nself-test (non-vacuity): inline fixtures + mutated inputs")
    ok = True

    # (A) FABRICATED affirmative-universal sentence the checker MUST REJECT.
    bad = "RNR outperforms every baseline on every input and always wins."
    ok &= line(is_affirmative_universal(bad),
               "inline negative control: fabricated affirmative-universal",
               "correctly REJECTED (flagged as affirmative-universal)"
               if is_affirmative_universal(bad) else "NOT FLAGGED -- inv2 is vacuous!")

    # (B) DISCLAIMER sentence the checker MUST ACCEPT (not flag).
    good = ("We do not claim universal dominance --- no lossless code can map "
            "every input of length n to a strictly shorter one.")
    ok &= line(not is_affirmative_universal(good),
               "inline positive control: real disclaimer sentence",
               "correctly ACCEPTED (not flagged)"
               if not is_affirmative_universal(good) else "FALSE POSITIVE -- disclaimer flagged!")

    # (B') a SECOND disclaimer phrased as an affirmative-looking universal but negated.
    good2 = "RNR does not assert that any fixed coder dominates any baseline on every input."
    ok &= line(not is_affirmative_universal(good2),
               "inline positive control: negated-universal disclaimer",
               "correctly ACCEPTED (not flagged)"
               if not is_affirmative_universal(good2) else "FALSE POSITIVE -- negated disclaimer flagged!")

    # (C) inv1 must fail if a cross-reference is stripped from section 1.1.
    mutated = text.replace("the conditionally-advantageous\npredicate of Definition 2.1", "XXX", 1)
    mutated = mutated.replace("Definition 2.1\n($\\mathbb{E}", "XXX\n($\\mathbb{E}", 1)
    # Brute-force: blank every "Definition 2.1" form in the section 1.1 region.
    regions_m = extract_regions(text)
    s11 = regions_m["section 1.1"]
    s11_stripped = re.sub(r"definition 2\.1", "XXX", s11, flags=re.IGNORECASE)
    regions_m["section 1.1"] = s11_stripped
    inv1_mut = check_inv1(regions_m, verbose=False)
    ok &= line(not inv1_mut, "inv1 detects a stripped cross-reference",
               "section 1.1 with Definition 2.1 removed correctly FAILS inv1"
               if not inv1_mut else "MUTATION NOT DETECTED -- inv1 is vacuous!")

    # (D) inv3 must fail if a cost term is deleted from section 10.2.1.
    regions_m2 = extract_regions(text)
    r = regions_m2["section 10.2.1"]
    r_del = re.sub(r"verification overhead", "XXX", r, flags=re.IGNORECASE)
    regions_m2["section 10.2.1"] = r_del
    inv3_mut = check_inv3(regions_m2, verbose=False)
    ok &= line(not inv3_mut, "inv3 detects a deleted cost term",
               "section 10.2.1 with 'verification overhead' removed correctly FAILS inv3"
               if not inv3_mut else "MUTATION NOT DETECTED -- inv3 is vacuous!")

    # (E) inv2 must fire if a real affirmative-universal sentence is injected.
    ok &= selftest_inv2_injection(text)

    # (F) inv4 must fail if a cross-referenced header is fabricated (removed).
    no_def = text.replace("\\textbf{Definition 2.1", "\\textbf{XXX 2.1", 1)
    inv4_mut = check_inv4(no_def, verbose=False)
    ok &= line(not inv4_mut, "inv4 detects a dangling (fabricated) cross-reference",
               "removing the Definition 2.1 header correctly FAILS inv4"
               if not inv4_mut else "MUTATION NOT DETECTED -- inv4 is vacuous!")

    return ok


def check_full_document(label, text):
    """Full battery (inv1-inv4 + self-test) for a document that carries the
    Abstract / section 1.1 / section 10.2.1 regions and the Definition 2.1 /
    Theorem 8.1 headers (monolith and Part I)."""
    print("\n" + "-" * 74)
    print(f"Document: {label}")
    regions = extract_regions(text)
    return {
        f"[{label}] inv1 conditional scoping + cross-refs (Abstract/1.1)":
            check_inv1(regions),
        f"[{label}] inv2 no affirmative universal-dominance sentence":
            check_inv2(text),
        f"[{label}] inv3 precondition-inequality term-completeness (10.2.1)":
            check_inv3(regions),
        f"[{label}] inv4 named cross-references resolve to headers":
            check_inv4(text),
        f"[{label}] self-test non-vacuity": selftest(text),
    }


def check_companion_part(label, text):
    """Companion-part battery (Parts II/III): the whole-paper inv2 scan follows
    the body text that moved out of the monolith, plus its injection
    self-test. inv1/inv3/inv4 anchors (Def 2.1, Thm 8.1, section 10.2.1)
    stayed in Part I and do not apply here."""
    print("\n" + "-" * 74)
    print(f"Document: {label}")
    inv2_ok = check_inv2(text)
    print("\nself-test (non-vacuity): injected affirmative-universal MUST fail")
    inj_ok = selftest_inv2_injection(text)
    return {
        f"[{label}] inv2 no affirmative universal-dominance sentence": inv2_ok,
        f"[{label}] self-test non-vacuity": inj_ok,
    }


def main():
    print("Verification: BUG-010 conditional-scope / no-universal-dominance invariants")
    print("=" * 74)
    print("DOCUMENT-INVARIANT GUARD -- asserts the central claim stays CONDITIONAL,")
    print("carries its Def 2.1 + Thm 8.1 + section 10.2 cross-references, makes NO")
    print("affirmative universal-dominance claim, and that the section 10.2.1")
    print("precondition inequality lists every cost term with a section anchor.")
    print("Checked in every document the invariant appears in after the split:")
    print("monolith (source of truth) and Part I carry all four invariants; the")
    print("whole-paper inv2 scan follows the moved body text into Parts II/III.")
    print("Does NOT measure any compression ratio (that is BUG-010-D, external).")

    results = {}
    results.update(check_full_document(
        "monolith tex/rnr_coding.tex", read(MAIN)))
    results.update(check_full_document(
        "Part I tex/papers/core/rnr_core.tex", read(CORE)))
    results.update(check_companion_part(
        "Part II tex/papers/random_access/rnr_random_access.tex", read(RA)))
    results.update(check_companion_part(
        "Part III tex/papers/dispersion/rnr_dispersion.tex", read(DISP)))

    print("\n" + "=" * 74)
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

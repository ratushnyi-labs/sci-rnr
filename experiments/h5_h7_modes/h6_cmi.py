#!/usr/bin/env python3
"""H6 -- Type-III-B mutual-information advantage I(H;L|C)
(exp-design section 2.5; main paper Theorems 6.3/6.3b).

Extends the preconditions campaign's measurement to the full corpus
set with the CM predictor as the strong context instrument.  Theorem
6.3's saving depends on the operative decoder-reproducible context C,
so three instruments are reported per corpus:

  prev1  C = previous byte (replication of the preconditions screen)
  cm     C = the CM predictor's argmax byte (integer context-mixing
         model, orders 0..3, decoder-reproducible function of the
         whole past -- the strong local instrument)
  plane  (record-structured corpora: telemetry_grid float32, stride 4;
         sqlite_synth page-structured is reported with stride 4 as a
         calibration column only)  C = (byte plane, previous
         same-plane byte); also decomposed per plane.

All estimates are Miller-Madow corrected plug-ins with chunk-bootstrap
percentile CIs (secondary, labeled).  Verdicts compare against the H6
bands: ASCII text 0.3-0.7 bpb, float32 1.0-2.0 bpb, falsification
floor 0.1 bpb; a class is ADVERSE only if every instrument sits below
the floor (the operative claim is existential over practical C).
"""

from __future__ import annotations

import numpy as np

import common
from common import loader

CORPORA = ["enwik8", "canterbury", "scripts_src", "calgary",
           "sqlite_synth", "telemetry_grid",
           "ctrl_urandom", "ctrl_zstd", "ctrl_base64"]

BANDS = {   # corpus -> (class label, band lo, band hi) or None
    "enwik8": ("ascii-text", 0.3, 0.7),
    "canterbury": ("ascii-text", 0.3, 0.7),
    "scripts_src": ("ascii-text", 0.3, 0.7),
    "calgary": ("mixed", None, None),
    "sqlite_synth": ("structured-binary", None, None),
    "telemetry_grid": ("float32", 1.0, 2.0),
    "ctrl_urandom": ("control", None, None),
    "ctrl_zstd": ("control", None, None),
    "ctrl_base64": ("control", None, None),
}
FLOOR = 0.1
PLANE_CORPORA = {"telemetry_grid": 4, "sqlite_synth": 4}


def _with_drift(c, s, c_card, seed, n_resamples, stride: int = 1) -> dict:
    """CI + two bias instruments.

    (a) conditional-permutation debias (common.cmi_perm_null): the MM
        statistic recomputed with low nibbles permuted within context
        strata -- exact null H (indep) L | C with the same support
        geometry; `debiased` = mm - null_mm is the reported point, and
        the CI is the mm CI shifted by the fixed null (the null's own
        sampling noise is not propagated -- labeled approximation).
    (b) split-half drift = mm(full) - mm(first half): positive drift
        means the estimate is still converging downward with n."""
    r = common.cmi_block_bootstrap(c, s, c_card, seed, n_resamples)
    r.update(common.cmi_mm(c, s, c_card))
    r["null_mm"] = common.cmi_perm_null(c, s, c_card, seed + 7777)
    r["debiased"] = r["mm"] - r["null_mm"]
    r["debiased_lo"] = r["lo"] - r["null_mm"]
    r["debiased_hi"] = r["hi"] - r["null_mm"]
    half = s.size // 2
    r["half_mm"] = common.cmi_mm(c[:half], s[:half], c_card)["mm"]
    r["half_drift"] = r["mm"] - r["half_mm"]
    return r


def instruments_for(name: str, data: bytes, data_cm: bytes,
                    n_resamples: int, seed: int) -> dict:
    """Cheap instruments (prev1, plane) run on `data` (MiB scale, to
    match the preconditions campaign's indicator slices); the CM
    instrument runs on `data_cm` (small slice; integer CM inference is
    ~15 KB/s)."""
    x = np.frombuffer(data, dtype=np.uint8).astype(np.int64)
    xc = np.frombuffer(data_cm, dtype=np.uint8).astype(np.int64)
    out = {}
    # MiB-scale instruments get a capped resample count (each resample
    # gathers the whole slice); their CIs are narrow at this n anyway.
    nr_cheap = min(n_resamples, 200)

    # prev1: C = previous byte
    out["prev1"] = _with_drift(x[:-1], x[1:], 256, seed, nr_cheap)

    # cm: C = CM argmax byte (position 0 predicts from Q_init; keep all)
    cm_ctx = common.cm_argmax_stream(data_cm, W=3)
    out["cm"] = _with_drift(cm_ctx, xc, 256, seed + 1, n_resamples)
    out["cm"]["slice_bytes"] = len(data_cm)

    # plane decomposition for record-structured corpora
    stride = PLANE_CORPORA.get(name)
    if stride:
        plane = np.arange(x.size) % stride
        c2 = plane[stride:] * 256 + x[:-stride]      # (plane, lag-k byte)
        s2 = x[stride:]
        out["plane_lag%d" % stride] = _with_drift(
            c2, s2, 256 * stride, seed + 2, nr_cheap, stride=stride)
        per_plane = {}
        plane_cm = np.arange(xc.size) % stride
        for p in range(stride):
            m = plane[stride:] == p                  # positions in plane p
            cp = x[:-stride][m]                      # previous same-plane byte
            sp = s2[m]
            rp = common.cmi_mm(cp, sp, 256)
            rp["null_mm"] = common.cmi_perm_null(cp, sp, 256, seed + 30 + p)
            rp["debiased"] = rp["mm"] - rp["null_mm"]
            mc = plane_cm[stride:] == p
            ccm, scm = cm_ctx[stride:][mc], xc[stride:][mc]
            rp_cm = common.cmi_mm(ccm, scm, 256)
            rp_cm["null_mm"] = common.cmi_perm_null(ccm, scm, 256,
                                                    seed + 40 + p)
            rp_cm["debiased"] = rp_cm["mm"] - rp_cm["null_mm"]
            per_plane[str(p)] = {"lag%d" % stride: rp, "cm": rp_cm}
        out["per_plane"] = per_plane
    return out


def verdict_for(name: str, inst: dict) -> dict:
    cls, lo, hi = BANDS[name]
    points = {k: v["debiased"] for k, v in inst.items() if k != "per_plane"}
    strongest = max(points, key=points.get)
    peak = points[strongest]
    v = {"class": cls, "points": points, "strongest": strongest,
         "peak_mm": peak}
    if lo is None:
        v["verdict"] = "reported (no H6 band for this class)"
        if cls == "control":
            v["verdict"] = ("anchor-consistent (expected ~0)"
                            if abs(peak) < 0.05 or name == "ctrl_base64"
                            else "UNEXPECTED nonzero")
        return v
    if peak < FLOOR:
        v["verdict"] = "ADVERSE (all instruments below 0.1 floor)"
    elif lo <= peak <= hi:
        v["verdict"] = "in predicted band [%g, %g]" % (lo, hi)
    elif peak > hi:
        v["verdict"] = "above predicted band [%g, %g] (favorable)" % (lo, hi)
    else:
        v["verdict"] = ("below predicted band [%g, %g] but above the "
                        "0.1 falsification floor" % (lo, hi))
    return v


def run_h6(slice_bytes: int, cm_slice_bytes: int, n_resamples: int) -> dict:
    res = {}
    for i, name in enumerate(CORPORA):
        n = min(slice_bytes, loader.corpus_bytes(name))
        n_cm = min(cm_slice_bytes, loader.corpus_bytes(name))
        data = loader.read_range(name, 0, n)
        data_cm = data[:n_cm]
        seed = common.derive_seed("h6-boot", i * 10)
        inst = instruments_for(name, data, data_cm, n_resamples, seed)
        res[name] = {
            "slice_bytes": n,
            "cm_slice_bytes": n_cm,
            "corpus_sha256": common.corpus_sha(name),
            "instruments": inst,
            "verdict": verdict_for(name, inst),
        }
    return res

#!/usr/bin/env python3
"""h1_measure.py -- H1 first measurement with the CM instrument.

Hypothesis H1 (exp-design section 2.1; Part I Theorems 4.1/4.2): on
English text and source code, an RNR Type-I coder achieves 1.5-2.0 bpb
with the shared 8-16 MB transformer predictor.  That predictor is an
external upgrade; THIS campaign measures the local hybrid instrument
(integer context-mixing predictor, impl_cm/cm_predictor.py) against
(a) the plain order-W n-gram reference predictor of impl/rnr1.py and
(b) the bench baseline compressors, on the text/code corpora of
data/MANIFEST.json, with the decomposition

    achieved bpb = model cross-entropy (predictor quality)
                 + coder redundancy   (arithmetic-coder quantization)
                 + container overhead (header + seek index)
                 [+ sync/state-reset cost, measured separately by a
                    single-sub-block control run]

Statistics per bench/METHODS.md (frozen): BCa bootstrap, B=10000,
95% CIs; paired per-block differences for coder comparisons; seeds via
metrics.derive_seed; deviations are listed in findings.md and flagged
in each record ("deviation" field).

Scales:
    check    -- reduced slices, gate use only (< 10 min with checks)
    campaign -- the recorded first-measurement campaign (default)
    full     -- exp-design-sized slices (parent's real campaign)

Usage:  python -u h1_measure.py [--scale campaign] [--jsonl PATH]
"""

import argparse
import hashlib
import json
import sys
import tarfile
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
sys.path.insert(0, str(_HERE.parent / "data"))
sys.path.insert(0, str(_HERE.parent / "bench"))
sys.path.insert(0, str(_HERE.parent / "impl"))

import loader  # noqa: E402
import metrics  # noqa: E402
import results as results_store  # noqa: E402
import baselines  # noqa: E402
import rnr1  # noqa: E402
from cm_coder import (  # noqa: E402
    TraceRecorder,
    archive_stats,
    pack_cm,
    pack_ngram,
    unpack_cm,
)

W_ORDER = 3          # context order for both rnr predictors
K_SYNC = 65536       # sub-block / sync spacing (block unit of METHODS 2)

CALGARY_TEXT = ["bib", "book1", "book2", "news", "paper1", "paper2",
                "progc", "progl", "progp", "trans"]
CANTERBURY_TEXT = ["alice29.txt", "asyoulik.txt", "cp.html", "fields.c",
                   "grammar.lsp", "lcet10.txt", "plrabn12.txt", "xargs.1"]

SCALES = {
    #            enwik8_dev  member-cap  scripts_src
    "check":     (128 * 1024,  16 * 1024,  64 * 1024),
    "campaign":  (1 << 20,    128 * 1024, 512 * 1024),
    "full":      (4 << 20,    None,       None),
}


def _tar_concat(corpus, members, cap):
    """Concatenate named tar members (each truncated to cap bytes)."""
    parts = []
    spec = []
    with tarfile.open(loader.payload_path(corpus)) as tf:
        index = {Path(m.name).name: m for m in tf.getmembers() if m.isfile()}
        for name in members:
            m = index[name]
            data = tf.extractfile(m).read()
            if cap is not None:
                data = data[:cap]
            parts.append(data)
            spec.append({"member": name, "bytes": len(data)})
    return b"".join(parts), spec


def build_corpora(scale):
    """Deterministic corpus slices for one scale: name -> (bytes, spec)."""
    dev_n, cap, scripts_n = SCALES[scale]
    out = {}
    dev = loader.read_range("enwik8", 1 << 20, dev_n)
    out["enwik8_dev"] = (dev, {"source": "enwik8", "offset": 1 << 20,
                               "bytes": len(dev)})
    cal, cal_spec = _tar_concat("calgary", CALGARY_TEXT, cap)
    out["calgary_text"] = (cal, {"source": "calgary", "members": cal_spec})
    can, can_spec = _tar_concat("canterbury", CANTERBURY_TEXT, cap)
    out["canterbury_text"] = (can, {"source": "canterbury",
                                    "members": can_spec})
    src = loader.payload_path("scripts_src").read_bytes()
    if scripts_n is not None:
        src = src[:scripts_n]
    out["scripts_src"] = (src, {"source": "scripts_src", "prefix_bytes":
                                len(src)})
    return out


def _block_bpbs(sizes_bytes, block_lens):
    return [8.0 * s / n for s, n in zip(sizes_bytes, block_lens)]


def measure_rnr(name, data, packer, unpacker, K=K_SYNC):
    """One rnr coder run: pack (traced), round-trip, decomposition."""
    rec = TraceRecorder()
    t0 = time.perf_counter()
    raw = packer(data, W=W_ORDER, K=K, trace=rec)
    enc_s = time.perf_counter() - t0
    t0 = time.perf_counter()
    back = unpacker(raw)
    dec_s = time.perf_counter() - t0
    roundtrip = back == data
    st = archive_stats(raw)
    n = len(data)
    ce_bits = sum(b["ce_bits"] for b in rec.blocks)
    block_lens = [b["n"] for b in rec.blocks]
    per_block_repair = st["per_block_coded_bytes"]
    measured = {
        "bpb": 8.0 * len(raw) / n,
        "compressed_size": len(raw),
        "original_size": n,
        "enc_s": enc_s,
        "dec_s": dec_s,
        "repair_bpb": 8.0 * st["repair_bytes"] / n,
        "model_ce_bpb": ce_bits / n,
        "coder_redundancy_bpb": (8.0 * st["repair_bytes"] - ce_bits) / n,
        "container_overhead_bpb": 8.0 * st["container_bytes"] / n,
        "roundtrip_bit_exact": bool(roundtrip),
        "raw_mode_blocks": int(sum(1 for m in st["modes"]
                                   if m == rnr1.MODE_RAW)),
    }
    detail = {
        "block_lens": block_lens,
        "per_block_repair_bpb": _block_bpbs(per_block_repair, block_lens),
        "per_block_ce_bpb": [b["ce_bits"] / b["n"] for b in rec.blocks],
    }
    print("  %-12s bpb=%.4f  modelCE=%.4f  redund=%.5f  cont=%.5f  "
          "rt=%s  enc=%.1fs dec=%.1fs"
          % (name, measured["bpb"], measured["model_ce_bpb"],
             measured["coder_redundancy_bpb"],
             measured["container_overhead_bpb"], roundtrip, enc_s, dec_s),
          flush=True)
    return measured, detail


def measure_baseline_stream(codec, data):
    r = baselines.run_baseline(codec, data)  # includes round-trip check
    return {
        "bpb": r.bits_per_byte,
        "compressed_size": r.compressed_size,
        "original_size": r.original_size,
        "enc_s": r.enc_s,
        "dec_s": r.dec_s,
        "peak_rss_estimate": r.peak_rss_estimate,
        "rss_mode": r.rss_mode,  # in-process-delta = lower bound, METHODS 3
    }, r.config


def measure_baseline_blocks(codec, data, K=K_SYNC):
    """Per-block compressed sizes (paired comparison unit, METHODS 2)."""
    cfg, comp, decomp = baselines._CODECS[codec]
    sizes, lens = [], []
    for off in range(0, len(data), K):
        block = data[off:off + K]
        blob = comp(block, cfg)
        if decomp(blob, cfg) != block:
            raise baselines.RoundTripError("%s block round-trip" % codec)
        sizes.append(len(blob))
        lens.append(len(block))
    return sizes, lens


def _ci_dict(ci):
    return {"point": ci.point, "lo": ci.lo, "hi": ci.hi,
            "method": ci.method, "n_resamples": ci.n_resamples}


def run_campaign(scale, jsonl_path=None, log=print):
    corpora = build_corpora(scale)
    codecs = baselines.available_codecs()
    versions = baselines.codec_versions()
    summary = {"scale": scale, "W": W_ORDER, "K": K_SYNC,
               "codecs": codecs, "versions": versions, "corpora": {}}
    deviation = ("timing single-run (repeats=1, no warmup) for wall-budget; "
                 "rates deterministic; see findings.md amendment note")

    def emit(corpus, coder, config, measured, sha, spec, **extra):
        if jsonl_path is None:
            return
        rec = results_store.make_record(
            corpus, coder, config,
            metrics.derive_seed("H1:%s:%s" % (corpus, coder), 0, tag="run"),
            measured,
            corpus_sha256=sha, slice_spec=spec,
            experiment="H1-first-measurement", scale=scale,
            deviation=deviation, versions=versions, **extra)
        results_store.append(jsonl_path, rec)

    for corpus, (data, spec) in corpora.items():
        sha = hashlib.sha256(data).hexdigest()
        n = len(data)
        log("[%s] %s: %d bytes sha=%s" % (scale, corpus, n, sha[:12]),
            flush=True)
        centry = {"bytes": n, "sha256": sha, "spec": spec, "coders": {}}

        cm_meas, cm_det = measure_rnr("rnr1-cm", data, pack_cm, unpack_cm)
        emit(corpus, "rnr1-cm", {"W": W_ORDER, "K": K_SYNC}, cm_meas, sha,
             spec, per_block=cm_det)
        centry["coders"]["rnr1-cm"] = {"measured": cm_meas, **cm_det}

        ng_meas, ng_det = measure_rnr("rnr1-ngram", data, pack_ngram,
                                      lambda raw: rnr1.unpack(raw))
        emit(corpus, "rnr1-ngram", {"W": W_ORDER, "K": K_SYNC}, ng_meas, sha,
             spec, per_block=ng_det)
        centry["coders"]["rnr1-ngram"] = {"measured": ng_meas, **ng_det}

        for codec in codecs:
            meas, cfg = measure_baseline_stream(codec, data)
            emit(corpus, codec, {**cfg, "mode": "stream"}, meas, sha, spec)
            sizes, lens = measure_baseline_blocks(codec, data)
            bmeas = {"bpb": 8.0 * sum(sizes) / n,
                     "compressed_size": sum(sizes), "original_size": n}
            emit(corpus, codec, {**cfg, "mode": "block"}, bmeas, sha, spec,
                 per_block={"block_lens": lens,
                            "per_block_bpb": _block_bpbs(sizes, lens)})
            centry["coders"][codec] = {
                "stream": meas, "block": bmeas,
                "per_block_bpb": _block_bpbs(sizes, lens)}
            log("  %-12s stream=%.4f bpb  block=%.4f bpb"
                % (codec, meas["bpb"], bmeas["bpb"]), flush=True)

        # Statistics: BCa CI on mean per-block repair bpb; paired diffs.
        stats = {}
        cm_blocks = cm_det["per_block_repair_bpb"]
        ng_blocks = ng_det["per_block_repair_bpb"]
        if len(cm_blocks) >= 3:
            stats["rnr1-cm_block_bpb_ci"] = _ci_dict(metrics.bootstrap_ci(
                cm_blocks,
                seed=metrics.derive_seed("H1:%s:rnr1-cm" % corpus, 0,
                                         tag="boot")))
            diffs = {"rnr1-ngram": ng_blocks}
            for codec in codecs:
                diffs[codec] = centry["coders"][codec]["per_block_bpb"]
            stats["paired_diff_vs"] = {}
            for other, blocks in diffs.items():
                ci = metrics.paired_bootstrap_diff(
                    cm_blocks, blocks,
                    seed=metrics.derive_seed(
                        "H1:%s:rnr1-cm-vs-%s" % (corpus, other), 0,
                        tag="boot"))
                stats["paired_diff_vs"][other] = {
                    **_ci_dict(ci),
                    "cohens_d": metrics.cohens_d_paired(cm_blocks, blocks)}
        # Paired CE comparison (predictor quality, per block).
        ce_ci = metrics.paired_bootstrap_diff(
            cm_det["per_block_ce_bpb"], ng_det["per_block_ce_bpb"],
            seed=metrics.derive_seed("H1:%s:ce-cm-vs-ngram" % corpus, 0,
                                     tag="boot")) \
            if len(cm_blocks) >= 2 else None
        if ce_ci is not None:
            stats["ce_diff_cm_vs_ngram"] = _ci_dict(ce_ci)
        centry["stats"] = stats
        summary["corpora"][corpus] = centry

    # Sync/state-reset cost control: single-sub-block CM run on enwik8_dev.
    if scale != "check":
        data, _ = corpora["enwik8_dev"]
        rec = TraceRecorder()
        t0 = time.perf_counter()
        raw = pack_cm(data, W=W_ORDER, K=len(data), trace=rec)
        enc_s = time.perf_counter() - t0
        ce = sum(b["ce_bits"] for b in rec.blocks) / len(data)
        base_ce = summary["corpora"]["enwik8_dev"]["coders"]["rnr1-cm"][
            "measured"]["model_ce_bpb"]
        summary["sync_reset_cost"] = {
            "corpus": "enwik8_dev", "K_control": len(data),
            "model_ce_bpb_no_reset": ce,
            "model_ce_bpb_K%d" % K_SYNC: base_ce,
            "state_reset_cost_bpb": base_ce - ce,
            "archive_bpb_no_reset": 8.0 * len(raw) / len(data),
            "enc_s": enc_s,
        }
        log("[control] no-reset CM on enwik8_dev: CE %.4f bpb "
            "(reset cost %.4f bpb)" % (ce, base_ce - ce), flush=True)
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--scale", choices=sorted(SCALES), default="campaign")
    ap.add_argument("--jsonl", default=None,
                    help="results store path (default results/h1_<scale>.jsonl)")
    ap.add_argument("--summary", default=None,
                    help="summary JSON path (default results/h1_<scale>_summary.json)")
    args = ap.parse_args(argv)
    jsonl = Path(args.jsonl) if args.jsonl else (
        _HERE / "results" / ("h1_%s.jsonl" % args.scale))
    summary_path = Path(args.summary) if args.summary else (
        _HERE / "results" / ("h1_%s_summary.json" % args.scale))
    t0 = time.perf_counter()
    summary = run_campaign(args.scale, jsonl_path=jsonl)
    summary["wall_s"] = time.perf_counter() - t0
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=1, sort_keys=True))
    print("campaign done in %.1f s; results -> %s, summary -> %s"
          % (summary["wall_s"], jsonl, summary_path), flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Paper-ready figures + summary tables from the Type-I scale campaign
JSONL stores (protocol.md section 6).

Deterministic: figures depend only on the JSONL contents (sorted
iteration everywhere; PDF CreationDate stripped), so re-running on the
same data reproduces byte-identical files.

Outputs (out/figures/ by default):
  fig_rate_vs_size.{pdf,png}      bpb vs size, one panel per type
  fig_enc_time_vs_size.{pdf,png}  encode wall s vs size, log-log
  fig_dec_time_vs_size.{pdf,png}  decode wall s vs size, log-log
  fig_pareto_<size>.{pdf,png}     bpb vs encode throughput, per type
  summary_<size>_<type>.csv       one table per (size, type) cell
  captions.md                     caption text for every figure/table

Design: colorblind-safe palette (validated: worst adjacent CVD
deltaE 25.0, >= 12 target); coder families carry the hue, variants the
linestyle, every coder a distinct marker (secondary encoding); no
in-figure titles (captions live in captions.md); every figure has a
table view (the CSVs).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent

TYPES = ("white_noise", "text", "video", "images",
         "archive_mix", "precompressed")

# Family hue (validated ordering), per-coder linestyle + marker.
FAMILY_COLOR = {
    "rnr1": "#2a78d6",     # blue
    "brotli": "#eb6834",   # orange
    "zstd": "#1baf7a",     # aqua
    "xz": "#4a3aa7",       # violet
    "gzip": "#008300",     # green
    "bzip2": "#e87ba4",    # magenta
}
CODER_STYLE = {
    # coder: (family, linestyle, marker)
    "rnr1-fast": ("rnr1", "-", "o"),
    "rnr1-cm-fast": ("rnr1", "--", "D"),
    "rnr1-ref": ("rnr1", ":", "s"),
    "zstd-19": ("zstd", "-", "^"),
    "zstd-22u": ("zstd", "--", "v"),
    "xz-6": ("xz", "-", "P"),
    "xz-9e": ("xz", "--", "X"),
    "brotli-q11": ("brotli", "-", "*"),
    "brotli-q9": ("brotli", "--", "p"),
    "gzip-9": ("gzip", "-", "d"),
    "bzip2-9": ("bzip2", "-", "h"),
}

RC = {
    "font.size": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#dddbd6",
    "grid.linewidth": 0.5,
    "axes.edgecolor": "#52514e",
    "axes.labelcolor": "#0b0b0b",
    "xtick.color": "#52514e",
    "ytick.color": "#52514e",
    "lines.linewidth": 1.4,
    "lines.markersize": 4.5,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
}

PDF_META = {"CreationDate": None}  # deterministic PDFs


def _style(coder: str):
    fam, ls, mk = CODER_STYLE.get(coder, ("rnr1", "-", "o"))
    return FAMILY_COLOR[fam], ls, mk


# ---------------------------------------------------------------------------
# Data loading / aggregation
# ---------------------------------------------------------------------------

def load_records(paths: list) -> list:
    recs = []
    for p in paths:
        with open(p) as f:
            for line in f:
                line = line.strip()
                if line:
                    recs.append(json.loads(line))
    # Drop superseded records (METHODS.md section 6).  Runset-store records
    # (per-run generated test sets, run_testset.py) have their own schema
    # and analysis pipeline -- not tier-campaign cells; skip them here.
    recs = [r for r in recs if "corpus" in r]
    superseded = {tuple(r["supersedes"]) for r in recs if r.get("supersedes")}
    out = []
    for r in recs:
        key = (r["corpus"], r["coder"],
               json.dumps(r["config"], sort_keys=True), r["seed"],
               r["timestamp"])
        if key not in superseded and list(key) not in [list(s) for s in
                                                       superseded]:
            out.append(r)
    return out


def size_mb(rec: dict) -> float:
    n = rec["measured"].get("original_size")
    if n is None:
        import corpora
        n = corpora.TIER_BYTES[rec["tier"]]
    return n / 1e6


def aggregate(recs: list) -> dict:
    """(ctype, size_mb, coder) -> dict of mean metrics + flags."""
    groups: dict = {}
    for r in recs:
        key = (r["ctype"], round(size_mb(r), 3), r["coder"])
        groups.setdefault(key, []).append(r)
    agg = {}
    for key, rs in sorted(groups.items()):
        ok = [r for r in rs if r.get("status") == "done"
              and r["measured"].get("bpb") is not None]
        censored = [r for r in rs if r.get("censored")]
        failed = [r for r in rs if r.get("status") == "failed"]

        def mean(field):
            vals = [r["measured"][field] for r in ok
                    if r["measured"].get(field) is not None]
            return sum(vals) / len(vals) if vals else None

        agg[key] = {
            "n_runs": len(ok),
            "n_censored": len(censored),
            "n_failed": len(failed),
            "bpb": mean("bpb"),
            "ratio": mean("ratio"),
            "enc_s": mean("enc_s"),
            "dec_s": mean("dec_s"),
            "enc_peak_rss": mean("enc_peak_rss"),
            "dec_peak_rss": mean("dec_peak_rss"),
            "roundtrip_ok": all(r["measured"].get("roundtrip_ok")
                                for r in ok) if ok else None,
            "censor_phase": (censored[0].get("censor_phase")
                             if censored else None),
        }
    return agg


def _present(agg, field="bpb"):
    """Sorted unique (types, sizes, coders) that actually have data."""
    types = sorted({k[0] for k in agg}, key=lambda t: TYPES.index(t)
                   if t in TYPES else 99)
    sizes = sorted({k[1] for k in agg})
    coders = sorted({k[2] for k in agg},
                    key=lambda c: list(CODER_STYLE).index(c)
                    if c in CODER_STYLE else 99)
    return types, sizes, coders


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _panel_grid(n):
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    return rows, cols


def _legend_handles(coders):
    from matplotlib.lines import Line2D
    handles = []
    for c in coders:
        col, ls, mk = _style(c)
        handles.append(Line2D([0], [0], color=col, linestyle=ls, marker=mk,
                              markersize=4.5, label=c))
    return handles


def _save(fig, out_dir: Path, stem: str, made: list):
    pdf, png = out_dir / f"{stem}.pdf", out_dir / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight", metadata=PDF_META)
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)
    made.extend([pdf, png])


def _metric_vs_size(agg, metric, ylabel, stem, out_dir, made,
                    logy=True, ylim_pad=None):
    types, sizes, coders = _present(agg)
    rows, cols = _panel_grid(len(types))
    fig, axes = plt.subplots(rows, cols,
                             figsize=(2.6 * cols + 0.6, 2.2 * rows + 0.8),
                             squeeze=False, sharex=True)
    for i, ctype in enumerate(types):
        ax = axes[i // cols][i % cols]
        for coder in coders:
            pts = [(s, agg[(ctype, s, coder)][metric]) for s in sizes
                   if (ctype, s, coder) in agg
                   and agg[(ctype, s, coder)][metric] is not None]
            if not pts:
                continue
            xs, ys = zip(*sorted(pts))
            col, ls, mk = _style(coder)
            ax.plot(xs, ys, color=col, linestyle=ls, marker=mk)
        ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
        ax.text(0.0, 1.02, ctype, transform=ax.transAxes, fontsize=7.5,
                va="bottom", color="#52514e")
        if i // cols == rows - 1:
            ax.set_xlabel("corpus size (MB)")
        if i % cols == 0:
            ax.set_ylabel(ylabel)
    for j in range(len(types), rows * cols):
        axes[j // cols][j % cols].set_visible(False)
    fig.legend(handles=_legend_handles(coders), loc="lower center",
               ncol=min(len(coders), 6), bbox_to_anchor=(0.5, -0.06),
               fontsize=7)
    fig.tight_layout()
    _save(fig, out_dir, stem, made)


def fig_rate_vs_size(agg, out_dir, made):
    _metric_vs_size(agg, "bpb", "rate (bits/byte)", "fig_rate_vs_size",
                    out_dir, made, logy=False)


def fig_enc_time(agg, out_dir, made):
    _metric_vs_size(agg, "enc_s", "encode wall time (s)",
                    "fig_enc_time_vs_size", out_dir, made, logy=True)


def fig_dec_time(agg, out_dir, made):
    _metric_vs_size(agg, "dec_s", "decode wall time (s)",
                    "fig_dec_time_vs_size", out_dir, made, logy=True)


def _pareto_front(points):
    """points: list of (enc_MBps, bpb).  Front = max throughput for any
    given rate: no other point has both higher throughput and lower bpb."""
    front = []
    for p in points:
        if not any(q[0] >= p[0] and q[1] <= p[1] and q != p for q in points):
            front.append(p)
    return sorted(front)


def fig_pareto(agg, out_dir, made, captions):
    types, sizes, coders = _present(agg)
    for s in sizes:
        stem = f"fig_pareto_{s:g}MB"
        present_types = [t for t in types
                         if any((t, s, c) in agg for c in coders)]
        if not present_types:
            continue
        rows, cols = _panel_grid(len(present_types))
        fig, axes = plt.subplots(rows, cols,
                                 figsize=(2.6 * cols + 0.6, 2.3 * rows + 0.8),
                                 squeeze=False)
        for i, ctype in enumerate(present_types):
            ax = axes[i // cols][i % cols]
            pts = []
            for coder in coders:
                a = agg.get((ctype, s, coder))
                if not a or a["bpb"] is None or not a["enc_s"]:
                    continue
                thr = (s * 1e6 / 1e6) / a["enc_s"]  # MB/s
                col, ls, mk = _style(coder)
                ax.scatter([thr], [a["bpb"]], color=col, marker=mk, s=26,
                           zorder=3)
                pts.append((thr, a["bpb"]))
            if len(pts) >= 2:
                fr = _pareto_front(pts)
                ax.plot([p[0] for p in fr], [p[1] for p in fr],
                        color="#52514e", linewidth=0.8, linestyle="-",
                        zorder=2, alpha=0.7)
            ax.set_xscale("log")
            ax.text(0.0, 1.02, ctype, transform=ax.transAxes,
                    fontsize=7.5, va="bottom", color="#52514e")
            if i // cols == rows - 1:
                ax.set_xlabel("encode throughput (MB/s)")
            if i % cols == 0:
                ax.set_ylabel("rate (bits/byte)")
        for j in range(len(present_types), rows * cols):
            axes[j // cols][j % cols].set_visible(False)
        fig.legend(handles=_legend_handles(coders), loc="lower center",
                   ncol=min(len(coders), 6), bbox_to_anchor=(0.5, -0.06),
                   fontsize=7)
        fig.tight_layout()
        _save(fig, out_dir, stem, made)
        captions.append(
            f"**{stem}**: Rate-versus-encode-throughput Pareto plot at "
            f"corpus size {s:g} MB, one panel per content type. Each point "
            f"is one coder (marker/hue as in the shared legend); the gray "
            f"staircase joins the Pareto-efficient coders (no coder is "
            f"both faster to encode and denser). Throughput axis is "
            f"logarithmic. Censored/failed cells are absent (see the "
            f"summary CSVs).")


# ---------------------------------------------------------------------------
# Summary CSVs
# ---------------------------------------------------------------------------

CSV_FIELDS = ["coder", "n_runs", "bpb", "ratio", "enc_s", "dec_s",
              "enc_MBps", "dec_MBps", "enc_peak_rss_MB", "dec_peak_rss_MB",
              "roundtrip_ok", "n_censored", "censor_phase", "n_failed"]


def write_csvs(agg, out_dir, made, captions):
    types, sizes, coders = _present(agg)
    for s in sizes:
        for ctype in types:
            rows = []
            for coder in coders:
                a = agg.get((ctype, s, coder))
                if a is None:
                    continue
                mbps = (lambda t: (s / t) if t else None)

                def fmt(x, nd=4):
                    return round(x, nd) if isinstance(x, float) else x

                rows.append({
                    "coder": coder,
                    "n_runs": a["n_runs"],
                    "bpb": fmt(a["bpb"]),
                    "ratio": fmt(a["ratio"], 5),
                    "enc_s": fmt(a["enc_s"], 3),
                    "dec_s": fmt(a["dec_s"], 3),
                    "enc_MBps": fmt(mbps(a["enc_s"]), 3),
                    "dec_MBps": fmt(mbps(a["dec_s"]), 3),
                    "enc_peak_rss_MB": fmt((a["enc_peak_rss"] or 0) / 1e6, 1)
                    if a["enc_peak_rss"] is not None else None,
                    "dec_peak_rss_MB": fmt((a["dec_peak_rss"] or 0) / 1e6, 1)
                    if a["dec_peak_rss"] is not None else None,
                    "roundtrip_ok": a["roundtrip_ok"],
                    "n_censored": a["n_censored"],
                    "censor_phase": a["censor_phase"] or "",
                    "n_failed": a["n_failed"],
                })
            if not rows:
                continue
            path = out_dir / f"summary_{s:g}MB_{ctype}.csv"
            with open(path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
                w.writeheader()
                w.writerows(rows)
            made.append(path)
    captions.append(
        "**summary_<size>MB_<type>.csv**: per-(size, type) table view of "
        "every metric -- mean over runs of rate (bpb), compression ratio, "
        "encode/decode wall seconds and derived MB/s, peak RSS per phase "
        "(/usr/bin/time -l child high-watermark, MB), round-trip "
        "verification, and censoring/failure counts. This is the "
        "accessible table view backing every figure.")


BASE_CAPTIONS = [
    "**fig_rate_vs_size**: Compression rate (bits per byte, lower is "
    "better) versus corpus size, one panel per content type; size axis "
    "logarithmic. Hue = coder family, linestyle = variant, marker = "
    "coder (colorblind-safe; validated worst adjacent CVD deltaE 25).",
    "**fig_enc_time_vs_size**: Encode wall-clock seconds versus corpus "
    "size, log-log, one panel per content type. Single-threaded CLI "
    "invocations measured around the whole child process.",
    "**fig_dec_time_vs_size**: Decode wall-clock seconds versus corpus "
    "size, log-log, one panel per content type; decode output is "
    "streamed to a hash, not written to disk (except RNR coders, see "
    "protocol 4.3).",
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", nargs="+", default=None,
                    help="JSONL stores (default: out/*/results.jsonl)")
    ap.add_argument("--out-dir", default=str(HERE / "out" / "figures"))
    args = ap.parse_args(argv)

    paths = ([Path(p) for p in args.results] if args.results
             else sorted((HERE / "out").glob("*/results.jsonl")))
    if not paths:
        print("no results.jsonl stores found", file=sys.stderr)
        return 2
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    recs = load_records(paths)
    agg = aggregate(recs)
    if not agg:
        print("no usable records", file=sys.stderr)
        return 2

    plt.rcParams.update(RC)
    made: list = []
    captions = list(BASE_CAPTIONS)
    fig_rate_vs_size(agg, out_dir, made)
    fig_enc_time(agg, out_dir, made)
    fig_dec_time(agg, out_dir, made)
    fig_pareto(agg, out_dir, made, captions)
    write_csvs(agg, out_dir, made, captions)

    cap_path = out_dir / "captions.md"
    with open(cap_path, "w") as f:
        f.write("# Figure and table captions (Type-I scale campaign)\n\n")
        f.write("Generated deterministically by make_figures.py from: "
                + ", ".join(str(p.relative_to(HERE)) if p.is_relative_to(HERE)
                            else str(p) for p in paths) + "\n\n")
        for c in captions:
            f.write(c + "\n\n")
    made.append(cap_path)

    for p in made:
        print(f"[figures] wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Append-only JSON-lines results store for the RNR empirical program.

Every measurement is one JSON object on one line, keyed by
(corpus, coder, config, seed, timestamp). Records are appended, never
edited or deleted (METHODS.md section 6): a correction is a NEW record
whose "supersedes" field names the key of the record it replaces.

summarize() renders the store as a Markdown table grouped by
(corpus, coder, config), aggregating over seeds/runs, with a bootstrap
95% CI on the mean bits-per-byte where at least 3 runs exist.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import metrics as _metrics

REQUIRED_KEYS = ("corpus", "coder", "config", "seed", "timestamp")


def make_record(
    corpus: str,
    coder: str,
    config: dict,
    seed: int,
    measured: dict,
    *,
    timestamp: str | None = None,
    **extra,
) -> dict:
    """Assemble a well-formed results record.

    `measured` holds the metric values, e.g. {"bpb": 2.31,
    "compressed_size": 302123, "original_size": 1048576,
    "enc_s": 0.91, "dec_s": 0.02, "peak_rss_estimate": 12345678}.
    Extra keyword fields (versions, hardware, corpus_sha256, ...) are
    stored verbatim.
    """
    record = {
        "corpus": corpus,
        "coder": coder,
        "config": dict(config),
        "seed": int(seed),
        "timestamp": timestamp
        or datetime.now(timezone.utc).isoformat(timespec="microseconds"),
        "measured": dict(measured),
    }
    record.update(extra)
    return record


def record_key(record: dict) -> tuple:
    """Identity key of a record: (corpus, coder, canonical-config, seed, timestamp)."""
    return (
        record["corpus"],
        record["coder"],
        json.dumps(record["config"], sort_keys=True),
        record["seed"],
        record["timestamp"],
    )


def append(path, record: dict) -> None:
    """Append one record as a single JSON line. Never rewrites the file."""
    missing = [k for k in REQUIRED_KEYS if k not in record]
    if missing:
        raise ValueError(f"record is missing required keys: {missing}")
    line = json.dumps(record, sort_keys=True, separators=(",", ":"))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", encoding="utf-8") as f:
        f.write(line + "\n")
        f.flush()


def load(path) -> list:
    """Read all records from a JSONL store; raises on malformed lines."""
    records = []
    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{lineno}: malformed JSON line") from exc
    return records


def _fmt(value, digits: int = 4) -> str:
    if value is None:
        return "--"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def summarize(records_or_path) -> str:
    """Markdown summary table grouped by (corpus, coder, config).

    Aggregates the "measured" dict over runs: mean bpb with a bootstrap
    95% CI (frozen methodology, METHODS.md; shown only for >= 3 runs),
    mean encode/decode seconds, and max peak-RSS estimate.
    """
    if isinstance(records_or_path, (str, Path)):
        records = load(records_or_path)
    else:
        records = list(records_or_path)

    groups: dict = {}
    for rec in records:
        key = (
            rec["corpus"],
            rec["coder"],
            json.dumps(rec["config"], sort_keys=True),
        )
        groups.setdefault(key, []).append(rec)

    header = (
        "| corpus | coder | config | runs | bpb mean | bpb 95% CI | "
        "enc s | dec s | peak RSS (MB) |\n"
        "|---|---|---|---|---|---|---|---|---|"
    )
    lines = [header]
    for (corpus, coder, config_json) in sorted(groups):
        recs = groups[(corpus, coder, config_json)]
        bpbs = [r["measured"].get("bpb") for r in recs]
        bpbs = [b for b in bpbs if b is not None]
        encs = [r["measured"].get("enc_s") for r in recs]
        encs = [e for e in encs if e is not None]
        decs = [r["measured"].get("dec_s") for r in recs]
        decs = [d for d in decs if d is not None]
        rsss = [r["measured"].get("peak_rss_estimate") for r in recs]
        rsss = [r for r in rsss if r is not None]

        bpb_mean = sum(bpbs) / len(bpbs) if bpbs else None
        if len(bpbs) >= 3 and min(bpbs) < max(bpbs):
            ci = _metrics.bootstrap_ci(bpbs, seed=0)
            ci_str = f"[{ci.lo:.4g}, {ci.hi:.4g}]"
        elif len(bpbs) >= 3:
            ci_str = "[degenerate]"
        else:
            ci_str = "--"
        enc_mean = sum(encs) / len(encs) if encs else None
        dec_mean = sum(decs) / len(decs) if decs else None
        rss_mb = max(rsss) / 1e6 if rsss else None

        lines.append(
            f"| {corpus} | {coder} | `{config_json}` | {len(recs)} "
            f"| {_fmt(bpb_mean)} | {ci_str} "
            f"| {_fmt(enc_mean, 3)} | {_fmt(dec_mean, 3)} | {_fmt(rss_mb, 3)} |"
        )
    if len(lines) == 1:
        lines.append("| (no records) | | | | | | | | |")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("usage: results.py <store.jsonl>  (prints the summary table)")
        raise SystemExit(2)
    print(summarize(sys.argv[1]))

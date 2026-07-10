#!/usr/bin/env python3
"""Full-factorial Type-I scale campaign runner (protocol.md).

Cells = tier x content-type x coder x run-index.  Each cell is one
encode -> size/hash -> decode-verify -> cleanup measurement (coders.py)
appended as one JSONL record (bench/results.py schema) to
out/<tier>/results.jsonl.

Per-cell checkpointing: every terminal cell (done / censored / failed)
appends one JSON line to out/<tier>/cells.journal; an interrupted
campaign re-run skips journaled cells, so progress is never lost and a
cell is never recorded twice.  The journal is written AFTER the results
record, so a crash between the two can at worst re-run one cell and
append a duplicate-keyed record (deduplicated by summaries on
(corpus, coder, config, seed); cell_id is embedded in each record).

Usage:
  run_campaign.py --tier smoke [--types a,b] [--coders x,y]
                  [--cap-min 90] [--out-dir DIR] [--runs N]

Time cap: a cell whose encode+decode exceeds the cap is SIGKILLed and
recorded as censored=true with the phase reached (protocol 5.3).
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "bench"))
sys.path.insert(0, str(HERE))

import corpora  # noqa: E402
import coders  # noqa: E402
import metrics as _metrics  # noqa: E402
import results as _results  # noqa: E402

# Pre-registered repetition policy (protocol 5.2; METHODS amendment note
# methods_amendment_2026-07-10.md): 3 timing runs at 500 MB, 1 at the
# larger tiers and at smoke.
RUNS_PER_TIER = {"smoke": 1, "smoke48": 1, "500": 3, "1500": 1, "5000": 1}
DEFAULT_CAP_MIN = 90.0


def cell_id(tier: str, ctype: str, coder: str, run_idx: int) -> str:
    return f"{tier}/{ctype}/{coder}/r{run_idx}"


def load_journal(path: Path) -> dict:
    """cell_id -> status for every terminal cell."""
    done = {}
    if path.exists():
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    rec = json.loads(line)
                    done[rec["cell"]] = rec["status"]
    return done


def append_journal(path: Path, cid: str, status: str, **extra) -> None:
    rec = {"cell": cid, "status": status,
           "ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    rec.update(extra)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as f:
        f.write(json.dumps(rec, sort_keys=True) + "\n")
        f.flush()


def host_info() -> dict:
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": sys.version.split()[0],
    }


def plan_cells(tier: str, types: list, coder_names: list, runs: int) -> list:
    cells = []
    for ctype in types:
        for coder in coder_names:
            for r in range(runs):
                cells.append((ctype, coder, r))
    return cells


def _clean_stale_scratch(scratch_root: Path) -> None:
    """Remove scratch subdirs left by dead campaign invocations only.
    (Live invocations own scratch/run-<pid>; never touch a live pid's.)"""
    if not scratch_root.exists():
        return
    for sub in scratch_root.iterdir():
        if not sub.is_dir() or not sub.name.startswith("run-"):
            continue
        try:
            pid = int(sub.name.split("-", 1)[1])
        except ValueError:
            continue
        try:
            os.kill(pid, 0)  # probe: raises if pid is gone
        except ProcessLookupError:
            shutil.rmtree(sub, ignore_errors=True)
        except PermissionError:
            pass  # alive under another uid: leave it


def _failed_record_key(results_path: Path, cid: str):
    """Key of the latest failed record for a cell (for `supersedes`)."""
    key = None
    if results_path.exists():
        with open(results_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rec = json.loads(line)
                if rec.get("cell") == cid and rec.get("status") == "failed":
                    key = list(_results.record_key(rec))
    return key


def run(tier: str, types: list, coder_names: list, out_dir: Path,
        cap_min: float, runs: int | None, keep_scratch: bool = False,
        retry_failed: bool = False) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    journal_path = out_dir / "cells.journal"
    results_path = out_dir / "results.jsonl"
    scratch_root = out_dir / "scratch"
    _clean_stale_scratch(scratch_root)
    scratch = scratch_root / f"run-{os.getpid()}"
    nruns = runs if runs is not None else RUNS_PER_TIER[tier]
    cap_s = cap_min * 60.0

    journal = load_journal(journal_path)
    absent = [t for t in types if not corpora.available(tier, t)]
    if absent:
        print(f"[campaign] tier={tier}: skipping undefined cells for types "
              f"{absent} (data/big manifest: text is capped at 1000 MB)")
        types = [t for t in types if t not in absent]
    cells = plan_cells(tier, types, coder_names, nruns)
    todo = [(c, k, r) for (c, k, r) in cells
            if cell_id(tier, c, k, r) not in journal
            or (retry_failed
                and journal[cell_id(tier, c, k, r)] == "failed")]
    print(f"[campaign] tier={tier} cells={len(cells)} "
          f"done={len(cells) - len(todo)} todo={len(todo)} cap={cap_min}min")

    # Materialize corpora once per type (idempotent) and re-hash once per
    # campaign invocation as an integrity check against the manifest pin.
    entries = {}
    for ctype in types:
        entry = corpora.materialize(tier, ctype)
        path = corpora.corpus_path(tier, ctype)
        live = corpora._sha256_file(path)
        if live != entry["sha256"]:
            raise RuntimeError(
                f"corpus {tier}/{ctype} sha256 drifted from manifest "
                f"({live[:16]} != {entry['sha256'][:16]}); delete and rebuild")
        entries[ctype] = entry

    counts = {"done": 0, "censored": 0, "failed": 0,
              "skipped": len(cells) - len(todo)}
    for ctype, coder, run_idx in todo:
        cid = cell_id(tier, ctype, coder, run_idx)
        entry = entries[ctype]
        seed = _metrics.derive_seed(
            f"type1_scale:{tier}:{ctype}:{coder}", run_idx, tag="run")
        corpus_name = f"{ctype}@{tier}"
        base_extra = dict(
            cell=cid, tier=tier, ctype=ctype, run_idx=run_idx,
            corpus_sha256=entry["sha256"],
            coder_version=coders.coder_version(coder),
            host=host_info(), rss_mode="time-l-child",
            cap_s=cap_s,
        )
        if retry_failed and journal.get(cid) == "failed":
            old_key = _failed_record_key(results_path, cid)
            if old_key is not None:
                base_extra["supersedes"] = old_key  # METHODS.md section 6
        print(f"[campaign] {cid} ...", flush=True)
        t0 = time.perf_counter()
        try:
            measured = coders.run_cell(
                coder, corpora.corpus_path(tier, ctype), entry["sha256"],
                scratch, cap_s)
            status = "done"
            if not measured["roundtrip_ok"]:
                status = "failed"  # verification failure is a hard failure
            rec = _results.make_record(
                corpus_name, coder, coders.REGISTRY[coder]["config"],
                seed, measured, censored=False, status=status, **base_extra)
        except coders.CellTimeout as exc:
            status = "censored"
            rec = _results.make_record(
                corpus_name, coder, coders.REGISTRY[coder]["config"],
                seed,
                {"enc_s": None, "dec_s": None, "bpb": None,
                 "roundtrip_ok": None},
                censored=True, status="censored",
                censor_phase=exc.phase, censor_elapsed_s=exc.elapsed_s,
                **base_extra)
        except Exception as exc:  # coder failure: record and continue
            status = "failed"
            rec = _results.make_record(
                corpus_name, coder, coders.REGISTRY[coder]["config"],
                seed,
                {"enc_s": None, "dec_s": None, "bpb": None,
                 "roundtrip_ok": False},
                censored=False, status="failed", error=str(exc)[:800],
                **base_extra)
        _results.append(results_path, rec)
        append_journal(journal_path, cid, status,
                       wall_s=round(time.perf_counter() - t0, 3))
        counts[status] = counts.get(status, 0) + 1
        print(f"[campaign] {cid} -> {status} "
              f"({time.perf_counter() - t0:.1f}s)", flush=True)

    if not keep_scratch and scratch.exists():
        shutil.rmtree(scratch, ignore_errors=True)
    print(f"[campaign] finished: {counts}")
    return counts


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tier", required=True,
                    choices=sorted(corpora.TIER_BYTES))
    ap.add_argument("--types", default=",".join(corpora.CONTENT_TYPES),
                    help="comma-separated content types")
    ap.add_argument("--coders", default=None,
                    help="comma-separated coder names "
                         "(default: all available at the tier)")
    ap.add_argument("--cap-min", type=float, default=DEFAULT_CAP_MIN,
                    help="per-cell time cap in minutes (default 90)")
    ap.add_argument("--runs", type=int, default=None,
                    help="override the tier repetition policy")
    ap.add_argument("--out-dir", default=None,
                    help="results directory (default out/<tier>)")
    ap.add_argument("--purge-corpora", action="store_true",
                    help="delete this tier's corpus files at the end")
    ap.add_argument("--retry-failed", action="store_true",
                    help="re-run journaled failed cells; the new record "
                         "supersedes the failed one (METHODS.md section 6)")
    args = ap.parse_args(argv)

    types = [t.strip() for t in args.types.split(",") if t.strip()]
    for t in types:
        if t not in corpora.CONTENT_TYPES:
            ap.error(f"unknown type {t!r}")
    coder_names = ([c.strip() for c in args.coders.split(",") if c.strip()]
                   if args.coders else coders.available_coders(args.tier))
    for c in coder_names:
        if c not in coders.REGISTRY:
            ap.error(f"unknown coder {c!r}; known: "
                     f"{sorted(coders.REGISTRY)}")
        allowed = coders.REGISTRY[c]["tiers"]
        if allowed is not None and args.tier not in allowed:
            ap.error(f"coder {c!r} is not scheduled at tier {args.tier} "
                     f"(protocol feasibility rule; tiers={sorted(allowed)})")

    out_dir = Path(args.out_dir) if args.out_dir else HERE / "out" / args.tier
    counts = run(args.tier, types, coder_names, out_dir, args.cap_min,
                 args.runs, retry_failed=args.retry_failed)
    if args.purge_corpora:
        for t in types:
            corpora.corpus_path(args.tier, t).unlink(missing_ok=True)
    return 0 if counts.get("failed", 0) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

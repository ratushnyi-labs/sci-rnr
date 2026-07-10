#!/usr/bin/env python3
"""PASS/FAIL gate for the Type-I scale campaign machinery (protocol.md
section 7).  Everything runs at smoke scale; exit code 0 iff all checks
pass.

  CK1  protocol + dated METHODS amendment note present and complete
  CK2  every matrix coder available at smoke, versions readable
  CK3  corpus determinism: seeded generators rebuild byte-identical
       corpora; the images smoke corpus equals the pinned data/ payload
  CK4  full smoke matrix (every smoke coder x {white_noise, text} at
       16 MiB): all metrics present, round-trips verified, rate sanity
  CK5  checkpoint-resume: a campaign killed mid-cell resumes without
       re-running or duplicating finished cells
  CK6  figures + summary CSVs generate from the smoke JSONL and are
       byte-identical on regeneration (determinism)
  CK7  record schema: required keys, seeds match the METHODS.md
       derivation, corpus hash pinned in every record

Run:  /Users/para/.venvs/rnr/bin/python -u run_checks.py [--skip-slow]
  --skip-slow drops rnr1-ref from CK4 (dev loop only; the full gate
  must run without it).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
PY = sys.executable
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "bench"))

import coders  # noqa: E402
import corpora  # noqa: E402
import run_campaign  # noqa: E402
import metrics as _metrics  # noqa: E402

CHECKS: list = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def _sha_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            c = f.read(1 << 23)
            if not c:
                return h.hexdigest()
            h.update(c)


# ---------------------------------------------------------------------------

@check("CK1 protocol and amendment note")
def ck1(ctx):
    proto = HERE / "protocol.md"
    amend = HERE / "methods_amendment_2026-07-10.md"
    assert proto.exists(), "protocol.md missing"
    text = proto.read_text()
    for marker in ("pre-registered", "500 MB", "1500 MB", "5000 MB",
                   "white_noise", "precompressed", "90 min", "censored",
                   "cmix", "brotli -q11"):
        assert marker in text, f"protocol.md lacks required marker {marker!r}"
    assert amend.exists(), "dated METHODS amendment note missing"
    atext = amend.read_text()
    for marker in ("2026-07-10", "3 timed runs", "warmup", "METHODS"):
        assert marker in atext, f"amendment note lacks marker {marker!r}"
    return "protocol.md + dated amendment note complete"


@check("CK2 coder availability")
def ck2(ctx):
    expected = {"zstd-19", "zstd-22u", "xz-6", "xz-9e", "brotli-q11",
                "brotli-q9", "gzip-9", "bzip2-9", "rnr1-ref"}
    avail = set(coders.available_coders("smoke"))
    missing = expected - avail
    assert not missing, f"coders missing at smoke tier: {sorted(missing)}"
    vers = {}
    for name in sorted(expected):
        v = coders.coder_version(name)
        assert v and v != "unavailable", f"no version for {name}"
        vers[name] = v
    fast = sorted(avail - expected)
    note = f"; impl_fast lane present: {fast}" if fast else \
        "; impl_fast lane absent (rnr1-ref covers smoke)"
    return f"{len(expected)} coders available{note}"


@check("CK3 corpus determinism")
def ck3(ctx):
    # Materialize the full smoke set (idempotent).
    entries = {t: corpora.materialize("smoke", t)
               for t in corpora.CONTENT_TYPES}
    # Seeded generators must reproduce byte-identical output.
    for ctype in ("white_noise", "video"):
        with tempfile.NamedTemporaryFile(delete=False) as f:
            tmp = Path(f.name)
        try:
            with open(tmp, "wb") as f:
                corpora._GENERATORS[ctype](
                    f, corpora.TIER_BYTES["smoke"], "smoke")
            got = _sha_file(tmp)
            want = entries[ctype]["sha256"]
            assert got == want, \
                f"{ctype} regeneration drifted: {got[:16]} != {want[:16]}"
        finally:
            tmp.unlink(missing_ok=True)
    # images smoke corpus is the pinned data/ payload.
    with open(REPO / "data" / "MANIFEST.json") as f:
        data_manifest = json.load(f)
    want = data_manifest["entries"]["telemetry_grid"]["sha256"]
    assert entries["images"]["sha256"] == want, \
        "images smoke corpus does not match pinned telemetry_grid"
    return "white_noise/video regenerate identically; images == pinned payload"


@check("CK4 smoke matrix end-to-end")
def ck4(ctx):
    names = coders.available_coders("smoke")
    if ctx["skip_slow"]:
        names = [n for n in names if n != "rnr1-ref"]
    types = ["white_noise", "text"]
    counts = run_campaign.run(
        "smoke", types, names, HERE / "out" / "smoke",
        cap_min=run_campaign.DEFAULT_CAP_MIN, runs=None)
    assert counts["failed"] == 0, f"failed cells: {counts}"
    assert counts["censored"] == 0, f"unexpected censoring at smoke: {counts}"

    recs = [json.loads(l) for l in
            open(HERE / "out" / "smoke" / "results.jsonl") if l.strip()]
    need = {"bpb", "ratio", "compressed_size", "original_size", "enc_s",
            "dec_s", "enc_peak_rss", "dec_peak_rss", "roundtrip_ok"}
    seen = set()
    for r in recs:
        if r["ctype"] not in types or r["coder"] not in names:
            continue
        if r.get("status") != "done":
            continue
        m = r["measured"]
        missing = {k for k in need if m.get(k) is None}
        assert not missing, f"{r['cell']}: missing metrics {missing}"
        assert m["roundtrip_ok"] is True, f"{r['cell']}: round-trip failed"
        assert m["original_size"] == corpora.TIER_BYTES["smoke"]
        assert m["enc_s"] > 0 and m["dec_s"] > 0
        assert m["enc_peak_rss"] > 1e6, f"{r['cell']}: implausible RSS"
        if r["ctype"] == "white_noise":
            assert 7.99 < m["bpb"] < 8.6, \
                f"{r['cell']}: noise bpb {m['bpb']:.3f} out of range"
        if r["ctype"] == "text":
            assert m["bpb"] < 5.0, \
                f"{r['cell']}: text bpb {m['bpb']:.3f} implausibly high"
        seen.add((r["ctype"], r["coder"]))
    want = {(t, c) for t in types for c in names}
    missing = want - seen
    assert not missing, f"matrix cells without done records: {missing}"
    return (f"{len(want)} matrix cells done: metrics complete, "
            f"round-trips verified, rates sane")


@check("CK5 kill-and-resume checkpointing")
def ck5(ctx):
    with tempfile.TemporaryDirectory(prefix="t1s_resume_") as td:
        out_dir = Path(td) / "out"
        argv = [PY, "-u", str(HERE / "run_campaign.py"), "--tier", "smoke",
                "--types", "text", "--coders", "gzip-9,bzip2-9,xz-6",
                "--out-dir", str(out_dir)]
        journal = out_dir / "cells.journal"
        proc = subprocess.Popen(argv, stdout=subprocess.DEVNULL,
                                stderr=subprocess.STDOUT,
                                start_new_session=True)
        try:
            deadline = time.time() + 300
            while time.time() < deadline:
                if journal.exists() and journal.read_text().strip():
                    break
                time.sleep(0.2)
            else:
                raise AssertionError("no cell finished within 300 s")
            os.killpg(proc.pid, signal.SIGKILL)  # kill mid-campaign
        finally:
            proc.wait()
        lines_before = [l for l in journal.read_text().splitlines() if l]
        n_before = len(lines_before)
        assert 1 <= n_before < 3, \
            f"kill landed outside a partial state ({n_before} cells)"
        recs_before = [json.loads(l) for l in
                       open(out_dir / "results.jsonl") if l.strip()]
        assert len(recs_before) == n_before, "journal/results out of sync"

        cp = subprocess.run(argv, capture_output=True, text=True,
                            timeout=600)
        assert cp.returncode == 0, f"resume run failed: {cp.stdout[-400:]}"
        assert f"done={n_before}" in cp.stdout and "todo=" in cp.stdout, \
            "resume did not skip journaled cells"
        lines_after = [l for l in journal.read_text().splitlines() if l]
        cells_after = [json.loads(l)["cell"] for l in lines_after]
        assert len(cells_after) == 3, f"expected 3 cells, {cells_after}"
        assert len(set(cells_after)) == 3, "duplicate cells after resume"
        recs = [json.loads(l) for l in open(out_dir / "results.jsonl")
                if l.strip()]
        assert len(recs) == 3, "duplicate/missing records after resume"
        # Finished-before-kill cells were NOT re-run.
        assert [json.loads(l)["cell"] for l in lines_before] == \
            cells_after[:n_before], "resume re-ordered finished cells"
    return (f"killed after {n_before}/3 cells; resume completed the rest "
            f"without re-running or duplicating")


@check("CK6 figures from smoke data, deterministic")
def ck6(ctx):
    live_stores = sorted((HERE / "out").glob("*/results.jsonl"))
    assert live_stores, "no results stores for figures"
    with tempfile.TemporaryDirectory(prefix="t1s_figs_") as td:
        # Snapshot the stores so a concurrently appending campaign cannot
        # make the two generations diverge.
        stores = []
        for i, s in enumerate(live_stores):
            snap = Path(td) / f"store_{i}.jsonl"
            snap.write_bytes(s.read_bytes())
            stores.append(snap)
        outs = []
        for sub in ("a", "b"):
            fd = Path(td) / sub
            cp = subprocess.run(
                [PY, "-u", str(HERE / "make_figures.py"),
                 "--results", *[str(s) for s in stores],
                 "--out-dir", str(fd)],
                capture_output=True, text=True, timeout=600)
            assert cp.returncode == 0, f"make_figures failed: {cp.stderr[-400:]}"
            outs.append(fd)
        a, b = outs
        files = sorted(p.name for p in a.iterdir())
        for stem in ("fig_rate_vs_size", "fig_enc_time_vs_size",
                     "fig_dec_time_vs_size"):
            for ext in (".pdf", ".png"):
                assert stem + ext in files, f"missing {stem}{ext}"
        assert any(f.startswith("fig_pareto_") and f.endswith(".pdf")
                   for f in files), "missing Pareto figure"
        csvs = [f for f in files if f.endswith(".csv")]
        assert any("white_noise" in f for f in csvs) and \
            any("text" in f for f in csvs), f"summary CSVs missing: {csvs}"
        assert "captions.md" in files, "captions.md missing"
        cap = (a / "captions.md").read_text()
        for stem in ("fig_rate_vs_size", "fig_enc_time_vs_size",
                     "fig_dec_time_vs_size", "summary_"):
            assert stem in cap, f"caption missing for {stem}"
        # Determinism: byte-identical regeneration.
        for f in files:
            if (a / f).read_bytes() != (b / f).read_bytes():
                raise AssertionError(f"non-deterministic output: {f}")
        # Copy the validated outputs into the standing figures dir.
        final = HERE / "out" / "figures"
        final.mkdir(parents=True, exist_ok=True)
        for f in files:
            (final / f).write_bytes((a / f).read_bytes())
    return (f"{len(files)} artifacts (PDF+PNG+CSV+captions) generated, "
            f"byte-identical on regeneration; copied to out/figures/")


@check("CK7 record schema and seed derivation")
def ck7(ctx):
    recs = [json.loads(l) for l in
            open(HERE / "out" / "smoke" / "results.jsonl") if l.strip()]
    assert recs, "no smoke records"
    for r in recs:
        for k in ("corpus", "coder", "config", "seed", "timestamp",
                  "measured", "cell", "tier", "ctype", "corpus_sha256",
                  "coder_version", "rss_mode", "cap_s"):
            assert k in r, f"record {r.get('cell')} missing key {k}"
        want = _metrics.derive_seed(
            f"type1_scale:{r['tier']}:{r['ctype']}:{r['coder']}",
            r["run_idx"], tag="run")
        assert r["seed"] == want, f"{r['cell']}: seed mismatch"
        assert len(r["corpus_sha256"]) == 64
        assert r["rss_mode"] == "time-l-child"
    return f"{len(recs)} records schema-complete with derived seeds"


# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-slow", action="store_true",
                    help="drop rnr1-ref from CK4 (dev only)")
    args = ap.parse_args(argv)
    ctx = {"skip_slow": args.skip_slow}

    print("Type-I scale campaign machinery gate "
          f"(smoke tier{', --skip-slow' if args.skip_slow else ''})")
    print("=" * 66)
    ok = True
    for name, fn in CHECKS:
        t0 = time.perf_counter()
        try:
            detail = fn(ctx)
            print(f"PASS  {name}  [{time.perf_counter() - t0:.1f}s]")
            print(f"      {detail}")
        except Exception as exc:
            ok = False
            print(f"FAIL  {name}  [{time.perf_counter() - t0:.1f}s]")
            print(f"      {exc}")
    print("=" * 66)
    print(f"OVERALL: {'PASS' if ok else 'FAIL'}")
    if args.skip_slow and ok:
        print("note: --skip-slow gate is NOT sufficient for paper tiers; "
              "run the full gate.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

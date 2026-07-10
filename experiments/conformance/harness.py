#!/usr/bin/env python3
"""Instrumented conformance harness for the reference RNR Type-I coder.

Runs the impl/rnr1.py coding path over a spec'd set of test blocks and
emits a platform-comparable JSON report containing, per block:

  - archive SHA-256 + length (H8: byte-identical archives),
  - per-sub-block *inference-trace* digests (H14): a SHA-256 over the
    canonical little-endian serialization of every integer the coding
    path derives from the predictor at every position -- the full
    integer distribution freq[256] (the "logits" of the count-based
    predictor), total, argmax byte, candidate list, class map, exact
    class partition sums, coded class, and the payload coding interval,
  - a per-sub-block adaptive-state hash ladder (H9): after each
    sub-block is coded, the predictor's full adaptive count state
    (all context tables, canonically ordered) is hashed; the ladder is
    chained so a single hex value summarizes the whole state evolution,
  - decoder-side trace + state ladder (encoder/decoder adaptation
    agreement on MODE_CODED sub-blocks),
  - round-trip verification, max arithmetic-coder total observed
    (overflow-guard headroom vs the 2^24 bound, R-3.2/R-5.2).

The traced encoder mirrors impl/rnr1.py's encode loop using rnr1's own
primitives; on every run the full archive it assembles is cross-checked
byte-for-byte against rnr1.pack() output, so the verdicts apply to the
actual reference implementation, not to a divergent copy.

The report deliberately contains ONLY platform-invariant fields in
report["blocks"]; environment details live in report["env"].  Two
conforming platforms must produce byte-identical report["blocks"].

Determinism note: this harness itself follows the integer-only rules on
the traced path; SHA-256 and canonical struct/numpy little-endian
serialization are platform-independent.

Usage (native or inside a container with the repo mounted at any path):

  python -u harness.py --spec WORK/blocks/spec.json \
      --out WORK/native_a/report.json \
      --archives-dir WORK/native_a/archives \
      [--double-encode 2] \
      [--cross-decode-archives DIR --cross-decode-manifest REPORT.json]
"""

import argparse
import hashlib
import json
import os
import platform
import struct
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "impl"))

import numpy as np  # noqa: E402

import rnr1  # noqa: E402


# ---------------------------------------------------------------------------
# Canonical serialization helpers (all little-endian, fixed widths)
# ---------------------------------------------------------------------------

def _le_bytes(arr, dtype):
    """Canonical little-endian bytes of a numpy array (explicit dtype)."""
    return np.ascontiguousarray(arr, dtype=dtype).tobytes()


def state_hash(pred):
    """SHA-256 of the predictor's full adaptive state, canonically ordered.

    Canonical order: context tables by order 0..W, entries by integer
    context key ascending (independent of dict insertion / hash order),
    each entry as (key <Q>, total <q>, counts 256x <i8>); then the
    current history window.  This captures the entire integer adaptive
    state of Definition 10.1 / Theorem 10.2's compatible-online path.
    """
    h = hashlib.sha256()
    h.update(struct.pack("<i", pred.W))
    for order, table in enumerate(pred.tables):
        h.update(struct.pack("<iq", order, len(table)))
        for key in sorted(table):
            ent = table[key]
            h.update(struct.pack("<Qq", key, int(ent[1])))
            h.update(_le_bytes(ent[0], "<i8"))
    h.update(bytes(pred.hist))
    return h.hexdigest()


class LadderChain:
    """Chained SHA-256 ladder: h_j = SHA256(h_{j-1} || item_j)."""

    def __init__(self, label):
        self._cur = hashlib.sha256(label.encode("ascii")).digest()
        self.items = []

    def add(self, hexdigest):
        self.items.append(hexdigest)
        self._cur = hashlib.sha256(self._cur + bytes.fromhex(hexdigest)).digest()

    def hexdigest(self):
        return self._cur.hex()


# ---------------------------------------------------------------------------
# Traced encode / decode (mirrors impl/rnr1.py loops exactly)
# ---------------------------------------------------------------------------

def _trace_position(h, pos, freq, total, bstar, cand, classid, class_freq,
                    cls, pay):
    """Fold one position's full integer inference record into hash h.

    Record layout (little-endian): pos <q>, total <q>, freq 256x <i8>,
    bstar <h>, cand 15x <i2>, classid 256x <i1>, class_freq 5x <i8>,
    cls <b>, payload interval (lo, hi, total) 3x <q> (zeros for class 0).
    """
    # Integer-path dtype audit (R-3.1): the coding path must be int64.
    assert freq.dtype == np.int64, freq.dtype
    assert class_freq.dtype == np.int64, class_freq.dtype
    h.update(struct.pack("<qq", pos, int(total)))
    h.update(_le_bytes(freq, "<i8"))
    h.update(struct.pack("<h", bstar))
    h.update(_le_bytes(cand, "<i2"))
    h.update(_le_bytes(classid, "<i1"))
    h.update(_le_bytes(class_freq, "<i8"))
    h.update(struct.pack("<bqqq", cls, *pay))


def encode_subblock_traced(block, W):
    """Mirror of rnr1.encode_subblock with a full inference trace.

    Returns (blob, trace_hex, state_hex, max_total).
    """
    pred = rnr1.NGramPredictor(W)
    enc = rnr1.ArithmeticEncoder()
    h = hashlib.sha256()
    max_total = 0
    for pos, x in enumerate(block):
        freq, total = pred.dist()
        bstar, cand, classid, class_freq = rnr1._partition(freq)
        cls = int(classid[x])
        ccum = np.concatenate(([0], np.cumsum(class_freq)))
        enc.encode(int(ccum[cls]), int(ccum[cls + 1]), int(total))
        if cls != rnr1.CLS_0:
            members = rnr1._class_members(cls, cand, classid)
            mcum = np.concatenate(([0], np.cumsum(freq[members])))
            idx = int(np.nonzero(members == x)[0][0])
            pay = (int(mcum[idx]), int(mcum[idx + 1]), int(mcum[-1]))
            enc.encode(*pay)
        else:
            pay = (0, 0, 0)
        _trace_position(h, pos, freq, total, bstar, cand, classid,
                        class_freq, cls, pay)
        max_total = max(max_total, int(total), pay[2])
        pred.update(x)
    return enc.finish(), h.hexdigest(), state_hash(pred), max_total


def decode_subblock_traced(blob, nbytes, W):
    """Mirror of rnr1.decode_subblock with the same trace record layout.

    Returns (data, trace_hex, state_hex).  On a conforming platform the
    trace equals the encoder-side trace of the same sub-block.
    """
    pred = rnr1.NGramPredictor(W)
    dec = rnr1.ArithmeticDecoder(blob)
    h = hashlib.sha256()
    out = bytearray(nbytes)
    for pos in range(nbytes):
        freq, total = pred.dist()
        bstar, cand, classid, class_freq = rnr1._partition(freq)
        ccum = np.concatenate(([0], np.cumsum(class_freq)))
        t = dec.target(int(total))
        cls = int(np.searchsorted(ccum, t, side="right")) - 1
        dec.consume(int(ccum[cls]), int(ccum[cls + 1]), int(total))
        if cls == rnr1.CLS_0:
            x = bstar
            pay = (0, 0, 0)
        else:
            members = rnr1._class_members(cls, cand, classid)
            mcum = np.concatenate(([0], np.cumsum(freq[members])))
            t = dec.target(int(mcum[-1]))
            idx = int(np.searchsorted(mcum, t, side="right")) - 1
            pay = (int(mcum[idx]), int(mcum[idx + 1]), int(mcum[-1]))
            dec.consume(*pay)
            x = int(members[idx])
        _trace_position(h, pos, freq, total, bstar, cand, classid,
                        class_freq, cls, pay)
        out[pos] = x
        pred.update(x)
    return bytes(out), h.hexdigest(), state_hash(pred)


def pack_traced(data, W, K):
    """Mirror of rnr1.pack that additionally returns the traces.

    Returns (archive, info) where info holds per-sub-block mode,
    inference-trace digests, state ladder, and max_total.  The archive
    is cross-checked byte-for-byte against rnr1.pack by the caller.
    """
    n = len(data)
    m = (n + K - 1) // K
    blobs, entries = [], []
    modes, traces = [], []
    ladder = LadderChain("rnr1-conformance-state-ladder-v1")
    trace_chain = LadderChain("rnr1-conformance-inference-trace-v1")
    repair_bit_offset = 0
    max_total = 0
    for j in range(m):
        block = data[j * K: (j + 1) * K]
        blob, trace_hex, state_hex, mt = encode_subblock_traced(block, W)
        max_total = max(max_total, mt)
        mode = rnr1.MODE_CODED
        if len(blob) >= len(block):
            blob = block
            mode = rnr1.MODE_RAW
        entries.append(struct.pack(
            rnr1.INDEX_ENTRY_FMT, j * K, repair_bit_offset, 0, mode,
            rnr1.h_v(block)))
        blobs.append(blob)
        repair_bit_offset += 8 * len(blob)
        modes.append(mode)
        traces.append(trace_hex)
        ladder.add(state_hex)
        trace_chain.add(trace_hex)
    header = struct.pack(
        rnr1.HEADER_FMT, rnr1.MAGIC, bytes(rnr1.VERSION), 0,
        rnr1.FLAG_SUBBLOCK_HASHES, W, K, n,
        rnr1.model_description_hash(W), rnr1.h_v(data), m)
    archive = header + b"".join(entries) + b"".join(blobs)
    info = {
        "modes": modes,
        "enc_traces": traces,
        "enc_trace_chain": trace_chain.hexdigest(),
        "enc_state_ladder": ladder.items,
        "enc_state_ladder_chain": ladder.hexdigest(),
        "max_total": max_total,
    }
    return archive, info


def unpack_traced(raw):
    """Traced full decode of an RNR1 archive (verifying, like unpack).

    Returns (data, dec_traces, dec_states): per-sub-block trace/state
    digests for MODE_CODED sub-blocks (None for MODE_RAW: the raw path
    never runs the predictor, per impl/rnr1.py decode_block).
    """
    arc = rnr1.Archive(raw)
    pieces, dec_traces, dec_states = [], [], []
    for j in range(arc.m):
        blk_len = arc._block_len(j)
        if arc.entries[j][3] == rnr1.MODE_RAW:
            out = arc._block_blob(j)[:blk_len]
            dec_traces.append(None)
            dec_states.append(None)
        else:
            out, trace_hex, state_hex = decode_subblock_traced(
                arc._block_blob(j), blk_len, arc.W)
            dec_traces.append(trace_hex)
            dec_states.append(state_hex)
        if rnr1.h_v(out) != arc.entries[j][4]:
            raise ValueError("sub-block %d hash mismatch" % j)
        pieces.append(out)
    data = b"".join(pieces)
    if rnr1.h_v(data) != arc.v:
        raise ValueError("archive verification failed: H_v mismatch")
    return data, dec_traces, dec_states


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------

def env_info():
    return {
        "machine": platform.machine(),
        "system": platform.system(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "byteorder": sys.byteorder,
        "pythonhashseed": os.environ.get("PYTHONHASHSEED", "<unset>"),
    }


def process_block(name, data, W, K, archives_dir, double_encode):
    t0 = time.perf_counter()
    archive, info = pack_traced(data, W, K)

    # Cross-check: the traced mirror must reproduce impl/rnr1.py exactly.
    impl_archive = rnr1.pack(data, W=W, K=K)
    if archive != impl_archive:
        raise AssertionError(
            "harness drift: traced archive != rnr1.pack archive for %r" % name)

    # Optional in-process repeatability probe (same-process double encode).
    double_ok = None
    if double_encode:
        double_ok = rnr1.pack(data, W=W, K=K) == impl_archive

    # Traced decode + round trip.
    decoded, dec_traces, dec_states = unpack_traced(archive)
    roundtrip_ok = decoded == data

    # Encoder/decoder adaptation agreement on coded sub-blocks (H9/H14).
    enc_dec_trace_ok = True
    enc_dec_state_ok = True
    n_coded = 0
    for j, mode in enumerate(info["modes"]):
        if mode != rnr1.MODE_CODED:
            continue
        n_coded += 1
        if dec_traces[j] != info["enc_traces"][j]:
            enc_dec_trace_ok = False
        if dec_states[j] != info["enc_state_ladder"][j]:
            enc_dec_state_ok = False

    if archives_dir:
        with open(os.path.join(archives_dir, name + ".rnr1"), "wb") as f:
            f.write(archive)

    block = {
        "name": name,
        "n": len(data),
        "W": W,
        "K": K,
        "input_sha256": hashlib.sha256(data).hexdigest(),
        "archive_sha256": hashlib.sha256(archive).hexdigest(),
        "archive_len": len(archive),
        "n_subblocks": len(info["modes"]),
        "n_raw_subblocks": sum(1 for m in info["modes"] if m == rnr1.MODE_RAW),
        "enc_traces": info["enc_traces"],
        "enc_trace_chain": info["enc_trace_chain"],
        "enc_state_ladder": info["enc_state_ladder"],
        "enc_state_ladder_chain": info["enc_state_ladder_chain"],
        "max_total": info["max_total"],
        "roundtrip_ok": roundtrip_ok,
        "enc_dec_trace_ok": enc_dec_trace_ok,
        "enc_dec_state_ok": enc_dec_state_ok,
        "n_coded_subblocks": n_coded,
    }
    if double_ok is not None:
        block["double_encode_ok"] = double_ok
    return block, time.perf_counter() - t0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--spec", default=None,
                    help="JSON spec: [{name,file,W,K,sha256}, ...]; file "
                         "paths are relative to the spec's directory; "
                         "omit for a decode-only (cross-decode) run")
    ap.add_argument("--out", required=True, help="report JSON path")
    ap.add_argument("--archives-dir", default=None,
                    help="write per-block archives here (for cross-decode)")
    ap.add_argument("--double-encode", type=int, default=0,
                    help="re-encode the first N blocks in-process and "
                         "compare (same-process repeatability probe)")
    ap.add_argument("--cross-decode-archives", default=None,
                    help="directory of archives produced elsewhere")
    ap.add_argument("--cross-decode-manifest", default=None,
                    help="the other platform's report.json (for expected "
                         "input hashes)")
    args = ap.parse_args(argv)

    if not args.spec and not (args.cross_decode_archives
                              and args.cross_decode_manifest):
        ap.error("need --spec and/or --cross-decode-archives+manifest")
    spec, spec_dir = [], "."
    if args.spec:
        with open(args.spec, "r", encoding="utf-8") as f:
            spec = json.load(f)
        spec_dir = os.path.dirname(os.path.abspath(args.spec))
    if args.archives_dir:
        os.makedirs(args.archives_dir, exist_ok=True)

    blocks = []
    total_s = 0.0
    for i, item in enumerate(spec):
        with open(os.path.join(spec_dir, item["file"]), "rb") as f:
            data = f.read()
        got = hashlib.sha256(data).hexdigest()
        if got != item["sha256"]:
            raise AssertionError("spec sha mismatch for %s" % item["name"])
        blk, dt = process_block(
            item["name"], data, item["W"], item["K"],
            args.archives_dir, double_encode=(i < args.double_encode))
        total_s += dt
        blocks.append(blk)
        print("block %-24s n=%-6d m=%-3d raw=%-2d %.2fs  arc=%s" % (
            item["name"], blk["n"], blk["n_subblocks"],
            blk["n_raw_subblocks"], dt, blk["archive_sha256"][:16]),
            flush=True)

    # Global chains for a one-line platform comparison.
    g_arc = LadderChain("rnr1-conformance-archives-v1")
    g_trace = LadderChain("rnr1-conformance-traces-v1")
    g_state = LadderChain("rnr1-conformance-states-v1")
    for blk in blocks:
        g_arc.add(blk["archive_sha256"])
        g_trace.add(blk["enc_trace_chain"])
        g_state.add(blk["enc_state_ladder_chain"])

    report = {
        "format": "rnr1-conformance-report-v1",
        "env": env_info(),
        "wall_s": round(total_s, 3),
        "blocks": blocks,
        "global": {
            "n_blocks": len(blocks),
            "archives_chain": g_arc.hexdigest(),
            "traces_chain": g_trace.hexdigest(),
            "states_chain": g_state.hexdigest(),
            "max_total": max([b["max_total"] for b in blocks] + [0]),
            "prob_total_bound": rnr1.PROB_TOTAL_BOUND,
        },
    }

    # Cross-decode leg: decode archives produced on another platform.
    if args.cross_decode_archives and args.cross_decode_manifest:
        with open(args.cross_decode_manifest, "r", encoding="utf-8") as f:
            other = json.load(f)
        xdec = []
        for oblk in other["blocks"]:
            path = os.path.join(args.cross_decode_archives,
                                oblk["name"] + ".rnr1")
            with open(path, "rb") as f:
                raw = f.read()
            out = rnr1.unpack(raw)  # verifying decode
            ok = hashlib.sha256(out).hexdigest() == oblk["input_sha256"]
            xdec.append({"name": oblk["name"], "ok": bool(ok)})
            print("cross-decode %-24s %s" % (oblk["name"],
                                             "ok" if ok else "MISMATCH"),
                  flush=True)
        report["cross_decode"] = {
            "n": len(xdec),
            "n_ok": sum(1 for x in xdec if x["ok"]),
            "results": xdec,
        }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1, sort_keys=True)
        f.write("\n")
    print("report: %s  (%d blocks, %.1fs)  env=%s" % (
        args.out, len(blocks), total_s, report["env"]), flush=True)


if __name__ == "__main__":
    main()

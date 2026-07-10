#!/usr/bin/env python3
"""Joint vs factorized coding arms with the reference Type-I coder.

Both arms use impl/rnr1.py's own NGramPredictor at MATCHED order W and
per-block cold starts (the sub-block/sync semantics of the container:
predictor reset to Q_init at each block boundary).

Bit accounting is by exact replay: the coder's class/payload chain costs
exactly -log2(freq[x]/total) per byte (chain rule over the exact integer
partition, see the rnr1.py module docstring), so accumulating that
quantity gives the code CONTENT length up to arithmetic-coder
quantization. validate_replay() cross-checks the replay total against a
real pack() archive's repair-stream size on the same bytes.

Arms
----
joint:       records concatenated in record order; the order-W context
             sees the preceding fields of the same record (and the tail
             of the previous record), i.e. joint conditioning.
factorized:  one stream per field (column), context = the same field's
             history across records only -- for exchangeable records this
             converges to the per-position marginal coder of
             Theorem 5.5(F). Per-block bits are summed over fields.

Per-block accounting on the SAME record partition in both arms yields
the paired per-block differences required by METHODS.md section 2.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "impl"))

import rnr1  # noqa: E402

LOG2 = (np.log2(np.arange(1, 1 << 20))).astype(np.float64)  # small cache


def _log2i(n: int) -> float:
    if n < LOG2.shape[0] + 1:
        return float(LOG2[n - 1])
    return math.log2(n)


def stream_bits(data: bytes | np.ndarray, W: int) -> float:
    """Exact replay code-content bits of one cold-started stream."""
    pred = rnr1.NGramPredictor(W)
    dist = pred.dist
    update = pred.update
    bits = 0.0
    if isinstance(data, np.ndarray):
        data = data.tobytes()
    for x in data:
        freq, total = dist()
        bits += _log2i(int(total)) - _log2i(int(freq[x]))
        update(x)
    return bits


def blockwise_bits(streams: list[np.ndarray], n_blocks: int,
                   W: int) -> np.ndarray:
    """Per-block replay bits, cold start per block, summed over streams.

    Each stream is an (M_records, L) uint8 array (L bytes contributed per
    record); the record axis is partitioned into n_blocks contiguous
    blocks and every stream restarts its predictor at each block start.
    Returns bits[n_blocks].
    """
    M = streams[0].shape[0]
    edges = np.linspace(0, M, n_blocks + 1).astype(int)
    out = np.zeros(n_blocks)
    for s in streams:
        for b in range(n_blocks):
            seg = np.ascontiguousarray(s[edges[b]:edges[b + 1]])
            out[b] += stream_bits(seg, W)
    return out


def block_sizes_bytes(M: int, R: int, n_blocks: int) -> np.ndarray:
    edges = np.linspace(0, M, n_blocks + 1).astype(int)
    return np.diff(edges) * R


def joint_streams(records: np.ndarray) -> list[np.ndarray]:
    """Single stream: records concatenated in order (M, R)."""
    return [records]


def factorized_streams(records: np.ndarray) -> list[np.ndarray]:
    """One stream per field: column f across records, shape (M, 1)."""
    return [records[:, f:f + 1] for f in range(records.shape[1])]


# ---------------------------------------------------------------------------
# Cross-validation of the replay accounting against real archives
# ---------------------------------------------------------------------------

def pack_content_bits(data: bytes, W: int, K: int) -> tuple[float, int, int]:
    """(repair-stream bits, archive bytes, n raw-mode blocks) of pack()."""
    raw = rnr1.pack(data, W=W, K=K)
    arc = rnr1.Archive(raw)
    n_raw = sum(1 for e in arc.entries if e[3] == rnr1.MODE_RAW)
    return 8.0 * len(arc.repair), len(raw), n_raw


def validate_replay(records: np.ndarray, W: int, cap_bytes: int = 98304,
                    ) -> dict:
    """Compare replay bits vs a real single-sub-block archive.

    Uses the joint stream truncated to cap_bytes, K = stream length
    (single sub-block, single cold start on both sides). Round-trips the
    archive as well. Only meaningful on compressible streams (no
    store-mode escape), which is checked and reported.
    """
    data = np.ascontiguousarray(records).tobytes()[:cap_bytes]
    K = max(256, len(data))
    ideal = stream_bits(np.frombuffer(data, dtype=np.uint8).reshape(1, -1), W)
    content, arc_bytes, n_raw = pack_content_bits(data, W, K)
    decoded = rnr1.unpack(rnr1.pack(data, W=W, K=K))
    return {
        "bytes": len(data),
        "replay_bits": ideal,
        "pack_content_bits": content,
        "archive_bytes": arc_bytes,
        "raw_mode_blocks": n_raw,
        "roundtrip_ok": decoded == data,
        "rel_err": abs(content - ideal) / max(ideal, 1.0),
    }

#!/usr/bin/env python3
"""H5 -- Type-III-A mode selection on heterogeneous archives
(exp-design section 2.4; main paper Theorems 6.1/6.2).

Harness: a seeded heterogeneous composite (text + sqlite + telemetry +
urandom segments, per the preconditions qualification table's
heterogeneity findings) is tiled; every tile is coded standalone by
every mode in common.MODES, giving the exact per-tile cost matrix
L_m(t) in bits.  Two selectors are evaluated against the per-tile
oracle min_m L_m(t) and against every fixed mode:

  prefix selector  U_m(t) = L_m(prefix of t) * |t|/|prefix|   -- the
      selector is trained on a prefix sample of each tile (the tile is
      the scheduling segment unit here; tiles never straddle segment
      boundaries).  Encoder-side selection, so an explicit mode flag of
      ceil(log2 |M|) bits per tile is charged to the selector.
  causal selector  U_m(t) = L_m(t-1) * |t|/|t-1| for tiles with a
      predecessor in the SAME segment (decoder-reproducible state, no
      flag); first tile of each segment uses a fixed declared default
      mode (zstd-22) and contributes no epsilon sample.

Epsilon instrumentation (documented per the task): because the full
cost matrix L_m(t) is measured for every mode, the per-tile prediction
error eps_t = max_m |U_m(t) - L_m(t)| is computed EXACTLY, not
estimated; the Theorem 6.1 bound L_sel(t) <= L_oracle(t) + 2 eps_t is
then checked tile by tile (it must hold identically -- a violation
means broken instrumentation), and the empirical content reported is
(a) the size of eps under each training rule and (b) the realized
regret vs the 5%/20% H5 thresholds.
"""

from __future__ import annotations

import math

import numpy as np

import common
from common import MODES, MODE_NAMES, MODE_FLAG_BITS, loader

SOURCES = ["enwik8", "sqlite_synth", "telemetry_grid", "ctrl_urandom"]
CAUSAL_DEFAULT = "zstd-22"
PREFIX_FRAC = 4          # prefix = tile_bytes / PREFIX_FRAC


def build_composition(seed: int, total_tiles: int, tile_bytes: int):
    """Seeded heterogeneous composite: segments of 1-3 tiles drawn from
    the four qualification-table sources at tile-aligned offsets; the
    first four segments cycle the sources so each class appears."""
    rng = np.random.default_rng(seed)
    segs = []
    tiles_left = total_tiles
    i = 0
    while tiles_left > 0:
        src = SOURCES[i % 4] if i < 4 else SOURCES[int(rng.integers(0, 4))]
        ln = min(int(rng.integers(1, 4)), tiles_left)
        nbytes = ln * tile_bytes
        max_tile_off = (loader.corpus_bytes(src) - nbytes) // tile_bytes
        off = int(rng.integers(0, max_tile_off)) * tile_bytes
        segs.append({"src": src, "offset": off, "tiles": ln})
        tiles_left -= ln
        i += 1
    tiles, tile_src, tile_pos_in_seg = [], [], []
    for s in segs:
        raw = loader.read_range(s["src"], s["offset"], s["tiles"] * tile_bytes)
        for j in range(s["tiles"]):
            tiles.append(raw[j * tile_bytes:(j + 1) * tile_bytes])
            tile_src.append(s["src"])
            tile_pos_in_seg.append(j)
    return tiles, tile_src, tile_pos_in_seg, segs


def cost_matrix(blocks) -> np.ndarray:
    """L[m, t] = standalone bit cost of block t under mode m."""
    L = np.zeros((len(MODE_NAMES), len(blocks)), dtype=np.int64)
    for mi, m in enumerate(MODE_NAMES):
        fn = MODES[m]
        for ti, blk in enumerate(blocks):
            L[mi, ti] = fn(blk)
    return L


def run_h5(seed_index: int, total_tiles: int, tile_bytes: int) -> dict:
    seed = common.derive_seed("h5-comp", seed_index)
    tiles, tile_src, pos_in_seg, segs = build_composition(
        seed, total_tiles, tile_bytes)
    T = len(tiles)
    prefix_bytes = tile_bytes // PREFIX_FRAC
    prefixes = [t[:prefix_bytes] for t in tiles]

    L = cost_matrix(tiles)                       # (M, T) bits
    Lp = cost_matrix(prefixes)                   # (M, T) bits on prefixes

    # --- prefix selector -------------------------------------------------
    U_prefix = Lp.astype(float) * (tile_bytes / prefix_bytes)
    sel_p = np.argmin(U_prefix, axis=0)
    eps_p = np.max(np.abs(U_prefix - L), axis=0)             # exact eps_t
    bits_sel_p = int(L[sel_p, np.arange(T)].sum()) + MODE_FLAG_BITS * T

    # --- causal selector --------------------------------------------------
    default_mi = MODE_NAMES.index(CAUSAL_DEFAULT)
    sel_c = np.empty(T, dtype=int)
    eps_c = np.full(T, np.nan)
    tile_len = np.array([len(t) for t in tiles], dtype=float)
    for t in range(T):
        if pos_in_seg[t] == 0:
            sel_c[t] = default_mi
        else:
            U = L[:, t - 1].astype(float) * (tile_len[t] / tile_len[t - 1])
            sel_c[t] = int(np.argmin(U))
            eps_c[t] = float(np.max(np.abs(U - L[:, t])))
    bits_sel_c = int(L[sel_c, np.arange(T)].sum())           # no flag bits

    # --- oracle and fixed modes -------------------------------------------
    orc = L.min(axis=0)
    bits_orc = int(orc.sum())
    n_bytes = int(sum(len(t) for t in tiles))
    fixed = {m: int(L[mi].sum()) for mi, m in enumerate(MODE_NAMES)}

    # --- Theorem 6.1 bound, tile by tile ----------------------------------
    lhs_p = L[sel_p, np.arange(T)].astype(float)
    bound_ok_p = bool(np.all(lhs_p <= orc + 2 * eps_p + 1e-9))
    ok_c = True
    for t in range(T):
        if pos_in_seg[t] > 0:
            if L[sel_c[t], t] > orc[t] + 2 * eps_c[t] + 1e-9:
                ok_c = False

    # --- regret ratios with paired-tile bootstrap CIs ----------------------
    diff_p = lhs_p + MODE_FLAG_BITS - orc                    # per-tile
    diff_c = L[sel_c, np.arange(T)].astype(float) - orc
    boot_seed = common.derive_seed("h5-boot", seed_index)
    reg_p = common.ratio_bootstrap(diff_p, orc.astype(float), boot_seed)
    reg_c = common.ratio_bootstrap(diff_c, orc.astype(float), boot_seed + 1)

    # --- per-source oracle-mode table --------------------------------------
    per_source = {}
    for src in SOURCES:
        idx = [t for t in range(T) if tile_src[t] == src]
        if not idx:
            continue
        wins = np.bincount(np.argmin(L[:, idx], axis=0),
                           minlength=len(MODE_NAMES))
        per_source[src] = {
            "tiles": len(idx),
            "oracle_mode_counts": {MODE_NAMES[i]: int(w)
                                   for i, w in enumerate(wins) if w},
            "oracle_bpb": float(orc[idx].sum() / (len(idx) * tile_bytes)),
        }

    eps_c_vals = eps_c[~np.isnan(eps_c)]
    return {
        "seed": int(seed), "seed_index": seed_index,
        "tile_bytes": tile_bytes, "tiles": T, "n_bytes": n_bytes,
        "segments": segs, "modes": MODE_NAMES,
        "mode_flag_bits": MODE_FLAG_BITS,
        "prefix_bytes": prefix_bytes,
        "cost_matrix_bits": L.tolist(),
        "tile_src": tile_src,
        "oracle_bits": bits_orc,
        "oracle_bpb": bits_orc / n_bytes,
        "fixed_mode_bits": fixed,
        "fixed_mode_bpb": {m: v / n_bytes for m, v in fixed.items()},
        "selector_prefix": {
            "bits": bits_sel_p, "bpb": bits_sel_p / n_bytes,
            "choices": [MODE_NAMES[i] for i in sel_p],
            "eps_mean_bpb": float(eps_p.mean() / tile_bytes),
            "eps_p95_bpb": float(np.percentile(eps_p, 95) / tile_bytes),
            "bound_holds_all_tiles": bound_ok_p,
            "regret_ratio": reg_p,
            "regret_le_2eps_aggregate":
                bool(diff_p.sum() <= 2 * eps_p.sum() + MODE_FLAG_BITS * T),
        },
        "selector_causal": {
            "bits": bits_sel_c, "bpb": bits_sel_c / n_bytes,
            "choices": [MODE_NAMES[i] for i in sel_c],
            "eps_mean_bpb": float(eps_c_vals.mean() / tile_bytes)
                if eps_c_vals.size else None,
            "eps_p95_bpb": float(np.percentile(eps_c_vals, 95) / tile_bytes)
                if eps_c_vals.size else None,
            "bound_holds_all_tiles": ok_c,
            "regret_ratio": reg_c,
            "default_mode": CAUSAL_DEFAULT,
            "eps_tiles": int(eps_c_vals.size),
        },
        "roundtrip_sampled_ok": bool(
            common.rnr_roundtrip_ok(tiles[0])
            and common.rnr_roundtrip_ok(tiles[-1], cm=True)),
    }

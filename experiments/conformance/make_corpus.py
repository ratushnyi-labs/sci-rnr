#!/usr/bin/env python3
"""Deterministic conformance test-corpus builder (exp-design section 8.1).

Materializes the test blocks as files under WORKDIR/blocks/ plus a
spec.json the harness consumes.  Fully deterministic: block offsets are
derived from the frozen master seed via bench/metrics.derive_seed (tag
"corpus"), so the block set is reproducible from the manifest'd corpora
alone and its SHA-256 hashes are published in the spec.

Full mode (exp-design section 8.1): 100 blocks of 1-16 KiB spanning all
nine materialized corpus classes (>= 10 per class), plus synthetic edge
blocks and the overflow-stress block.  Fast mode: 1 block per corpus
plus the same edge/stress blocks.

Runs natively only (needs data/payloads); the container legs receive
the materialized block files, never the raw corpora.
"""

import argparse
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "data"))
sys.path.insert(0, os.path.join(REPO, "bench"))

import loader  # noqa: E402
import metrics  # noqa: E402

# Cycles are fixed here (frozen with the suite): sizes cover 1-16 KiB,
# K always >= 256 (format constraint) and < size for most blocks so the
# H9 state ladder has multiple rungs; W covers the [1, 8] envelope ends
# plus the default.
SIZES = [1024, 2048, 3072, 4096, 5120, 6144, 8192, 10240, 12288, 14336, 16384]
KS = [512, 1024, 2048, 512, 1024, 2048, 4096, 2048, 4096, 2048, 4096]
WS = [3, 3, 2, 3, 4, 3, 1, 3, 3, 5, 3]

EXPERIMENT_ID = "conformance-h8-h9-h14"


def corpus_blocks(per_corpus):
    """Deterministic (name, corpus, offset, size, W, K) tuples."""
    items = []
    for cname in sorted(loader.corpus_names()):
        nbytes = loader.corpus_bytes(cname)
        for i in range(per_corpus):
            size = SIZES[i % len(SIZES)]
            K = KS[i % len(KS)]
            W = WS[i % len(WS)]
            seed = metrics.derive_seed(
                "%s:%s" % (EXPERIMENT_ID, cname), i, tag="corpus")
            offset = seed % max(1, nbytes - size)
            items.append(("%s_%02d" % (cname, i), cname, offset, size, W, K))
    return items


def synthetic_blocks():
    """Edge-case blocks (deterministic content, no RNG needed)."""
    per = (b"AB" * 2048)
    ramp = bytes(range(256)) * 9  # 2304 bytes, crosses K=2048 boundary
    return [
        # name, data, W, K
        ("edge_empty", b"", 3, 4096),
        ("edge_one_byte", b"\x00", 3, 256),
        ("edge_zeros_4k", b"\x00" * 4096, 3, 1024),
        ("edge_periodic_4k", per, 2, 1024),
        ("edge_exact_k", b"x" * 2048 + b"", 3, 2048),      # n == K
        ("edge_k_plus_1", ramp[:2049], 3, 2048),           # 2nd sub-block = 1B
        ("edge_min_k", bytes(ramp[:700]), 3, 256),         # minimum legal K
        # Overflow-stress: > COUNT_CAP repeats of one byte in ONE sub-block
        # drives every context total to the cap; the maximum total seen by
        # dist() must be exactly 256 + 32*(COUNT_CAP-1) = 2097376, i.e. the
        # documented 8x headroom below the 2^24 arithmetic-coder bound.
        ("stress_cap_repeat", b"\xaa" * (66 * 1024), 3, 128 * 1024),
    ]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--full", action="store_true",
                    help="100-block campaign corpus (exp-design 8.1); "
                         "default is the fast reduced-scale corpus")
    args = ap.parse_args(argv)

    blocks_dir = os.path.join(args.workdir, "blocks")
    os.makedirs(blocks_dir, exist_ok=True)

    per_corpus = 11 if args.full else 1
    items = corpus_blocks(per_corpus)
    if args.full:
        # 9 corpora x 11 = 99; trim to exactly 100 with the stress block
        # counted in, keeping >= 10 per class: drop the 8 largest-index
        # extras evenly is NOT needed -- 99 + 8 synthetic = 107 blocks
        # total, of which 100 are the exp-design 8.1 census (99 corpus
        # blocks + stress) and 7 are boundary edges kept as extras.
        pass

    spec = []
    for (name, cname, offset, size, W, K) in items:
        data = loader.read_range(cname, offset, size)
        assert len(data) == size
        _write(blocks_dir, spec, name, data, W, K,
               source="%s@%d" % (cname, offset))
    for (name, data, W, K) in synthetic_blocks():
        _write(blocks_dir, spec, name, data, W, K, source="synthetic")

    spec_path = os.path.join(blocks_dir, "spec.json")
    with open(spec_path, "w", encoding="utf-8") as f:
        json.dump(spec, f, indent=1, sort_keys=True)
        f.write("\n")
    total = sum(s["n"] for s in spec)
    print("corpus: %d blocks, %d bytes total -> %s"
          % (len(spec), total, spec_path))


def _write(blocks_dir, spec, name, data, W, K, source):
    fname = name + ".bin"
    with open(os.path.join(blocks_dir, fname), "wb") as f:
        f.write(data)
    spec.append({
        "name": name,
        "file": fname,
        "n": len(data),
        "W": W,
        "K": K,
        "sha256": hashlib.sha256(data).hexdigest(),
        "source": source,
    })


if __name__ == "__main__":
    main()

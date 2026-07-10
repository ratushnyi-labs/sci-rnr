#!/usr/bin/env python3
"""Batch-produce reference archives for the whole non-deferred corpus at
W in {2,3}, K=64KiB, using ref_driver.parallel_pack (normative rnr1 code
paths, process-pool over independent sub-blocks).  Writes to the directory
given as argv[1].  Logs one line per archive with timing."""

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "data"))
from ref_driver import parallel_pack  # noqa: E402
import loader  # noqa: E402

OUT = sys.argv[1]
K = 65536


def main():
    os.makedirs(OUT, exist_ok=True)
    names = loader.corpus_names()  # non-deferred entries only
    # Largest first so the tail of the run is short files.
    names.sort(key=lambda n: -loader.corpus_bytes(n))
    for W in (2, 3):
        for name in names:
            dst = os.path.join(OUT, "%s.W%d.rnr1" % (name, W))
            if os.path.exists(dst):
                print("skip %s" % dst, flush=True)
                continue
            path = loader.payload_path(name)
            with open(path, "rb") as f:
                data = f.read()
            t0 = time.time()
            raw = parallel_pack(data, W, K, jobs=8)
            dt = time.time() - t0
            with open(dst + ".tmp", "wb") as f:
                f.write(raw)
            os.rename(dst + ".tmp", dst)
            print("ref-pack %-14s W=%d: %9d -> %9d bytes  %7.1fs  %.1f KB/s"
                  % (name, W, len(data), len(raw), dt,
                     len(data) / 1024.0 / dt), flush=True)
    print("REF BATCH DONE", flush=True)


if __name__ == "__main__":
    main()

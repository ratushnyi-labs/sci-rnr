# impl_fast — fast C Type-I engine (bit-identical to impl/rnr1.py)

C port of the normative reference coder `impl/rnr1.py` (archive format:
`impl/FORMAT.md`).  The Python reference stays normative; this engine is
accepted only through byte-identity gates (see `run_checks.py` and
`NOTES.md`).

## Files

| file | purpose |
|------|---------|
| `rnr1_fast.c` | single-file C implementation (predictor, arithmetic coder, container, store mode, SHA-256 via CommonCrypto) with optional multi-threaded pack/unpack over independent sub-blocks |
| `rnr1fast.py` | thin ctypes wrapper: `pack(data, W, K, threads)` / `unpack(raw, verify, threads)` + CLI |
| `build.sh` | `clang -O3` build of `librnr1fast.dylib` |
| `run_checks.py` | PASS/FAIL acceptance gate (build, bit-identity on smoke classes + full corpus, cross-decode both directions, error paths, determinism, throughput) |
| `ref_driver.py` | parallel driver over the normative reference (process pool over independent sub-blocks; used by the gate to make corpus-scale reference runs feasible) |
| `ref_batch.py` | batch producer of cached reference archives for the gate |
| `NOTES.md` | engineering notes, measured throughput/memory, honest caveats, dated methodology amendment |

## Quick start

```
sh impl_fast/build.sh
/Users/para/.venvs/rnr/bin/python impl_fast/rnr1fast.py pack   IN OUT.rnr --W 3 --K 65536 --threads 8
/Users/para/.venvs/rnr/bin/python impl_fast/rnr1fast.py unpack OUT.rnr RESTORED
```

## Acceptance gate

```
/Users/para/.venvs/rnr/bin/python impl_fast/run_checks.py [--refdir DIR] [--skip-corpus]
```

All checks print PASS/FAIL and an OVERALL line.  The corpus matrix
(F3/F4) regenerates cached reference archives on a cold run (~25 min on
8 cores); `--skip-corpus` runs the fast subset.

Measured on this host (enwik8, W=3, K=64 KiB, gate run 2026-07-10):
9.07 MB/s encode single thread, 45.4 MB/s at 8 threads (byte-identical
output); decode 9.9 / 50.5 MB/s.  Full gate: OVERALL PASS 8/8
(`run_checks_2026-07-10.log`).  Details and caveats in `NOTES.md`.

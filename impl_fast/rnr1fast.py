#!/usr/bin/env python3
"""Thin ctypes wrapper around librnr1fast.dylib (rnr1_fast.c).

API mirrors the normative reference impl/rnr1.py:

    import rnr1fast
    raw = rnr1fast.pack(data, W=3, K=65536)     # bytes -> archive bytes
    data = rnr1fast.unpack(raw, verify=True)    # archive bytes -> bytes

CLI mirrors the reference CLI (pack/unpack subcommands):

    rnr1fast.py pack   <input> <output> [--W 3] [--K 65536]
    rnr1fast.py unpack <input> <output> [--no-verify]

Build the library first:  sh impl_fast/build.sh
"""

import argparse
import ctypes
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
LIB_PATH = os.path.join(HERE, "librnr1fast.dylib")

_lib = None


def _load():
    global _lib
    if _lib is None:
        if not os.path.exists(LIB_PATH):
            raise RuntimeError("librnr1fast.dylib not built; run impl_fast/build.sh")
        lib = ctypes.CDLL(LIB_PATH)
        lib.rnr1_pack.restype = ctypes.c_int64
        lib.rnr1_pack.argtypes = [
            ctypes.c_char_p, ctypes.c_uint64, ctypes.c_uint32, ctypes.c_uint32,
            ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_uint64)]
        lib.rnr1_unpack.restype = ctypes.c_int64
        lib.rnr1_unpack.argtypes = [
            ctypes.c_char_p, ctypes.c_uint64, ctypes.c_int,
            ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_uint64)]
        lib.rnr1_pack_mt.restype = ctypes.c_int64
        lib.rnr1_pack_mt.argtypes = [
            ctypes.c_char_p, ctypes.c_uint64, ctypes.c_uint32, ctypes.c_uint32,
            ctypes.c_uint32,
            ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_uint64)]
        lib.rnr1_unpack_mt.restype = ctypes.c_int64
        lib.rnr1_unpack_mt.argtypes = [
            ctypes.c_char_p, ctypes.c_uint64, ctypes.c_int, ctypes.c_uint32,
            ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_uint64)]
        lib.rnr1_free.restype = None
        lib.rnr1_free.argtypes = [ctypes.c_void_p]
        lib.rnr1_strerror.restype = ctypes.c_char_p
        lib.rnr1_strerror.argtypes = [ctypes.c_int64]
        _lib = lib
    return _lib


def _take(lib, rc, outp, outlen):
    if rc != 0:
        raise ValueError(lib.rnr1_strerror(rc).decode("ascii"))
    try:
        return ctypes.string_at(outp.value, outlen.value)
    finally:
        lib.rnr1_free(outp.value)


def pack(data: bytes, W: int = 3, K: int = 65536, threads: int = 1) -> bytes:
    """Encode; threads > 1 parallelizes over independent sub-blocks
    (byte-identical output, verified by run_checks F3b)."""
    lib = _load()
    outp = ctypes.c_void_p()
    outlen = ctypes.c_uint64(0)
    if threads > 1:
        rc = lib.rnr1_pack_mt(data, len(data), W, K, threads,
                              ctypes.byref(outp), ctypes.byref(outlen))
    else:
        rc = lib.rnr1_pack(data, len(data), W, K,
                           ctypes.byref(outp), ctypes.byref(outlen))
    return _take(lib, rc, outp, outlen)


def unpack(raw: bytes, verify: bool = True, threads: int = 1) -> bytes:
    lib = _load()
    outp = ctypes.c_void_p()
    outlen = ctypes.c_uint64(0)
    if threads > 1:
        rc = lib.rnr1_unpack_mt(raw, len(raw), 1 if verify else 0, threads,
                                ctypes.byref(outp), ctypes.byref(outlen))
    else:
        rc = lib.rnr1_unpack(raw, len(raw), 1 if verify else 0,
                             ctypes.byref(outp), ctypes.byref(outlen))
    return _take(lib, rc, outp, outlen)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="rnr1fast")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pack")
    p.add_argument("input"); p.add_argument("output")
    p.add_argument("--W", type=int, default=3)
    p.add_argument("--K", type=int, default=65536)
    p.add_argument("--threads", type=int, default=1)
    p = sub.add_parser("unpack")
    p.add_argument("input"); p.add_argument("output")
    p.add_argument("--no-verify", action="store_true")
    p.add_argument("--threads", type=int, default=1)
    args = ap.parse_args(argv)
    with open(args.input, "rb") as f:
        data = f.read()
    t0 = time.time()
    if args.cmd == "pack":
        out = pack(data, W=args.W, K=args.K, threads=args.threads)
    else:
        out = unpack(data, verify=not args.no_verify, threads=args.threads)
    dt = time.time() - t0
    with open(args.output, "wb") as f:
        f.write(out)
    src = len(data) if args.cmd == "pack" else len(out)
    print("%s %d -> %d bytes in %.3fs (%.2f MB/s of source)"
          % (args.cmd, len(data), len(out), dt, src / 1e6 / max(dt, 1e-9)),
          file=sys.stderr)


if __name__ == "__main__":
    main()

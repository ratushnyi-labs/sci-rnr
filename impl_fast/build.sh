#!/bin/sh
# Build librnr1fast.dylib from rnr1_fast.c (macOS; CommonCrypto is in
# libSystem, no extra link flags needed).
set -e
cd "$(dirname "$0")"
clang -O3 -std=c11 -Wall -Wextra -Wno-deprecated-declarations \
      -fPIC -shared -o librnr1fast.dylib rnr1_fast.c
echo "built librnr1fast.dylib"

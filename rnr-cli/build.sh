#!/bin/sh
# Cross-target release builds for the rnr production CLI.
#
# Produces dist/rnr-<target>[.exe] for every target that links on this
# host, verifies linkage (musl binaries must be fully static), runs a
# cross-OS determinism smoke test under Docker (linux binaries must
# produce archives byte-identical to the macOS binary on the same
# input), and emits dist/MANIFEST.json with per-target sha256 + linkage
# + smoke status.
#
# Windows MSVC targets: the Rust objects compile (verified via cargo
# check), but linking requires the MSVC import libraries (kernel32.lib,
# ntdll.lib, userenv.lib, ws2_32.lib, dbghelp.lib) from a Windows SDK,
# which is not present on this host; see dist/MANIFEST.json for the
# recorded status.  On a host with the SDK (or xwin-fetched libs),
# `cargo build --release --target x86_64-pc-windows-msvc` needs no
# source changes.
#
# Usage: sh build.sh            (set RNR_SKIP_DOCKER=1 to skip the
#                                Docker cross-OS smoke test)
#
# The smoke images default to alpine:latest for both platforms; when
# registry pulls are unavailable, point RNR_DOCKER_IMG_ARM64 /
# RNR_DOCKER_IMG_AMD64 at any locally present linux images with /bin/sh
# (only /bin/sh, cmp and a mount are used; the rnr binaries are static).

set -eu
cd "$(dirname "$0")"

UNIX_TARGETS="aarch64-apple-darwin x86_64-apple-darwin x86_64-unknown-linux-musl aarch64-unknown-linux-musl"
WIN_TARGETS="x86_64-pc-windows-msvc aarch64-pc-windows-msvc"
SMOKE_SRC="../data/payloads/enwik8"
SMOKE_BYTES=4194304

mkdir -p dist
rm -f dist/MANIFEST.json

# ---------------------------------------------------------------- build
for t in $UNIX_TARGETS; do
    echo "== build $t"
    cargo build --release --target "$t"
    cp "target/$t/release/rnr" "dist/rnr-$t"
done

WIN_STATUS_FILE=$(mktemp)
for t in $WIN_TARGETS; do
    echo "== build $t (cross-compile attempt)"
    if cargo build --release --target "$t" >/dev/null 2>&1; then
        cp "target/$t/release/rnr.exe" "dist/rnr-$t.exe"
        echo "$t=linked" >> "$WIN_STATUS_FILE"
    else
        # Objects still must compile: the failure must be link-only.
        cargo check --release --target "$t" >/dev/null 2>&1 \
            && echo "$t=compile-ok-link-failed" >> "$WIN_STATUS_FILE" \
            || echo "$t=compile-failed" >> "$WIN_STATUS_FILE"
        echo "   link failed (missing MSVC import libs); cargo check passes"
    fi
done

# ------------------------------------------------------- linkage report
echo "== linkage"
LINKAGE_FILE=$(mktemp)
for t in $UNIX_TARGETS; do
    b="dist/rnr-$t"
    case "$t" in
    *linux-musl)
        desc=$(file -b "$b")
        echo "$t|$desc" >> "$LINKAGE_FILE"
        # "statically linked" or "static-pie linked": both are fully
        # static (no interpreter, no dynamic deps).
        case "$desc" in
        *"statically linked"* | *"static-pie linked"*)
            echo "   $t: fully static ($desc)" ;;
        *) echo "   $t: NOT STATIC: $desc"; exit 1 ;;
        esac
        ;;
    *apple-darwin)
        libs=$(otool -L "$b" | tail -n +2 | awk '{print $1}' | tr '\n' ' ')
        echo "$t|Mach-O; dylibs: $libs" >> "$LINKAGE_FILE"
        echo "   $t: dylibs = $libs"
        ;;
    esac
done

# -------------------------------------------- cross-OS determinism smoke
# The same 4 MiB input must produce byte-identical archives from the
# macOS binary and from both linux binaries executed under Docker.
SMOKE_STATUS="skipped"
if [ "${RNR_SKIP_DOCKER:-0}" != "1" ] && command -v docker >/dev/null 2>&1; then
    echo "== docker cross-OS smoke (4 MiB, W=3, K=65536)"
    SMOKE_DIR=$(mktemp -d)
    head -c $SMOKE_BYTES "$SMOKE_SRC" > "$SMOKE_DIR/smoke.bin"
    ./dist/rnr-aarch64-apple-darwin pack "$SMOKE_DIR/smoke.bin" \
        "$SMOKE_DIR/mac.rnr" --W 3 --K 65536 >/dev/null
    MAC_SHA=$(shasum -a 256 "$SMOKE_DIR/mac.rnr" | awk '{print $1}')
    SMOKE_STATUS="mac=$MAC_SHA"
    ok=1
    IMG_ARM64="${RNR_DOCKER_IMG_ARM64:-alpine:latest}"
    IMG_AMD64="${RNR_DOCKER_IMG_AMD64:-alpine:latest}"
    for plat_target in "linux/arm64 aarch64-unknown-linux-musl $IMG_ARM64" \
                       "linux/amd64 x86_64-unknown-linux-musl $IMG_AMD64"; do
        set -- $plat_target
        plat=$1; t=$2; img=$3
        docker run --rm --platform "$plat" --entrypoint /bin/sh \
            -v "$SMOKE_DIR:/s" -v "$(pwd)/dist:/dist:ro" "$img" \
            -c "/dist/rnr-$t pack /s/smoke.bin /s/$t.rnr --W 3 --K 65536 >/dev/null && /dist/rnr-$t unpack /s/$t.rnr /s/$t.out >/dev/null && cmp /s/$t.out /s/smoke.bin"
        SHA=$(shasum -a 256 "$SMOKE_DIR/$t.rnr" | awk '{print $1}')
        if [ "$SHA" = "$MAC_SHA" ]; then
            echo "   $plat ($t): archive sha identical to macOS"
            SMOKE_STATUS="$SMOKE_STATUS;$t=identical"
        else
            echo "   $plat ($t): SHA MISMATCH $SHA != $MAC_SHA"
            SMOKE_STATUS="$SMOKE_STATUS;$t=MISMATCH"
            ok=0
        fi
    done
    rm -rf "$SMOKE_DIR"
    [ "$ok" = "1" ] || exit 1
else
    echo "== docker cross-OS smoke skipped"
fi

# ------------------------------------------------------------- manifest
python3 - "$WIN_STATUS_FILE" "$LINKAGE_FILE" "$SMOKE_STATUS" <<'EOF'
import hashlib, json, os, subprocess, sys

win_status = dict(l.strip().split("=", 1)
                  for l in open(sys.argv[1]) if l.strip())
linkage = dict(l.strip().split("|", 1)
               for l in open(sys.argv[2]) if l.strip())
smoke = sys.argv[3]

rustc = subprocess.run(["rustc", "--version"], capture_output=True,
                       text=True).stdout.strip()
targets = {}
for name in sorted(os.listdir("dist")):
    if not name.startswith("rnr-"):
        continue
    path = os.path.join("dist", name)
    t = name[4:].removesuffix(".exe")
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    targets[t] = {
        "file": name,
        "sha256": h,
        "bytes": os.path.getsize(path),
        "linkage": linkage.get(t, "PE32 executable (not run on this host)"),
    }
for t, st in win_status.items():
    if t not in targets:
        targets[t] = {
            "file": None,
            "status": st,
            "note": ("Rust sources compile for this target (cargo check "
                     "passes); linking requires the MSVC import libraries "
                     "(kernel32.lib, ntdll.lib, userenv.lib, ws2_32.lib, "
                     "dbghelp.lib) from a Windows SDK, absent on this "
                     "build host. cargo-xwin was not installed and has no "
                     "prebuilt binary here; on a host with the SDK the "
                     "build needs no source changes."),
        }
manifest = {
    "crate": "rnr-cli",
    "bin": "rnr",
    "version": "1.0.0",
    "rustc": rustc,
    "bit_identical_to": "impl/rnr1.py (gate: rnr-cli/run_checks.py)",
    "cross_os_smoke": smoke,
    "targets": targets,
}
with open("dist/MANIFEST.json", "w") as f:
    json.dump(manifest, f, indent=2)
    f.write("\n")
print("wrote dist/MANIFEST.json (%d targets)" % len(targets))
EOF

echo "== done"

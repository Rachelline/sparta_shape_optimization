#!/usr/bin/env python3
"""
ad_convert.py -- mechanical double -> sfloat conversion of the SPARTA physics
core, so that a -DSPARTA_AD build carries forward-mode derivatives through the
whole solver (sfloat = Sacado SFad in the AD build, == double in the stock
build; see src/sfloat.h).

WHAT IT DOES (per top-level src/*.cpp,*.h file, excluding the passive set):
  - token `double`     -> `sfloat`
  - token `MPI_DOUBLE` -> `MPI_SFLOAT`   (macro resolves per build)
  - `union ubuf { double d; ... }` reverted to keep `double d`: it bit-puns 8
    bytes with an int64, so the pun target must stay a real double.
  - prepends a marker comment so converted files are identifiable / idempotent.

Only top-level src/ is processed (no recursion): src/KOKKOS, src/FFT,
src/STUBS, src/PYTHON are never touched -- consistent with the AD build being
serial/stub, non-package.

STOCK build is a functional no-op (sfloat==double), so it must still compile
and reproduce gold logs. Hazards only bite the AD build and are caught by the
mpi-stubs-ad-sacado gold-log CI job for any path the suite exercises.

EXCLUDE LIST -- passive by design (kept double):
  sfloat.h/spatype.h (the AD seam), math_extra (hand-templated), random_*
  (RNG must not carry derivatives -- frozen-stream semantics), timer
  (wallclock), memory (byte allocator; gets its own calloc-under-AD edit),
  library (C API stays double-facing), main.cpp, spawindows.h.

KNOWN LATENT HAZARDS (audited 2026-07-21, deliberately NOT handled in this
first pass -- shape-opt never exercises these paths; harden later if needed):
  - Binary I/O (fwrite/fread of doubles) in write_restart/read_restart, dump_*,
    *_custom, read_isurf/write_isurf: under -DSPARTA_AD these serialize
    sizeof(sfloat) (=N+1 doubles) instead of 8 bytes, so AD-build restart/dump
    files are self-consistent but NOT interchangeable with stock files, and
    embed meaningless derivative data. The right fix is per-site spval() on
    write / read-into-double; deferred.
  - memcpy/memset sized by sizeof(double) in ~30 files (cut3d, fix_ablate, ...):
    correct only where the buffer truly holds doubles; audit before trusting
    those paths under AD.
These are documented in docs/PLAN.md ("double->sfloat conversion").

USAGE:
  python3 tools/ad_convert.py --dry-run   # report scope, write nothing
  python3 tools/ad_convert.py             # convert in place (idempotent)
  python3 tools/ad_convert.py --revert     # strip conversion (restore stock)
"""

import os
import re
import sys

SRC = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "src"))

MARKER = "/* AD-CONVERTED: double->sfloat by tools/ad_convert.py (see sfloat.h) */\n"

EXCLUDE_FILES = {
    "sfloat.h", "spatype.h",
    "math_extra.h", "math_extra.cpp",        # hand-templated
    "random_mars.h", "random_mars.cpp",      # RNG: passive (frozen stream)
    "random_knuth.h", "random_knuth.cpp",
    "timer.h", "timer.cpp",                  # wallclock: passive
    "memory.h", "memory.cpp",                # byte allocator (calloc-under-AD)
    "library.h", "library.cpp",              # C API stays double-facing
    "main.cpp",
    "spawindows.h",
}

UBUF_RE = re.compile(r"(union\s+ubuf\s*\{[^}]*\})", re.S)

# Low-level headers that USE sfloat but don't (transitively) include spatype.h,
# so they must include sfloat.h directly. spatype.h itself is the main bootstrap
# (excluded, edited by hand); every other converted file gets sfloat via it.
BOOTSTRAP_INCLUDE = {"geometry.h", "math_eigen.h"}
INCLUDE_LINE = '#include "sfloat.h"\n'
GUARD_RE = re.compile(r'(#define\s+SPARTA_\w+_H[ \t]*\n)')


def convert_text(text, fname=None):
    text = re.sub(r"\bMPI_DOUBLE\b", "MPI_SFLOAT", text)
    text = re.sub(r"\bdouble\b", "sfloat", text)
    text = UBUF_RE.sub(lambda m: m.group(1).replace("sfloat", "double"), text)
    if fname in BOOTSTRAP_INCLUDE and INCLUDE_LINE not in text:
        # insert right after the include guard's #define
        text = GUARD_RE.sub(lambda m: m.group(1) + "\n" + INCLUDE_LINE, text, count=1)
    return text


def iter_src_files():
    for fname in sorted(os.listdir(SRC)):
        path = os.path.join(SRC, fname)
        if os.path.isdir(path):
            continue
        if not (fname.endswith(".cpp") or fname.endswith(".h")):
            continue
        if fname in EXCLUDE_FILES:
            continue
        yield fname, path


def main():
    dry = "--dry-run" in sys.argv
    revert = "--revert" in sys.argv
    converted = skipped = reverted = 0
    ndoubles = 0

    for fname, path in iter_src_files():
        with open(path) as f:
            text = f.read()
        has_marker = text.startswith(MARKER)

        if revert:
            if has_marker:
                body = text[len(MARKER):]
                # restore: sfloat->double, MPI_SFLOAT->MPI_DOUBLE, drop the
                # bootstrap include; protect the "sfloat.h" filename token.
                body = body.replace(INCLUDE_LINE, "")
                body = re.sub(r"\bMPI_SFLOAT\b", "MPI_DOUBLE", body)
                body = re.sub(r"\bsfloat\b(?!\.h)", "double", body)
                if not dry:
                    with open(path, "w") as f:
                        f.write(body)
                reverted += 1
            continue

        if has_marker:                      # idempotent
            skipped += 1
            continue

        ndoubles += len(re.findall(r"\bdouble\b", text))
        if not dry:
            with open(path, "w") as f:
                f.write(MARKER + convert_text(text, fname))
        converted += 1

    if revert:
        print(f"{'[dry] ' if dry else ''}reverted {reverted} files")
    else:
        verb = "would convert" if dry else "converted"
        print(f"{verb} {converted} files ({ndoubles} `double` tokens), "
              f"{skipped} already-marked, {len(EXCLUDE_FILES)} on the exclude list")


if __name__ == "__main__":
    main()

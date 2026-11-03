"""Choice of the AMPL executable for the tests: a sanitizer build
(e.g., the first 'ampl' on PATH in a development setup) prints its
diagnostics into every solver log and runs slower, so a release AMPL
can be given as a fallback."""
import os
import shutil

SANITIZER_MARKERS = (b"__asan_init", b"__ubsan_handle", b"__tsan_init")


def is_sanitizer_build(path: str) -> bool:
    """Whether the executable at path contains sanitizer runtime symbols"""
    try:
        with open(path, "rb") as f:
            data = f.read()
    except OSError:
        return False
    return any(m in data for m in SANITIZER_MARKERS)


def resolve(ampl: str, fallback: str = "") -> str:
    """The AMPL executable to use: ampl, or fallback if ampl is not
    found or is a sanitizer build (and fallback is found)"""
    path = shutil.which(ampl) or ampl
    reason = ""
    if not os.path.isfile(path):
        reason = "not found"
    elif is_sanitizer_build(path):
        reason = "a sanitizer build"
    if not reason:
        return path
    fb = shutil.which(fallback) if fallback else None
    if fb:
        print("AMPL executable '{}' is {}: using the fallback '{}'".format(
            path, reason, fb))
        return fb
    if fallback:
        print("WARNING: AMPL executable '{}' is {}, fallback '{}' "
              "not found".format(path, reason, fallback))
    return path

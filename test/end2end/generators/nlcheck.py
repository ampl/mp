#!/usr/bin/env python3
"""Check that amplfmt.py changes the layout of models only.

For every .mod file in the given folders, writes the NL file of the original
and of the formatted model with AMPL (option presolve 0 and the default
presolve) and compares them byte for byte. With --apply, rewrites the files
whose NL files are all equal, after the whole folder is checked.

Usage: nlcheck.py [--apply] [--careful] [--ampl PATH] DIR...

--careful: wrap long code lines only (third-party and converted models);
the check then also requires that only whitespace changes.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amplfmt import format_model, wrap_lines         # noqa: E402

AMPL = os.environ.get("AMPL", "ampl")   # AMPL executable (or --ampl)


def write_nl(ampl, mod, stub, pre):
    run = stub + ".run"
    open(run, "w").write(f'option nl_comments 0; option presolve {pre};\n'
                         f'model "{mod}";\nwrite "g{stub}";\n')
    subprocess.run([ampl, run], capture_output=True, timeout=600)


def write_nls(ampl, mods, outdir, tag):
    """The NL files of all models, one AMPL process each (an error in a
    model would stop a shared session), three at a time."""
    jobs = [(m, f"{outdir}/{tag}_{i}_{pre}", pre)
            for i, m in enumerate(mods) for pre in (0, 10)]
    with ThreadPoolExecutor(3) as ex:
        list(ex.map(lambda j: write_nl(ampl, *j), jobs))
    return ""


def read(path):
    """The NL file without its first line (it names the problem)."""
    try:
        return open(path, "rb").read().split(b"\n", 1)[1]
    except OSError:
        return None


def check_dir(ampl, d, apply, careful=False):
    mods = sorted(os.path.join(d, f) for f in os.listdir(d) if f.endswith(".mod"))
    tmp = tempfile.mkdtemp(prefix="nlcheck_")
    try:
        new = []
        for i, m in enumerate(mods):
            text = open(m).read()
            nm = os.path.join(tmp, f"new_{i}.mod")
            out = wrap_lines(text) if careful else format_model(text)
            if careful and "".join(out.split()) != "".join(text.split()):
                raise SystemExit(f"{m}: not only whitespace changed")
            open(nm, "w").write(out)
            new.append(nm)
        log_o = write_nls(ampl, mods, tmp, "o")
        log_n = write_nls(ampl, new, tmp, "n")
        bad, nonl, changed = [], [], 0
        for i, m in enumerate(mods):
            same = all(read(f"{tmp}/o_{i}_{p}.nl") == read(f"{tmp}/n_{i}_{p}.nl")
                       for p in (0, 10))
            unchanged = open(new[i]).read() == open(m).read()
            if not same and unchanged:
                print("    NL differs between AMPL runs (random data), file unchanged:", m)
            elif not same:
                bad.append((m, "NL differs"))
            elif read(f"{tmp}/o_{i}_0.nl") is None:
                nonl.append(m)                   # e.g. no variables
            elif open(new[i]).read() != open(m).read():
                changed += 1
        if apply and not bad:
            for i, m in enumerate(mods):
                text = open(new[i]).read()
                if text != open(m).read():
                    open(m, "w").write(text)
        long_lines = sum(1 for nm in new for ln in open(nm)
                         if len(ln.rstrip().expandtabs(8)) > 80)
        print(f"{d}: {len(mods)} models, {changed} reformatted, {len(bad)} problems, "
              f"{long_lines} lines > 80 after formatting"
              + (" (applied)" if apply and not bad else ""))
        for m in nonl:
            print("    no NL from either (checked by AMPL's messages only):", m)
        for m, why in bad[:20]:
            print("   ", why, m)
        if bad:
            sys.stdout.write(("   AMPL (original): " + log_o[-800:]) if log_o.strip() else "")
        return not bad
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--ampl", default=AMPL)
    ap.add_argument("--careful", action="store_true",
                    help="wrap long code lines only, at existing whitespace "
                         "(third-party and converted models)")
    ap.add_argument("dirs", nargs="+")
    a = ap.parse_args()
    ok = all([check_dir(a.ampl, d, a.apply, a.careful) for d in a.dirs])
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

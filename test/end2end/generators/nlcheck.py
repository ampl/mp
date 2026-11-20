#!/usr/bin/env python3
"""Check that amplfmt.py changes the layout of models only.

For every .mod file in the given folders, writes the NL file of the original
and of the formatted model with AMPL (option presolve 0 and the default
presolve) and compares them byte for byte. Only what reaches the NL file is
compared: a change in an unused declaration or a comment is not seen.

With --apply, the reformatted files of a folder are rewritten after the whole
folder is checked, and only if no model of it has a problem. A problem is:
AMPL fails or writes no NL file for the original or the formatted model
(also for a valid model without variables, or one that presolve finds
infeasible), the NL files differ, or the source changed during the check.
One such model leaves the whole folder as it is.

The rewriting itself is done file by file. If a source changes between the
check and its rewriting, that file is kept and reported, the files rewritten
before it stay rewritten, the summary line lacks "(applied)" and the exit
code is 1.

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
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from amplfmt import format_model, wrap_lines         # noqa: E402

AMPL = os.environ.get("AMPL", "ampl")   # AMPL executable (or --ampl)


def write_nl(ampl, mod, stub, pre):
    run = stub + ".run"
    Path(run).write_text(f'option nl_comments 0; option presolve {pre};\n'
                         f'model "{mod}";\nwrite "g{stub}";\n')
    return subprocess.run([ampl, run], capture_output=True, text=True,
                          errors="replace", timeout=600)


def write_nls(ampl, mods, outdir, tag):
    """The NL files of all models, one AMPL process each (an error in a
    model would stop a shared session), three at a time."""
    jobs = [(m, f"{outdir}/{tag}_{i}_{pre}", pre)
            for i, m in enumerate(mods) for pre in (0, 10)]
    with ThreadPoolExecutor(3) as ex:
        return list(ex.map(lambda j: write_nl(ampl, *j), jobs))


def read(path):
    """The NL file without its first line (it names the problem)."""
    try:
        with open(path, "rb") as source:
            _, _, body = source.read().partition(b"\n")
        return body or None
    except OSError:
        return None


def apply_model(model, original, formatted):
    """Replace an unchanged source without truncating it in place."""
    path = Path(model).resolve()
    with tempfile.TemporaryDirectory(prefix=".nlcheck-", dir=path.parent) as tmp:
        candidate = Path(tmp) / path.name
        candidate.write_text(formatted)
        shutil.copymode(path, candidate)
        if path.read_text() != original:
            return False
        candidate.replace(path)
    return True


def check_dir(ampl, d, apply, careful=False):
    mods = sorted(os.path.join(d, f) for f in os.listdir(d) if f.endswith(".mod"))
    tmp = tempfile.mkdtemp(prefix="nlcheck_")
    try:
        inputs, new = [], []
        long_lines = 0
        for i, m in enumerate(mods):
            text = Path(m).read_text()
            inputs.append(text)
            nm = os.path.join(tmp, f"new_{i}.mod")
            out = wrap_lines(text) if careful else format_model(text)
            if careful and "".join(out.split()) != "".join(text.split()):
                raise SystemExit(f"{m}: not only whitespace changed")
            Path(nm).write_text(out)
            long_lines += sum(len(ln.rstrip().expandtabs(8)) > 80
                              for ln in out.splitlines())
            new.append(nm)
        original_runs = write_nls(ampl, mods, tmp, "o")
        formatted_runs = write_nls(ampl, new, tmp, "n")
        bad, changed = [], []
        for i, m in enumerate(mods):
            pairs = [(read(f"{tmp}/o_{i}_{p}.nl"), read(f"{tmp}/n_{i}_{p}.nl"))
                     for p in (0, 10)]
            runs = original_runs[2*i:2*i+2] + formatted_runs[2*i:2*i+2]
            if any(run.returncode for run in runs) or any(
                    a is None or b is None for a, b in pairs):
                messages = "\n".join((run.stderr or run.stdout).strip()
                                     for run in runs if run.stderr or run.stdout)
                bad.append((m, "AMPL failed or did not produce NL", messages))
                continue
            same = all(a == b for a, b in pairs)
            formatted = Path(new[i]).read_text()
            unchanged = formatted == inputs[i]
            if Path(m).read_text() != inputs[i]:
                bad.append((m, "Source changed during validation", ""))
                continue
            if not same and unchanged:
                print("    NL differs between AMPL runs (random data), file unchanged:", m)
            elif not same:
                bad.append((m, "NL differs", ""))
            elif not unchanged:
                changed.append((m, inputs[i], formatted))
        if apply and not bad:
            for m, original, formatted in changed:
                if not apply_model(m, original, formatted):
                    bad.append((m, "Source changed before applying", ""))
        print(f"{d}: {len(mods)} models, {len(changed)} reformatted, {len(bad)} problems, "
              f"{long_lines} lines > 80 after formatting"
              + (" (applied)" if apply and not bad else ""))
        for m, why, messages in bad[:20]:
            print("   ", why, m)
            if messages:
                print(messages[-800:])
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

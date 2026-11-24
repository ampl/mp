"""AMPLS regression: a failed AMPLSSolve() returns a code and a message.

Usage: python ampls-solve-error-test.py DRIVER_LIB [--solver NAME]
       [--libdir DIR ...] [--ampl AMPL]

DRIVER_LIB is a shared MP driver library exporting AMPLSOpen_<NAME>
(default NAME: highs, e.g., highsmp-lib.dll; HiGHS needs no license).
The solver must take multiple objectives through MP's emulator, which
does not support a repeated solve. --libdir adds a folder with the
solver's runtime libraries to the search path.

Checks, on one model with two objectives:
- with obj:multi=1, the first AMPLSSolve() returns 0 and is reported as
  solved; the second one fails without an exception leaving it: it
  returns nonzero and its message is among AMPLSGetMessages(). The test
  then reports the failure as a caller does, with AMPLSReportError(),
  and checks the .sol file (solve result and message);
- with one objective, two AMPLSSolve() calls both succeed;
- AMPLSGetSolveResult(): NOT_SET after loading, the reported result after a
  solve, the error's code (or 500) after a failed call; a corrupt NL file
  fails to load with -1 and gives 500.
"""

import argparse
import ctypes as ct
import os
from pathlib import Path
import re
import subprocess
import tempfile


PTR, INT, STR = ct.c_void_p, ct.c_int, ct.c_char_p
MODEL = """var x >= 0, <= 10;
var y >= 0, <= 10;
s.t. cap: x + y <= 12;
maximize f1: x;
maximize f2: y;
"""
REPEAT_MESSAGE = b"A repeated solve is not supported"
# sol::NOT_SET (mp/common.h): no solve result yet
NOT_SET = -200


def bind(library, name, result, *arguments):
    function = getattr(library, name)
    function.restype = result
    function.argtypes = arguments
    return function


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def solve_result(sol_path):
    """The solve result code written in a .sol file."""
    found = re.findall(rb"^objno -?\d+ (-?\d+)", sol_path.read_bytes(),
                       re.MULTILINE)
    require(found, f"No solve result in {sol_path}")
    return int(found[-1])


class Ampls:
    def __init__(self, driver, name):
        self.name = name
        self.open = bind(driver, f"AMPLSOpen_{name}", PTR, INT,
                         ct.POINTER(STR))
        self.close = bind(driver, f"AMPLSClose_{name}", None, PTR)
        self.load = bind(driver, "AMPLSLoadNLModel", INT, PTR, STR,
                         ct.POINTER(STR))
        self.solve = bind(driver, "AMPLSSolve", INT, PTR)
        self.report = bind(driver, "AMPLSReportResults", INT, PTR, STR)
        self.report_error = bind(driver, "AMPLSReportError", INT, PTR, INT,
                                 STR)
        self.messages = bind(driver, "AMPLSGetMessages", ct.POINTER(STR),
                             PTR)
        self.solve_result = bind(driver, "AMPLSGetSolveResult", INT, PTR)

    def all_messages(self, solver):
        array, text, i = self.messages(solver), b"", 0
        while array[i] is not None:
            text += array[i]
            i += 1
        return text

    def run(self, nl, options, check, load_ok=True):
        argv = (STR * 2)(self.name.encode(), None)  # must outlive the solver
        solver = self.open(1, argv)
        require(solver, "AMPLS could not open the solver")
        try:
            opts = (STR * (len(options) + 1))(*options, None)
            loaded = self.load(solver, os.fsencode(nl), opts)
            if load_ok:
                require(loaded == 0, "AMPLS could not load the NL file")
                require(self.solve_result(solver) == NOT_SET,
                        f"Solve result before solving: "
                        f"{self.solve_result(solver)}, expected {NOT_SET}")
            check(solver) if load_ok else check(solver, loaded)
        finally:
            self.close(solver)


def check_failed_repeat(ampls, nl):
    def check(solver):
        first = nl.with_name("first.sol")
        require(ampls.solve(solver) == 0, "The first solve failed")
        require(ampls.report(solver, os.fsencode(first)) == 0,
                "Reporting the first solve failed")
        require(solve_result(first) < 200,
                f"The first solve is not reported solved: {solve_result(first)}")
        require(ampls.solve_result(solver) == solve_result(first),
                "AMPLSGetSolveResult() differs from the reported result")
        code = ampls.solve(solver)
        require(code != 0, "The repeated multiobjective solve did not fail")
        require(code == -1 or code > 0,
                f"Unexpected return code of the failed solve: {code}")
        messages = ampls.all_messages(solver)
        require(REPEAT_MESSAGE in messages,
                "The failed solve added no message")
        # Report the failure as a caller does (e.g., the MP2NL wrapper):
        # AMPLSReportError() writes the .sol file of the last report
        result = code if code > 0 else 500
        require(ampls.solve_result(solver) == result,
                f"AMPLSGetSolveResult() after the failed solve: "
                f"{ampls.solve_result(solver)}, expected {result}")
        require(ampls.report_error(solver, result, messages) == 0,
                "AMPLSReportError failed")
        require(ampls.solve_result(solver) == result,
                "AMPLSReportError() changed AMPLSGetSolveResult()")
        require(solve_result(first) == result,
                f"The failure is reported with solve result "
                f"{solve_result(first)}, expected {result}")
        require(REPEAT_MESSAGE in first.read_bytes(),
                "The .sol file of the failure lacks its message")
    ampls.run(nl, [b"obj:multi=1"], check)


def check_repeat_one_objective(ampls, nl):
    def check(solver):
        for i in (1, 2):
            require(ampls.solve(solver) == 0, f"Solve {i} failed")
            before_report = ampls.solve_result(solver)
            sol = nl.with_name(f"single{i}.sol")
            require(ampls.report(solver, os.fsencode(sol)) == 0,
                    f"Reporting solve {i} failed")
            require(solve_result(sol) < 200,
                    f"Solve {i} is not reported solved: {solve_result(sol)}")
            require(ampls.solve_result(solver) == solve_result(sol),
                    "AMPLSGetSolveResult() differs from the reported result")
            require(before_report == solve_result(sol),
                    f"AMPLSGetSolveResult() before the report of solve {i}: "
                    f"{before_report}, expected {solve_result(sol)}")
    ampls.run(nl, [], check)


def check_corrupt_nl(ampls, nl):
    """An NL file that cannot be read: the load returns -1 (no solve
    result code: formatted read errors once carried 1, i.e., 'solved'),
    and AMPLSGetSolveResult() gives 500"""
    bad = nl.with_name("corrupt.nl")
    bad.write_bytes(b"g3 1 1 0\n not an NL header\n")
    def check(solver, loaded):
        require(loaded == -1, f"Loading a corrupt NL file returned {loaded}")
        require(ampls.solve_result(solver) == 500,
                f"AMPLSGetSolveResult() after the failed load: "
                f"{ampls.solve_result(solver)}, expected 500")
    ampls.run(bad, [], check, load_ok=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("driver", type=Path)
    parser.add_argument("--solver", default="highs")
    parser.add_argument("--libdir", type=Path, action="append", default=[])
    parser.add_argument("--ampl", default="ampl")
    args = parser.parse_args()
    dirs = [args.driver.resolve().parent] + [d.resolve() for d in args.libdir]
    # Solver libraries loaded at run time (e.g., HiGHS) use the PATH
    os.environ["PATH"] = os.pathsep.join(
        [str(d) for d in dirs] + [os.environ.get("PATH", "")])
    handles = [os.add_dll_directory(str(d)) for d in dirs] \
        if hasattr(os, "add_dll_directory") else []
    try:
        ampls = Ampls(ct.CDLL(str(args.driver.resolve())), args.solver)
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            (directory / "case.mod").write_text(MODEL)
            (directory / "write.run").write_text(
                "option presolve 0; model case.mod; write gcase;\n")
            subprocess.run([args.ampl, "write.run"], cwd=directory,
                           check=True)
            nl = directory / "case.nl"
            check_failed_repeat(ampls, nl)
            print("failed repeated solve: passed")
            check_repeat_one_objective(ampls, nl)
            print("repeated solve, one objective: passed")
            check_corrupt_nl(ampls, nl)
            print("corrupt NL file: passed")
    finally:
        for handle in handles:
            handle.close()


if __name__ == "__main__":
    main()

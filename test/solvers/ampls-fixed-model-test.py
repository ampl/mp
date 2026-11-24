"""AMPLS regression for reporting and re-solving the same model.

Usage: python ampls-fixed-model-test.py SOLVER DRIVER_DLL SDK_DLL [--ampl AMPL]
SOLVER is cplex, gurobi or xpress. The library paths (DLL on Windows,
dylib or so elsewhere) refer to the build under test.
"""

import argparse
import contextlib
import ctypes as ct
import os
from pathlib import Path
import re
import subprocess
import tempfile


PTR, INT, STR, REAL = ct.c_void_p, ct.c_int, ct.c_char_p, ct.c_double
MODEL_DIR = (Path(__file__).resolve().parents[1] /
             "end2end/cases/categorized/fast")


def bind(library, name, result, *arguments):
    function = getattr(library, name)
    function.restype = result
    function.argtypes = arguments
    return function


def require(condition, message):
    if not condition:
        raise AssertionError(message)


# Every message of the fixed solve (success, skipped) starts with this text.
FIXED_SOLVE_PREFIX = b"Fixed MIP for mip:basis"
# The line reporting success, complete with its line end.
FIXED_SOLVE_SUCCESS = {
    "cplex": FIXED_SOLVE_PREFIX
             + rb": [0-9]+ (barrier|simplex) iterations?[\r\n]",
    "xpress": FIXED_SOLVE_PREFIX + rb": nonlinear solve[^\r\n]*[\r\n]",
    # Gurobi, and Xpress for a linear MIP, report the fixed solve only
    # when it took iterations
    "gurobi": None,
    "xpress-linear": None,
}


def check_fixed_solve_message(name, text):
    """Check that the report text describes one successful fixed solve.

    The text must be that of this report alone: a message kept from an
    earlier report of the same model would hide a failure of this one.
    """
    require(b" failed in Do" not in text, "The fixed solve failed")
    require(b"WARNING" not in text, "The report carries a warning")
    success = FIXED_SOLVE_SUCCESS[name]
    if success is None:
        require(text.count(FIXED_SOLVE_PREFIX) <= 1,
                "The report repeats the message of an earlier fixed solve")
        return
    require(len(re.findall(success, text)) == 1,
            "The fixed solve did not report success exactly once")
    require(text.count(FIXED_SOLVE_PREFIX) == 1,
            "The report repeats the message of an earlier fixed solve")


def cplex_checks(driver, sdk, solver, model, binary):
    env = bind(driver, "AMPLSGetEnv_cplex", PTR, PTR)(solver)
    problem_type = bind(sdk, "CPXgetprobtype", INT, PTR, PTR)
    get_primal = bind(sdk, "CPXgetx", INT, PTR, PTR, ct.POINTER(REAL), INT, INT)
    set_sense = bind(sdk, "CPXchgobjsen", INT, PTR, PTR, INT)
    value = REAL()

    def incumbent():
        require(get_primal(env, model, ct.byref(value), binary, binary) == 0,
                "CPLEX did not return the binary variable value")
        return value.value

    return lambda: problem_type(env, model), incumbent, lambda: set_sense(env, model, 1)


def gurobi_checks(sdk, model, binary):
    get_int = bind(sdk, "GRBgetintattr", INT, PTR, STR, ct.POINTER(INT))
    get_value = bind(sdk, "GRBgetdblattrelement", INT, PTR, STR, INT,
                     ct.POINTER(REAL))
    set_int = bind(sdk, "GRBsetintattr", INT, PTR, STR, INT)
    is_mip, value = INT(), REAL()

    def state():
        require(get_int(model, b"IsMIP", ct.byref(is_mip)) == 0,
                "Gurobi did not return IsMIP")
        return is_mip.value

    def incumbent():
        require(get_value(model, b"X", binary, ct.byref(value)) == 0,
                "Gurobi did not return the binary variable value")
        return value.value

    return state, incumbent, lambda: set_int(model, b"ModelSense", 1)


def sol_duals(path):
    """The dual values in the SOL file at path (text format)."""
    lines = path.read_bytes().decode("latin-1").splitlines()
    i = lines.index("Options") + 1
    i += int(lines[i]) + 1                  # the option values
    ncon, nduals = int(lines[i]), int(lines[i + 1])
    return [float(v) for v in lines[i + 4:i + 4 + nduals]] if ncon else []


def check_continuous_resolve(sdk, solver, model, binary, report, solve, path):
    """Gurobi: make the model continuous natively, solve and report.

    The duals must be those of the LP relaxation (dual of Lin 1), not those
    of the fixed model of the previous solve (minimization, all 0).
    """
    set_char = bind(sdk, "GRBsetcharattrelement", INT, PTR, STR, INT,
                    ct.c_char)
    set_int = bind(sdk, "GRBsetintattr", INT, PTR, STR, INT)
    require(set_char(model, b"VType", binary, b"C") == 0
            and set_int(model, b"ModelSense", -1) == 0,
            "Could not make the model continuous")
    solve(solver)
    require(report(solver, os.fsencode(path)) == 0, "AMPLS reporting failed")
    duals = sol_duals(path)
    require(duals and max(abs(d) for d in duals) > 0.5,
            f"The duals are those of an earlier fixed model: {duals}")


def xpress_checks(sdk, model, binary, first_binary):
    column_type = bind(sdk, "XPRSgetcoltype", INT, PTR, STR, INT, INT)
    get_primal = bind(sdk, "XPRSgetsolution", INT, PTR, ct.POINTER(INT),
                      ct.POINTER(REAL), INT, INT)
    set_bounds = bind(sdk, "XPRSchgbounds", INT, PTR, INT, ct.POINTER(INT),
                      STR, ct.POINTER(REAL))
    set_sense = bind(sdk, "XPRSchgobjsense", INT, PTR, INT)
    kind, value = ct.create_string_buffer(1), REAL()

    def state():
        require(column_type(model, kind, binary, binary) == 0,
                "Xpress did not return the binary column type")
        return kind.raw

    def incumbent():
        require(get_primal(model, None, ct.byref(value), binary, binary) == 0,
                "Xpress did not return the binary variable value")
        return value.value

    def reverse():
        """Force the other value of the binary.

        Nonlinear case (first value 0): lower bound 1. The objective stays
        as loaded: its exponential terms reach Xpress as an inequality
        valid for minimization only, and the solution check uses the
        objective of the NL file. Linear case (first value 1): minimize.
        """
        if first_binary == 1:
            return set_sense(model, 1)      # XPRS_OBJ_MINIMIZE
        return set_bounds(model, 1, ct.byref(INT(binary)), b"L",
                          ct.byref(REAL(1.0)))

    return state, incumbent, reverse


def check_repeated_solve(name, driver, sdk, nl, binary, first_binary,
                         message):
    open_solver = bind(driver, f"AMPLSOpen_{name}", PTR, INT, ct.POINTER(STR))
    close = bind(driver, f"AMPLSClose_{name}", None, PTR)
    native = bind(driver, f"AMPLSGetModel_{name}", PTR, PTR)
    load = bind(driver, "AMPLSLoadNLModel", INT, PTR, STR, ct.POINTER(STR))
    extras = bind(driver, "AMPLSReadExtras", None, PTR)
    solve = bind(driver, "AMPLSSolve", None, PTR)
    report = bind(driver, "AMPLSReportResults", INT, PTR, STR)
    argv = (STR * 2)(name.encode(), None)   # a driver may keep the pointer:
    solver = open_solver(1, argv)           # it must outlive the solver
    require(solver, "AMPLS could not open the solver")
    try:
        options = (STR * 3)(b"mip:basis=1", b"threads=1", None)
        require(load(solver, os.fsencode(nl), options) == 0,
                "AMPLS could not load the NL fixture")
        extras(solver)
        model = native(solver)
        if name == "cplex":
            state, incumbent, reverse = cplex_checks(driver, sdk, solver, model, binary)
        elif name == "gurobi":
            state, incumbent, reverse = gurobi_checks(sdk, model, binary)
        else:
            state, incumbent, reverse = xpress_checks(sdk, model, binary,
                                                      first_binary)
        solve(solver)
        original_type = state()
        require(abs(incumbent() - first_binary) < 1e-5,
                "The first solve returned an unexpected incumbent")
        sol = nl.with_suffix(".sol")

        def check_report(path, binary_value):
            require(report(solver, os.fsencode(path)) == 0, "AMPLS reporting failed")
            check_fixed_solve_message(message, path.read_bytes())
            require(native(solver) == model, "Reporting replaced the native model")
            require(state() == original_type, "Reporting changed the native problem type")
            require(abs(incumbent() - binary_value) < 1e-5,
                    "Reporting changed the original incumbent")

        # A report that fails (the solution file cannot be written) must
        # still leave the user's problem, not the fixed one
        bad = sol.with_name("no such directory") / "bad.sol"
        require(report(solver, os.fsencode(bad)) != 0,
                "Reporting to an unwritable file did not fail")
        require(state() == original_type,
                "A failed report left the native problem fixed")
        check_report(sol, first_binary)
        # A repeated report of the same solve carries the same notes:
        # the fixed solve is not repeated, its message is not lost.
        check_report(sol.with_name("repeat.sol"), first_binary)
        require(reverse() == 0, "Could not reverse the objective sense")
        solve(solver)
        second_binary = 1 - first_binary
        require(abs(incumbent() - second_binary) < 1e-5,
                "The second solve did not find the new incumbent")
        check_report(sol.with_name("second.sol"), second_binary)
        if name == "gurobi":
            check_continuous_resolve(sdk, solver, model, binary, report,
                                     solve, sol.with_name("third.sol"))
    finally:
        close(solver)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", choices=("cplex", "gurobi", "xpress"))
    parser.add_argument("driver", type=Path)
    parser.add_argument("sdk", type=Path)
    parser.add_argument("--ampl", default="ampl")
    args = parser.parse_args()
    # CPLEX: an MIQCP (fixed on a copy) and a linear MIP (fixed in the
    # native problem, which must be restored after each report).
    # Xpress: a nonlinear and a linear MIP (both fixed on a copy).
    # Each source with the binary's value in the first solve and the key
    # of the expected fixed-solve message.
    sources = {"cplex": (("mip_basis/fixed_miqcp_cplex.mod", 1, "cplex"),
                         ("mip_basis/fixed_milp_cplex.mod", 1, "cplex")),
               "gurobi": (("mip_basis/fixed_milp_cplex.mod", 1, "gurobi"),),
               "xpress": (("suf_std/nl_dual_01.mod", 0, "xpress"),
                          ("mip_basis/fixed_milp_cplex.mod", 1,
                           "xpress-linear"))}[args.solver]
    dll_dirs = contextlib.ExitStack()
    if hasattr(os, "add_dll_directory"):    # Windows
        for library in (args.driver, args.sdk):
            dll_dirs.enter_context(
                os.add_dll_directory(str(library.resolve().parent)))
    with dll_dirs:
        driver = ct.CDLL(str(args.driver.resolve()))
        sdk = ct.CDLL(str(args.sdk.resolve()))
        for source, first_binary, message in sources:
            with tempfile.TemporaryDirectory() as temp:
                directory = Path(temp)
                (directory / "case.mod").write_bytes(
                    (MODEL_DIR / source).read_bytes())
                (directory / "write.run").write_text(
                    "option presolve 0; option auxfiles rc; "
                    "model case.mod; write gcase;\n")
                subprocess.run([args.ampl, "write.run"], cwd=directory,
                               check=True)
                binary = (directory / "case.col").read_text().splitlines(
                    ).index("b")
                check_repeated_solve(args.solver, driver, sdk,
                                     directory / "case.nl", binary,
                                     first_binary, message)
            print(f"{source}: passed")


if __name__ == "__main__":
    main()

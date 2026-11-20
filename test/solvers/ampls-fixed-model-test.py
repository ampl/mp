"""Windows AMPLS regression for reporting and re-solving the same model.

Usage: python ampls-fixed-model-test.py SOLVER DRIVER_DLL SDK_DLL [--ampl AMPL]
SOLVER is cplex or xpress. The DLL paths refer to the build under test.
"""

import argparse
import ctypes as ct
import os
from pathlib import Path
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


def xpress_checks(sdk, model, binary):
    column_type = bind(sdk, "XPRSgetcoltype", INT, PTR, STR, INT, INT)
    get_primal = bind(sdk, "XPRSgetsolution", INT, PTR, ct.POINTER(INT),
                      ct.POINTER(REAL), INT, INT)
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

    return state, incumbent, lambda: set_sense(model, -1)


def check_repeated_solve(name, driver, sdk, nl, binary):
    open_solver = bind(driver, f"AMPLSOpen_{name}", PTR, INT, ct.POINTER(STR))
    close = bind(driver, f"AMPLSClose_{name}", None, PTR)
    native = bind(driver, f"AMPLSGetModel_{name}", PTR, PTR)
    load = bind(driver, "AMPLSLoadNLModel", INT, PTR, STR, ct.POINTER(STR))
    extras = bind(driver, "AMPLSReadExtras", None, PTR)
    solve = bind(driver, "AMPLSSolve", None, PTR)
    report = bind(driver, "AMPLSReportResults", INT, PTR, STR)
    solver = open_solver(1, (STR * 2)(name.encode(), None))
    require(solver, "AMPLS could not open the solver")
    try:
        options = (STR * 3)(b"mip:basis=1", b"threads=1", None)
        require(load(solver, os.fsencode(nl), options) == 0,
                "AMPLS could not load the NL fixture")
        extras(solver)
        model = native(solver)
        if name == "cplex":
            state, incumbent, reverse = cplex_checks(driver, sdk, solver, model, binary)
        else:
            state, incumbent, reverse = xpress_checks(sdk, model, binary)
        solve(solver)
        original_type = state()
        first_binary = 1 if name == "cplex" else 0
        require(abs(incumbent() - first_binary) < 1e-5,
                "The first solve returned an unexpected incumbent")
        sol = nl.with_suffix(".sol")
        expected = (b"barrier iteration" if name == "cplex" else b"nonlinear solve")

        def check_report(path, binary_value):
            require(report(solver, os.fsencode(path)) == 0, "AMPLS reporting failed")
            require(expected in path.read_bytes(), "The fixed solve did not report success")
            require(native(solver) == model, "Reporting replaced the native model")
            require(state() == original_type, "Reporting changed the native problem type")
            require(abs(incumbent() - binary_value) < 1e-5,
                    "Reporting changed the original incumbent")

        check_report(sol, first_binary)
        require(reverse() == 0, "Could not reverse the objective sense")
        solve(solver)
        second_binary = 1 - first_binary
        require(abs(incumbent() - second_binary) < 1e-5,
                "The second solve did not find the new incumbent")
        check_report(sol.with_name("second.sol"), second_binary)
    finally:
        close(solver)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", choices=("cplex", "xpress"))
    parser.add_argument("driver", type=Path)
    parser.add_argument("sdk", type=Path)
    parser.add_argument("--ampl", default="ampl")
    args = parser.parse_args()
    source = ("mip_basis/fixed_miqcp_cplex.mod" if args.solver == "cplex"
              else "suf_std/nl_dual_01.mod")
    with (os.add_dll_directory(str(args.driver.resolve().parent)),
          os.add_dll_directory(str(args.sdk.resolve().parent)),
          tempfile.TemporaryDirectory() as temp):
        directory = Path(temp)
        (directory / "case.mod").write_bytes((MODEL_DIR / source).read_bytes())
        (directory / "write.run").write_text(
            "option presolve 0; option auxfiles rc; model case.mod; write gcase;\n")
        subprocess.run([args.ampl, "write.run"], cwd=directory, check=True)
        binary = (directory / "case.col").read_text().splitlines().index("b")
        check_repeated_solve(args.solver, ct.CDLL(str(args.driver.resolve())),
                             ct.CDLL(str(args.sdk.resolve())),
                             directory / "case.nl", binary)


if __name__ == "__main__":
    main()

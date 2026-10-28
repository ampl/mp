"""Table of the constraints and expressions that MP-based drivers receive natively.

Runs `<solver> -v` and `<solver> -c` for each driver in DRIVERS that is found on
PATH (drivers in SOURCES that are not on PATH are read from their source code),
saves the results in doc/source/tables/acceptance.json, and writes
doc/source/tables/acceptance.rst, which the "Full constraint list" section of
modeling-tools.rst includes. A driver is updated only when it is on PATH and
its driver date is not older than the saved one, so a run without solvers
(e.g., a CI documentation build) leaves the saved data and the table as they are.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys

LOCAL_PATH = os.path.dirname(os.path.abspath(__file__))
TABLES = os.path.join(LOCAL_PATH, "..", "source", "tables")
DATA = os.path.join(TABLES, "acceptance.json")
TABLE = os.path.join(TABLES, "acceptance.rst")

# Driver executable -> name in the table header
DRIVERS = {
    "gurobi": "Gurobi", "xpress": "Xpress", "cplex": "CPLEX", "copt": "COPT",
    "mosek": "MOSEK", "highs": "HiGHS", "scip": "SCIP", "gcg": "GCG",
    "cbc": "CBC", "cuopt": "cuOpt", "cpsat": "CP-SAT", "baronmp": "BARONMP",
}

# Groups of rows: (title, rows); a row is (label, AMPL form, [MP constraint names]).
# Several names in one row form a family; a cell then shows the members'
# common acceptance, or '*' if they differ.
GROUPS = [
    ("Logical and conditional", [
        ("and", "``and``, ``forall``", ["AndConstraint"]),
        ("or", "``or``, ``exists``", ["OrConstraint"]),
        ("not", "``not``", ["NotConstraint"]),
        ("implication", "``==>`` ... ``else``", ["ImplicationConstraint"]),
        ("equivalence", "``<==>``", ["EquivalenceConstraint"]),
        ("if-then-else", "``if`` ... ``then`` ... ``else``", ["IfThenConstraint"]),
        ("indicator, linear", "binary ``==>`` linear", [
            "IndicatorConstraintLinEQ", "IndicatorConstraintLinGE",
            "IndicatorConstraintLinLE"]),
        ("indicator, quadratic", "binary ``==>`` quadratic", [
            "IndicatorConstraintQuadEQ", "IndicatorConstraintQuadGE",
            "IndicatorConstraintQuadLE"]),
    ]),
    ("Counting and constraint programming", [
        ("alldiff", "``alldiff``", ["AllDiffConstraint"]),
        ("count", "``count``, ``atmost``, ``atleast``, ``exactly``",
         ["CountConstraint"]),
        ("numberof", "``numberof``", [
            "NumberofConstConstraint", "NumberofVarConstraint"]),
    ]),
    ("Piecewise-linear, min, max, abs, SOS", [
        ("piecewise-linear", "``<<breakpoints; slopes>>``", ["PLConstraint"]),
        ("max", "``max``", ["MaxConstraint"]),
        ("min", "``min``", ["MinConstraint"]),
        ("abs", "``abs``", ["AbsConstraint"]),
        ("SOS1, SOS2", "suffixes ``.sosno``, ``.ref``",
         ["SOS1Constraint", "SOS2Constraint"]),
    ]),
    ("Integer arithmetic and rounding", [
        ("product", "``x*y``", ["MulConstraint"]),
        ("div", "``div``", ["IntDivConstraint"]),
        ("mod", "``mod``", ["ModConstraint"]),
        ("floor", "``floor``", ["FloorConstraint"]),
        ("ceil", "``ceil``", ["CeilConstraint"]),
        ("round", "``round``", ["RoundConstraint"]),
    ]),
    ("Quadratic and conic", [
        ("quadratic constraint", "quadratic ``<=``, ``=``, ``>=``",
         ["QuadConLE", "QuadConEQ", "QuadConGE"]),
        ("quadratic range", "``lb <= quadratic <= ub``", ["QuadConRange"]),
        ("quadratic expression", "quadratic term in an expression",
         ["QuadraticFunctionalConstraint"]),
        ("second-order cone", "recognized from algebraic form",
         ["QuadraticConeConstraint"]),
        ("rotated second-order cone", "recognized from algebraic form",
         ["RotatedQuadraticConeConstraint"]),
        ("exponential cone", "recognized from algebraic form",
         ["ExponentialConeConstraint"]),
        ("power cone", "", ["PowerConeConstraint"]),
        ("geometric cone", "", ["GeometricConeConstraint"]),
    ]),
    ("Nonlinear functions", [
        ("exp", "``exp``", ["ExpConstraint"]),
        ("a^x", "constant ``^`` expression", ["ExpAConstraint"]),
        ("log", "``log``", ["LogConstraint"]),
        ("log, base a", "``log10``", ["LogAConstraint"]),
        ("x^a", "expression ``^`` constant", ["PowConstExpConstraint"]),
        ("x^y", "expression ``^`` expression", ["PowConstraint"]),
        ("division", "``/``", ["DivConstraint"]),
        ("sin", "``sin``", ["SinConstraint"]),
        ("cos", "``cos``", ["CosConstraint"]),
        ("tan", "``tan``", ["TanConstraint"]),
        ("asin", "``asin``", ["AsinConstraint"]),
        ("acos", "``acos``", ["AcosConstraint"]),
        ("atan", "``atan``", ["AtanConstraint"]),
        ("sinh", "``sinh``", ["SinhConstraint"]),
        ("cosh", "``cosh``", ["CoshConstraint"]),
        ("tanh", "``tanh``", ["TanhConstraint"]),
        ("asinh", "``asinh``", ["AsinhConstraint"]),
        ("acosh", "``acosh``", ["AcoshConstraint"]),
        ("atanh", "``atanh``", ["AtanhConstraint"]),
        ("logistic", "``1/(1+exp(-x))``", ["LogisticConstraint"]),
        ("signed power", "``sign(x)*abs(x)^a``", ["SignpowConstExpConstraint"]),
        ("user functions", "imported functions", ["CallConstraint"]),
    ]),
    ("Complementarity", [
        ("linear", "``complements``, linear", ["ComplementarityLinear"]),
        ("quadratic", "``complements``, quadratic",
         ["ComplementarityQuadratic"]),
        ("nonlinear", "``complements``, nonlinear", ["NLComplementarity"]),
    ]),
    ("Linear and internal (technical)", [
        ("linear", "linear ``<=``, ``=``, ``>=``",
         ["LinConLE", "LinConEQ", "LinConGE"]),
        ("linear range", "``lb <= linear <= ub``", ["LinConRange"]),
        ("linear expression", "", ["LinearFunctionalConstraint"]),
        ("conditional linear", "comparison used as a logical value", [
            "CondLinConEQ", "CondLinConGE", "CondLinConGT", "CondLinConLE",
            "CondLinConLT"]),
        ("conditional quadratic", "comparison used as a logical value", [
            "CondQuadConEQ", "CondQuadConGE", "CondQuadConGT", "CondQuadConLE",
            "CondQuadConLT"]),
        ("nonlinear constraint", "expression tree", ["NLConstraint"]),
        ("nonlinear assignment", "", ["NLAssignEQ", "NLAssignGE", "NLAssignLE"]),
        ("nonlinear logical", "", ["NLLogical"]),
        ("nonlinear reified", "", ["NLReifEquiv", "NLReifImpl", "NLReifRimpl"]),
        ("unary encoding", "", ["UnaryEncodingConstraint"]),
        ("SDP dot product", "", ["SDPDotProdConstraint"]),
    ]),
]

LEVEL = {"NativeRecommended": 2, "NativeAcceptedButNotRecommended": 1,
         "NotAccepted": 0}

# Drivers that can be read from their source code when not on PATH
# (e.g., they need hardware not available here):
# driver -> (model API header with ACCEPT_CONSTRAINT/ACCEPT_EXPRESSION
#            declarations, driver on PATH whose -c gives convertibility)
SOURCES = {"cuopt": ("cuoptmp/cuoptmpmodelapi.h", "highs")}
SOURCE_LEVEL = {"Recommended": 2, "AcceptedButNotRecommended": 1,
                "NotAccepted": 0}
MP_ROOT = os.path.join(LOCAL_PATH, "..", "..")


def from_source(driver, data):
    """Acceptance of a driver from its source code, convertibility from the
    reference driver's saved data; date: latest MP release in CHANGES.mp.md."""
    header, ref = SOURCES[driver]
    if ref not in data:
        return None, None
    text = open(os.path.join(MP_ROOT, "solvers", header)).read()
    with open(os.path.join(MP_ROOT, "CHANGES.mp.md")) as f:
        m = re.search(r"^## (\d{8})", f.read(), re.M)
    date = m.group(1) if m else None
    cons = {}
    for name, e in data[ref]["cons"].items():
        cons[name] = [e[0], 0, 0, e[3]]
    for kind, name, level in re.findall(
            r"ACCEPT_(CONSTRAINT|EXPRESSION)\(\s*(\w+)\s*,\s*(\w+)", text):
        if name in cons:
            cons[name][1 if kind == "CONSTRAINT" else 2] = SOURCE_LEVEL.get(level, 0)
    return date, cons


def driver_date(exe):
    try:
        out = subprocess.run([exe, "-v"], capture_output=True, text=True,
                             timeout=60, stdin=subprocess.DEVNULL).stdout
    except (OSError, subprocess.TimeoutExpired):
        return None, None
    m = re.search(r"driver\((\d{8})\)", out)
    return (m.group(1) if m else None), out.strip().splitlines()[0] if out else ""


def acceptance(exe):
    """{constraint name: [convertible, as constraint, as expression, options]}"""
    try:
        out = subprocess.run([exe, "-c"], capture_output=True, text=True,
                             timeout=60, stdin=subprocess.DEVNULL).stdout
    except (OSError, subprocess.TimeoutExpired):
        return None
    result = {}
    for line in out.splitlines():
        m = re.match(r"\s*\d+\.\s+(\w+)\s*:\s*(\w+);\s*(\w+);\s*(\w+);?\s*(.*)$",
                     line)
        if m:
            name, conv, con, expr, opts = m.groups()
            result[name] = [conv == "Convertible", LEVEL.get(con, 0),
                            LEVEL.get(expr, 0), opts.split()]
    return result or None


def cell(entries):
    """The level reported by `<solver> -c`: C/E - accepted natively as a
    constraint/expression, lower case if MP's reformulation is preferred
    by default; '-' - not accepted, MP converts it; 'x' - neither accepted
    nor convertible; '' - type unknown to this (older) driver;
    '*' - the members of the row differ."""
    def one(e):
        if e is None:
            return ""
        conv, con, expr, _ = e
        s = ("C" if con == 2 else "c" if con == 1 else "") + \
            ("E" if expr == 2 else "e" if expr == 1 else "")
        return s or ("–" if conv else "×")
    vals = {one(e) for e in entries}
    return vals.pop() if len(vals) == 1 else "\\*"  # escaped for reST


def write_table(data):
    solvers = [s for s in DRIVERS if s in data]
    head = [f"{DRIVERS[s]} {data[s]['date']}" + (" (source)" if data[s].get("source")
                                                 else "") for s in solvers]
    lines = [
        ".. This file is written by doc/scripts/acceptanceTable.py; do not edit.",
        "",
    ]
    for title, rows in GROUPS:
        lines += [f".. dropdown:: {title}", "",
                  "   .. csv-table::", "      :class: solver-support",
                  "      :header-rows: 1", ""]
        lines.append("      " + ",".join(
            ['"Construct (option)"'] + [f'"{h}"' for h in head]))
        for label, form, names in rows:
            cells = []
            for s in solvers:
                cells.append(cell([data[s]["cons"].get(n) for n in names]))
            opt = ""
            for s in solvers:
                e = data[s]["cons"].get(names[0])
                if e and e[3]:
                    opt = "``" + e[3][0] + "``"
                    break
            name = form if "``" in form else label
            row = [name + (f" ({opt})" if opt else "")] + cells
            lines.append("      " + ",".join(
                '"' + c.replace('"', '""') + '"' for c in row))
        lines.append("")
    with open(TABLE, "w") as f:
        f.write("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python acceptanceTable.py              # all drivers on PATH\n"
               "  PATH=~/mp/build/bin:$PATH python acceptanceTable.py --solvers gurobi,highs")
    parser.add_argument(
        "--solvers",
        help="comma-separated driver names to update (default: all in DRIVERS); "
             "the others keep their saved data. Example: --solvers gurobi,cpsat")
    args = parser.parse_args()
    selected = args.solvers.split(",") if args.solvers else list(DRIVERS)

    data = {}
    if os.path.exists(DATA):
        with open(DATA) as f:
            data = json.load(f)
    changed = False
    for s in selected:
        exe = shutil.which(s)
        if not exe and s in SOURCES:
            date, cons = from_source(s, data)
            if not date or not cons:
                print(f"{s}: not on PATH, no source data, kept")
                continue
            if s in data and date < data[s]["date"]:
                print(f"{s}: source {date} older than saved {data[s]['date']}, kept")
                continue
            entry = {"date": date, "banner": "from source " + SOURCES[s][0],
                     "source": True, "cons": cons}
            if data.get(s) != entry:
                data[s] = entry
                changed = True
                print(f"{s}: updated from source ({date})")
            continue
        if not exe:
            print(f"{s}: not on PATH, kept")
            continue
        date, banner = driver_date(exe)
        if not date:
            print(f"{s}: no driver date in '{s} -v', kept")
            continue
        if s in data and date < data[s]["date"]:
            print(f"{s}: driver {date} older than saved {data[s]['date']}, kept")
            continue
        cons = acceptance(exe)
        if cons is None:
            print(f"{s}: '{s} -c' gave no table, kept")
            continue
        if data.get(s) != {"date": date, "banner": banner, "cons": cons}:
            data[s] = {"date": date, "banner": banner, "cons": cons}
            changed = True
            print(f"{s}: updated from {exe} ({date})")
    if changed:
        with open(DATA, "w") as f:
            json.dump({s: data[s] for s in DRIVERS if s in data}, f, indent=1,
                      sort_keys=False)
    known = {n for _, rows in GROUPS for _, _, names in rows for n in names}
    for s in data:
        missing = sorted(set(data[s]["cons"]) - known)
        if missing:
            print(f"{s}: not in GROUPS, not shown: {', '.join(missing)}")
    if data:
        write_table(data)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:  # never fail a documentation build
        print(f"acceptanceTable.py: {e}", file=sys.stderr)

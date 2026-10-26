import enum
import re
import subprocess


class SolverCategory(enum.Enum):
    """Which logical and combinatorial constraints a driver accepts
    natively, i.e., which ones its model converter does not reformulate.

    Converter statistics (e.g., 'Initial.flat1__U__indeq' with debug=1)
    depend on it, so a test case can give expected values per category:
    'values_per_solver_category' in modellist.json.
    Nonlinear functions do not matter for the category."""
    MILP = "MILP"                           # linear constraints only
    MIPWithIndicators = "MIPWithIndicators" # + indicator constraints
    MILPWithLogicals = "MILPWithLogicals"   # + and, or
    CP = "CP"                       # + count, numberof, alldiff, if-then...
    CPLegacy = "CPLegacy"           # old MP framework: no model converter,
                                    # no converter statistics

    @staticmethod
    def fromAccepted(acc: set):
        """Category of a driver accepting the constraints 'acc'
        (names of its acc:... options without the prefix)"""
        if acc & {"alldiff", "count", "numberofconst", "numberofvar"}:
            return SolverCategory.CP
        if "indle" in acc and {"and", "or"} <= acc:
            return SolverCategory.MILPWithLogicals
        if "indle" in acc:
            return SolverCategory.MIPWithIndicators
        return SolverCategory.MILP

    @staticmethod
    def ofDriver(exePath: str):
        """Category implied by the driver's acc:... options
        ('<driver> -=acc:'), or None if it lists none (old MP framework,
        ASL drivers) or cannot be run"""
        try:
            out = subprocess.run([exePath, "-=acc:"], capture_output=True,
                                 text=True, timeout=30).stdout
        except (OSError, subprocess.SubprocessError):
            return None
        acc = set(re.findall(r"^acc:(\w+)", out, re.MULTILINE)) - {"_all"}
        return SolverCategory.fromAccepted(acc) if acc else None

from runModels import runModels, writeNLFiles
from sys import platform
import Solver

# Set some parameters
if platform == "win32":
    bin = "/build/vs64/bin/debug"
else:
    bin = "/../../mp/build/bin"
solversbin = "../.." + bin
timeout = 5
nthreads = 8
modDir = ""

# Write NL files
# writeNLFiles("../../test/models/lindo")

# Create two solver objects
# s = Solver.LindoSolver(solversbin + "/lindoglobal-timebound", timeout, nthreads)
# o = Solver.OcteractSolver("C:/Program Files (x86)/Octeract/bin/octeract-engine", timeout, nthreads)
g = Solver.GurobiDirectSolver("gurobidirect", timeout, nthreads)
# b = Solver.BaronSolver(solversbin +"/baron-timebound", timeout, nthreads)
# Execute a comparison exporting to CSV in the current directory

runModels(modDir, g, recursive=True)


# Example for defining a new solver to be used
class UserDefSolver(Solver.AMPLSolver):
    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None):
        super().__init__(exeName, timeout, nthreads, otherOptions)

    def _setTimeLimit(self, seconds):
        # syntax to set the timelimit in the solver
        return "time={}".format(seconds)

    def _setNThreads(self, threads):
        # syntax to set the number of threads in the solver
        return "threads={}".format(threads)

    def _doParseSolution(self, st, stdout=None):
        # This function should set the following values in the self.stats dictionaty:
        # self.stats["outmsg"] to the solver message
        # self.stats["timelimit"] to True or False indicating whether the time limit was reached before optimal
        # self.stats["solution"] to the objective value
        if not st:
            self.stats["outmsg"] = "Solution file empty"
            self.stats["timelimit"] = False
            return None
        self.stats["outmsg"] = st[0]
        self.stats["timelimit"] = "time limit" in st[0]
        tag = "objective "
        if tag in st[0]:
            n = st[0][st[0].index(tag) + len(tag):]
            try:
                self.stats["solution"] = float(n)
            except:
                print("No solution, string: {}".format(n))
                self.stats["solution"] = None

# Creates the user defined solver instance, pointing to the executable
# and setting some solver option (passed to the solver in addition to the
# timeout and nthreads options

# uds = UserDefSolver(solversbin + "/cplex", timeout, nthreads, "mipalg=4")
# runModels(modDir, uds)

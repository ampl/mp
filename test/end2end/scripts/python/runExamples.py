from sys import platform
import argparse

from runModels import runModels
import Solver
import SolverCollection


class Tester:
    def __init__(self):
        self._parser = argparse.ArgumentParser(description='AMPL solver testing script.')

    def runTests(self):
        self.parseOptions()
        self.initSolvers()
        if self._args.printsolvers:
            self.printSolvers()
            return
        self.collectAndRunCases()

    def parseOptions(self):
        self._parser.add_argument('solvers', metavar='solver', type=str, nargs='+',
                            help='a solver to test')

        self._parser.add_argument('--lpmethod', type=str, metavar='', default="",
                            help='algorithm for lp: SIMPLEX or BARRIER')

        self._parser.add_argument('--options', type=str, metavar='', default="",
                            help='extra solver options')
        self._parser.add_argument('--binPath', type=str, metavar='', default="",
                            help='default path to look for solver executables')
        self._parser.add_argument('--reportstub', type=str, metavar='', default="report",
                            help='stub for CSV test report filename, e.g., /tmp/report, default: report')
        self._parser.add_argument('--printsolvers', action="store_true",
                            help='print available solvers and exit')
        self._parser.add_argument('--timeout', type=int, metavar='T', default=1200,
                        help='timeout per instance, seconds')
        self._parser.add_argument('--nthreads', type=int, metavar='N', default=8,
                        help='number of threads in a solver')
        self._parser.add_argument('--dir', type=str, metavar='path', default="",
                        help='path to the test case folder')
        self._parser.add_argument('--nonrecursive', action="store_true",
                        help='non-recursive case collection')
        self._parser.add_argument('--allfiles', action="store_true",
                        help='collect all .mod, .nl files; otherwise local modellist.json only. ' +
                                  'If modellist.json is to be used for comparison data, each case name' +
                             '\'s first word should match the file stem unless "files" are specified')
        self._parser.add_argument('--preferNL', action="store_true",
                        help='prefer NL models if both AMPL and NL are present')
        self._parser.add_argument('--justNL', action="store_true",
                        help='only run NL models. Useful when AMPL executable is not available')
        self._parser.add_argument('-k', '--keepLogs', action="store_true",
                        help='keep log files for each run, filenames: {model}.{solver}.log')

        self._args = self._parser.parse_args()

    def initSolvers(self):
        self._solvers = SolverCollection.SolverCollection()
        SolverCollection.addStdSolvers(self._solvers, self._args.binPath)
        self.setSolverParameters()

    def setSolverParameters(self):
        for name, slv in self._solvers.getSolvers():
            slv.setTimeout(self._args.timeout)
            slv.setNThreads(self._args.nthreads)
            if self._args.lpmethod:
                slv.setLPMethod(self._args.lpmethod)

    def printSolvers(self):
        print("Available solvers:\n  * ", end='')
        print(*(self._solvers.getSolverNames()), sep="\n  * ")

    def collectAndRunCases(self):
        runModels(self._args.dir,
                  self._solvers.getSolversByNames(self._args.solvers),
                  solverOptions=self._args.options,
                  exportFile=self._args.reportstub,
                  recursive=not self._args.nonrecursive,
                  modellist=not self._args.allfiles,
                  preferAMPLModels=not self._args.preferNL,
                  justNL=self._args.justNL,
                  keepLogs=self._args.keepLogs)


def runTester():
    tester = Tester()
    tester.runTests()

if __name__ == "__main__":
    runTester()

# Write NL files
# writeNLFiles("../../test/models/lindo")

# Create two solver objects
# s = Solver.LindoSolver(solversbin + "/lindoglobal-timebound", timeout, nthreads)
# o = Solver.OcteractSolver("C:/Program Files (x86)/Octeract/bin/octeract-engine", timeout, nthreads)
# g = Solver.GurobiDirectSolver("gurobidirect", timeout, nthreads)
# b = Solver.BaronSolver(solversbin +"/baron-timebound", timeout, nthreads)
# Execute a comparison exporting to CSV in the current directory


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

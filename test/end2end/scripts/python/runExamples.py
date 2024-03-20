from sys import platform
import argparse

from runModels import runModels
import Exporter
import Solver
import SolverCollection


class Tester:
    def runTestsAPI(self, solvers: list, lpmethod: str = None, nlpmethod: str = None,
                    options: str = None, bin_path: str= "", reportstub: str=None,printsolvers:bool = False,
                    timeout: int = 2400, nthreads: int = 8, dir: str="", 
                    benchmark: bool = False, junit: bool=False, nonrecursive: bool=False,
                    allfiles: bool=False, prefer_nl: bool = False, export_lp: bool = False,
                    just_nl: bool = False, keep_logs: bool = False, verbose: bool = False,
                    exporter = None):
        self.initSolvers(timeout, nthreads, bin_path, lpmethod, nlpmethod, export_lp)
        if printsolvers:
            self.printSolvers()
            return
        self.collectAndRunCases(solvers, lpmethod, nlpmethod, options, bin_path, reportstub, 
                                timeout, nthreads, dir,
                                benchmark, junit, nonrecursive,
                                allfiles, prefer_nl, export_lp, 
                                just_nl, keep_logs, verbose, exporter)
    def run_from_command_line(self):
        args = self.parseOptions()
        args= vars(args)
        self.runTestsAPI(args["solvers"],
                         args["lpmethod"],
                         args["nlpmethod"],
                         args["options"],
                         args["bin_path"],
                         args["reportstub"],
                         args["printsolvers"],
                         args["timeout"],
                         args["nthreads"],
                         args["dir"],
                         args["benchmark"],
                         args["junit"],
                         args["nonrecursive"],
                         args["allfiles"],
                         args["prefer_nl"],
                         args["export_lp"],
                         args["just_nl"],
                         args["keep_logs"],
                         args["verbose"])
        
    def parseOptions(self):
        parser = argparse.ArgumentParser(description='AMPL solver testing script.')
        parser.add_argument("solvers", metavar="solver", type=str, nargs="+",
                            help="a solver to test")

        parser.add_argument("--lpmethod", type=str, metavar="", default="",
                            help="algorithm for lp: SIMPLEX or BARRIER")
        parser.add_argument("--nlpmethod", type=str, metavar="", default="",
                            help="nl support: REFORMULATION, NATIVE or NATIVEPL")
        parser.add_argument("--options", type=str, metavar="", default="",
                            help="extra solver options")
        parser.add_argument("--bin_path", type=str, metavar="", default="",
                            help="default path to look for solver executables")
        parser.add_argument("--reportstub", type=str, metavar="", default="report",
                            help="stub for CSV test report filename, e.g., /tmp/report, default: report")
        parser.add_argument("--printsolvers", action="store_true",
                            help="print available solvers and exit")
        parser.add_argument("--timeout", type=int, metavar="T", default=2400,
                        help="timeout per instance, seconds")
        parser.add_argument("--nthreads", type=int, metavar="N", default=8,
                        help="number of threads in a solver")
        parser.add_argument("--dir", type=str, metavar="path", default="",
                        help="path to the test case folder")
        parser.add_argument("--benchmark", action="store_true",
                        help="benchmark file output")
        parser.add_argument("--junit", action="store_true",
                help="Junit test file output")
        parser.add_argument("--nonrecursive", action="store_true",
                        help="non-recursive case collection")
        parser.add_argument("--allfiles", action="store_true",
                        help="collect all .mod, .nl files; otherwise local modellist.json only. " +
                                  "If modellist.json is to be used for comparison data, each case name" +
                             "'s first word should match the file stem unless 'files' are specified")
        parser.add_argument("--prefer_nl", action="store_true",
                        help="prefer NL models if both AMPL and NL are present")
        parser.add_argument("--export_lp", action="store_true",
                                  help="Export solver-level LP file")
        parser.add_argument("--just_nl", action="store_true",
                        help="only run NL models. Useful when AMPL executable is not available")
        parser.add_argument("-k", "--keep_logs", action="store_true",
                        help="keep log files for each run, filenames: {model}.{solver}.log")
        parser.add_argument("-v", "--verbose", action="store_true",
                        help="print solver output")

        return parser.parse_args()

    def initSolvers(self,  timeout: int, nthreads: int, bin_path:str = None, 
                    lpmethod: str = None, nlpmethod: str = None,
                    export_lp: bool = False):
        self._solvers = SolverCollection.SolverCollection()
        SolverCollection.addStdSolvers(self._solvers, bin_path)
        for name, slv in self._solvers.getSolvers():
            slv.setTimeout(timeout)
            slv.setNThreads(nthreads)
            if lpmethod:
                slv.setLPMethod(lpmethod)
            if nlpmethod:
                slv.setNLPMethod(nlpmethod)
            slv.setExportLP(export_lp)
            
    def printSolvers(self):
        print("Available solvers:\n  * ", end='')
        print(*(self._solvers.getSolverNames()), sep="\n  * ")
        
    def create_report_file_suffix(self, lpmethod, nlpmethod):
        suffixes = []
        if lpmethod:
            suffixes.append(lpmethod)
        if nlpmethod:
            suffixes.append(nlpmethod) 
        if len(suffixes)==0:
            return ""
        else: return "-".join(suffixes)
        
    def collectAndRunCases(self,
                    solvers: list, lpmethod: str = None, nlpmethod: str = None,
                    options: str = None, binPath: str= None, reportstub: str=None,
                    timeout: int = 2400, nthreads: int = 8, dir: str="", 
                    benchmark: bool = False, junit: bool=False, nonrecursive: bool=False,
                    allfiles: bool=False, prefer_nl: bool = False, exportLP: bool = False,
                    just_nl: bool = False, keep_logs: bool = False, verbose: bool = False,
                    exporter: Exporter.Exporter=None):
        
        if not exporter:
            if benchmark:
                from BenchmarkExporter import BenchmarkExporter
                exporter = BenchmarkExporter()
            elif junit:
                from JunitExporter import JunitExporter
                exporter = JunitExporter()
            else:
                exporter = Exporter.CSVTestExporter()

        runModels(dir,
                  self._solvers.getSolversByNames(solvers),
                  solverOptions=options,
                  exportFile=reportstub,
                  exporter=exporter,
                  recursive=not nonrecursive,
                  modellist=not allfiles,
                  preferAMPLModels=not prefer_nl,
                  justNL=just_nl,
                  keepLogs=keep_logs,
                  verbose=verbose, 
                  defaultReportSuffix=self.create_report_file_suffix(lpmethod, nlpmethod))

def runTester():
    tester = Tester()
    tester.run_from_command_line()
    

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

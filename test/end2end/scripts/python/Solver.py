import sys
import math
from pathlib import PurePath
import subprocess

from Model import Model, ModelTags
from TimeMe import TimeMe


class Solver(object):

    """description of class"""
    @staticmethod
    def getExecutableName(name):
        path = PurePath(name)
        if sys.platform == "win32":
            if path.suffix:
                return name
            else:
                return str(path.with_suffix(".exe"))
        else:
            return name

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None,
                 writeSolverName=False, unsupportedTags=None):
        self._exePath = Solver.getExecutableName(exeName)
        self._timeout = timeout
        self._nthreads = nthreads
        self._otherOptions = otherOptions
        self._writeSolverName = writeSolverName
        self._unsupportedTags = unsupportedTags


    def _doRun(self,  model: Model):
        """Method to be overriden for solvers implementations"""
        raise Exception("Not implemented")

    def _doParseSolution(self, st, stdout):
        """Method to be overriden for solvers implementations"""
        raise Exception("Not implemented")

    def _getSolution(self, model, stdout=None):
        self._stats["outmsg"] = "No solution file"
        self._stats["timelimit"] = False
        self._stats["solution"] = None
        sol = model.getSolFilePath()
        st = None
        if sol.exists():
            st = sol.read_text().splitlines()
        return self._doParseSolution(st, stdout)

    def _evaluateRun(self, model: Model):
        expsol = model.getExpectedSolution()
        if expsol is not None:
            self._stats["eval_done"] = True
            self._assertAndRecord(expsol, self._stats["solution"],
                                  "objective")

    def _assertAndRecord(self, expval, val, msg):
        b1 = isinstance(expval, (int, float))
        b2 = isinstance(val, (int, float))
        uneq = not math.isclose(expval, val, rel_tol=1e-6) if \
            b1 and b2 else expval != val
        if uneq:
            self._stats["eval_fail_msg"] = msg + \
                ": value " + str(val) + \
                ", expected " + str(expval)

    def runAndEvaluate(self, model: Model):
        t = TimeMe()
        self._stats = { "solver": self.getExecutable() }
        sol = model.getSolFilePath()
        if sol.exists():
            sol.unlink()
        with t:
            stdout = self._doRun(model)
        self._stats["solutionTime"] = t.interval
        self._getSolution(model, stdout)
        self._evaluateRun(model)

    def getName(self):
        path = PurePath(self._exePath)
        return path.stem

    def getExecutable(self):
        return self._exePath

    def getSolutionStats(self):
        return self._stats

    def setNThreads(self, nt):
        self._nthreads = nt

    def getNThreads(self):
        return self._nthreads

    def setTimeout(self, t):
        self._timeout = t

    def getTimeout(self):
        return self._timeout

    def getUnsupportedTags(self):
        return self._unsupportedTags

class AMPLSolver(Solver):
    """ Main class to inherit from when adding a new solver.

    The functions setTimeLimit and setNThreads must be implemented in all cases.
    The function _doParseSolution only if the solver is executed directly on the nl file,
    because if the models are executed through AMPL, the latter takes care of reporting 
    the statistics back.
    """

    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None, unsupportedTags=None):
        super().__init__(exeName, timeout, nthreads, otherOptions, unsupportedTags=unsupportedTags)

    def _setTimeLimit(self, seconds):
        raise Exception("Not implemented in base class")

    def _setNThreads(self, nthreads):
        raise Exception("Not implemented in base class")

    def _getAMPLOptionsName(self):
        raise Exception("Not implemented in base class")

    def _doParseSolution(self, st, stdout=None):
        raise Exception("Not implemented in base class")

    def _doRun(self,  model: Model,):
        toption = ""
        if self._timeout:
            try:
                toption = self._setTimeLimit(self._timeout)
            except:
                pass
        if self._nthreads:
            toption = "{} {}".format(toption, self._setNThreads(self._nthreads))
        if self._otherOptions:
            toption = "{} {}".format(toption, self._otherOptions)
        try:
            if toption:
                return subprocess.check_output([self._exePath, model.getFilePath(), "-AMPL", toption])
            else:
                return subprocess.check_output([self._exePath, model.getFilePath(), "-AMPL"],
                                               timeout=self._timeout)
        except subprocess.TimeoutExpired:
            pass
        except subprocess.CalledProcessError as e:
            print(str(e))

    def getAMPLOptions(self):
        name = "{}_options".format(self._getAMPLOptionsName())
        value = ""
        if self._timeout:
            value += self._setTimeLimit(self._timeout)
        if self._nthreads:
            value += " "
            value += self._setNThreads(self._nthreads)
        if self._otherOptions:
            value += " "
            value += self._otherOptions
        if value:
            return (name, value)



class LindoSolver(AMPLSolver):
    def _setTimeLimit(self, seconds):
        return "maxtime={}".format(seconds)

    def _setNThreads(self, threads):
        return "threads={}".format(threads)

    def _getAMPLOptionsName(self):
        return "lindoglobal"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        super().__init__(exeName, timeout, nthreads, otherOptions)

    def _doParseSolution(self, st, stdout=None):
        if st:
            tag = "OBJECTIVE VALUE:"
            for line in st:
                if "LOCAL" in line:
                    self._stats["timelimit"] = True
                if line.startswith(tag):
                    n = line[len(tag):]
                    self._stats["outmsg"] = prev
                    self._stats["solution"] = float(n)
                    return
                prev = line
        self._stats["outmsg"] = stdout


class GurobiSolver(AMPLSolver):

    def _setTimeLimit(self, seconds):
        return "timelim={}".format(seconds)

    def _setNThreads(self, threads):
        return "threads={}".format(threads)

    def _getAMPLOptionsName(self):
        return "gurobi"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        utags = [ModelTags.nonlinear, ModelTags.complementarity]
        super().__init__(exeName, timeout, nthreads, otherOptions,utags)

    def _doParseSolution(self, st, stdout=None):
        if not st:
            self._stats["outmsg"] = "Solution file empty"
            self._stats["timelimit"] = False
            return None
        self._stats["outmsg"] = st[0]
        self._stats["timelimit"] = "time limit" in st[0]
        tag = "objective "
        if tag in st[0]:
            n = st[0][st[0].index(tag) + len(tag):]
            try:
                self._stats["solution"] = float(n)
            except:
                print("No solution, string: {}".format(n))
                self._stats["solution"] = None


class GurobiDirectSolver(AMPLSolver):

    def _setTimeLimit(self, seconds):
        return "timelim={}".format(seconds)

    def _setNThreads(self, threads):
        return "threads={}".format(threads)

    def _getAMPLOptionsName(self):
        return "gurobidirect"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        utags = [ModelTags.nonlinear, ModelTags.complementarity]
        super().__init__(exeName, timeout, nthreads, otherOptions,utags)

    def _doParseSolution(self, st, stdout=None):
        if not st:
            self._stats["outmsg"] = "Solution file empty"
            self._stats["timelimit"] = False
            return None
        self._stats["outmsg"] = st[0]
        self._stats["timelimit"] = "time limit" in st[0]
        tag = "objective "
        if tag in st[0]:
            n = st[0][st[0].index(tag) + len(tag):]
            try:
                self._stats["solution"] = float(n)
            except:
                print("No solution, string: {}".format(n))
                self._stats["solution"] = None


class CPLEXDirectSolver(GurobiDirectSolver):
    def _getAMPLOptionsName(self):
        return "cplexdirect"


class CPLEXSolver(AMPLSolver):
    def _setTimeLimit(self, seconds):
        return "time={}".format(seconds)

    def _setNThreads(self, threads):
        return "threads={}".format(threads)

    def _getAMPLOptionsName(self):
        return "cplex"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        utags = [ModelTags.nonlinear, ModelTags.complementarity]
        super().__init__(exeName, timeout, nthreads, otherOptions,utags)

    def _doParseSolution(self, st, stdout=None):
        if not st:
            self._stats["outmsg"] = "Solution file empty"
            self._stats["timelimit"] = False
            return None
        self._stats["outmsg"] = st[0]
        self._stats["timelimit"] = "time limit" in st[0]
        tag = "objective "
        if tag in st[0]:
            n = st[0][st[0].index(tag) + len(tag):]
            try:
                self._stats["solution"] = float(n)
            except:
                print("No solution, string: {}".format(n))
                self._stats["solution"] = None


class BaronSolver(AMPLSolver):
    def _setTimeLimit(self, seconds):
        return "maxtime={}".format(seconds)

    def _setNThreads(self, threads):
        return "threads={}".format(threads)

    def _getAMPLOptionsName(self):
        return "baron"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        super().__init__(exeName, timeout, nthreads, otherOptions)

    def _doParseSolution(self, st, stdout=None):
        if not st:
            self._stats["outmsg"] = "Solution file empty"
            self._stats["timelimit"] = False
            return None
        self._stats["outmsg"] = st[0]
        self._stats["timelimit"] = "time limit" in st[0]
        tag = "Objective "
        if tag in st[0]:
            n = st[0][st[0].index(tag) + len(tag):]
            try:
              self._stats["solution"] = float(n)
            except:
              print("No solution, string: {}".format(n))
              self._stats["solution"] = None

class OcteractSolver(AMPLSolver):
    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        super().__init__(exeName, timeout, nthreads, otherOptions)

    def _doRun(self,  model: Model):
        optionFile = model.getSolFilePath().parent.joinpath("octeract.opt")
        if optionFile.exists():
            optionFile.unlink()
        if self._timeout:
            self._writeOptionFile(str(optionFile), self._timeout)
        if not self._nthreads:
            self._nthreads = 1
        try:
            return subprocess.check_output([self._exePath, str(self._nthreads), model.getFilePath(), "-AMPL",
                                            "-o", str(optionFile)], encoding="utf-8",  stderr=subprocess.STDOUT)
        except:
            pass
    def _getAMPLOptionsName(self):
        return "octeract"

    def _setTimeLimit(self, seconds):
        return "max_time={}".format(seconds)

    def _setNThreads(self, threads):
        return "num_cores={}".format(threads)


    def _writeOptionFile(self, file, timeout):
        with open(file, "w") as f:
            f.write("MAX_SOLVER_TIME={}\n".format(timeout))

    def _doParseSolution(self, st, stdout):
        if not st:
            return
        self._stats["outmsg"] = st[0]

        if not stdout:
            return
        for l in stdout.splitlines():
            if "Objective value" in l:
                tag = "Objective value at best solution found:"
                self._stats["timelimit"] = tag in l
                if not self._stats["timelimit"]:
                    tag = "Objective value at global solution:"
                n = l[l.index(tag) + len(tag):]
                self._stats["solution"] = float(n)
                return

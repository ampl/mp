import sys
import math
from pathlib import PurePath
import subprocess

from Model import Model, ModelTags
from TimeMe import TimeMe
import psutil
import time
from threading import Timer  

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
                 writeSolverName=False,
                 supportedTags=None,
                 unsupportedTags=None,
                 lpmethod = None):
        self._exePath = Solver.getExecutableName(exeName)
        self._timeout = timeout
        self._nthreads = nthreads
        if lpmethod:
            if not lpmethod in ["BARRIER", "SIMPLEX"]:
                raise Exception("Valid methods to be forced are BARRIER and SIMPLEX")
        self._lpmethod = lpmethod
        self._otherOptions = otherOptions
        self._writeSolverName = writeSolverName
        self._supportedTags = supportedTags
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
        self._stats["objective"] = None
        sol = model.getSolFilePath()
        st = None
        if sol.exists():
            st = sol.read_text().splitlines()
        return self._doParseSolution(st, stdout)

    def _evaluateRun(self, model: Model):
        expsol = model.getExpectedObjective()
        if expsol is not None:
            self._stats["eval_done"] = True
            self._assertAndRecord(expsol, self._stats["objective"],
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

    def runAndEvaluate(self, model: Model, logFile : str = None):
        t = TimeMe()
        self._stats = { "solver": self.getName() }
        sol = model.getSolFilePath()
        if sol.exists():
            sol.unlink()
        with t:
            stdout = self._doRun(model, logFile = logFile)
            if logFile is not None:
                with open(logFile, 'w') as f:
                    f.write(stdout)
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

    def setLPMethod(self, lpmethod):
        if not lpmethod in ["BARRIER", "SIMPLEX"]:
                raise Exception("Valid methods to be forced are BARRIER and SIMPLEX")
        self._lpmethod = lpmethod

    def getNThreads(self):
        return self._nthreads

    def setTimeout(self, t):
        self._timeout = t

    def getTimeout(self):
        return self._timeout

    def getSupportedTags(self):
        return self._supportedTags

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
                 otherOptions=None,
                 supportedTags=None, unsupportedTags=None):
        super().__init__(exeName, timeout, nthreads, otherOptions,
                         supportedTags=supportedTags,
                         unsupportedTags=unsupportedTags)

    def _setLPMethod(self, method : str):
        raise Exception("Not implemented in base class")

    def _setTimeLimit(self, seconds):
        raise Exception("Not implemented in base class")

    def _setNThreads(self, nthreads):
        raise Exception("Not implemented in base class")

    def _getAMPLOptionsName(self):
        raise Exception("Not implemented in base class")

    def setLogFile(self, name):
        return None

    def _doParseSolution(self, st, stdout=None):
        raise Exception("Not implemented in base class")

    def _doRun(self,  model: Model, logFile : str):
        toption = ""
        if self._timeout:
            try:
                toption = self._setTimeLimit(self._timeout)
            except:
                pass
        if self._nthreads:
            toption = "{} {}".format(toption, self._setNThreads(self._nthreads))
        if self._lpmethod:
            toption = "{} {}".format(toption, self._setLPMethod(self._lpmethod))
        if logFile is not None:
            toption = "{} {} outlev=1".format(toption, self.setLogFile(logFile))
        if self._otherOptions:
            toption = "{} {}".format(toption, self._otherOptions)
        try:
            if toption:
                return self._runProcess([self._exePath, model.getFilePath(), "-AMPL", toption],
                                        timeout=self._timeout, logFile = logFile)
                
            else:
                return self._runProcess([self._exePath, model.getFilePath(), "-AMPL"],
                                               timeout=self._timeout, logFile = logFile)
        except subprocess.TimeoutExpired:
            pass
        except subprocess.CalledProcessError as e:
            print(str(e))

    def stopProcess(p):
        # forcefully terminate solvers that do not terminate automatically 
        # after timeout
        p.terminate()

    def _runProcess(self, args : list, vestigial=False, timeout=None, logFile = None):
      # ritorna stdout
      # throws if not successfull
         if vestigial:
           if timeout:
              return subprocess.check_output(args, timeout=self._timeout)
           else:
              return subprocess.check_output(args)
         else:
              resultTable = []
              SLICE_IN_SECONDS = 1
              p = subprocess.Popen(args, universal_newlines=True, stdout=subprocess.PIPE)
              if self._timeout:
                  t = Timer(self._timeout+5, AMPLSolver.stopProcess, [p])
                  t.start()
              try:
                ps = psutil.Process(p.pid)
              except:
                ps = None
              while p.poll() == None:
                if ps:
                  try:
                    resultTable.append(ps.memory_info())
                  except:
                    pass
                  time.sleep(SLICE_IN_SECONDS)
              (out,err) = p.communicate()
              if logFile is not None:
                if self.setLogFile(logFile) is None: # if not handled via option
                    print(out, flush=True)
                    print(err, flush=True)
              rss = max([p.rss for p in resultTable])
              vms = max([p.vms for p in resultTable])
              self._stats["rss"]= rss
              self._stats["vms"]= vms
              return out


    def getAMPLOptions(self):
        name = "{}_options".format(self._getAMPLOptionsName())
        value = ""
        if self._timeout:
            value += self._setTimeLimit(self._timeout)
        if self._nthreads:
            value += " "
            value += self._setNThreads(self._nthreads)
        if self._lpmethod:
            value += " "
            value += self._setLPMethod(self._lpmethod)
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

    def _setLPMethod(self, method : str):
        return ""

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        super().__init__(exeName, timeout, nthreads, otherOptions)

    def _doParseSolution(self, st, stdout=None):
        if st:
            tag = "OBJECTIVE VALUE:"
            prev = ""
            for line in st:
                if "LOCAL" in line:
                    self._stats["timelimit"] = True
                if line.startswith(tag):
                    n = line[len(tag):]
                    self._stats["outmsg"] = prev
                    self._stats["objective"] = float(n)
                    return
                prev = line
        self._stats["outmsg"] = stdout


class GurobiSolver(AMPLSolver):

    def setLogFile(self, name):
        return f"logfile=\"{name}\""

    def _setTimeLimit(self, seconds):
        return "timelim={}".format(seconds)

    def _setNThreads(self, threads):
        return "threads={}".format(threads)

    def _setLPMethod(self, method : str):
        m  = "1" if method == "SIMPLEX" else "2"
        return f"lpmethod {m}"

    def _getAMPLOptionsName(self):
        return "gurobi"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {
                 ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.linear, ModelTags.quadratic}
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

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
                self._stats["objective"] = float(n)
            except:
                print("No solution, string: {}".format(n))
                self._stats["objective"] = None


# Mp Direct / FlatConverter drivers
class MPDirectSolver(AMPLSolver):
    def _setLPMethod(self, method : str):
        # typically have to reimplement this
        m  = "0" if method == "SIMPLEX" else "2"
        return f"alg:method {m}"

    def _setTimeLimit(self, seconds):
        return "timelim={}".format(seconds)

    def _setNThreads(self, threads):
        return "threads={}".format(threads)
    
    def setLogFile(self, name):
        return f"logfile=\"{name}\""

    def _getAMPLOptionsName(self):
        raise Exception("Not implemented in base class")

    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None, stags=None):
        sDefault = {ModelTags.continuous,
                    ModelTags.linear,
                    ModelTags.script}
        if stags is None:
            stags = sDefault
        else:
            stags = stags | sDefault
            # Direct/FlatConverter drivers with integer vars:
            if ModelTags.integer in stags or ModelTags.binary in stags:
                stags = stags | {ModelTags.complementarity, ModelTags.logical}
            # Direct/FlatConverter drivers with non-convex quadratics:
            if ModelTags.quadraticnonconvex in stags:
                stags = stags | {ModelTags.polynomial}
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

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
                self._stats["objective"] = float(n)
            except:
                print("No solution, string: {}".format(n))
                self._stats["objective"] = None




class CPLEXSolver(AMPLSolver):
    def _setLPMethod(self, method: str):
        return "" if method == "SIMPLEX" else "baropt"

    def _setTimeLimit(self, seconds):
        return "time={}".format(seconds)

    def _setNThreads(self, threads):
        return "threads={}".format(threads)

    def setLogFile(self, name):
        return f"logfile=\"{name}\""

    def _getAMPLOptionsName(self):
        return "cplex"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.quadratic}
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

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
                self._stats["objective"] = float(n)
            except:
                print("No solution, string: {}".format(n))
                self._stats["objective"] = None

class BaronSolver(AMPLSolver):

    def _setLPMethod(self, method : str):
        return "" 

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
              self._stats["objective"] = float(n)
            except:
              print("No solution, string: {}".format(n))
              self._stats["objective"] = None

class OcteractSolver(AMPLSolver):
    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        super().__init__(exeName, timeout, nthreads, otherOptions)

    def _setLPMethod(self, method : str):
        return "" 

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
                self._stats["objective"] = float(n)
                return

class MindoptSolver(AMPLSolver):

    def _setLPMethod(self, method : str):
        m  = "1" if method == "SIMPLEX" else "2"
        return f"method={m}"

    def _setTimeLimit(self, seconds):
        return "max_time={}".format(seconds)

    def _setNThreads(self, threads):
        return "num_threads={}".format(threads)

    def _getAMPLOptionsName(self):
        return "mindopt"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        super().__init__(exeName, timeout, nthreads, otherOptions)

    def _doParseSolution(self, st, stdout=None):
        if not st:
            self._stats["outmsg"] = "Solution file empty"
            self._stats["timelimit"] = False
            return None
        self._stats["outmsg"] = st[1]
        self._stats["timelimit"] = "time limit" in st[0]
        tag = "objective"
        if tag in st[1]:
            n = st[1][st[1].index(tag) + len(tag):]
            try:
              n = n.split(",")[0]
              self._stats["objective"] = float(n)
            except:
              print("No solution, string: {}".format(n))
              self._stats["objective"] = None

class XpressSolver(AMPLSolver):

    def _setLPMethod(self, method : str):
        return "" if method == "SIMPLEX" else "barrier"

    def _setTimeLimit(self, seconds):
        return "maxtime={}".format(seconds)

    def setLogFile(self, name):
        return f"logfile=\"{name}\""

    def _setNThreads(self, threads):
        return "threads={}".format(threads)

    def _getAMPLOptionsName(self):
        return "xpress"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.quadratic}
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

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
                self._stats["objective"] = float(n)
            except:
                print("No solution, string: {}".format(n))
                self._stats["objective"] = None

class GurobiDirectSolver(MPDirectSolver):

    def _setLPMethod(self, method: str):
        m = "0" if method == "SIMPLEX" else "2"
        return f"alg:method {m}"

    def _getAMPLOptionsName(self):
        return "gurobi"

    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None):
        stags = {
                 ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.plinear,
                 ModelTags.quadratic, ModelTags.quadraticnonconvex,
                 ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                 ModelTags.unbdd, ModelTags.return_mipgap,
                 ModelTags.sos, ModelTags.presosenc, ModelTags.sens,
                 ModelTags.lazy_user_cuts, ModelTags.funcpieces,

                 ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                 ModelTags.multiobj, ModelTags.obj_priority,
                 ModelTags.multisol, ModelTags.sstatus,
                 ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

                 ModelTags.gurobi_cloud, ModelTags.gurobi_server,
                }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class CPLEXDirectSolver(MPDirectSolver):
    def _getAMPLOptionsName(self):
        return "cplexdirect"


    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None):
        stags = {
                 ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 # ModelTags.plinear,
                 # ModelTags.quadratic, ModelTags.quadraticnonconvex,
                 # ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric

                 ModelTags.relax,
                 ModelTags.multiobj
                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class XPRESSDirectSolver(MPDirectSolver):
    def _getAMPLOptionsName(self):
        return "xpress"


    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None):
        stags = {
                 ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.sos,  ModelTags.quadratic,

                 ModelTags.return_mipgap,

                 ModelTags.warmstart, ModelTags.mipstart,

                 ModelTags.multisol, ModelTags.sstatus,
                 ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax
                 
                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

class HighsSolver(MPDirectSolver):
    def _setLPMethod(self, method : str):
        m  = "simplex" if method == "SIMPLEX" else "ipm"
        return f"alg:method {m}"

    def _getAMPLOptionsName(self):
        return "highs"

    def _setNThreads(self, threads):
        parallel="parallel=on" if threads!=1 else ""
        return f"threads={threads} {parallel}"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.quadratic}
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class COPTSolver(MPDirectSolver):
    def _setLPMethod(self, method : str):
        m  = "1" if method == "SIMPLEX" else "2"
        return f"lp:method {m}"

    def _getAMPLOptionsName(self):
        return "copt"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.quadratic}
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class MosekSolver(MPDirectSolver):
    def _setLPMethod(self, method : str):
        m  = "1" if method == "SIMPLEX" else "2"
        return f"lp:method {m}"

    def _getAMPLOptionsName(self):
        return "mosek"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.quadratic
                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

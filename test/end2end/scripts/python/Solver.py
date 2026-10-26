import sys
import math
from pathlib import PurePath
import subprocess
import re 

from Model import Model, ModelTags
from SolverCategory import SolverCategory
from TimeMe import TimeMe
import psutil
import time
from threading import Timer  

class Solver(object):
    LPMETHODS = ["BARRIER", "SIMPLEX", "PDLP"]
    NLPMETHODS = ["REFORMULATION", "NATIVE", "NATIVEPL"]
    
    @staticmethod
    def resolveExecutable(name: str) -> str:
        """Path of the executable that will actually run for 'name':
           a name with a directory part is used as given (absolute),
           a bare name is looked up through amplpy modules and PATH"""
        import os, shutil
        if os.path.dirname(name):
            p = os.path.abspath(name)
            return p if os.path.isfile(p) else p + ": not found"
        try:
            from amplpy import modules
            modules.find(name)     # extends the search path with installed modules
        except Exception:
            pass
        found = shutil.which(name)
        return found if found else name + ": not found"

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
    
    def support_times(self) -> bool:
        """Return true to tell that the solver driver supports the option report_times"""
        return False
    
    def getExportLPFileName(self, modelfilename):
        fn = PurePath(modelfilename)
        ext = fn.with_stem(f"{fn.stem}_{self.getName()}").with_suffix('.lp')
        return ext.as_posix()

    def __init__(self, exeName, 
                 timeout=None, nthreads=None, otherOptions=None,
                 writeSolverName=False,
                 supportedTags=None,
                 unsupportedTags=None,
                 lpmethod = None,
                 nlpmethod=None,
                 name = None,):
        
        self._exePath = Solver.getExecutableName(exeName)
        self._timeout = timeout
        self._nthreads = nthreads
        if lpmethod:
            
            if not lpmethod in Solver.LPMETHODS :
                raise Exception("Valid methods to be forced are: {}").format(".".join(Solver.LPMETHODS))
        self._lpmethod = lpmethod
        self._otherOptions = otherOptions
        self._writeSolverName = writeSolverName
        self._supportedTags = supportedTags
        self._unsupportedTags = unsupportedTags
        self._nlpmethod = nlpmethod
        self._version = None
        self._name = name

    def get_version(self):
        return self._version
    def set_version(self,value):
        self._version = value
    
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

    def _evaluateRun(self, model: Model, output):
        expsol = model.getExpectedObjective()
        if expsol is not None:
            self._stats["eval_done"] = True
            self._assertAndRecord(expsol, self._stats["objective"],
                                  "objective")
        if model.hasExpectedOutput():
            self._getSolution(model)
            output += self._stats["outmsg"]
            output += self._getSolMessage(model)  # can span several lines
            for oitem in model.getExpectedOutput():
                self._stats["eval_done"] = True
                if str(oitem) not in output:
                   self._stats.setdefault("eval_fail_msg", []).append( \
                    "string  '" + str(oitem) + \
                    "'  is not found in the output")

    @staticmethod
    def _getSolMessage(model: Model) -> str:
        """Solver message of the .sol file: the lines before 'Options'"""
        sol = model.getSolFilePath()
        if not sol.exists():
            return ""
        lines = sol.read_text().splitlines()
        if "Options" in lines:
            lines = lines[:lines.index("Options")]
        return "\n".join(lines)

    def _assertAndRecord(self, expval, val, msg):
        b1 = isinstance(expval, (int, float))
        b2 = isinstance(val, (int, float))
        uneq = not math.isclose(expval, val, rel_tol=1e-6) if \
            b1 and b2 else expval != val
        if uneq:
            self._stats.setdefault("eval_fail_msg", []).append(msg + \
                ": value " + str(val) + \
                ", expected " + str(expval))

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
        self._output = stdout
        self._getSolution(model, stdout)
        self._evaluateRun(model, stdout)

    def get_output(self) -> str:
        """The solver's console output of the last run (NL runs)"""
        return getattr(self, "_output", "")

    def getName(self):
        if self._name:
            return self._name
        path = PurePath(self._exePath)
        return path.stem

    def named(self, name: str):
        """Set a solver name different from the executable's,
           e.g. to run one executable as two solvers"""
        self._name = name
        return self

    def acceptsMPOptions(self) -> bool:
        """Whether the solver accepts MP options and cases
           needing an MP driver (Model.needsMPDriver())"""
        return False

    def mapCaseOptionName(self, name: str) -> str:
        """AMPL option receiving a case's option 'name'"""
        return name

    def getFixedAMPLOptions(self) -> list:
        """Further (name, value) AMPL options set for every case"""
        return []

    ## SolverCategory of the driver, selecting the case's expected
    ## 'values_per_solver_category'; None: plain 'values' only
    CATEGORY = None

    def getCategory(self):
        return self.CATEGORY

    def checkCategory(self):
        """Warning if the constraints the driver accepts natively
        imply a category other than the declared one"""
        if self.CATEGORY is None or self.CATEGORY == SolverCategory.CPLegacy:
            return None
        actual = SolverCategory.ofDriver(self.getExecutable())
        if actual is None or actual == self.CATEGORY:
            return None
        return "WARNING: solver '{}' is declared {}, but its acc: options " \
            "imply {}; update {}.CATEGORY".format(self.getName(),
                self.CATEGORY.value, actual.value, type(self).__name__)

    def getExecutable(self):
        return self._exePath

    def getSolutionStats(self):
        return self._stats

    def setNThreads(self, nt):
        self._nthreads = nt
        
    def setNLPMethod(self, nlpmethod: str): # REFORMULATION, NATIVE, NATIVEPL 

        if not nlpmethod in Solver.NLPMETHODS:
                raise RuntimeError("Valid methods are {}".format(",".join(Solver.NLPMETHODS)))
        self._nlpmethod = nlpmethod
        
    def setLPMethod(self, lpmethod):
        if not lpmethod in Solver.LPMETHODS:
                raise RuntimeError("Valid methods are {}".format(",".join(Solver.LPMETHODS)))
        self._lpmethod = lpmethod

    def getNThreads(self):
        return self._nthreads
    def setExportLP(self, value : bool):
        self._exportLP = value
    def setTimeout(self, t):
        self._timeout = t

    def getTimeout(self):
        return self._timeout

    def getSupportedTags(self):
        return self._supportedTags

    def getUnsupportedTags(self):
        return self._unsupportedTags


def on_terminate(proc):
    print("Process {} terminated with exit code {}".format(proc, proc.returncode))


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
        sDefault = {
            ModelTags.script,
            ModelTags.cd}
        if supportedTags is None:
            supportedTags = sDefault
        else:
            supportedTags = supportedTags | sDefault

        super().__init__(exeName, timeout, nthreads, otherOptions,
                         supportedTags=supportedTags,
                         unsupportedTags=unsupportedTags)

    def _setExportLP(self, name: str):
        raise Exception("Not implemented in base class")

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

    
    def _getNLRunArgs(self) -> list:
        """Further command-line arguments when running an NL file
        directly (no AMPL options then, see getFixedAMPLOptions())"""
        return []

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
        if self._exportLP:
            toption = "{} {}".format(toption, self._setExportLP(self.getExportLPFileName(model.getFilePath())))

        try:
            if toption:
                return self._runProcess([self._exePath, model.getFilePath(), "-AMPL", toption]
                                        + self._getNLRunArgs(),
                                        timeout=self._timeout, logFile = logFile)
                
            else:
                return self._runProcess([self._exePath, model.getFilePath(), "-AMPL"]
                                        + self._getNLRunArgs(),
                                        timeout=self._timeout, logFile = logFile)
        except subprocess.TimeoutExpired:
            return "subprocess timeout expired"
        except subprocess.CalledProcessError as e:
            print(str(e))
            return str(e)
        return ""

    def stopProcess(p):
        # forcefully terminate solvers that do not terminate automatically 
        # after timeout
        p.kill()

    def _runProcess(self, args : list, vestigial=True, timeout=None, logFile = None):
      # ritorna stdout
      # throws if not successfull
         exit_code = 0
         if vestigial:
           if timeout:   ## https://psutil.readthedocs.io/en/latest/#terminate-my-children
             proc = psutil.Popen(args, universal_newlines=True,
             stdout=subprocess.PIPE,        ## PIPE still from subprocess?
             stderr=subprocess.PIPE)
             out = ""
             err = ""
             try:
               out, err = proc.communicate(timeout=timeout)
             except subprocess.TimeoutExpired:
               print("Timeout expired...")
               procs = proc.children()
               for p in procs:
                   print("Terminate ", p.pid)
                   p.terminate()
               gone, alive = psutil.wait_procs(procs, timeout=3, callback=on_terminate)
               for p in alive:
                   print("Kill ", p.pid)
                   p.kill()
               out, err = proc.communicate()
             exit_code = proc.returncode
             self._stats["exit_code"] = exit_code
             return out + '\n\nSTDERR:\n' + err
              ## This does not timeout with MP2NL:
              ## return subprocess.check_output(args, text=True,
              ##       stderr=subprocess.PIPE, timeout=self._timeout)
           else:
              return subprocess.check_output(args, text=True,
                     stderr=subprocess.PIPE)
         else:  # Since adding an NL-file case, this requires Ctrl-C after all tests
              resultTable = []
              SLICE_IN_SECONDS = 1
              p = subprocess.Popen(args, universal_newlines=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
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
                  # time.sleep(SLICE_IN_SECONDS)
                  try:    # Don't wait the whole second:
                      out, err = p.communicate(timeout=SLICE_IN_SECONDS)
                  except subprocess.TimeoutExpired:
                      pass
              out, err = p.communicate()
              exit_code = p.returncode
              self._stats["exit_code"] = exit_code
              if logFile is not None:
                if self.setLogFile(logFile) is None: # if not handled via option
                    print(out, flush=True)
                    print(err, flush=True)
              rss = max([p.rss for p in resultTable])
              vms = max([p.vms for p in resultTable])
              self._stats["rss"]= rss
              self._stats["vms"]= vms
              return out + '\n' + str(err)

    def _setNLPMethod(self, converter: str):
        return ""

    def getAMPLOptions(self, model) -> tuple:
        name = "{}_options".format(self._getAMPLOptionsName())
        value = ""
        if self._timeout:
            value += self._setTimeLimit(self._timeout)
        if self._nthreads:
            value += f" {self._setNThreads(self._nthreads)}"
        if self._lpmethod:
            value += f" {self._setLPMethod(self._lpmethod)}"
        if self._nlpmethod:
            value += f" {self._setNLPMethod(self._nlpmethod)}"
        if self._otherOptions:
            value += f" {self._otherOptions}"
        if self._exportLP:
            value += f" {self._setExportLP(self.getExportLPFileName(model.getFilePath()))}"
        if value:
            return (name, value)

class LindoSolver(AMPLSolver):
    def _setTimeLimit(self, seconds):
        return "maxtime={}".format(int(seconds))

    def _setNThreads(self, threads):
        return "threads={}".format(threads)

    def _getAMPLOptionsName(self):
        return "lindoglobal"

    def _setLPMethod(self, method : str):
        return ""

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
            ModelTags.linear,
            ModelTags.plinear,
            ModelTags.quadratic,
            ModelTags.quadratic_obj,
            ModelTags.quadraticnonconvex,

            ModelTags.socp,      
            ModelTags.socp_hard_to_recognize,
            ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric}
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

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

class IPOptSolver(AMPLSolver):
    def _setTimeLimit(self, seconds):
        return "max_cpu_time={}".format(seconds)

    def _setNThreads(self, threads):
        return ""

    def _getAMPLOptionsName(self):
        return "ipopt"

    def _setLPMethod(self, method : str):
        return ""

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
       stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.linear,
                 ModelTags.plinear,
                 ModelTags.quadratic,
                 ModelTags.quadratic_obj,
                 ModelTags.quadraticnonconvex,
                 ModelTags.formulas,

                 # ModelTags.complementarity,  -- not supported without MP2NL.

                 ModelTags.socp,      
                 ModelTags.socp_hard_to_recognize,
                 ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric}

       super().__init__(exeName, timeout, nthreads, otherOptions, stags)
    
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

        
class CouenneSolver(AMPLSolver):
    def _setTimeLimit(self, seconds):
        return "max_cpu_time={}".format(seconds)

    def _setNThreads(self, threads):
        return ""

    def _setLPMethod(self, method : str):
        return ""

    def _getAMPLOptionsName(self):
        return "couenne"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                  ModelTags.linear,
                  ModelTags.plinear,
                  ModelTags.quadratic,
                  ModelTags.quadratic_obj,
                  ModelTags.quadraticnonconvex,

                  ModelTags.socp,
                  ModelTags.socp_hard_to_recognize,
                  ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric}
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


class BonminSolver(AMPLSolver):
    def _setTimeLimit(self, seconds):
        return "max_cpu_time={}".format(seconds)

    def _setNThreads(self, threads):
        return ""

    def _setLPMethod(self, method : str):
        return ""

    def _getAMPLOptionsName(self):
        return "bonmin"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                  ModelTags.linear,
                  ModelTags.plinear,
                  ModelTags.quadratic,
                  ModelTags.quadratic_obj,
                  ModelTags.quadraticnonconvex,

                  ModelTags.socp,
                  ModelTags.socp_hard_to_recognize,
                  ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric}
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


class LoqoSolver(AMPLSolver):
    def _setTimeLimit(self, seconds):
        return "timlim={}".format(seconds)

    def _setNThreads(self, threads):
        return ""

    def _setLPMethod(self, method : str):
        return ""

    def _getAMPLOptionsName(self):
        return "loqo"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                  ModelTags.linear,
                  ModelTags.plinear,
                  ModelTags.quadratic,
                  ModelTags.quadratic_obj,
                  ModelTags.quadraticnonconvex,

                  ModelTags.socp,
                  ModelTags.socp_hard_to_recognize,
                  ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric}
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


class LgoSolver(AMPLSolver):
    def _setTimeLimit(self, seconds):
        return "timelim={}".format(seconds)

    def _setNThreads(self, threads):
        return ""

    def _getAMPLOptionsName(self):
        return "lgo"

    def _setLPMethod(self, method : str):
        return ""

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {
                 ModelTags.continuous, #ModelTags.integer, ModelTags.binary,
                 ModelTags.linear, ModelTags.quadratic,
                 ModelTags.quadraticnonconvex,
                 ModelTags.nonlinear, ModelTags.complementarity}
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

    def _doParseSolution(self, st, stdout=None):
        if not st:
            self._stats["outmsg"] = "Solution file empty"
            self._stats["timelimit"] = False
            return None
        self._stats["outmsg"] = st[0]
        self._stats["timelimit"] = "time limit" in st[1]
        
        pattern = r"Objective ([\\d.]*)\\n"
        for i,s in enumerate(st):
            match = re.search(pattern, st)
            if match:
                # Extract the matched number
                self._stats["objective"]  = float(match.group(1))
                break
            if i==4: # parse max the first 4 lines
                print("No solution, string: {}".format( st[1]))
                self._stats["objective"] = None
                break


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
                  ModelTags.linear, ModelTags.quadratic, ModelTags.sos,
                 ModelTags.return_mipgap,
                 ModelTags.sstatus,
                 ModelTags.multisol,
                 ModelTags.multiobj,
                 ModelTags.iis,
                 ModelTags.feasrelax, ModelTags.feasrelax_modes, ModelTags.feasrelax_squared}
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
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.linear,
                 ModelTags.plinear,
                 ModelTags.quadratic,
                 ModelTags.quadratic_obj,
                 ModelTags.quadraticnonconvex,

                 ModelTags.socp,      
                 ModelTags.socp_hard_to_recognize,
                 ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric}

        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

    def _doParseSolution(self, st, stdout=None):
        if not st:
            self._stats["outmsg"] = "Solution file empty"
            self._stats["timelimit"] = False
            return None
        self._stats["outmsg"] = st[0]
        self._stats["timelimit"] = "time limit" in st[0]
        tag = "Objective "
        if tag in st[1]:
            n = st[1][st[1].index(tag) + len(tag):]
            try:
              self._stats["objective"] = float(n)
            except:
              print("No solution, string: {}".format(n))
              self._stats["objective"] = None

class KnitroSolver(AMPLSolver):

    def _setLPMethod(self, method : str):
        return "act_lpalg={}".format(1 if method == "SIMPLEX" else 3)

    def _setTimeLimit(self, seconds):
        return f"maxtime_real={seconds} ma_maxtime_real={seconds}"

    def _setNThreads(self, threads):
        return f"threads={threads}"

    def _getAMPLOptionsName(self):
        return "knitro"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.linear,
                 ModelTags.quadratic,
                 ModelTags.quadratic_obj,
                 ModelTags.quadraticnonconvex,

                 ModelTags.socp,      
                 ModelTags.socp_hard_to_recognize,
                 ModelTags.expcones,
                 ModelTags.polynomial, ModelTags.polynomial_4,
                 ModelTags.nonlinear,
                 ModelTags.complementarity,
                 ModelTags.log,
                 ModelTags.trigonometric, ModelTags.htrigonometric}

        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

    def _doParseSolution(self, st, stdout=None):
        if not st:
            self._stats["outmsg"] = "Solution file empty"
            self._stats["timelimit"] = False
            return None
        self._stats["outmsg"] = st[0]
        self._stats["timelimit"] = "Unknown" in st[0] or "Time limit" in st[0]
        tag = "objective "
        pattern = r'\bobjective\s+(\d+\.?\d*(?:e[+-]?\d+)?)\b'
        match = re.search(pattern, st[1])
        if match:
            # Extract the matched number
            self._stats["objective"]  = float(match.group(1))
        else:
            print("No solution, string: {}".format( st[1]))
            self._stats["objective"] = None
              

class ConoptSolver(AMPLSolver):

    def _setLPMethod(self, method : str):
        return "" 

    def _setTimeLimit(self, seconds):
        return "maxftime={}".format(seconds)

    def _setNThreads(self, threads):
        return ""

    def _getAMPLOptionsName(self):
        return "conopt"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.linear,
                 ModelTags.plinear,
                 ModelTags.quadratic,
                 ModelTags.quadratic_obj,
                 ModelTags.quadraticnonconvex,

                 ModelTags.socp,      
                 ModelTags.socp_hard_to_recognize,
                 ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric}

        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

    def _doParseSolution(self, st, stdout=None):
        if not st:
            self._stats["outmsg"] = "Solution file empty"
            self._stats["timelimit"] = False
            return None
        self._stats["outmsg"] = st[0]
        self._stats["timelimit"] = "time limit" in st[0]
        tag = "Objective "
        if tag in st[1]:
            n = st[1][st[1].index(tag) + len(tag):]
            try:
              self._stats["objective"] = float(n)
            except:
              print("No solution, string: {}".format(n))
              self._stats["objective"] = None

class SnoptSolver(AMPLSolver):
    """SNOPT ASL drivers snopt, snopt777 (local NLP solver, no integer variables)."""

    def _setLPMethod(self, method : str):
        return ""

    def _setTimeLimit(self, seconds):
        return "lim:time={}".format(seconds)

    def _setNThreads(self, threads):
        return ""

    def _getAMPLOptionsName(self):
        return self.getName()           # snopt_options, snopt777_options

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous,
                 ModelTags.linear,
                 ModelTags.plinear,
                 ModelTags.quadratic,
                 ModelTags.quadratic_obj,
                 ModelTags.quadraticnonconvex,

                 ModelTags.socp,
                 ModelTags.socp_hard_to_recognize,
                 ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                 ModelTags.warmstart, ModelTags.sstatus}

        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

    def _doParseSolution(self, st, stdout=None):
        if not st:
            self._stats["outmsg"] = "Solution file empty"
            self._stats["timelimit"] = False
            return None
        self._stats["outmsg"] = st[0]
        self._stats["timelimit"] = "limit" in st[0]
        tag = "objective "
        if tag in st[0]:
            n = st[0][st[0].index(tag) + len(tag):].split(";")[0]
            try:
                self._stats["objective"] = float(n)
            except:
                print("No solution, string: {}".format(n))
                self._stats["objective"] = None

class OcteractSolver(AMPLSolver):
    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.plinear,
                 ModelTags.quadratic,
                 ModelTags.quadratic_obj,
                 ModelTags.quadraticnonconvex,

                 ModelTags.socp,      
                 ModelTags.socp_hard_to_recognize,
                 ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric}

        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

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
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                  ModelTags.linear, ModelTags.quadratic, ModelTags.sos}
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

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
        # The ASL driver's maxtime is an integer keyword; "maxtime=5.0"
        # leaves ".0" as an unknown keyword and aborts the solve
        return "maxtime={}".format(int(seconds))

    def setLogFile(self, name):
        return f"logfile=\"{name}\""

    def _setNThreads(self, threads):
        return "threads={}".format(threads)

    def _getAMPLOptionsName(self):
        return "xpress"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.linear,
                 ModelTags.integer, ModelTags.binary,
                 ModelTags.quadratic,
                 ModelTags.quadraticnonconvex,
                 ModelTags.quadratic_obj_nonconvex,
                 ModelTags.socp,  # ModelTags.socp_hard_to_recognize,
                 ModelTags.plinear}
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


# MP Direct / FlatConverter drivers
class MPDirectSolver(AMPLSolver):
    
    def support_times(self) -> bool:
        return True
    
    def _supported_nl(self) ->list:
        """Return all constraint types supported by MP"""
        # This function should now return all constraint types supported by MP
        # because solvers are not rejecting the option starting from the solvers
        # released on 20240320
        return ["acc:abs", "acc:acos", "acc:acosh", "acc:and", "acc:asin", 
            "acc:asinh", "acc:atan", "acc:atanh", "acc:cos", "acc:cosh", 
            "acc:div", "acc:exp", "acc:expa", "acc:indeq", "acc:indge", 
            "acc:indle", "acc:log", "acc:max", "acc:min", "acc:or", 
            "acc:pow", "acc:quadeq", "acc:quadge", "acc:quadle", 
            "acc:sin", "acc:sinh", "acc:sos2", "acc:tan", "acc:tanh"]
     
    def _setNLPMethod(self, method: str): # reformulation, native, nativepl 
        if method == "REFORMULATION":
            return " ".join([f"{opt}=1" for opt in self._supported_nl()])
        if method == "NATIVE":
            return " ".join([f"{opt}=?" for opt in self._supported_nl()])
        else:
            raise RuntimeError(f"{method} not supported in this solver")
        
    def _setLPMethod(self, method : str):
        # typically have to reimplement this
        m  = "0" if method == "SIMPLEX" else "2"
        return f"alg:method {m}"

    def _setExportLP(self, name: str):
        return f"writeprob={name}"

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
                    ModelTags.mpdriver
                    }
        if stags is None:
            stags = sDefault
        else:
            stags = stags | sDefault
            # Direct/FlatConverter drivers with integer vars:
            if ModelTags.integer in stags or ModelTags.binary in stags:
                stags = stags | {
                    ModelTags.complementarity,
                    ModelTags.logical,
                    ModelTags.plinear,
                    ModelTags.nonlinear,
                    ModelTags.prod_with_bin_var,
                    ModelTags.log,
                    ModelTags.multiobj,
                    ModelTags.int_ops
                }
            # Direct/FlatConverter drivers with non-convex quadratics:
            if ModelTags.quadraticnonconvex in stags:
                stags = stags | {ModelTags.polynomial, ModelTags.prod_int}
        
        opts = self._get_timing_options()
        if otherOptions:
            opts  = f"{opts} {otherOptions}"
        super().__init__(exeName, timeout, nthreads, opts, stags)

    def _get_timing_options(self):
        return "timing=1"

    def acceptsMPOptions(self) -> bool:
        return True

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
    CATEGORY = SolverCategory.MILPWithLogicals

    
    def _setNLPMethod(self, method: str): # reformulation, native, nativepl 
        if method == "REFORMULATION":
            return " ".join([f"{opt}=1" for opt in self._supported_nl()])
        if method == "NATIVE":
            param = " ".join([f"{opt}=?" for opt in self._supported_nl()])
            return param
        if method == "NATIVEPL":
            param = " ".join([f"{opt}=2" for opt in self._supported_nl()])
            param += " pre:funcnonlinear=-1"
            return param
        
    def _setLPMethod(self, method: str):
        m = "0" if method == "SIMPLEX" else "2"
        return f"alg:method {m}"

    def _getAMPLOptionsName(self):
        return "gurobi"

    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None):
        stags = {
                 ModelTags.gurobi,

                 ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.plinear,
                 ModelTags.quadratic,
                 ModelTags.quadratic_obj,
                 ModelTags.quadraticnonconvex,
                 ModelTags.formulas,

                 ModelTags.socp,      ## MP transforms cones to quadratics
                 ModelTags.socp_hard_to_recognize,
                 ModelTags.expcones,  ## Although might not have the expconic solver

                 ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,
                 ModelTags.nlobj_native,

                 ModelTags.signpow,
                 ModelTags.logistic,

                 ModelTags.unbdd,
                 ModelTags.qcpdual,
                 ModelTags.return_mipgap,
                 ModelTags.sos, ModelTags.presosenc, ModelTags.sens,
                 ModelTags.lazy_user_cuts, ModelTags.funcpieces,
                 ModelTags.integralityfocus,

                 ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                 ModelTags.multiobj, ModelTags.obj_priority,
                 ModelTags.multisol, ModelTags.sstatus, ModelTags.fixmodel,
                 ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax, ModelTags.feasrelax_modes, ModelTags.feasrelax_squared,

                 ModelTags.writelp, ModelTags.writesol,

                 ModelTags.check_pl_approx_exp,
                 ModelTags.check_pl_approx_expA,
                 ModelTags.check_pl_approx_log,

                 ModelTags.check_pl_approx_sin,
                 ModelTags.check_pl_approx_cos,
                 ModelTags.check_pl_approx_tan,

                 ModelTags.check_pl_approx_pow,
                 ModelTags.check_qp_approx_pow,

                 ModelTags.check_pl2sos2,

                 ModelTags.check_sos2_from_pl,

                 ModelTags.gurobi_cloud, ModelTags.gurobi_server, ModelTags.gurobi_wls
                }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

class IlogCPSolver(MPDirectSolver):
    """IBM ILOG CP Optimizer driver (ilogcp), old MP framework;
    uses CPLEX for problems without logical/nonlinear constraints."""
    CATEGORY = SolverCategory.CPLegacy

    def _getAMPLOptionsName(self):
        return "ilogcp"

    def _setLPMethod(self, method : str):
        return ""

    def _setNLPMethod(self, method: str):
        return ""

    def _setTimeLimit(self, seconds):
        return "timelimit={}".format(seconds)

    def _setNThreads(self, threads):
        return "workers={}".format(threads)

    def setLogFile(self, name):
        return ""

    def _setExportLP(self, name: str):
        return f"exportfile={name}"

    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None):
        stags = {
                 ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.linear, ModelTags.logical,
                 ModelTags.quadratic, ModelTags.quadratic_obj,
                 ModelTags.multisol, ModelTags.multiobj,
                 ModelTags.iis,
                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)
        self._supportedTags = self._supportedTags - {ModelTags.int_ops}  # old MP framework

class CPLEXDirectSolver(MPDirectSolver):
    CATEGORY = SolverCategory.MIPWithIndicators

    def _getAMPLOptionsName(self):
        return "cplex"


    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None):

        stags = {
                 ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.sos,
                 
                 ModelTags.quadratic, ModelTags.quadratic_obj, ModelTags.quadratic_obj_nonconvex,
                 ModelTags.plinear,
                 
                 ModelTags.socp,      ## MP transforms cones to quadratics
                 ModelTags.socp_hard_to_recognize,

                 ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,
                 
                 ModelTags.unbdd,
                 ModelTags.return_mipgap,
                 ModelTags.sens,
                 
                 ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,
                 ModelTags.multiobj, ModelTags.obj_priority,
                 ModelTags.multisol, ModelTags.sstatus, 
                 ModelTags.iis, ModelTags.feasrelax, ModelTags.feasrelax_modes, ModelTags.feasrelax_squared, ModelTags.fixmodel,
                 ModelTags.writelp, ModelTags.writesol,

                 ModelTags.check_pl2sos2,
                 ModelTags.check_sos2_from_pl

                 
                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

class XPRESSDirectSolver(MPDirectSolver):
    CATEGORY = SolverCategory.MILPWithLogicals

       
    def _getAMPLOptionsName(self):
        return "xpress"

    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None):
        stags = {
                 ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.sos,
                 
                 ModelTags.quadratic, ModelTags.quadratic_obj, ModelTags.quadraticnonconvex,
                 ModelTags.plinear,
                 
                 ModelTags.socp,
                 ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,
                 ModelTags.formulas,
                 
                 ModelTags.return_mipgap,

                 ModelTags.warmstart, ModelTags.mipstart,
                 ModelTags.multiobj,
                 ModelTags.multisol, ModelTags.sstatus,
                 ModelTags.iis, ModelTags.fixmodel,
                 
                 ModelTags.writelp, ModelTags.writesol,

                 ModelTags.option_flag_barrier,

                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)



class CbcMPSolver(MPDirectSolver):
    def _setLPMethod(self, method : str):
        m  = "simplex" if method == "SIMPLEX" else "ipm"
        return f"alg:method {m}"

    def _getAMPLOptionsName(self):
        return "cbc"

    def _setNThreads(self, threads):
        return ""

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.quadratic_obj, ModelTags.sos,

                 ModelTags.writelp, ModelTags.writesol,

                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

class MP2NLSolver(MPDirectSolver):
            def _getAMPLOptionsName(self):
                return "mp2nl"

            def _setNThreads(self, threads):
                return "solver=gurobi"  ## For auto-testing

            def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
                stags = {
                ModelTags.mp2nl,

                ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                ModelTags.plinear,
                ModelTags.quadratic,
                ModelTags.quadratic_obj,
                ModelTags.quadraticnonconvex,

                ModelTags.socp,      ## MP transforms cones to quadratics
                ModelTags.socp_hard_to_recognize,
                ModelTags.expcones,  ## Although might not have the expconic solver

                ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                ModelTags.signpow,
                ModelTags.logistic,

                ModelTags.unbdd,
                ModelTags.qcpdual,
#                ModelTags.return_mipgap,
                ModelTags.sos, #ModelTags.presosenc,
                ModelTags.sens,
                ModelTags.lazy_user_cuts,
#                ModelTags.funcpieces,
#                ModelTags.integralityfocus,

                ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                ModelTags.multiobj,
                ModelTags.obj_priority,
#                ModelTags.multisol,
                ModelTags.sstatus,
#                ModelTags.fixmodel,
 #               ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

#                ModelTags.writelp, ModelTags.writesol,

                         }
                super().__init__(exeName, timeout, nthreads, otherOptions, stags)

class BaronMPSolver(MPDirectSolver):
    CATEGORY = SolverCategory.MILP

    def _getAMPLOptionsName(self):
        return "baron"
    
    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary, 
                 ModelTags.quadratic,
                 ModelTags.quadratic_obj,
                 ModelTags.quadratic_obj_nonconvex,
                 ModelTags.quadraticnonconvex,
                 ModelTags.formulas,
                 ModelTags.socp, ModelTags.socp_hard_to_recognize,
                 ModelTags.sos,
                 ModelTags.nonlinear,
                 ModelTags.log,
                 #ModelTags.trigonometric,
                 ModelTags.return_mipgap,
                 ModelTags.warmstart,  ModelTags.mipstart,
                 ModelTags.writelp, ModelTags.writesol,

                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class SCIPSolver(MPDirectSolver):
    CATEGORY = SolverCategory.MILPWithLogicals

    def _getAMPLOptionsName(self):
        return "scip"
     
    def _setNThreads(self, threads):
        return ""

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.fixmodel, ModelTags.continuous, ModelTags.linear,
                 ModelTags.integer, ModelTags.binary,
                 ModelTags.quadratic,
                 ModelTags.quadratic_obj,
                 ModelTags.quadratic_obj_nonconvex,
                 ModelTags.quadraticnonconvex,
                 ModelTags.socp, ModelTags.socp_hard_to_recognize,
                 ModelTags.sos,
                 ModelTags.nonlinear, ModelTags.formulas_scip,
                 ModelTags.log,
                 ModelTags.trigonometric,
                 ModelTags.return_mipgap,
                 ModelTags.multisol, ModelTags.iis,
                 ModelTags.warmstart, ModelTags.mipstart,

                 ModelTags.writelp, ModelTags.writesol,

                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class GCGSolver(MPDirectSolver):
    CATEGORY = SolverCategory.MILP

    def _getAMPLOptionsName(self):
        return "gcg"

    def _setNThreads(self, threads):
        return ""

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.linear,
                 ModelTags.integer, ModelTags.binary,
                 ModelTags.return_mipgap,

                 ModelTags.writelp, ModelTags.writesol,

                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class HighsSolver(MPDirectSolver):
    CATEGORY = SolverCategory.MILP

    def _setLPMethod(self, method : str):
        mymap = { "SIMPLEX" : "simplex",
                  "BARRIER" : "ipm",
                  "PDLP"    : "pdlp"}
        return f"lpmethod {mymap[method]}"

    def _getAMPLOptionsName(self):
        return "highs"

    def _setNThreads(self, threads):
        parallel="parallel=on" if threads!=1 else ""
        return f"threads={threads} {parallel}"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.fixmodel, ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.quadratic_obj,

                 ModelTags.sstatus, 
                 ModelTags.warmstart, ModelTags.return_mipgap, ModelTags.iis,
                 ModelTags.sens, ModelTags.feasrelax, ModelTags.mipstart,

                 ModelTags.writelp, ModelTags.writesol,

                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

        
class KnitroMPSolver(MPDirectSolver):
    def _setLPMethod(self, method : str):
        mymap = { "SIMPLEX" : 1,
                  "BARRIER" : 2
                  }
        return f"lp:method {mymap[method]}"

    def _getAMPLOptionsName(self):
        return "knitro"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.quadratic, ModelTags.quadratic_obj,
                 ModelTags.quadraticnonconvex,
                 ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,
                 ModelTags.quadratic_obj_nonconvex,

                 ModelTags.socp,      ## MP transforms cones to quadratics
                 ModelTags.socp_hard_to_recognize,

                 ModelTags.writelp,

                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)




class COPTSolver(MPDirectSolver):
    CATEGORY = SolverCategory.MIPWithIndicators

    def _setLPMethod(self, method : str):
        mymap = { "SIMPLEX" : 1,
                  "BARRIER" : 2,
                  "PDLP"    : 6}
        return f"lp:method {mymap[method]}"

    def _getAMPLOptionsName(self):
        return "copt"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.fixmodel, ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.sstatus, ModelTags.return_mipgap, ModelTags.sens,
                 ModelTags.unbdd,     ## unbdd_01 fails: COPT returns no primal point, and its
                                      ## simplex/barrier rays do not improve the objective
                 ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,
                 ModelTags.sos, ModelTags.lazy_user_cuts,
                 ModelTags.iis, ModelTags.feasrelax, ModelTags.feasrelax_modes,
                 ModelTags.writesol,
                 ModelTags.quadratic, ModelTags.quadratic_obj,
                 ModelTags.quadraticnonconvex,
                 ModelTags.formulas,
                 ModelTags.nlobj_native,

                 ModelTags.socp,      ## MP transforms cones to quadratics
                 ModelTags.socp_hard_to_recognize,

                 ModelTags.writelp,

                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

class cuOptSolver(MPDirectSolver):
    def _setLPMethod(self, method : str):
        mymap = { "SIMPLEX" : 2,
                  "PDLP"    : 1}
        return f"alg:method {mymap[method]}"


    def _getAMPLOptionsName(self):
        return "cuoptmp"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.continuous, ModelTags.linear,
                 ModelTags.integer, ModelTags.binary,
                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

class CpsatSolver(MPDirectSolver):
    """OR-Tools CP-SAT: integer variables and integer coefficients only."""
    CATEGORY = SolverCategory.CP

    def _setLPMethod(self, method : str):
        return ""                               # No LP methods

    def _getAMPLOptionsName(self):
        return "cpsat"

    def setLogFile(self, name):
        return ""                               # No log file option, only outlev

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.integer, ModelTags.binary,
                 ModelTags.return_mipgap,
                 ModelTags.prod_int,     # products, integer powers natively
                 ModelTags.multisol,     # improving solutions found
                 ModelTags.mipstart,     # solution hints
                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)
        # MPDirectSolver adds tags for drivers with continuous variables
        # (and for nonlinear and complementarity models reformulated to
        # them). CP-SAT does not have them.
        self._supportedTags = self._supportedTags - {
            ModelTags.continuous,
            ModelTags.complementarity,
            ModelTags.plinear,
            ModelTags.nonlinear,
            ModelTags.prod_with_bin_var,
            ModelTags.log
        }
class GecodeSolver(MPDirectSolver):
    """Gecode, old MP framework: integer variables and constraint
    programming (logical and global constraints) only."""
    CATEGORY = SolverCategory.CPLegacy

    def _setLPMethod(self, method : str):
        return ""                               # No LP methods

    def _setNLPMethod(self, method: str):
        return ""

    def _getAMPLOptionsName(self):
        return "gecode"

    def setLogFile(self, name):
        return ""                               # No log file option

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.integer, ModelTags.binary,
                 ModelTags.logical,
                 ModelTags.multisol,
                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)
        # As for CP-SAT: no continuous variables, so none of the tags
        # MPDirectSolver adds for drivers with them
        self._supportedTags = self._supportedTags - {
            ModelTags.continuous,
            ModelTags.complementarity,
            ModelTags.plinear,
            ModelTags.nonlinear,
            ModelTags.prod_with_bin_var,
            ModelTags.log,
            ModelTags.int_ops   # old MP framework
        }

class MosekSolver(MPDirectSolver):
    CATEGORY = SolverCategory.MIPWithIndicators

    def _setLPMethod(self, method : str):
        m  = "1" if method == "SIMPLEX" else "4"
        return f"lp:method {m}"

    def _getAMPLOptionsName(self):
        return "mosek"

    def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
        stags = {ModelTags.fixmodel, ModelTags.continuous, ModelTags.linear, ModelTags.integer,
                 ModelTags.binary,

                 ModelTags.quadratic, ModelTags.quadratic_obj,
                 ModelTags.socp, ModelTags.socp_hard_to_recognize,
                 ModelTags.expcones,

                 ModelTags.warmstart, ModelTags.mipstart,
                 ModelTags.return_mipgap, ModelTags.sens, ModelTags.sstatus,

                 ModelTags.writelp, ModelTags.writesol,
                 }
        # ModelTags.quadratic
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class CPLEXODHSolver(MPDirectSolver):
    def _getAMPLOptionsName(self):
        return "cplexodh"

    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None):
        stags = {
                 ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.plinear, ModelTags.sos,
                 # ModelTags.quadratic, ModelTags.quadraticnonconvex,
                 ModelTags.nonlinear,
                 ModelTags.log,
                 # ModelTags.trigonometric
                 ModelTags.relax,
                 ModelTags.multiobj,
                 ModelTags.sstatus,
                 ModelTags.return_mipgap,

                 ModelTags.writelp, ModelTags.writesol,

                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)

class GUROBIODHSolver(MPDirectSolver):
    def _getAMPLOptionsName(self):
        return "gurobiodh"

    def __init__(self, exeName, timeout=None, nthreads=None,
                 otherOptions=None):
        stags = {
                 ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                 ModelTags.plinear, ModelTags.sos,
                 # ModelTags.quadratic, ModelTags.quadraticnonconvex,
                 ModelTags.nonlinear,
                 ModelTags.log,
                 # ModelTags.trigonometric
                 ModelTags.relax,
                 ModelTags.multiobj,
                 ModelTags.sstatus,
                 ModelTags.return_mipgap,

                 ModelTags.writelp, ModelTags.writesol,

                 }
        super().__init__(exeName, timeout, nthreads, otherOptions, stags)


### ASL+MP2NL solvers (to run with mp2nl=1)
class ASLMP2NLSolver(MPDirectSolver):
    """An ASL driver built with the MP2NL wrapper, run through MP2NL:
       mp2nl=1 in '<config>_options' makes the wrapper call MP2NL, which
       reformulates the model with the case's MP options (mp2nl_options)
       and runs the driver on the new NL file. The driver reads the time
       limit, threads and LP method from its own '<engine>_options'."""
    def _getMP2NLConfig(self):
        """MP2NL configuration of the wrapper (MP2NL_SOLVERNAME)"""
        return self._getAMPLOptionsName()

    def _getEngineOptionsName(self):
        """Options variable read by the ASL driver itself"""
        return self._getAMPLOptionsName()

    def mapCaseOptionName(self, name: str) -> str:
        if name == f"{self._getAMPLOptionsName()}_options":
            return "mp2nl_options"     # written for the MP driver
        return name

    def getAMPLOptions(self, model) -> tuple:
        value = self._otherOptions if self._otherOptions else ""
        if self._exportLP:
            value += f" {self._setExportLP(self.getExportLPFileName(model.getFilePath()))}"
        return ("mp2nl_options", value)

    def getFixedAMPLOptions(self) -> list:
        native = []
        if self._timeout:
            native.append(self._setTimeLimit(self._timeout))
        if self._nthreads:
            native.append(self._setNThreads(self._nthreads))
        if self._lpmethod:
            native.append(self._setLPMethod(self._lpmethod))
        native = " ".join([o for o in native if o])
        cfg = f"{self._getMP2NLConfig()}_options"
        eng = f"{self._getEngineOptionsName()}_options"
        if cfg == eng:
            return [(cfg, f"mp2nl=1 {native}".strip())]
        return [(cfg, "mp2nl=1"), (eng, native)]

    def _getNLRunArgs(self) -> list:
        # The wrapper recognizes mp2nl=1 only as a separate argument
        return ["mp2nl=1"]


class GurobiASLMP2NLSolver(ASLMP2NLSolver):
            def _getAMPLOptionsName(self):
                return "gurobi"

            def _getMP2NLConfig(self):
                return "gurobiasl"

            def _get_timing_options(self):
                return ""

            def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
                stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                ModelTags.plinear,
                ModelTags.quadratic,
                ModelTags.quadratic_obj,
                ModelTags.quadraticnonconvex,

                ModelTags.socp,      ## MP transforms cones to quadratics
                ModelTags.socp_hard_to_recognize,
                ModelTags.expcones,  ## Although might not have the expconic solver

                ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                ModelTags.unbdd,
                ModelTags.qcpdual,
#                ModelTags.return_mipgap,
                ModelTags.sos, #ModelTags.presosenc,
                ModelTags.sens,
                ModelTags.lazy_user_cuts,
#                ModelTags.funcpieces,
#                ModelTags.integralityfocus,

                ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                ModelTags.multiobj,
                ModelTags.obj_priority,
#                ModelTags.multisol,
                ModelTags.sstatus,
#                ModelTags.fixmodel,
 #               ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

#                ModelTags.writelp, ModelTags.writesol,

                         }
                super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class XpressASLMP2NLSolver(ASLMP2NLSolver):
            def _getAMPLOptionsName(self):
                return "xpress"

            def _getMP2NLConfig(self):
                return "xpressasl"

            def _setTimeLimit(self, seconds):
                return "maxtime={}".format(int(seconds))   # integer keyword

            def _get_timing_options(self):
                return ""

            def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
                stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                ModelTags.plinear,
                ModelTags.quadratic,
                ModelTags.quadratic_obj,
                ModelTags.quadraticnonconvex,

                ModelTags.socp,      ## MP transforms cones to quadratics
                ModelTags.socp_hard_to_recognize,
                ModelTags.expcones,  ## Although might not have the expconic solver

                ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                ModelTags.unbdd,
                ModelTags.qcpdual,
#                ModelTags.return_mipgap,
                ModelTags.sos, #ModelTags.presosenc,
                ModelTags.sens,
                ModelTags.lazy_user_cuts,
#                ModelTags.funcpieces,
#                ModelTags.integralityfocus,

                ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                ModelTags.multiobj,
                ModelTags.obj_priority,
#                ModelTags.multisol,
                ModelTags.sstatus,
#                ModelTags.fixmodel,
 #               ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

#                ModelTags.writelp, ModelTags.writesol,

                         }
                super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class CPLEXASLMP2NLSolver(ASLMP2NLSolver):
            def _getAMPLOptionsName(self):
                return "cplex"

            def _getMP2NLConfig(self):
                return "cplexasl"

            def _setTimeLimit(self, seconds):
                return "time={}".format(seconds)

            def _get_timing_options(self):
                return ""

            def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
                stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                ModelTags.plinear,
                ModelTags.quadratic,
                ModelTags.quadratic_obj,
                ModelTags.quadraticnonconvex,

                ModelTags.socp,      ## MP transforms cones to quadratics
                ModelTags.socp_hard_to_recognize,
                ModelTags.expcones,  ## Although might not have the expconic solver

                ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                ModelTags.unbdd,
                ModelTags.qcpdual,
#                ModelTags.return_mipgap,
                ModelTags.sos, #ModelTags.presosenc,
                ModelTags.sens,
                ModelTags.lazy_user_cuts,
#                ModelTags.funcpieces,
#                ModelTags.integralityfocus,

                ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                ModelTags.multiobj,
                ModelTags.obj_priority,
#                ModelTags.multisol,
                ModelTags.sstatus,
#                ModelTags.fixmodel,
 #               ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

#                ModelTags.writelp, ModelTags.writesol,

                         }
                super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class ConoptMP2NLSolver(ASLMP2NLSolver):
            def _getAMPLOptionsName(self):
                return "conopt"

            def _setTimeLimit(self, seconds):
                return "maxftime={}".format(seconds)

            def _setNThreads(self, threads):
                return ""

            def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
                stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                ModelTags.plinear,
                ModelTags.quadratic,
                ModelTags.quadratic_obj,
                ModelTags.quadraticnonconvex,

                ModelTags.socp,      ## MP transforms cones to quadratics
                ModelTags.socp_hard_to_recognize,
                ModelTags.expcones,  ## Although might not have the expconic solver

                ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                ModelTags.unbdd,
                ModelTags.qcpdual,
#                ModelTags.return_mipgap,
                ModelTags.sos, #ModelTags.presosenc,
                ModelTags.sens,
                ModelTags.lazy_user_cuts,
#                ModelTags.funcpieces,
#                ModelTags.integralityfocus,

                ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                ModelTags.multiobj,
                ModelTags.obj_priority,
#                ModelTags.multisol,
                ModelTags.sstatus,
#                ModelTags.fixmodel,
 #               ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

#                ModelTags.writelp, ModelTags.writesol,

                         }
                super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class BonminMP2NLSolver(ASLMP2NLSolver):
            def _getAMPLOptionsName(self):
                return "bonmin"

            def _setNThreads(self, threads):
                return ""

            def _get_timing_options(self):
                return ""

            def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
                stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                ModelTags.plinear,
                ModelTags.quadratic,
                ModelTags.quadratic_obj,
                ModelTags.quadraticnonconvex,

                ModelTags.socp,      ## MP transforms cones to quadratics
                ModelTags.socp_hard_to_recognize,
                ModelTags.expcones,  ## Although might not have the expconic solver

                ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                ModelTags.unbdd,
                ModelTags.qcpdual,
#                ModelTags.return_mipgap,
                ModelTags.sos, #ModelTags.presosenc,
                ModelTags.sens,
                ModelTags.lazy_user_cuts,
#                ModelTags.funcpieces,
#                ModelTags.integralityfocus,

                ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                ModelTags.multiobj,
                ModelTags.obj_priority,
#                ModelTags.multisol,
                ModelTags.sstatus,
#                ModelTags.fixmodel,
 #               ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

#                ModelTags.writelp, ModelTags.writesol,

                         }
                super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class CouenneMP2NLSolver(ASLMP2NLSolver):
            def _getAMPLOptionsName(self):
                return "couenne"

            def _setTimeLimit(self, seconds):
                return "max_cpu_time={}".format(seconds)

            def _setNThreads(self, threads):
                return ""

            def _get_timing_options(self):
                return ""

            def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
                stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                ModelTags.plinear,
                ModelTags.quadratic,
                ModelTags.quadratic_obj,
                ModelTags.quadraticnonconvex,

                ModelTags.socp,      ## MP transforms cones to quadratics
                ModelTags.socp_hard_to_recognize,
                ModelTags.expcones,  ## Although might not have the expconic solver

                ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                ModelTags.unbdd,
                ModelTags.qcpdual,
#                ModelTags.return_mipgap,
                ModelTags.sos, #ModelTags.presosenc,
                ModelTags.sens,
                ModelTags.lazy_user_cuts,
#                ModelTags.funcpieces,
#                ModelTags.integralityfocus,

                ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                ModelTags.multiobj,
                ModelTags.obj_priority,
#                ModelTags.multisol,
                ModelTags.sstatus,
#                ModelTags.fixmodel,
 #               ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

#                ModelTags.writelp, ModelTags.writesol,

                         }
                super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class LindoglobalMP2NLSolver(ASLMP2NLSolver):
            def _getAMPLOptionsName(self):
                return "lindoglobal"

            def _setTimeLimit(self, seconds):
                return "maxtime={}".format(int(seconds))

            def _setNThreads(self, threads):
                return "threads={}".format(threads)

            def _get_timing_options(self):
                return ""

            def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
                stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                ModelTags.plinear,
                ModelTags.quadratic,
                ModelTags.quadratic_obj,
                ModelTags.quadraticnonconvex,
                ModelTags.formulas,

                ModelTags.socp,      ## MP transforms cones to quadratics
                ModelTags.socp_hard_to_recognize,
                ModelTags.expcones,  ## Although might not have the expconic solver

                ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                ModelTags.unbdd,
                ModelTags.qcpdual,
#                ModelTags.return_mipgap,
                ModelTags.sos, #ModelTags.presosenc,
                ModelTags.sens,
                ModelTags.lazy_user_cuts,
#                ModelTags.funcpieces,
#                ModelTags.integralityfocus,

                ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                ModelTags.multiobj,
                ModelTags.obj_priority,
#                ModelTags.multisol,
                ModelTags.sstatus,
#                ModelTags.fixmodel,
 #               ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

#                ModelTags.writelp, ModelTags.writesol,

                         }
                super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class BaronMP2NLSolver(ASLMP2NLSolver):
            def _getAMPLOptionsName(self):
                return "baron"

            def _setTimeLimit(self, seconds):
                return "maxtime={}".format(seconds)

            def _setNThreads(self, threads):
                return "threads={}".format(threads)

            def _get_timing_options(self):
                return ""

            def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
                stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                ModelTags.plinear,
                ModelTags.quadratic,
                ModelTags.quadratic_obj,
                ModelTags.quadraticnonconvex,

                ModelTags.socp,      ## MP transforms cones to quadratics
                ModelTags.socp_hard_to_recognize,
                ModelTags.expcones,  ## Although might not have the expconic solver

                ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                ModelTags.unbdd,
                ModelTags.qcpdual,
#                ModelTags.return_mipgap,
                ModelTags.sos, #ModelTags.presosenc,
                ModelTags.sens,
                ModelTags.lazy_user_cuts,
#                ModelTags.funcpieces,
#                ModelTags.integralityfocus,

                ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                ModelTags.multiobj,
                ModelTags.obj_priority,
#                ModelTags.multisol,
                ModelTags.sstatus,
#                ModelTags.fixmodel,
 #               ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

#                ModelTags.writelp, ModelTags.writesol,

                         }
                super().__init__(exeName, timeout, nthreads, otherOptions, stags)


class KnitroMP2NLSolver(ASLMP2NLSolver):
            def _getAMPLOptionsName(self):
                return "knitro"

            def _setLPMethod(self, method : str):
                return "act_lpalg={}".format(1 if method == "SIMPLEX" else 3)

            def _setTimeLimit(self, seconds):
                return f"maxtime_real={seconds} ma_maxtime_real={seconds}"

            def _setNThreads(self, threads):
                return f"threads={threads}"

            def __init__(self, exeName, timeout=None, nthreads=None, otherOptions=None):
                stags = {ModelTags.continuous, ModelTags.integer, ModelTags.binary,
                ModelTags.plinear,
                ModelTags.quadratic,
                ModelTags.quadratic_obj,
                ModelTags.quadraticnonconvex,
                ModelTags.formulas,

                ModelTags.socp,      ## MP transforms cones to quadratics
                ModelTags.socp_hard_to_recognize,
                ModelTags.expcones,  ## Although might not have the expconic solver

                ModelTags.nonlinear, ModelTags.log, ModelTags.trigonometric,

                ModelTags.unbdd,
                ModelTags.qcpdual,
#                ModelTags.return_mipgap,
                ModelTags.sos, #ModelTags.presosenc,
                ModelTags.sens,
                ModelTags.lazy_user_cuts,
#                ModelTags.funcpieces,
#                ModelTags.integralityfocus,

                ModelTags.relax, ModelTags.warmstart, ModelTags.mipstart,

                ModelTags.multiobj,
                ModelTags.obj_priority,
#                ModelTags.multisol,
                ModelTags.sstatus,
#                ModelTags.fixmodel,
 #               ModelTags.iis, ModelTags.iisforce, ModelTags.feasrelax,

#                ModelTags.writelp, ModelTags.writesol,

                         }
                super().__init__(exeName, timeout, nthreads, otherOptions, stags)

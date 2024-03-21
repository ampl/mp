from threading import Lock
import math
from pathlib import Path
from shutil import which
from token import LBRACE

from Solver import Solver
from amplpy import AMPL, Kind, OutputHandler, ErrorHandler, Environment
from Model import Model
import time
from TimeMe import TimeMe

class InnerOutputHandler(OutputHandler):
    def isInvalidOption(msg):
        for error in ["Unknown option", "invalid key"]:
            if error in msg:
                return True
        return False
    def __init__(self, appendError, printOutput = False, storeOutput=False):
      self._printOutput = printOutput
      self._appendError = appendError
      self.file = None
      self._storeOutput = storeOutput
      if storeOutput:
        self._msgs = list()
        self._kinds = list()
    def clear_output(self):
        self._msgs = []
        self._kinds = []
    def get_output(self) -> str:
       return "".join(self._msgs)
    def append_error_to_list(self, msg: str):
       self._msgs.append(msg)
       self._kinds.append("ERROR")
    def __del__(self):
       self.close_log()
    def is_recording_output(self) -> bool:
        return self._storeOutput
    
    def set_log_file(self, filename: str):
       try:
        self.file = open(filename , 'w')
       except IOError:
            print("Error: Could not open file.")
    def close_log(self):
       if self.file:
           self.file.close()
           self.file = None
    def output(self, kind, msg):
        # Check if the message contains errors coming from the solver, 
        # that are passed as normal messages from AMPL but could help
        # with better error detection
      if InnerOutputHandler.isInvalidOption(msg):
         self._appendError(msg)
      if self.file:
         self.file.write(msg)
      if self._storeOutput:
        self._kinds.append(kind)
        self._msgs.append(msg)
      if self._printOutput:
         print(msg, flush=True)
      pass
    def getMessages(self):
      return self._msgs
    def getKinds(self):
      return self._kinds

class InnerErrorHandler(ErrorHandler):
  def __init__(self, appendError):
    self._appendError=appendError
  def error(self, exc):
    self._appendError(exc)
    pass
  def warning(self, exc): 
    pass

class AMPLRunner(object):

    def __init__(self, solver=None, optionsExtra=None, writeSolverName = False,
                 printOutput=False,
                 storeOutput = False):
        self.isBenchmark = False
        if solver:
            self.setSolver(solver)
        else:
            self._solver = None
        self._writeSolverName = writeSolverName
        self._amplInitialized = False
        self._storeOutput = storeOutput
        self._printOutput = printOutput
        self._optionsExtra = optionsExtra
        self._logFile = None
        self._version = None
        
    def _initAMPL(self, model):
        if self._amplInitialized:
          self._ampl.reset()
          self._ampl.eval("reset options;")
          self._setSolverInAMPL(model)
          self._outputHandler.clear_output()
          return
        if self.isBenchmark: # Issues with non-server licenses
            time.sleep(.5)   # so wait until the license is released
        # To do: revert back to x-ampl when issue 162 in escrow is fixed
        e = Environment(binary_name="ampl")
        self._ampl = AMPL(e)
        self._outputHandler = InnerOutputHandler(self.appendError, printOutput=self._printOutput, 
                                                 storeOutput=self._storeOutput)
        if self._logFile is not None:
            self._outputHandler.set_log_file(self._logFile)
        self._ampl.setOutputHandler(self._outputHandler)
        self._ampl.setErrorHandler(InnerErrorHandler(self.appendError))
        solver_msg = self._logFile is not None or self._printOutput or self._storeOutput
        self._ampl.setOption("solver_msg", solver_msg)
        if self._solver:
          self._setSolverInAMPL(model)
        self._amplInitialized = True
   
    def _terminateAMPL(self):
      self._ampl.close()
      self._amplInitialized = False

    def appendError(self, exception):
      self._lastError = exception
      if self._outputHandler.is_recording_output():
         self._outputHandler.append_error_to_list(str(exception))

    def readModel(self, model: Model):
        mp = Path(model.getFilePath())
        if not mp.exists():
            raise Exception("Model {} not found".format(model.getFilePath()))
        if model.isScript():
          self._ampl.eval("include '{}';".format(str(mp.absolute().resolve())))
        else:
          self._ampl.read(str(mp.absolute().resolve()))
          files = model.getAdditionalFiles()
          if files:
              for f in files:
                  fp = Path(f)
                  if fp.suffix == ".dat":
                      self._ampl.readData(str(f.absolute().resolve()))
                  else:
                      self._ampl.read(str(f.absolute().resolve()))
        return mp

    def writeNL(self, model, outdir=None):
        """ Write an NL file corresponding to the specified model. 
            By default it writes in the model directory, unless outdir is specified"""
        print(f"Opening {model.getName()}.", end=" ", flush=True)
        self.doInit(model)
        self.setupOptions(model)
        mp = self.doReadModel(model)
        if outdir:
            dir = str(Path(outdir).absolute().resolve())
        else:
            dir = str(mp.parent.resolve())

        self._ampl.cd(dir)
        nlname = model.getName()
        print(f"Writing NL file {nlname}.nl.", end=" ", flush=True)
        self._ampl.eval("write 'g{}';".format(nlname))
        print("Done.", flush=True)
        self._terminateAMPL()

    def setSolver(self, solver: Solver):
        self._solver = solver
        sp = which(solver.getExecutable())
        if sp is None:
            raise Exception("Solver '{}' not found.".
                format(solver.getExecutable()))

    def _setSolverInAMPL(self, model):
        sp = self._solver.getExecutable()
        
        self._ampl.setOption("solver", sp)
        if not self._solver.get_version():
           self._solver.set_version(self._getSolverVersion(sp))
        (name, value) = self._solver.getAMPLOptions(model)
        if self._logFile is not None:
            logOption = self._solver.setLogFile(self._logFile)
            if logOption:
                value += f" {logOption}"
        self._ampl.setOption(name, value)
        
    def _getSolverVersion(self, solver) -> tuple:
        """Returns driver version (e.g. 20240316), version string"""
        command = f"shell \"{solver} -v\";"
        print(f"\nExecuting {command}")
        v=self._ampl.get_output(command)
        print(f"\nGotten output: {v}")
        driver_version = None
        import re
        for l in v.splitlines():
            driver_match = re.search(r'driver\((\d+)\)', l)
            if driver_match:
                driver_version = driver_match.group(1)
                print(f"Matched: {driver_version} in \n{l}")
                all_string=l
                print(f"Returning \ndriver_version={driver_version}\nall_string={all_string}")
                break
        if not driver_version:
           driver_version=v
           all_string=v
        return (driver_version, all_string)


    def tryGetObjective(self):
      try: 
          return self._ampl.getCurrentObjective().value()
      except:
          try:
             # Get the first objective
             objs = self._ampl.getObjectives()
             obj = next(objs)[1]
             try:
               # Try as scalar objective
               return obj.value()
             except Exception as e: 
               pass
             for index, instance in obj:
                 # Return the first instance if indexed
                 return instance.value()
          except Exception as e:
            pass
      return None

    def doInit(self, model: Model):
        self.stats = { "solver": self._solver.getName() if  self._solver else "" }
        self._initAMPL(model)
        self._lastError = None

    def printInitMessage(self, model: Model):
      msg = "Generating and solving {}".format(model.getName())
      if self._solver:
        if  self._writeSolverName:
            msg += " with {}".format(self._solver.getName())
        if self._solver.getNThreads():
            msg += " using {} threads".format(
                self._solver.getNThreads())
      print( "{0: <80}".format(msg), end="", flush=True)

    def doReadModel(self, model:Model):
        mp = str(Path(model.getFilePath()).parent.absolute().resolve())
        if model.isScript() or model.doCd():
            self._ampl.cd(mp)
        return self.readModel(model)
    
    def get_output(self) -> str:
     return self._outputHandler.get_output()
    
    def runAndEvaluate(self, model: Model, logFile : str):
        self._run(model, logFile)
        self._evaluateRun(model)
        # Todo: check bug in the API (or in amplpy) by which old objectives are 
        # reported after reset. Terminating AMPL each time takes care of the problem
        # but it's hardly efficient
        # self._terminateAMPL()

    def _run(self, model: Model,  logFile : str = None):
      self._logFile = logFile
      self.doInit(model)
      self.setupOptions(model)
      if self.isBenchmark:
         print("\n\t\t{0: <20}: Reading... ".format(self._solver.getName()), flush=True, end="")
      amplStats = { "AMPLreadTime" : "-",
                   "AMPLgenerationTime" : "-",
                   "AMPLsolveTime" : "-"
                   }
      self.stats["AMPLstats"] = amplStats
      t = TimeMe()
      t.tick()
      mp = self.doReadModel(model)
      amplStats["AMPLreadTime"]= t.toc()
      
      t.tick()
      try:
          if self.isBenchmark:
                print("Generating... ", flush=True, end="")
          ncontvars = self._ampl.getValue("_nvars - _snbvars - _snivars")
          nintvars = self._ampl.getValue("_snbvars  + _snivars")
          nconstr =  self._ampl.getValue("_ncons")
          nnz  = self._ampl.getValue("_snzcons")
      except:
          ncontvars = 0
          nintvars = 0
          nconstr =  0
          nnz  = 0
      amplStats["AMPLgenerationTime"] = t.toc()
      
      self.stats["modelStats"] = {"nvars" : int(ncontvars),
                                   "nintvars" : int(nintvars),
                                   "nconstr" : int(nconstr),
                                   "nnz" : int(nnz)}
      solve_result = None
      if not model.isScript():
          if self.isBenchmark:
            print("Solving... ", end="")
          t.tick()
          self._ampl.solve()
          solve_result = self._ampl.get_value("solve_result")
          if solve_result != "solved":
              print("WARNING: not solved (solve_result: {})".format(solve_result))
          amplStats["AMPLsolveTime"]= t.toc()
      # Close log before any further amplapi operations
      self._outputHandler.close_log()
      
      # Gather results
      self.stats["AMPLstats"] = amplStats
      if not solve_result:
         solve_result = self._ampl.getValue("solve_result")
      self.stats["solutionTime"] = self._ampl.getValue("_solve_elapsed_time")
      if self._solver: # Get times suffix if supported
        if self._solver.support_times():
           try:
              solver_time =self._ampl.get_value("Initial.time_solver")
              setup_time=self._ampl.get_value("Initial.time_setup")
              self.stats["times"]= {
                 "setup" : setup_time,
                 "solver": solver_time
                  }
           except:
              pass
      self.stats["objective"] = self.tryGetObjective()
      if self._lastError:
        self.stats["outmsg"] = str(self._lastError)
        self.stats["timelimit"] = solve_result
        self.stats["errormsg"] = self._lastError
      else:
        self.stats["outmsg"] = self._ampl.getValue("solve_message")
        self.stats["timelimit"] = solve_result
      return

    def setupOptions(self, model: Model):
        if not self._solver: # useful for writeNL
            return
        (slvname, slvval) = self._solver.getAMPLOptions(model)         # slvname = 'gurobi_options' e.g.
        if model.hasOptions():
            optmap = model.getOptions()
            for name, val in optmap.items():
                if name.endswith("SOLVER_options"):               # Any-solver option
                    if not slvname in optmap:                     # When no 'gurobi_options'
                        name = slvname
                    else:
                        continue                                  # Skip as solver-specific given
                if slvname==name:
                    slvval = slvval + ' ' + val                   # Prepend 'default' options like nthreads
                else:
                    self._ampl.setOption(name, val)
        if self._optionsExtra:
            slvval = slvval + ' ' + self._optionsExtra
        if (len(slvval)>0):
            self._ampl.setOption(slvname, slvval)

    def _evaluateRun(self, model: Model):
        expsol = model.getExpectedObjective()
        if expsol is not None:
            self.stats["eval_done"] = True
            self._assertAndRecord(expsol, self.stats["objective"],
                                  "objective")
        if model.hasExpectedValues():
            for name, ev in model.getExpectedValues().items():
                self.stats["eval_done"] = True
                try:
                    val = self._ampl.getValue(name)
                    self._assertAndRecord(ev, val,
                        "value of entity '{}'".format(name))
                except:
                    self.stats["eval_fail_msg"] = "error retrieving '{}'".format(name)

    def _assertAndRecord(self, expval, val, msg):
        b1 = isinstance(expval, (int, float))
        b2 = isinstance(val, (int, float))
        uneq = not \
            math.isclose(expval, val, rel_tol=1e-5, abs_tol=1e-5) \
            if b1 and b2 else expval != val
        if uneq:
            self.stats["eval_fail_msg"] = msg + \
                ": value " + str(val) + \
                ", expected " + str(expval)

    def getName(self):
        return "ampl-" + self._solver.getName()

    def getExecutable(self):
        return self._solver._exePath

    def getSolver(self):
        return self._solver

    def getSolutionStats(self):
        return self.stats

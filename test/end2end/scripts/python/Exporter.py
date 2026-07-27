from ModelRunner import ModelRunner
import math
import platform
from Solver import Solver


class Exporter(object):
    """Base class to export ModelRunner results"""
    def __init__(self):
        self.collector_=None
        

    def get_last_progress(self) -> list:
        return None
    
    def export(self):
        """To be called when finished running models to finalize"""
        if self.collector_: self.collector_.export(
            Exporter.get_platform_string())
        self._export()
        

    def initialize_collector(self, mr: ModelRunner):
        runs=mr.getRuns()
        runners = mr.getRunners()
        
        solvers_runs = zip(runners, runs)

        for (runner, solver_run) in solvers_runs:
            last_run=solver_run[-1]
            solver=last_run["solver"]
            
            if isinstance(solver, Solver):
                solvername=solver.getName()
            else:
                solvername=solver
                solver=runner
                
            self.collector_.add_solver_def(solvername, solver.get_version())
             
    def exportInstanceResults(self, mr: ModelRunner):
         
        i = len( mr.getRuns()[0] )
        if self.collector_ and i == 1:
            self.initialize_collector(mr)
        
        m = mr.getModels()[i-1]
        for r in mr.getRuns():
            last_run=r[-1]
            solver=last_run["solver"]
            if "Skipped" in last_run["outmsg"]:
                rescollector="Skipped"
            elif "script failure" in last_run["outmsg"]:
                rescollector="AMPL Failure"
            elif "eval_fail_msg" in last_run:
                rescollector="Failure"
            else:
                rescollector = "OK"
            last_run["parsed_result"]=rescollector
            time= {"solutionTime" : r[-1].get("solutionTime", 0)}

            if "times" in r[-1]:
                for p in ["setup", "solver"]:  
                    time[p] = r[-1]["times"].get(p, 0)

            if isinstance(solver, Solver):
                solver=solver.getName()
            if self.collector_:
                self.collector_.add_results(solver, m.getName(), rescollector, time)
        self._exportInstanceResults(mr)
    
    def _exportInstanceResults(self, mr: ModelRunner):
        raise NotImplementedError("Not implemented in base class")
    
    def _export(self):
        pass
    

    # Return False if failed
    def printStatus(self, model, run):
        if "eval_done" in run:                # Evaluation happened
            if "eval_fail_msg" not in run:    # Evaluation ok
                print("  Ok.", end='', flush=True)
            else:
                print("\n  FAILED({}):\n    - {}\n{: <68}".format(
                    run["solver"], "\n    - ".join(run["eval_fail_msg"]), ' '),
                    end='', flush=True)
                return False
        else:
            print("    <no cmp data>", end='')
        return True

    def assignFile(self, filename: str):
        self._fileName=filename
        self.__init__(filename)
        
    def add_statistic_collector(self, collector):
        self.collector_=collector
        
    def get_platform_string():
        system = platform.system()
        machine = platform.machine()

        if system == "Windows":
            if "64" in machine:
                return "win64"
            else:
                return "win32"
        elif system == "Linux":
            if "aarch64" in machine:
                return "linuxaarch64"
            elif "64" in machine:
                return "linux64"
        elif system == "Darwin":  # MacOS
            return "osx64"
        return "unknown"

class CSVTestExporter(Exporter):
    def __init__(self, fileName=""):
        super().__init__()
        self._fileName = fileName

    def sanifyString(self, s):
        try:
            s = s.replace("\n", " - ")
            s = s.replace(",", ";")
            return s
        except:
          return s

    def _getHeader(self, mr: ModelRunner):
        hasModelStats= False
        hasDriverStats=False
        for r in mr.getRuns():
            if "modelStats" in r[-1]:
                hasModelStats=True
                break
        for r in mr.getRuns():
            if "AMPLstats" in r[-1]:
                hasDriverStats=True
                break
        hdr = "Name,\tExpected_Obj,\tVariables,\tInt_Variables,\tConstraints,\tNnz"
        for (i,r) in enumerate(mr.getRuns()):
            sname = r[-1]["solver"]
            hdr += f",\t{sname}-Obj,\t{sname}-Time,\t{sname}-TimeLimit,\t{sname}-rss,\t{sname}-vms"
            if hasDriverStats:
                hdr += f",\t{sname}-AMPLreadtime,\t{sname}-AMPLgenerationtime,\t{sname}-AMPLsolvetime"
        for (i,r) in enumerate(mr.getRuns()):
            hdr += ",\t{}-SolverMsg".format(r[-1]["solver"])
        return hdr
    def getModelsStats(run):
        if not "modelStats" in run[-1]:
            return None
        stats = run[-1]["modelStats"]
        res = ",\t{},\t{},\t{},\t{}".format(
            stats["nvars"], stats["nintvars"], stats["nconstr"], stats["nnz"])
        return res
    def getAMPLStats(run):
        if not "AMPLstats" in run[-1]:
            return ",\t-,\t-,\t-"
        stats = run[-1]["AMPLstats"]
        res = ",\t{},\t{},\t{}".format(
            stats["AMPLreadTime"], stats["AMPLgenerationTime"], stats["AMPLsolveTime"])
        return res

    def _getLastResultLine(self, mr: ModelRunner):
        i = len( mr.getRuns()[0] )
        m = mr._models[i-1]
        res = "{},\t{}".format(m.getName(), m.getExpectedObjective())
        # try writing stats - works only for AMPL runners
        writtenStats = False
        for (i,r) in enumerate(mr.getRuns()):
            stats = CSVTestExporter.getModelsStats(r)
            if stats != None:
                writtenStats = True
                res += stats
                break
        if not writtenStats:
            res += ",\t-,\t-,\t-,\t-"
        for (i,r) in enumerate(mr.getRuns()):
            res += ",\t{},\t{},\t{},\t{},\t{}".format(
              self._getDictMemberOrMissingStr(r[-1], "objective"),
              self._getDictMemberOrMissingStr(r[-1], "solutionTime"),
              self._getDictMemberOrMissingStr(r[-1], "timelimit"),
              self._getDictMemberOrMissingStr(r[-1], "rss"),
              self._getDictMemberOrMissingStr(r[-1], "vms"))
            res += CSVTestExporter.getAMPLStats(r)

        for (i,r) in enumerate(mr.getRuns()):
            res += ",\t{}".format(
              self.sanifyString(
                self._getDictMemberOrMissingStr(r[-1], "outmsg")))
        return res
                
    def _getDictMemberOrMissingStr(self, dct, key):
        try:
            return dct[key]
        except:
            return "-"

    def _exportInstanceResults(self, mr: ModelRunner):
        i = len( mr.getRuns()[0] )
        filemode = "w" if i == 1 else "a+"
        with open(self._fileName, filemode) as file:
                if i == 1:
                    file.write( "{}\n".format( self._getHeader(mr) ) )
                file.write( "{}\n".format( self._getLastResultLine(mr) ) )


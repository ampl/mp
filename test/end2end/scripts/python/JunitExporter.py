from junitparser import TestCase, TestSuite, JUnitXml,  Failure, Skipped
from Exporter import Exporter
import Solver
import ModelRunner
from os.path import splitext
import xml.sax.saxutils
import platform

class JunitExporter(Exporter):
    def __init__(self, fileName=""):
        self._fileName=fileName
        self._test_suites={}
        self._versions={}
        self.collector_=None
        
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

    def get_file_name(self, solver: str):
        base_name, _= splitext(self._fileName)
        return f"{solver}-{base_name}.xml"
    
    def exportInstanceResults(self, mr: ModelRunner):
        i = len( mr.getRuns()[0] )
        if i == 1:
            self.initialize(mr)
        self.append_last_results(mr)
        self.save()
        
    def initialize(self, mr: ModelRunner):
        runs=mr.getRuns()
        runners = mr.getRunners()
        
        solvers_runs = zip(runners, runs)

        for (runner, solver_run) in solvers_runs:
            last_run=solver_run[-1]
            solver=last_run["solver"]
            
            if isinstance(solver, Solver.Solver):
                solvername=solver.getName()
            else:
                solvername=solver
                solver=runner
                
            if self.collector_:
                self.collector_.add_solver_def(solvername, solver.get_version())
             
            name = f"{solvername}-{JunitExporter.get_platform_string()}"
            self._test_suites[solvername]= TestSuite(name)
            
        
    def append_last_results(self, mr: ModelRunner):
      i = len( mr.getRuns()[0] )
      m = mr.getModels()[i-1]
      res = [m.getName(), m.getExpectedObjective()]
      for r in mr.getRuns():
            last_run=r[-1]
            solver=last_run["solver"]
            if isinstance(solver, Solver.Solver):
                solver=solver.getName()

            time= {"solutionTime" : r[-1].get("solutionTime", 0)}
            if "times" in r[-1]:
                for p in ["setup", "solver"]:  
                 time[p] = r[-1]["times"].get(p, 0)

            tc=TestCase(res[0], time=time["solutionTime"])
            if "Skipped" in last_run["outmsg"]:
                tc.result=[Skipped(last_run["outmsg"])]
                res="Skipped"
            elif "eval_fail_msg" in last_run:
                safe_string = xml.sax.saxutils.escape(last_run["output"])
                safe_string =safe_string.replace('\b', '')
                tc.system_out=safe_string
                tc.result=[Failure(last_run["eval_fail_msg"])]
                res="Failure"
            else:
                res = "OK"
                
            self._test_suites[solver].add_testcase(tc)
            
            if self.collector_:
               self.collector_.add_results(solver, res, time) 
                   
    def save(self):
      for solver,ts in self._test_suites.items():
        xml = JUnitXml()
        xml.add_testsuite(ts)
        xml.write(self.get_file_name(solver), pretty=True)
    def export(self):
        if self.collector_: self.collector_.export(
            JunitExporter.get_platform_string())
            
    def add_statistic_collector(self, collector):
        self.collector_=collector
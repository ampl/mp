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
            if "64" in machine:
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
        for solver_run in runs:
            last_run=solver_run[-1]
            solver=last_run["solver"]
            if isinstance(solver, Solver.Solver):
                solver=solver.getName() 
            name = f"{solver}-{JunitExporter.get_platform_string()}"
            self._test_suites[solver]= TestSuite(solver)
        
    def append_last_results(self, mr: ModelRunner):
      i = len( mr.getRuns()[0] )
      m = mr.getModels()[i-1]
      res = [m.getName(), m.getExpectedObjective()]
      for r in mr.getRuns():
            last_run=r[-1]
            solver=last_run["solver"]
            if isinstance(solver, Solver.Solver):
                solver=solver.getName()
            tc=TestCase(res[0], time= r[-1].get("solutionTime", 0))
            if "Skipped" in last_run["outmsg"]:
                tc.result=[Skipped(last_run["outmsg"])]
            if "eval_fail_msg" in last_run:
                safe_string = xml.sax.saxutils.escape(last_run["output"])
                safe_string =safe_string.replace('\b', '')
                tc.system_out=safe_string
                tc.result=[Failure(last_run["eval_fail_msg"])]
            
            
            self._test_suites[solver].add_testcase(tc)
                   
    def save(self):
      
      for solver,ts in self._test_suites.items():
        xml = JUnitXml()
        xml.add_testsuite(ts)
        xml.write(self.get_file_name(solver), pretty=True)
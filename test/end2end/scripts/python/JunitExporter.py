from junitparser import TestCase, TestSuite, JUnitXml,  Failure, Skipped
from Exporter import Exporter
import Solver
import ModelRunner
from os.path import splitext
import xml.sax.saxutils

class JunitExporter(Exporter):
    def __init__(self, fileName=""):
        self._fileName=fileName
        self._test_suites={}
        

    def get_file_name(self):
        base_name, _= splitext(self._fileName)
        return base_name + '.xml'
    
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
                safe_string =safe_string .replace('\b', '')
                tc.system_out=safe_string
                tc.result=[Failure(last_run["eval_fail_msg"])]
            
            
            self._test_suites[solver].add_testcase(tc)
                   
    def save(self):
      xml = JUnitXml()
      for ts in self._test_suites.values():
        xml.add_testsuite(ts)
      xml.write(self.get_file_name(), pretty=True)

      # for res in results:
      #   n = res[0]
      #   r = res[1]
      #   tc=TestCase(n, time= r[1])
      #   if r[0] != 0:
      #     outputlog = r[2].replace("\n", "")
      #     failuremessage = f"Exit code: {r[0]}"
      #     if len(outputlog) > 0:
      #       failuremessage += f"\n{outputlog}"
      #     else:
      #       failuremessage += f"\nstdout and stderr empty."
      #     tc.result=[Failure(failuremessage, "Runtime failure")]
      #   ts.add_testcase(tc)
      # xml = JUnitXml()
      # xml.add_testsuite(ts)
      # xml.write(self.get_file_name(), pretty=True)
  
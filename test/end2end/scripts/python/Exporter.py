from ModelRunner import ModelRunner
import math


class Exporter(object):
    """Base class to export ModelRunner results"""

    def export(self, mr: ModelRunner):
        raise Exception("Not defined for base class")

    # Return False if failed
    def printStatus(self, model, run):
        if "eval_done" in run:                # Evaluation happened
            if "eval_fail_msg" not in run:    # Evaluation ok
                print("  Ok.", end='')
            else:
                print("          FAILED({}): {}\n{: <80}".format(
                    run["solver"], run["eval_fail_msg"], ' '),
                    end='', flush=True)
                return False
        else:
            print("    <no cmp data>", end='')
        return True


class CSVTestExporter(Exporter):
    def __init__(self, fileName):
        self._fileName = fileName

    def sanifyString(self, s):
        try:
            s = s.replace("\n", " - ")
            s = s.replace(",", ";")
            return s
        except:
          return s

    def _getHeader(self, mr: ModelRunner):
        hdr = "Name,\tExpected_Obj"
        for (i,r) in enumerate(mr.getRuns()):
            hdr += ",\t{}-Obj,\t{}-Time,\t{}-TimeLimit".format(
              r[-1]["solver"], r[-1]["solver"], r[-1]["solver"])
        for (i,r) in enumerate(mr.getRuns()):
            hdr += ",\t{}-SolverMsg".format(r[-1]["solver"])
        return hdr
            
    def _getLastResultLine(self, mr: ModelRunner):
        i = len( mr.getRuns()[0] )
        m = mr._models[i-1]
        res = "{},\t{}".format(m.getName(), m.getExpectedObjective())
        for (i,r) in enumerate(mr.getRuns()):
            res += ",\t{},\t{},\t{}".format(
              self._getDictMemberOrMissingStr(r[-1], "objective"),
              self._getDictMemberOrMissingStr(r[-1], "solutionTime"),
              self._getDictMemberOrMissingStr(r[-1], "timelimit"))
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

    def exportInstanceResults(self, mr: ModelRunner):
        i = len( mr.getRuns()[0] )
        filemode = "w" if i == 1 else "a+"
        with open(self._fileName, filemode) as file:
                if i == 1:
                    file.write( "{}\n".format( self._getHeader(mr) ) )
                file.write( "{}\n".format( self._getLastResultLine(mr) ) )


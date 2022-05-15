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
                print("  Ok.", end='', flush=True)
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

    def exportInstanceResults(self, mr: ModelRunner):
        i = len( mr.getRuns()[0] )
        filemode = "w" if i == 1 else "a+"
        with open(self._fileName, filemode) as file:
                if i == 1:
                    file.write( "{}\n".format( self._getHeader(mr) ) )
                file.write( "{}\n".format( self._getLastResultLine(mr) ) )


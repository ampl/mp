from Solver import Solver
from AMPLRunner import AMPLRunner
from TimeMe import TimeMe
from Model import Model

class ModelRunner(object):
    """Class to run a set of models and capture their outputs"""

    def __init__(self, runners, optionsExtra=None):
        self._runners = runners
        self._amplRunners = None
        self._runs = [ list() for r in self._runners ]
        self._optionsExtra = optionsExtra

    def getRuns(self):
        return self._runs

    def getModels(self) -> list:
        return self._models
    def getRunners(self)->list:
        return self._runners
    def getLogFileName(m : Model, s : Solver):
        return f"{m.getName()}.{s.getName()}.log"

    def run(self, modelList: list, exporter=None, keepLogs = False, verbose=False):
        """Run the models in this instance. If exporter != None, it exports the results as it goes"""
        self._models = modelList
        n = 0
        nFailed = 0
        nSkipped = 0
        
        # Set to true for junit to add the solver output to the test result
        keep_output=False
        try:
        # Attempt to access the class
            __import__("JunitExporter")
            import JunitExporter
            if isinstance(exporter, JunitExporter.JunitExporter):
                keep_output=True
        except:
            pass
        
        for m in modelList:
            n += 1
            if m.isNL():
                cr = self._runners
                msg = "{}. Solving as NL: '{}'".format(n, m.getName())
            else:
                if not self._amplRunners:
                    self._amplRunners = [
                        AMPLRunner(r, self._optionsExtra, 
                                   printOutput=verbose, storeOutput=keep_output)
                        for r in self._runners ]
                cr = self._amplRunners
                msg = "{}. Solving with AMPL: '{}'".format(n, m.getName())
            print("{0: <80}".format(msg), end="", flush=True)
            t = TimeMe()
            failedSome = False
            skippedSome = False
            with t:
                for (i,r) in enumerate(cr):
                    r.isBenchmark = len(cr) > 1
                    self._runs[i].append({})
                    if isinstance(r, AMPLRunner):
                        ss = r.getSolver()
                    else:
                        ss = r
                    if not m.isSubsetOfTags(ss.getSupportedTags()) or \
                            m.hasAnyTag(ss.getUnsupportedTags()):
                        self._runs[i][-1]["outmsg"] = "Skipped, unsupported tags"
                        self._runs[i][-1]["solver"] = ss
                        if r.isBenchmark:
                            print("\n\t\tSkipped due to unsupported tags", flush=True)
                        else:
                            print("Skipped due to unsupported tags", flush=True)
                            skippedSome = True
                        continue
                    if keepLogs:
                        r.runAndEvaluate(m, logFile=ModelRunner.getLogFileName(m, ss))
                    else:
                        r.runAndEvaluate(m, logFile=None)
                    
                    stats = r.getSolutionStats()
                    EFM ="eval_fail_msg"
                    EM ="errormsg"
                    if EFM in stats:
                        if EM in stats:
                            stats[EFM]= stats[EM]
                    if keep_output:
                         stats["output"]=r.get_output()
                    self._runs[i][-1] = stats
                    if exporter:
                        if not exporter.printStatus(m, stats):
                            failedSome = True
            nFailed += failedSome
            nSkipped += skippedSome
            if exporter:
                self.export(exporter)
            print("  (%.4fs, %d failed, %d skipped)" % (t.interval, nFailed, nSkipped))
        exporter.export()

    def export(self, exporter):
        exporter.exportInstanceResults(self)


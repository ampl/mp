from Solver import Solver
from AMPLRunner import AMPLRunner
from TimeMe import TimeMe


class ModelRunner(object):
    """Class to run a set of models and capture their outputs"""

    def __init__(self, runners, optionsExtra=None):
        self._runners = runners
        self._amplRunners = None
        self._runs = [ list() for r in self._runners ]
        self._optionsExtra = optionsExtra
        
    def getRuns(self):
        return self._runs

    def run(self, modelList: list, exporter=None):
        """Run the models in this instance. If exporter != None, it exports the results as it goes"""
        self._models = modelList
        n = 0
        nFailed = 0
        for m in modelList:
            n += 1
            if m.isNL():
                cr = self._runners
                msg = "{}. Solving as NL: '{}'".format(n, m.getName())
            else:
                if not self._amplRunners:
                    self._amplRunners = [
                        AMPLRunner(r, self._optionsExtra) for r in self._runners ]
                cr = self._amplRunners
                msg = "{}. Solving with AMPL: '{}'".format(n, m.getName())
            print("{0: <80}".format(msg), end="", flush=True)
            t = TimeMe()
            failedSome = False
            lastModel = None
            with t:
                for (i,r) in enumerate(cr):
                    self._runs[i].append({})
                    if isinstance(r, AMPLRunner):
                        ss = r.getSolver()
                    else:
                        ss = r
                    if m.hasAnyTag(ss.getUnsupportedTags()):
                        self._runs[i][-1]["outmsg"] = "Skipped, unsupported tags"
                        self._runs[i][-1]["solver"] = ss
                        print("Skipped due to unsupported tags")
                        continue
                    r.runAndEvaluate(m)
                    stats = r.getSolutionStats()
                    self._runs[i][-1] = stats
                    if exporter:
                        if not exporter.printStatus(m, stats):
                            failedSome = True
            nFailed += failedSome
            if exporter:
                self.export(exporter)
            print("  (%.4fs, %d failed)" % (t.interval, nFailed))

    def export(self, exporter):
        exporter.exportInstanceResults(self)


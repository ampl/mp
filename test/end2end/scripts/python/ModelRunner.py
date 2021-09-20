from Solver import Solver
from AMPLRunner import AMPLRunner
from TimeMe import TimeMe


class ModelRunner(object):
    """Class to run a set of models and capture their outputs"""

    def __init__(self, runners):
        self._runners = runners
        self._amplRunners = None

    def run(self, modelList: list, exporter=None):
        """Run the models in this instance. If exporter != None, it exports the results as it goes"""
        self._models = modelList
        self._runs = [ list() for r in self._runners ]
        n = 0
        nFailed = 0
        for m in modelList:
            n += 1
            if m.isNL():
                cr = self._runners
                msg = "{}. Solving as NL: '{}'".format(n, m.getName())

            else:
                if not self._amplRunners:
                    self._amplRunners = [ AMPLRunner(r) for r in self._runners ]
                cr = self._amplRunners
                msg = "{}. Solving with AMPL: '{}'".format(n, m.getName())
            print("{0: <80}".format(msg), end="", flush=True)
            t = TimeMe()
            failedSome = False
            with t:
                for (i,r) in enumerate(cr):
                    if isinstance(r, AMPLRunner):
                      ss = r.getSolver()
                    else:
                      ss = r
                    if m.hasAnyTag(ss.getUnsupportedTags()):
                      print("Skipped due to unsupported tags")
                      continue
                    r.runAndEvaluate(m)
                    stats = r.getSolutionStats()
                    self._runs[i].append(stats)
                    if exporter:
                        self.export(exporter)
                        if not exporter.printStatus(m, stats):
                            failedSome = True
            nFailed += failedSome
            print("  (%.4fs, %d failed)" % (t.interval, nFailed))

    def export(self, exporter):
        exporter.export(self)


class ModelComparer(ModelRunner):
    """Class to run a set of models and capture their outputs"""

    def __init__(self, runner1, runner2):
        super().__init__(runner1)
        self._runner2 = runner2
        self.amplRunner2 = None

    def run(self, modelList: list, exporter=None):
        self._models = modelList
        self._runs = list()
        self._runs2 = list()
        self._runners._writeSolverName=True
        self._runner2._writeSolverName=True
        for m in modelList:
            if m.isNL():
              cr1 = self._runners
              cr2 = self._runner2
            else:
                if not self._amplRunners:
                  self._amplRunners = AMPLRunner(writeSolverName = True)
                  self._amplRunners.setSolver(self._runners)
                  self._amplRunner2 = AMPLRunner(writeSolverName = True)
                  self._amplRunner2.setSolver(self._runner2)
                cr1 = self._amplRunners
                cr2 = self._amplRunner2
            cr1.runAndEvaluate(m)
            stat1=cr1.getSolutionStats()
            if not exporter.printStatus(m, stat1):
                            failedSome = True
            cr2.runAndEvaluate(m)
            stat2=cr2.getSolutionStats()
            if not exporter.printStatus(m, stat2):
                            failedSome = True
            self._runs.append(stat1)
            self._runs2.append(stat2)
            exporter.exportOne(self)
           

    def getRunnerNames(self):
        return (self._runners.getName(), self._runner2.getName())

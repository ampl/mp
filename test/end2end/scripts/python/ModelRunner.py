from Solver import Solver
from AMPLRunner import AMPLRunner


class ModelRunner(object):
    """Class to run a set of models and capture their outputs"""

    def __init__(self, runner):
        self._runner = runner
        self._amplRunner = None

    def run(self, modelList: list, exporter=None):
        """Run the models in this instance. If exporter != None, it exports the results as it goes"""
        self._models = modelList
        self._runs = list()
        n = 0
        for m in modelList:
            if m.isNL():
                cr = self._runner

            else:
                if not self._amplRunner:
                    self._amplRunner = AMPLRunner()
                    self._amplRunner.setSolver(self._runner)
                cr = self._amplRunner
            n += 1
            print("{}. ".format(n), end='')
            cr.runAndEvaluate(m)
            self._runs.append(cr.getSolutionStats())
            if exporter:
                exporter.exportOne(self)
            print("")

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
        self._runner._writeSolverName=True
        self._runner2._writeSolverName=True
        for m in modelList:
            if m.isNL():
              cr1 = self._runner
              cr2 = self._runner2
            else:
                if not self._amplRunner:
                  self._amplRunner = AMPLRunner(writeSolverName = True)
                  self._amplRunner.setSolver(self._runner)
                  self._amplRunner2 = AMPLRunner(writeSolverName = True)
                  self._amplRunner2.setSolver(self._runner2)
                cr1 = self._amplRunner
                cr2 = self._amplRunner2
            cr1.run(m)
            cr2.run(m)
            self._runs.append(cr1.getSolutionStats())
            self._runs2.append(cr2.getSolutionStats())
            exporter.exportOne(self)

    def getRunnerNames(self):
        return (self._runner.getName(), self._runner2.getName())

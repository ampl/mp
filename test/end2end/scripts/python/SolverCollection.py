# This Python file uses the following encoding: utf-8

import Solver


class SolverCollection:
    def __init__(self):
        self._solvers = {}

    def addSolver(self, solver: Solver):
        name = solver.getExecutable()
        if name in self._solvers:
            raise "Solver '{}' already defined".format(name)
        self._solvers[name] = solver

    def getSolvers(self):
        return self._solvers.items()

    def getSolverNames(self):
        return self._solvers.keys()

    def getSolver(self, name):
        return self._solvers.get(name)


def addStdSolvers(solvers: SolverCollection):
    solvers.addSolver(Solver.LindoSolver("lindoglobal-timebound"))
    solvers.addSolver(Solver.OcteractSolver("octeract-engine"))
    solvers.addSolver(Solver.GurobiSolver("gurobi"))
    solvers.addSolver(Solver.GurobiDirectSolver("gurobidirect"))
    solvers.addSolver(Solver.CPLEXSolver("cplex"))
    solvers.addSolver(Solver.BaronSolver("baron-timebound"))


# if __name__ == "__main__":
#     pass

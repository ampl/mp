# This Python file uses the following encoding: utf-8

import Solver
from os import path

class SolverCollection:
    def __init__(self):
        self._solvers = {}

    def addSolver(self, solver: Solver):
        name = solver.getName()
        if name in self._solvers:
            raise "Solver '{}' already defined".format(name)
        self._solvers[name] = solver

    def getSolvers(self):
        return self._solvers.items()

    def getSolversByNames(self, names):
        return [self._solvers[x] for x in names]

    def getSolverNames(self):
        return self._solvers.keys()

    def getSolver(self, name):
        return self._solvers.get(name)

def addStdSolvers(solvers: SolverCollection, binPath=""):
    solvers.addSolver(Solver.LindoSolver(path.join(binPath, "lindoglobal")))
    solvers.addSolver(Solver.OcteractSolver(path.join(binPath, "octeract-engine")))
    solvers.addSolver(Solver.GurobiSolver(path.join(binPath,"gurobi")))
    solvers.addSolver(Solver.GurobiDirectSolver(path.join(binPath,"x-gurobi")))
    solvers.addSolver(Solver.CPLEXSolver(path.join(binPath,"cplex")))
    solvers.addSolver(Solver.CPLEXDirectSolver(path.join(binPath,"cplexdirect")))
    
    solvers.addSolver(Solver.BaronSolver(path.join(binPath,"baron")))
    solvers.addSolver(Solver.COPTSolver(path.join(binPath,"copt")))
    solvers.addSolver(Solver.COPTSolver(path.join(binPath,"coptdirect")))
# if __name__ == "__main__":
#     pass

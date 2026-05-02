# This Python file uses the following encoding: utf-8

import Solver
from os import path

class SolverCollection:
    def __init__(self):
        self._solvers = {}
        self._aliases = {} # alias -> solvername map

    def addSolver(self, solver: Solver, aliases: list = None):
        name = solver.getName()
        if name in self._solvers:
            raise "Solver '{}' already defined".format(name)
        self._solvers[name] = solver
        self._aliases[name]=name
        if aliases:
            for a in aliases:
                self._aliases[a]=name

    def getSolvers(self):
        return self._solvers.items()

    def getSolversByNames(self, names):
        return [self.getSolver(x) for x in names]

    def getSolverNames(self):
        return self._aliases.keys()

    def getSolver(self, nameoralias: str):
        return self._solvers.get(self._aliases[nameoralias])

def addStdSolvers(solvers: SolverCollection, binPath=""):
    ### ASL solvers
    solvers.addSolver(Solver.BaronSolver(path.join(binPath,"baronasl")))
    solvers.addSolver(Solver.LindoSolver(path.join(binPath, "lindoglobalasl")))
    solvers.addSolver(Solver.KnitroSolver(path.join(binPath,"knitroasl")))
    solvers.addSolver(Solver.CouenneSolver(path.join(binPath, "couenneasl")))
    solvers.addSolver(Solver.BonminSolver(path.join(binPath, "bonminasl")))
    solvers.addSolver(Solver.ConoptSolver(path.join(binPath,"conoptasl")))
    solvers.addSolver(Solver.ConoptSolver(path.join(binPath,"conopt4")))
    solvers.addSolver(Solver.OcteractSolver(path.join(binPath, "octeract-engine")))
    solvers.addSolver(Solver.MindoptSolver(path.join(binPath,"mindoptampl")))
    solvers.addSolver(Solver.LgoSolver(path.join(binPath, "lgo")))
    solvers.addSolver(Solver.IPOptSolver(path.join(binPath, "ipopt")))
    solvers.addSolver(Solver.LoqoSolver(path.join(binPath, "loqo")))
    solvers.addSolver(Solver.XpressSolver(path.join(binPath,"xpressasl")))
    solvers.addSolver(Solver.GurobiSolver(path.join(binPath,"gurobiasl")))
    solvers.addSolver(Solver.CPLEXSolver(path.join(binPath,"cplexasl")))


    ### ASL+MP2NL solvers (run with mp2nl=1 until becomes default)
    solvers.addSolver(Solver.BaronMP2NLSolver(path.join(binPath,"baron")))
    solvers.addSolver(Solver.LindoglobalMP2NLSolver(path.join(binPath, "lindoglobal")))
    solvers.addSolver(Solver.KnitroMP2NLSolver(path.join(binPath,"knitro")))
    solvers.addSolver(Solver.CouenneMP2NLSolver(path.join(binPath, "couenne")))
    solvers.addSolver(Solver.BonminMP2NLSolver(path.join(binPath, "bonmin")))
    solvers.addSolver(Solver.ConoptMP2NLSolver(path.join(binPath,"conopt")))
    solvers.addSolver(Solver.GurobiASLMP2NLSolver(path.join(binPath,"gurobiaslmp2nl")))


    ### MP solvers
    solvers.addSolver(Solver.GurobiDirectSolver(path.join(binPath,"gurobi")))
    solvers.addSolver(Solver.KnitroMPSolver(path.join(binPath,"knitromp")))
    solvers.addSolver(Solver.CPLEXDirectSolver(path.join(binPath,"cplex"))) 
    solvers.addSolver(Solver.BaronMPSolver(path.join(binPath,"baronmp")))
    solvers.addSolver(Solver.COPTSolver(path.join(binPath,"copt")))
    solvers.addSolver(Solver.cuOptSolver(path.join(binPath,"cuoptmp")))
    solvers.addSolver(Solver.HighsSolver(path.join(binPath,"highs")), aliases=["highsmp"])
    solvers.addSolver(Solver.MP2NLSolver(path.join(binPath,"mp2nl")))
    solvers.addSolver(Solver.XPRESSDirectSolver(path.join(binPath,"xpress")))
    solvers.addSolver(Solver.MosekSolver(path.join(binPath,"mosek")))
    solvers.addSolver(Solver.CbcMPSolver(path.join(binPath, "cbc")), aliases=["cbcmp"])
    solvers.addSolver(Solver.GCGSolver(path.join(binPath, "gcg")), aliases=["gcgmp"])
    solvers.addSolver(Solver.SCIPSolver(path.join(binPath, "scip")), aliases=["scipmp"])
    solvers.addSolver(Solver.CPLEXODHSolver(path.join(binPath, "cplexodh")))
    solvers.addSolver(Solver.GUROBIODHSolver(path.join(binPath, "gurobiodh")))


# if __name__ == "__main__":
#     pass

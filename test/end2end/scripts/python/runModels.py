from ModelsDiscovery import ModelsDiscovery
from ModelRunner import ModelRunner, ModelComparer
from Exporter import CSVTestExporter
from Solver import Solver, LindoSolver, GurobiSolver, OcteractSolver, CPLEXSolver
from pathlib import Path
from sys import platform
from AMPLRunner import AMPLRunner


def writeNLFiles(directory, recursive=False):
    m = ModelsDiscovery()
    modelList = m.FindModelsGeneral(directory, recursive=recursive)
    amplRunner = AMPLRunner()
    toGenerate = filter(lambda m: not m.isNL(), modelList)
    for m in toGenerate:
        amplRunner.writeNL(m)


def runModels(directory, solver: Solver, solver2: Solver = None,
              exporter=None, exportFile=None, justNL=False,
              recursive=False, preferAMPLModels=False):
    """Convenient wrapper function for testing.

          With no optional argument specified, runs as:
          runModels(dir, solver) - executes all models in the specified directory with the specified solver,
            exporting the results to CSV in the cwd in a file with name run-dir-platform-solvername.csv

          Parameters
          ----------

          directory : str - Directory to find the models in
          solver : Solver - Solver to use, choose from the module Solver or implenent your own
          solver2 : Solver, optional - A second solver object. If specified, the models will be ran on both solvers
                                       and the report will contain all the statistics side by side
          exporter : Exporter, optional - An exporter object that overrides the default CSV one. In case this is specified, 
                                          the parameter "exportDir" is ignored
          exportFile: str, optional - Override the output file name for the default exporter
          justNL: bool     - If True, only considers the NL files in the expored directories. Useful when no AMPL license
                                   or amplpy is available
          recursive : bool - If True, finds models in the subdirectories also
          preferAMPLModels:bool - If True, executes the AMPL version of a model if both NL and AMPL versions are present.
    """
    if not exporter:
        if not exportFile:
            ename = solver.getName() if not solver2 else "{}-vs-{}".format(solver.getName(), solver2.getName())
            exportFile = "run-{}-{}-{}.csv".format(Path(directory).stem, platform, ename)
        exporter = CSVTestExporter(exportFile)
    if solver2:
        runner = ModelComparer(solver, solver2)
    else:
        runner = ModelRunner(solver)

    m = ModelsDiscovery()
    modelList = m.FindModelsGeneral(directory, recursive=recursive,
                                    preferAMPLModels=preferAMPLModels,
                                    justNL=justNL)
    runner.run(modelList, exporter)


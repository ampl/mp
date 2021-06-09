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


def runModels(directory, solvers,
              exporter=None, exportFile=None, justNL=False,
              recursive=False, preferAMPLModels=False):
    """Convenient wrapper function for testing.

          With no optional argument specified, runs as:
          runModels(dir, solver) - executes all models in the specified directory with the specified solver,
            exporting the results to CSV in the cwd in a file with name run-dir-platform-solvername.csv

          Parameters
          ----------

          directory : str - Directory to find the models in
          solver : list(Solver) - Solvers to use, choose from the module Solver or implement your own
          exporter : Exporter, optional - An exporter object that overrides the default CSV one. In case this is specified, 
                                          the parameter "exportDir" is ignored
          exportFile: str, optional - Override the output file name for the default exporter
          justNL: bool     - If True, only considers the NL files in the expored directories. Useful when no AMPL license
                                   or amplpy is available
          recursive : bool - If True, finds models in the subdirectories also
          preferAMPLModels:bool - If True, executes the AMPL version of a model if both NL and AMPL versions are present.
    """
    solvernames = [slv.getExecutable() for slv in solvers]
    if not exporter:
        if not exportFile:
            ename = "-".join(solvernames)
            exportFile = "run-{}-{}-{}.csv".format(Path(directory).stem, platform, ename)
        exporter = CSVTestExporter(exportFile)
    runner = ModelRunner(solvers)

    m = ModelsDiscovery()
    modelList = m.FindModelsGeneral(directory, recursive=recursive,
                                    preferAMPLModels=preferAMPLModels,
                                    justNL=justNL)
    if not modelList:
        print("No models or case descriptions found.")
    else:
        msg = "Running {} test cases with solvers {}".format(len(modelList), solvernames)
        if solvers:
            if solvers[0].getNThreads():
                msg += " using {} threads".format( solvers[0].getNThreads() )
        msg += '.'
        print(msg)
        runner.run(modelList, exporter)


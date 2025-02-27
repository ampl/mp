from cmath import e
from ModelsDiscovery import ModelsDiscovery
from ModelRunner import ModelRunner
from Exporter import CSVTestExporter
from Solver import Solver, LindoSolver, GurobiSolver, OcteractSolver, CPLEXSolver
from pathlib import Path
from sys import platform
from AMPLRunner import AMPLRunner
from Model import ModelTags

def writeModels(ampl, directory,modelList=True, justNL=False, recursive=False,preferAMPLModels=False, writeMPS=False):
    m = ModelsDiscovery()
    modelList = m.FindModelsGeneral(directory, recursive=recursive, modellist=modelList,
                                    preferAMPLModels=preferAMPLModels,
                                    justNL=justNL)
    amplRunner = AMPLRunner(ampl)
    if not writeMPS:
        toGenerate = filter(lambda m: not m.isNL(), modelList)
    else:
        toGenerate = modelList
    for m in toGenerate:
        amplRunner.writeModel(m, writeMPS=writeMPS)


def runModels(directory, ampl: str, solvers: list,
              solverOptions=None,
              exporter=None, exportFile=None, modellist=True, justNL=False,
              recursive=False, preferAMPLModels=False, keepLogs = False, verbose=False,
              defaultReportSuffix:str = None):
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
          keepLogs: bool - If True, store the logs of all executions; useful for debugging or for benchmarking
          verbose: bool - If True, print solver output
          defaultReportSuffix: str - A suffix to be added to the report file name, used by Tester to differentiate between
                                     executions with different params (e.g. simplex vs IPM or native vs reformulation NL support)
    """
    solvernames = [Path(slv.getExecutable()).stem for slv in solvers]
    if not exportFile:
        exportFile = "run"
    ename = "-".join(solvernames)
    suffix = f"-{defaultReportSuffix}" if defaultReportSuffix else ""
        
    exportFile += "-{}-{}-{}{}.csv".format(Path(directory).stem, platform, ename, suffix)
    if not exporter:
        exporter = CSVTestExporter(exportFile)

    exporter.assignFile(exportFile)

    runner = ModelRunner(ampl, solvers, solverOptions)

    m = ModelsDiscovery()
    
    found_models = exporter.get_last_progress()
    
    modelList = m.FindModelsGeneral(directory, recursive=recursive,
                                    modellist=modellist,
                                    preferAMPLModels=preferAMPLModels,
                                    justNL=justNL)
    if not modelList:
        print("No models or case descriptions found.")
    else:
        msg = "Running {} test cases with solvers {}".format(len(modelList), solvernames)
        if found_models is not None:
            end = "s" if len(found_models)>1 else ""
            msg +=f"\nContinuing previous run, discarding {len(found_models)} model{end}."
            modelList = [m for m in modelList if m.getName() not in found_models]

            mnames =  [m.getName() for m in modelList]
            notfound=[m for m in found_models if m not in mnames]
            if len(notfound) > 0:
                print(f"{len(notfound)} models found in previous export and not found in current list.")
                print(",".join(notfound))
            msg += "\nActually running {} test cases with solvers {}".format(len(modelList), solvernames)
        if solvers:
            if solvers[0].getNThreads():
                msg += " using {} threads".format( solvers[0].getNThreads() )
        msg += '.'
        print(msg)
        runner.run(modelList, exporter, keepLogs, verbose)

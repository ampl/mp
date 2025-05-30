import sys
import pathlib

import multiprocessing

if __name__ == '__main__':
    multiprocessing.freeze_support()

    ## Access to module files
    scriptpath = pathlib.Path(__file__).parent
    libpath = scriptpath.joinpath('scripts').joinpath('python')

    sys.path.insert(1, str(libpath))

    ## Import the tester app module
    from runExamples import runTester
    runTester()
    #from runModels import writeNLFiles
    #writeNLFiles("E:\\OneDrive\\mm-lpsimp-mps")


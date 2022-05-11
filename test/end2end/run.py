import sys
import pathlib

scriptpath = pathlib.Path(__file__).parent
libpath = scriptpath.joinpath('scripts').joinpath('python')

sys.path.insert(1, str(libpath))

from runExamples import runTester
runTester()
#from runModels import writeNLFiles
#writeNLFiles("E:\\OneDrive\\Documents\\Benchmarks\\models\\lpmodels-mps\\non-mm-mps")

# DLR2 mindopt
#Interior point method terminated. Time : 3307.930s
#
#OPTIMAL; objective 11601762.31
#2275457 simplex iterations

import sys
import pathlib

scriptpath = pathlib.Path(__file__).parent
libpath = scriptpath.joinpath('scripts').joinpath('python')

sys.path.insert(1, str(libpath))

from runExamples import runTester

runTester()

import sys
import pathlib

## Access to module files
scriptpath = pathlib.Path(__file__).parent
libpath = scriptpath.joinpath('scripts').joinpath('python')

sys.path.insert(1, str(libpath))

## Import & run the app module
from vizCvtGraphApp import runVizCvtGraphApp
runVizCvtGraphApp()


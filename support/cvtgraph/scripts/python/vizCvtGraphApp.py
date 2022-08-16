import argparse

from readCvtGraphJSONL import GraphReaderJSONL
from dataCvtGraph import CvtGraph
from drawGraph import GraphDrawer


# The Model Conversion Graph Vizualization App
class VizCvtGraphApp:
    # Constructor
    def __init__(self):
        self._parser = argparse.ArgumentParser(
            description='Model Conversion Graph Vizualizer')

    # Class entry point
    def run(self):
        self.parseOptions()
        reader = GraphReaderJSONL()
        graph = CvtGraph()
        reader.read(self._args.graphFile, graph)
        drawer = GraphDrawer(self._args.title)
        graph.drawGraph(drawer, self._args.W, self._args.H)

    # Parse cmdline options
    def parseOptions(self):
        self._parser.add_argument('graphFile', metavar='graph.jsonl', type=str,
                                  help='input file ' +
                                  '(JSON Lines, https://jsonlines.org/)')
        self._parser.add_argument('--title', metavar='<title>', type=str,
                                  help='drawing title')
        self._parser.add_argument('-W', '--width', type=float, metavar='<w>',
                                  dest='W', default=1000.0,
                                  help='drawing width')
        self._parser.add_argument('-H', '--height', type=float, metavar='<h>',
                                  dest='H', default=800.0,
                                  help='drawing height')
        self._args = self._parser.parse_args()


# Global entry point
def runVizCvtGraphApp():
    app = VizCvtGraphApp()
    app.run()


# Check if main()
if __name__ == "__main__":
    runVizCvtGraphApp()

# This Python file uses the following encoding: utf-8
import json

from dataCvtGraph import CvtGraph


# Graph file reader: JSON Lines format.
# Try to be incremental, assume more lines arrive later
class GraphReaderJSONL:
    # Read a file
    def read(self, filename: str, graph: CvtGraph):
        self._graph = graph
        with open(filename, 'r') as f:
            for line in f:
                # removing the new line characters
                self._processLine(line.rstrip())

    # Process next line
    def _processLine(self, line: str):
        values = json.loads(line)
        self._graph.addDataChunk(values)

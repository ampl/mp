# This Python file uses the following encoding: utf-8


# Record describing a range of nodes of certain type
class NodeRange:
    # Constructor
    def __init__(self):
        self.nodeType = None
        self.start = 0          # start index
        self.end = 0            # last+1 index

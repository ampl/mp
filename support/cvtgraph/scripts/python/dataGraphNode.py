# This Python file uses the following encoding: utf-8


# Graph node, internal representation
class GraphNode:
    # Constructor
    def __init__(self):
        self.nItems = 0
        self.rankMin = 1000000000
        self.rankMax = 0
        self.hasIn = False
        self.hasOut = False

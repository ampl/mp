# This Python file uses the following encoding: utf-8
from dataGraphNodeRange import NodeRange
from dataGraphNode import GraphNode
from dataGraphArc import GraphArc


# Model Conversion Graph representation
class CvtGraph:
    # Constructor
    def __init__(self):
        self._nChunksAdded = 0     # number of data chunks added
        self._nodes = dict()       # nodes
        self._arcs = dict()        # current link's dest nodes

    # Add data chunk described as a dictionary
    def addDataChunk(self, chunk: dict):
        self._nChunksAdded += 1
        if "link_type" in chunk:
            self._addLink(chunk)
        elif "node_type" in chunk:
            self._addNode(chunk)
        else:
            raise RuntimeError(
                  "Unknown graph data in chunk {}: {}".
                  format(self._nChunksAdded, chunk))

    # Draw the graph.
    # @param drawer: the drawing object
    def drawGraph(self, drawer, W, H):
        self._computeLayout(W, H)
        self._exportLayout(drawer)

    # Add link
    def _addLink(self, data: dict):
        lType = data["link_type"]
        srcNodes = self._compileLinkTerminalNodes(data["source_nodes"])
        destNodes = self._compileLinkTerminalNodes(data["dest_nodes"])
        self._addArcs_deduceNodes(lType, srcNodes, destNodes)

    # Add node
    def _addNode(self, data: dict):
        raise("Node input not implemented")

    # Compile sources or targets of a link.
    # Input: list of node ranges
    # @return compiled list
    def _compileLinkTerminalNodes(self, data: list):
        nodes = []
        for node in data:
            nodes.append(self._compileLinkTerminalNodeRange(node))
        return nodes

    # Compile single link terminal node range
    # @return the compiled struct
    def _compileLinkTerminalNodeRange(self, data: dict):
        range = NodeRange()
        assert 1 == len(data)            # 1 range per list item
        for key, rng in data.items():    # each item: { type, range }
            range.nodeType = key
            if isinstance(rng, int):
                range.start, range.end = rng, rng+1
            else:
                range.start, range.end = rng[0], rng[1]+1
        return range

    # Add graph arcs.
    # Currenlty this also saves nodes (node types),
    # as we don't expect node records from the file.
    # @param src_nodes: source node ranges
    # @param dest_nodes: destination node ranges
    def _addArcs_deduceNodes(self,
                             lType: str, srcNodes: list, destNodes: list):
        for src in srcNodes:
            srcRepr = self._nodes.setdefault(src.nodeType, GraphNode())
            srcRepr.rankMin = min(srcRepr.rankMax, srcRepr.rankMin)  # very src
            srcRepr.nItems = max(srcRepr.nItems, src.end)
            srcRepr.hasOut = True
            for dest in destNodes:
                destRepr = self._nodes.setdefault(dest.nodeType, GraphNode())
                # Update dest node's infos
                destRepr.nItems = max(destRepr.nItems, dest.end)
                destRepr.rankMin = min(destRepr.rankMin, srcRepr.rankMin+1)
                destRepr.rankMax = max(destRepr.rankMax, srcRepr.rankMax+1)
                destRepr.hasIn = True
                # Create/add meta-arc
                arcRepr = self._arcs.setdefault((src.nodeType, dest.nodeType),
                                                GraphArc())
                arcRepr.nItems += src.end-src.start

    # Compute layout
    def _computeLayout(self, W, H):
        # Node groups
        self._collectNodeTypes()
        self._sortNodeList(self._nodesSrc)
        self._sortNodeList(self._nodesDest)
        # Coords
        wGap = W/10    # margins
        hGap = H/10
        self._assignCoordsTopDown(self._nodesSrc, wGap, H-hGap,
                                  0, 0.8*H/(len(self._nodesSrc)-1))
        self._assignCoordsTopDown(self._nodesMid, 2.5*wGap, H-2*hGap,
                                  0.5*W/(self._rankMidMax-1),
                                  0.6*H/(len(self._nodesMid)-1), True)
        self._assignCoordsTopDown(self._nodesDest, W-wGap, H-hGap,
                                  0, 0.8*H/(len(self._nodesDest)-1))

    # Collect node types (source, middle, dest).
    # Relying on Python to keep dict() in the original order.
    def _collectNodeTypes(self):
        self._nodesSrc = []
        self._nodesMid = []
        self._nodesDest = []
        self._rankMidMax = 0
        for ndName, nd in self._nodes.items():
            if nd.hasIn:
                if nd.hasOut:
                    self._nodesMid.append(ndName)
                    self._rankMidMax = max(self._rankMidMax, nd.rankMax)
                else:
                    self._nodesDest.append(ndName)
            else:
                assert nd.hasOut
                self._nodesSrc.append(ndName)

    # Sort source or target node list
    def _sortNodeList(self, names: list):
        names.sort(key=lambda x: 1 if x.find('_vars') >= 0 else
                   2 if x.find('_cons') >= 0 else
                   3 if x.find('_objs') >= 0 else 4)

    # Assign coordinates top-down
    # dX is per rank increase, if fRank
    def _assignCoordsTopDown(self, ndNames: list,
                             xLeft, yTop, dX, dY, fRank=False):
        y = yTop
        for name in ndNames:
            node = self._nodes[name]
            node.x = xLeft
            if fRank:
                assert 1 <= node.rankMin     # only middle nodes with fRank
                node.x += (node.rankMax-1) * dX
            node.y = y
            y -= dY

    # Export to the drawer object
    def _exportLayout(self, drawer):
        # Nodes
        for name, node in self._nodes.items():
            drawer.addNode(node.x, node.y, node.nItems,
                           name, "# of items: "+str(node.nItems))
        # Edges
        for name, edge in self._arcs.items():
            nd0 = self._nodes[name[0]]
            nd1 = self._nodes[name[1]]
            drawer.addArc(nd0.x, nd0.y, nd1.x, nd1.y)
        # Draw
        drawer.draw()


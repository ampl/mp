# This Python file uses the following encoding: utf-8

# if __name__ == "__main__":
#     pass

class DiGraph:
  """
  A simple digraph or a wrapper around some graph library
  """

  def __init__(self):
    self._nodes = []
    self._arcs = []

  def AddNode(self, data=None):
    self._nodes.append(data)
    return len(self._nodes)-1

  def GetNode(self, idx):
    return self._nodes[idx]

  def ToText(self):
    return str(self._nodes)

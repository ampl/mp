# This Python file uses the following encoding: utf-8

# if __name__ == "__main__":
#     pass

class ModelView:
  """
  A view of a (sub) model
  """

  def __init__(self):
    self._data = None

    self._vars = {"Variables": []}
    self._cons = {"Constraints": []}
    self._objs = {"Objectives": []}

  def SetData(self, data):
    self._data = data

  def GetData(self):
    return self._data

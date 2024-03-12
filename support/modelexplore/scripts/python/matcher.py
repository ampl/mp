# This Python file uses the following encoding: utf-8

# if __name__ == "__main__":
#     pass

from scripts.python.model import Model
from scripts.python.modelview import ModelView

class Matcher:
  """
  Selects a submodel
  """

  def __init__(self):
    self.data = None


def MatchSubmodel(m: Model, patt: str, fwd: bool, bwd: bool):
  """
  Match a submodel containg the \a pattern,
  optionally extended by forward/backward
  reformulation graph search
  """
  mv1 = ModelView()
  mv2 = ModelView()
  mv1.SetData(m.MatchOrigModel(patt))
  mv2.SetData(m.MatchFinalModel(patt))
  return mv1, mv2

# This Python file uses the following encoding: utf-8

# if __name__ == "__main__":
#     pass

import streamlit as st

from scripts.python.graph import DiGraph

class Model:
  """
  An optimization model with conversion
  graph
  """

  def __init__(self):
    self._graph = DiGraph()    ## Underlyng graph

    self._vars = []            ## Pointers to various parts of the graph
    self._dvars = []
    self._cons_NL_all = []
    self._cons_NL = {        ## NL + SOS
      "All" : [],
      "Nonlinear" : [],
      "Linear" : [],
      "Logical": [],
      "SOS1": [],
      "SOS2": []}
    self._cons_Flat = {}
    self._cons_Flat_Group = {}
    self._objs_NL = []
    self._objs = []

  def UpdateVar(self, idx, data):
    self._updateNodeData(self._vars, idx, data)

  def UpdateDefVar(self, idx, data):
    self._updateNodeData(self._dvars, idx, data)

  def UpdateNLObj(self, idx, data):
    self._updateNodeData(self._objs_NL, idx, data)

  def UpdateFlatObj(self, idx, data):
    self._updateNodeData(self._objs, idx, data)

  def UpdateNLCon(self, type, idx, data):
    if "nonlin"==type:
      self._updateNodeData(self._cons_NL["Nonlinear"],
      len(self._cons_NL["Nonlinear"]),   ## these just 1x
      data)
    elif "lin"==type:
      self._updateNodeData(self._cons_NL["Linear"],
      len(self._cons_NL["Linear"]), data)
    elif "logical"==type:
      self._updateNodeData(self._cons_NL["Logical"],
      len(self._cons_NL["Logical"]), data)
    elif "_sos1"==type:
      self._updateNodeData(self._cons_NL["SOS1"],
      len(self._cons_NL["SOS1"]), data)
    elif "_sos2"==type:
      self._updateNodeData(self._cons_NL["SOS2"],
      len(self._cons_NL["SOS2"]), data)
    else:
      raise Exception("Unknown NL constraint type: "+type)

  def UpdateFlatConGroup(self, type, data):
    self._cons_Flat_Group[type] = data

  def UpdateFlatCon(self, type, idx, data):
    if type not in self._cons_Flat:
      self._cons_Flat[type] = []
    self._updateNodeData(self._cons_Flat[type], idx, data)
    if 0==data["depth"] \
        and type.startswith('_sos') \
        and "printed" in data:   ## we need the final status
      self.UpdateNLCon(type, 0, data)

  def _updateNodeData(self, specnodecnt, idx, data):
    data1, upd = self._updateItemData(specnodecnt, idx, data)
    if (not upd):
      idx = self._graph.AddNode(data1)
      data1["node_index"] = idx

  def _updateItemData(self, specnodecnt, idx, data):
    if len(specnodecnt)<=idx:
      specnodecnt.insert(idx, {})
    if (specnodecnt[idx] is None):    ## No such item
      specnodecnt[idx] = {}
    ifEmpty = 0==len(specnodecnt[idx])
    self._updateMap(specnodecnt[idx], data)
    return specnodecnt[idx], ifEmpty

  def _updateMap(self, data1, data2):
    data1.update(data2)

  # Match keyword to the original model
  def MatchOrigModel(self, keyw):
    result = {}
    result["NL Variables"] = self._matchRecords(self._vars, keyw, "is_from_nl")
    result["NL Defined Variables"] = self._matchRecords(self._dvars, keyw)
    result["NL Objectives"] = self._matchRecords(self._objs_NL, keyw)
#    result["NL Constraints"] \
#      = self._matchRecords(self._cons_NL.get("All"), keyw)
    result["NL Nonlinear Constraints"] \
      = self._matchRecords(self._cons_NL.get("Nonlinear"), keyw)
    result["NL Linear Constraints"] \
      = self._matchRecords(self._cons_NL.get("Linear"), keyw)
    result["NL Logical Constraints"] \
      = self._matchRecords(self._cons_NL.get("Logical"), keyw)
    result["NL SOS1 Constraints"] \
      = self._matchRecords(self._cons_NL.get("SOS1"), keyw)
    result["NL SOS2 Constraints"] \
      = self._matchRecords(self._cons_NL.get("SOS2"), keyw)
    return result

  # Match keyword to the final model
  def MatchFinalModel(self, keyw):
    result = {}
    result["Variables"] = self._matchRecords(self._vars, keyw)
    result["Objectives"] = self._matchRecords(self._objs, keyw)
    for ct, cv in self._cons_Flat.items():
      result["Constraints '" + ct + "'"] \
        = self._matchRecords(self._cons_Flat[ct], keyw)
    return result

  # Add records containing keyword
  # @return array of strings
  def _matchRecords(self, cnt, keyw, keyNeed1=None):
    result = ""
    if cnt is None:
      return result
    for i in cnt:
      pr = str(i)         ## TODO printed form
      if "printed" in i:
        pr = i["printed"]
      assert len(pr)
      if ';'!=pr[-1]:
        pr = pr + ';'
      if (""==keyw or keyw in pr) \
      and (keyNeed1==None \
        or (keyNeed1 in i and 1==i[keyNeed1])):
        result = result + "  \n" + pr        ## Markdown: 2x spaces + EOL
    return result

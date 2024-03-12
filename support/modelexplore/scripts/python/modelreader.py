# This Python file uses the following encoding: utf-8

# if __name__ == "__main__":
#     pass

import streamlit as st

import json
from scripts.python.model import Model

class ModelReader:
  """
  Model reader
  """

  def __init__(self):
    self._model = Model()

  def ReadModel(self, uploader):
    for line in uploader:
      # removing the new line characters
      self._processLine(line.rstrip())
    return self._model

  # Process next line
  def _processLine(self, line: str):
    values = json.loads(line)
    self._addDataChunk(values)

  # Add data chunk as a JSON-like object
  def _addDataChunk(self, chunk):
    if "VAR_index" in chunk:
      self._model.UpdateVar(chunk["VAR_index"], chunk)
    elif "NL_OBJECTIVE_index" in chunk:
      self._model.UpdateNLObj(chunk["NL_OBJECTIVE_index"], chunk)
    elif "NL_CON_TYPE" in chunk:
      self._model.UpdateNLCon(chunk["NL_CON_TYPE"], chunk["index"], chunk)
    elif "OBJECTIVE_index" in chunk:
      self._model.UpdateFlatObj(chunk["OBJECTIVE_index"], chunk)
    elif "CON_GROUP" in chunk:
        self._model.UpdateFlatConGroup(chunk["CON_TYPE"], chunk)
    elif "CON_TYPE" in chunk:
      self._model.UpdateFlatCon(chunk["CON_TYPE"], chunk["index"], chunk)


def ReadExplorerModel(uploader):
  mr = ModelReader()
  return mr.ReadModel(uploader)

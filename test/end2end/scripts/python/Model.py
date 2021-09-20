import enum
from pathlib import Path
from collections.abc import Iterable


class ModelTags(enum.Enum):
    NA = 0,

    linear = 1
    quadratic = 2
    socp = 3
    nonlinear = 4
    quadraticnonsdp = 5,
    complementarity = 6,
    arc = 7,
    plinear = 8

    continuous = 10,
    integer = 11,
    binary = 12,

    logical = 13,

    trigonometric = 100,
    htrigonometric = 101,
    log = 102,

    script = 1000,


    @staticmethod
    def fromString(l):
        ret = list()
        for s in l:
            ret.append(ModelTags[s])
        return ret


class Model(object):

    """Represents a model"""

    def __init__(self, filename, expsolution, tags, otherFiles=None, overrideName=None,description=None):
        if isinstance(filename, str):
            filename = Path(filename)
        self._filename = filename
        self._expsolution = expsolution
        self._tags = tags
        if otherFiles:
            self._otherFiles = otherFiles if type(
                otherFiles) is list else [otherFiles]
        else:
            self._otherFiles = None
        self._name = overrideName if overrideName else filename.stem
        self._description = description

    def getSolvers(self):
        return self._description["solvers"]

    def hasOptions(self):
        return (self._description is not None) and \
            ("options" in self._description)

    def getOptions(self):
        return self._description["options"]

    def hasExpectedValues(self):
        return (self._description is not None) and \
            ("values" in self._description)

    def getExpectedValues(self):
        return self._description["values"]

    def getExpectedObjective(self):
        return self._expsolution

    def getName(self):
        return self._name

    def isNL(self):
        return self._filename.suffix.lower() == ".nl"

    def getFilePath(self):
        return str(self._filename)

    def getSolFilePath(self):
        return self._filename.with_suffix(".sol")

    def getAdditionalFiles(self):
        return self._otherFiles

    def hasTag(self, tag : ModelTags):
      if isinstance(self._tags , Iterable):
        return tag in self._tags
      else:
        return tag == self._tags

    def hasAnyTag(self, tags : list):
      if tags is None:
        return False
      for tag in tags:
        if self.hasTag(tag):
          return True
      return False

    def isScript(self):
        return self.hasTag(ModelTags.script)

    @staticmethod
    def LINEAR(filename, expsolution, vars: ModelTags = ModelTags.continuous):
        return Model(filename, expsolution, [ModelTags.linear, vars])

    @staticmethod
    def QUADRATIC(filename, expsolution, vars: ModelTags = ModelTags.continuous):
        return Model(filename, expsolution, [ModelTags.quadratic, vars])

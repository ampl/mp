from Model import Model, ModelTags
import sys
import json
from pathlib import Path, PurePath
from itertools import repeat

class ModelsDiscovery(object):
    """Find all models in the specified directory"""

    def _ReadModelsDescription(self, directory, filename="modellist.json"):
        p = PurePath(directory)
        f = str(p.joinpath(filename))
        try:
            with open(f) as json_file:
                try:
                    return json.load(json_file)
                except:
                    print("Model list file '{}': bad format.".format(f),
                          sys.exc_info())
        except:
            return None

    def _CreateModel(self, f, desc):
        n = f.stem
        if desc:
            md = next((item for item in desc if item["name"] == n), None)
            if md:
              solution = md["solution"] if "solution" in md else None
              if "files" in md:
                files = md["files"]
                f = f.parent.joinpath(files[0])
                files = [f.parent.joinpath(ff) for ff in files[1:]]
                return Model(f, solution, ModelTags.fromString(md["tags"]),
                  otherFiles=files, description=md)
              else:
                return Model(f, solution, ModelTags.fromString(md["tags"]),
                  description=md)
        return Model(f, None, ModelTags.NA)

    def GetStem(self):
        return self._stem

    def _findModelsFunction(self, directory : Path, preferAMPLModels = False, justNL= False):
        models = []
        NLmodels = self.FindNLModels(directory)
        if justNL:
          return NLmodels

        if str(directory.absolute()).endswith("-mps"):
          AMPLmodels = self.FindAMPLMPSModels(directory)
        else:
          AMPLmodels = self.FindAMPLModels(directory)

        if not AMPLmodels:
          return NLmodels
        if not NLmodels:
          return AMPLmodels

        # If we get here, we have some models available in both AMPL and NL format
        onlyAMPL = []
        onlyNL = []
        both = []
        for n in NLmodels:
          try:
            a = next(filter(lambda m : m.getName() == n.getName(), AMPLmodels))
            both.append((a,n))
          except StopIteration:
            onlyNL.append(n)
        for a in AMPLmodels:
          try:
            n = next(filter(lambda m : m.getName() == n.getName(), NLmodels))
          except StopIteration:
            onlyAMPL.append(a)

        # If we have only one version, append it
        models.extend(onlyNL)
        models.extend(onlyAMPL)
        # If we have both NL and AMPL, use the AMPL version if "regenerate" is specified,
        # otherwise use the NL version
        if preferAMPLModels:
          models.extend(a for (a, n) in both)
        else:
          models.extend(n for (a, n) in both)

        return models

    def FindModelsGeneral(self, directory, recursive = True, preferAMPLModels = False, justNL= False):
        """Find models in the specified directory, recursively unless otherwise specified. 
           If both AMPL and NL versions are found, the parameter preferAMPLModels chooses what to choose.
           If justNL is specified, it ignores all AMPL models (useful for system without AMPL)
           If the directory has the suffix "-mps" it tries to find one single model file and expects one dat file each,
           useful for models exported from GAMS/MPS as many found in NETLIB"""
        self._desc = None
        p = Path(directory)
        models = self._findModelsFunction(p, preferAMPLModels, justNL)
        if not models:
            print("No models or case descriptions found.")
        if not recursive:
          return models
        
        # Recurse baby!
        
        dirs = list(filter(lambda x: x.is_dir(),  p.iterdir()))
        if recursive and len(dirs) > 0:
          for d in dirs:
            models.extend(self.FindModelsGeneral(d, True, preferAMPLModels, justNL))
        return models


    def FindAMPLMPSModels(self, directory):
        """Finds NETLIB-like models: one base model file and a different datafile for each instance"""
        self._desc = None
        p = Path(directory)
        self._stem = p.stem
        _desc  = self._ReadModelsDescription(directory)
        base = [f for f in p.iterdir() if not f.is_dir()
                and f.suffix == ".mod"][0]
        if not base:
            raise Exception("Base model not found")
        files = [f for f in p.iterdir() if not f.is_dir()
                 and f.suffix == ".dat"]
        models = []
        for f in files:
            n = f.stem
            if _desc:
                md = next(
                    (item for item in self._desc if item["name"] == n), None)
                if md:
                    m = Model(base, md["solution"], ModelTags.fromString(
                        md["tags"]), otherFiles=f, overrideName=n, description=md)
                else:
                    m = Model(base, 0, ModelTags.NA,
                              otherFiles=f, overrideName=n)
            else:
                m = Model(base, 0, ModelTags.NA, otherFiles=f, overrideName=n)
            models.append(m)
        return models

    def FindAMPLModels(self, directory):
        """Find all AMPL models. Creates one entry for each .mod file and for each entry in the modellist.json file.
        Tries to match the first with the latter to get the model properties"""
        return self.FindModels(directory, ".mod")

    def FindNLModels(self, directory):
        """Find all NL models"""
        return self.FindModels(directory, ".nl")

    def FindModels(self, directory, extension):
        self._desc = None
        p = Path(directory)
        self._stem = p.stem
        desc = self._ReadModelsDescription(directory)
        if desc:
            print("  Path '{}': loaded model description with {} items... ".format(directory, len(desc)))
            files = [f for f in p.iterdir() if not f.is_dir()
                     and f.suffix.lower() == extension]
            return list(map(self._CreateModel, files, repeat(desc)))
        return list()

from Model import Model, ModelTags
import sys
import json
from pathlib import Path, PurePath
from itertools import repeat

class ModelsDiscovery(object):
    """Find all models in the specified directory"""

    
    def FilterByTag(models : list, excludeTags : list):
        return filter(lambda x : x.hasAnyTag(excludeTags) != True, models)

    def FindModelsGeneral(self, directory, recursive = True,
                          modellist=True, preferAMPLModels = False, justNL= False):
        """Find models in the specified directory, recursively unless otherwise specified.
           If modellist=True, strictly follow modellist.json.
           Otherwise collect .mod, .nl files, and:
           If both AMPL and NL versions are found, the parameter preferAMPLModels chooses what to choose.
           If justNL is specified, it ignores all AMPL models (useful for system without AMPL)
           If the directory has the suffix "-mps" it tries to find one single model file and expects one dat file each,
           useful for models exported from GAMS/MPS as many found in NETLIB"""
        self._desc = None
        p = Path(directory)
        if modellist==True:
            models = self._parseModelList(p, preferAMPLModels)
        else:
            models = self._findModelsFunction(p, preferAMPLModels, justNL)
        if not recursive:
          return models
        
        # Recurse baby!
        
        dirs = list(filter(lambda x: x.is_dir(),  p.iterdir()))
        if recursive and len(dirs) > 0:
            for d in dirs:
                newmodels = self.FindModelsGeneral(d, True, modellist, preferAMPLModels, justNL)
                models.extend(newmodels)
        return models

    def _parseModelList(self, dir: Path, preferAMPLModels):
        desc = self._ReadModelsDescription(dir)
        if desc:
            print("  Path '{}': loaded model description with {} items... ".format(dir, len(desc)))
            result = list()
            for md in desc:
                f, files = self._determineModelFile(dir, md, preferAMPLModels)
                if not f:
                    print("None of the files found: ", files)
                    exit(-1)
                result.append(self._CreateModelFromListEntryWithKnownFile(f, md))
            return result
        return list()

    def _determineModelFile(self, dir: Path, md, preferAMPLModels):
        if "files" in md:
            files = [dir.joinpath( md["files"][0] )]
        else:
            name = md["name"].split()[0]
            files = [dir.joinpath(name + suf) for suf in [".mod", ".nl"]]
            if not preferAMPLModels:
                files.reverse()
        for f in files:
            if f.is_file():
                return f, files
        return None, files

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

    def FindAMPLMPSModels(self, directory):
        """Finds NETLIB-like models: one base model file and a different datafile for each instance"""
        self._desc = None
        p = Path(directory)
        self._stem = p.stem
        self._desc  = self._ReadModelsDescription(directory)
        base = [f for f in p.iterdir() if not f.is_dir()
                and f.suffix == ".mod"][0]
        if not base:
            raise Exception("Base model not found")
        files = [f for f in p.iterdir() if not f.is_dir()
                 and f.suffix == ".dat"]
        models = []
        for f in files:
            n = f.stem
            if self._desc:
                md = next(
                    (item for item in self._desc if item["name"] == n), None)
                if md:
                    m = Model(base, md["objective"], ModelTags.fromString(
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
            return list(map(self._CreateModelFindingDescription, files, repeat(desc)))
        return list()

    def _ReadModelsDescription(self, directory, filename="modellist.json"):
        p = PurePath(directory)
        f = str(p.joinpath(filename))
        if not Path(f).is_file():
            return None
        with open(f) as json_file:
            try:
                return json.load(json_file)
            except:
                print("Model list file '{}': bad format.".format(f),
                      sys.exc_info())
                raise

    def _CreateModelFindingDescription(self, f, desc):
        n = f.stem
        if desc:
            md = next((item for item in desc if item["name"].split()[0] == n), None)
            if md:
                return self._CreateModelFromListEntryWithKnownFile(f, md)
        return Model(f, None, ModelTags.NA)

    def _CreateModelFromListEntryWithKnownFile(self, f, md):
        solution = md["objective"] if "objective" in md else None
        if "files" in md:
            files = md["files"]
            f = f.parent.joinpath(files[0])
            files = [f.parent.joinpath(ff) for ff in files[1:]]
            return Model(f, solution, ModelTags.fromString(md["tags"]),
                         otherFiles=files, description=md, overrideName=md["name"])
        else:
            return Model(f, solution, ModelTags.fromString(md["tags"]),
                         description=md, overrideName=md["name"])

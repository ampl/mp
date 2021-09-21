from ModelRunner import ModelRunner, ModelComparer
import math


class Exporter(object):
    """Base class to export ModelRunner results"""

    def export(self, mr: ModelRunner):
        raise Exception("Not defined for base class")

    # Return False if failed
    def printStatus(self, model, run):
        if "eval_done" in run:                # Evaluation happened
            if "eval_fail_msg" not in run:    # Evaluation ok
                print("  Ok.", end='')
            else:
                print("          FAILED({}): {}\n{: <80}".format(
                    run["solver"], run["eval_fail_msg"], ' '),
                    end='', flush=True)
                return False
        else:
            print("    <no cmp data>", end='')
        return True


class StringTestExporter(Exporter):

    def printString(self, m, r):
        print(m.getName(), m.getExpectedSolution(),
              r["solution"], r["solutionTime"])

    def printStringCompare(self, m, r, r2):
        print(m.getName(), m.getExpectedSolution(),
              r["solution"], r["solutionTime"], r2["solution"],
              r2["solutionTime"])

    def export(self, mr: ModelRunner):
        if isinstance(mr, ModelComparer):
            for (m, r, r2) in zip(mr._models, mr._runs, mr._runs2):
                self.printStringCompare(m, r, r2)
        else:
            for (m, r) in zip(mr._models, mr._runs):
                self.printString(m, r)

    def exportOne(self, mr: ModelRunner):
        i = len(mr._runs)
        m = mr._models[i-1]
        r = mr._runs[i-1]
        if isinstance(mr, ModelComparer):
            r2 = mr._runs2[i-1]
            self.printStringCompare(m, r, r2)
        else:
            self.printString(m, r)


class CSVTestExporter(Exporter):
    def __init__(self, fileName):
        self._fileName = fileName

    def sanifyString(self, s):
        try:
            s = s.replace("\n", " - ")
            s = s.replace(",", ";")
            return s
        except:
          return s

    def getHeaderSingleRun(self):
        return "Name, Expected Solution, Solution, Time, Time Limit, SolverMsg\n"

    def getStringSingleRun(self, m, r):
        return "{}, {}, {}, {}, {}, {}\n".format(m.getName(), m.getExpectedObjective(), r["objective"], r["solutionTime"], r["timelimit"], self.sanifyString(r["outmsg"]))

    def getHeaderCompareRun(self, mc: ModelComparer):
        (r1, r2) = mc.getRunnerNames()
        return "Name, Expected Solution, {}-Solution, {}-Time, {}-TimeLimit, {}-Solution, {}-Time, {}-TimeLimit, {}-SolverMsg, {}-SolverMsg\n".format(r1, r1, r1, r2, r2, r2, r1, r2)

    def getStringCompareRun(self, m, r1, r2):
        return "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n".format(m.getName(), m.getExpectedObjective(),
                                                                 r1["objective"], r1["solutionTime"], r1["timelimit"],
                                                                 r2["objective"], r2["solutionTime"], r2["timelimit"],
                                                                 self.sanifyString(r1["outmsg"]), self.sanifyString(r2["outmsg"]))

    def export(self, mr: ModelRunner):
        if isinstance(mr, ModelComparer):
            with open(self._fileName, "w") as file:
                file.write(self.getHeaderCompareRun(mr))
                for (m, r, r2) in zip(mr._models, mr._runs, mr._runs2):
                    file.write(self.getStringCompareRun(m, r, r2))
        else:
            with open(self._fileName, "w") as file:
                file.write(self.getHeaderSingleRun())
                for (m, r) in zip(mr._models, mr._runs):
                    file.write(self.getStringSingleRun(m, r[-1]))

    def exportOne(self, mr: ModelRunner):
        i = len(mr._runs)
        m = mr._models[i-1]
        r = mr._runs[i-1]

        filemode = "w" if i == 1 else "a+"
        if isinstance(mr, ModelComparer):
            r2 = mr._runs2[i-1]
            with open(self._fileName, filemode) as file:
                if i == 1:
                    file.write(self.getHeaderCompareRun(mr))
                file.write(self.getStringCompareRun(m, r, r2))
        else:
            with open(self._fileName, filemode) as file:
                if i == 1:
                    file.write(self.getHeaderSingleRun())
                file.write(self.getStringSingleRun(m, r))
            self.printStatus(m, r)

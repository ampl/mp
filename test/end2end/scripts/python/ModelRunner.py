from Solver import Solver
from AMPLRunner import AMPLRunner
from TimeMe import TimeMe
from Model import Model

class ModelRunner(object):
    """Class to run a set of models and capture their outputs"""

    def __init__(self, ampl, runners, optionsExtra=None):
        self._ampl = ampl
        self._runners = runners
        self._amplRunners = None
        self._runs = [ list() for r in self._runners ]
        self._optionsExtra = optionsExtra

    def getRuns(self):
        return self._runs

    def getModels(self) -> list:
        return self._models
    def getRunners(self)->list:
        return self._runners
    def getLogFileName(m : Model, s : Solver):
        return f"{m.getName()}.{s.getName()}.log"

    @staticmethod
    def getAMPLVersion(amplPath: str) -> str:
        """First line of 'ampl -v', or a note when it cannot be run"""
        import subprocess
        try:
            out = subprocess.run([amplPath, "-v"], capture_output=True, text=True,
                                 stdin=subprocess.DEVNULL, timeout=30)
            first = (out.stdout or out.stderr).strip().splitlines()
            return first[0] if first else "unavailable (no output)"
        except Exception as exc:
            return "unavailable ({})".format(exc)

    STARS = "*" * 60

    def _folderSummaryLines(self, folder, nCases, counters, start):
        """Summary of one folder: counters are cumulative, start is their
           snapshot when the folder began"""
        nFailedSolver, nAbortedSolver, nFailedScriptOrAMPL, nSkipped = counters
        lines = []
        for i, r in enumerate(self.getRunners()):
            slv = f" {r.getName()}:" if len(self.getRunners()) > 1 else ""
            lines.append("Folder '%s':%s %d cases, %d fail (%d abrt) slv, "
                         "%d fail AMPL/PY/scrpt, %d skip" % (
                folder, slv, nCases,
                nFailedSolver[i] - start[0][i], nAbortedSolver[i] - start[1][i],
                nFailedScriptOrAMPL[i] - start[2][i], nSkipped[i] - start[3][i]))
        return lines

    def _printFolderSummary(self, lines):
        for l in lines:
            print(l)
        print(ModelRunner.STARS, flush=True)

    def _printFolderBanner(self, folder):
        print("")
        print(ModelRunner.STARS)
        print(folder)
        print(ModelRunner.STARS, flush=True)

    def run(self, modelList: list, exporter=None, keepLogs = False, verbose=False):
        """Run the models in this instance. If exporter != None, it exports the results as it goes"""
        self._models = modelList
        n = 0
        nFailedSolver = [0 for r in self.getRunners()]
        nAbortedSolver = [0 for r in self.getRunners()]
        nFailedScriptOrAMPL = [0 for r in self.getRunners()]
        nSkipped = [0 for r in self.getRunners()]
        counters = (nFailedSolver, nAbortedSolver, nFailedScriptOrAMPL, nSkipped)
        folderSummaries = []     # (folder, summary lines) in run order
        curFolder = None
        folderStart = None
        nFolderCases = 0
        EFM ="eval_fail_msg"
        EM ="errormsg"
        # Set to true for junit to add the solver output to the test result
        keep_output=False
        try:
        # Attempt to access the class
            __import__("JunitExporter")
            import JunitExporter
            if isinstance(exporter, JunitExporter.JunitExporter):
                keep_output=True
        except:
            pass
        
        for m in modelList:
            folder = m.getCaseDir() or ""
            if folder != curFolder:
                if curFolder is not None:
                    lines = self._folderSummaryLines(curFolder, nFolderCases, counters, folderStart)
                    self._printFolderSummary(lines)
                    folderSummaries.append(lines)
                curFolder = folder
                folderStart = tuple(list(c) for c in counters)
                nFolderCases = 0
                self._printFolderBanner(folder)
            n += 1
            nFolderCases += 1
            if m.isNL():
                cr = self._runners
                msg = "{}. As NL: '{}'".format(n, m.getName())
            else:
                if not self._amplRunners:
                    amplPath = Solver.resolveExecutable(self._ampl)
                    print("AMPL executable: '{}'".format(amplPath))
                    print("AMPL version: {}".format(ModelRunner.getAMPLVersion(amplPath)), flush=True)
                    self._amplRunners = [
                        AMPLRunner(self._ampl, r, self._optionsExtra,
                                   printOutput=verbose, storeOutput=keep_output,
                                   timeout=r.getTimeout())
                        for r in self._runners ]
                cr = self._amplRunners
                msg = "{}. Via AMPL: '{}'".format(n, m.getName())
            print("{0: <60}".format(msg), end="", flush=True)
            for (i,r) in enumerate(cr):
                t = TimeMe()
                with t:
                  try:
                    r.isBenchmark = len(cr) > 1
                    self._runs[i].append({})
                    if isinstance(r, AMPLRunner):
                        ss = r.getSolver()
                    else:
                        ss = r
                    if not m.isSubsetOfTags(ss.getSupportedTags()) or \
                            m.hasAnyTag(ss.getUnsupportedTags()) or \
                            (m.hasSolvers() and ss.getName() not in m.getSolvers()):
                        self._runs[i][-1]["outmsg"] = "Skipped, unsupported tags or solver restriction"
                        self._runs[i][-1]["solver"] = ss
                        if r.isBenchmark:
                            print("\n\t\t{0: <20}: Skipped due to unsupported tags or solver restriction".
                              format(ss.getName()), flush=True)
                        else:
                            print("Skipped due to unsupported tags or solver restriction", flush=True)
                        nSkipped[i] += 1
                    else:
                      if keepLogs:
                        r.runAndEvaluate(m, logFile=ModelRunner.getLogFileName(m, ss))
                      else:
                        r.runAndEvaluate(m, logFile=None)
                    
                      stats = r.getSolutionStats()
                      
                      if EFM in stats:
                        if EM in stats:
                            stats.setdefault(EFM, []).append( stats[EM] )
                      if keep_output:
                         stats["output"]=r.get_output()
                      self._runs[i][-1] = stats
                      if exporter:
                        if not exporter.printStatus(m, stats):
                            nFailedSolver[i] += 1
                      if 0 != stats["exit_code"]:
                          nAbortedSolver[i] += 1
                  except Exception as exc:
                    self._runs[i][-1]["outmsg"] = "AMPL(PY)/script failure"
                    self._runs[i][-1]["solver"] = ss
                    self._runs[i][-1].setdefault(EFM, []).append(str(exc))
                    print("   EXCEPTION: ", exc)
                    nFailedScriptOrAMPL[i] += 1
                print("  (%.2fs, %d fail (%d abrt) slv, %d fail AMPL/PY/scrpt, %d skip)" %
                  (t.interval, nFailedSolver[i], nAbortedSolver[i], nFailedScriptOrAMPL[i], nSkipped[i]),
                  end="", flush=True)
            if exporter:
                self.export(exporter)
            print("")
        if curFolder is not None:
            lines = self._folderSummaryLines(curFolder, nFolderCases, counters, folderStart)
            self._printFolderSummary(lines)
            folderSummaries.append(lines)
        if len(folderSummaries) > 1:
            print("")
            print(ModelRunner.STARS)
            print("Per-folder summary")
            print(ModelRunner.STARS)
            for lines in folderSummaries:
                for l in lines:
                    print(l)
        if curFolder is not None:
            print("")
            print(ModelRunner.STARS)
            for i, r in enumerate(self.getRunners()):
                slv = f" {r.getName()}:" if len(self.getRunners()) > 1 else ""
                print("Overall:%s %d cases, %d fail (%d abrt) slv, "
                      "%d fail AMPL/PY/scrpt, %d skip" % (
                    slv, n, nFailedSolver[i], nAbortedSolver[i],
                    nFailedScriptOrAMPL[i], nSkipped[i]))
            print(ModelRunner.STARS, flush=True)
        exporter.export()

    def export(self, exporter):
        exporter.exportInstanceResults(self)


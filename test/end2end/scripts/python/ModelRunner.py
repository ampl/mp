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
    SUMMARY_HEADER = ("Folder", "Solver", "Cases", "Fail slv", "Abrt", "Fail AMPL/PY/scrpt", "Skip")

    def _folderSummaryRows(self, folder, nCases, counters, start):
        """Summary rows (one per solver) of one folder: counters are cumulative,
           start is their snapshot when the folder began"""
        nFailedSolver, nAbortedSolver, nFailedScriptOrAMPL, nSkipped = counters
        return [(folder, r.getName(), nCases,
                 nFailedSolver[i] - start[0][i], nAbortedSolver[i] - start[1][i],
                 nFailedScriptOrAMPL[i] - start[2][i], nSkipped[i] - start[3][i])
                for i, r in enumerate(self.getRunners())]

    def _summaryLine(self, row):
        """One-line form of a summary row, used right after a folder"""
        folder, slv, nCases, nFail, nAbrt, nScrpt, nSkip = row
        slv = f" {slv}:" if len(self.getRunners()) > 1 else ""
        return "Folder '%s':%s %d cases, %d fail (%d abrt) slv, %d fail AMPL/PY/scrpt, %d skip" % (
            folder, slv, nCases, nFail, nAbrt, nScrpt, nSkip)

    def _printFolderSummary(self, rows):
        for row in rows:
            print(self._summaryLine(row))
        print(ModelRunner.STARS, flush=True)

    @staticmethod
    def _formatTable(header, rows, separatorBefore=None):
        """Plain-text table: columns padded to their widest cell, text columns
           left-aligned, numbers right-aligned; a separator line before the
           row index separatorBefore"""
        cells = [[str(c) for c in row] for row in rows]
        widths = [max(len(h), *(len(r[i]) for r in cells)) for i, h in enumerate(header)]
        def fmt(row):
            return "  ".join(str(c).rjust(w) if isinstance(c, (int, float)) else str(c).ljust(w)
                             for c, w in zip(row, widths))
        separator = "  ".join("-" * w for w in widths)
        lines = [fmt(header), separator]
        for k, row in enumerate(rows):
            if k == separatorBefore:
                lines.append(separator)
            lines.append(fmt(row))
        return lines

    def _printFinalSummary(self, folderRows, nCases, counters):
        """Per-folder recap (multi-folder runs) and overall summary as a table"""
        nFailedSolver, nAbortedSolver, nFailedScriptOrAMPL, nSkipped = counters
        overall = [("Overall", r.getName(), nCases, nFailedSolver[i], nAbortedSolver[i],
                    nFailedScriptOrAMPL[i], nSkipped[i]) for i, r in enumerate(self.getRunners())]
        print("")
        print(ModelRunner.STARS)
        if len(folderRows) > len(overall):      # more than one folder
            print("Summary per folder and overall")
            rows = folderRows + overall
            sep = len(folderRows)
        else:
            print("Overall summary")
            rows = overall
            sep = None
        print(ModelRunner.STARS)
        for l in ModelRunner._formatTable(ModelRunner.SUMMARY_HEADER, rows, sep):
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
        folderRows = []          # summary rows (folder, solver, counts) in run order
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
                    rows = self._folderSummaryRows(curFolder, nFolderCases, counters, folderStart)
                    self._printFolderSummary(rows)
                    folderRows.extend(rows)
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
                            (m.hasSolvers() and ss.getName() not in m.getSolvers()) or \
                            (m.hasSolversExcluded() and ss.getName() in m.getSolversExcluded()):
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
            rows = self._folderSummaryRows(curFolder, nFolderCases, counters, folderStart)
            self._printFolderSummary(rows)
            folderRows.extend(rows)
            self._printFinalSummary(folderRows, n, counters)
        exporter.export()

    def export(self, exporter):
        exporter.exportInstanceResults(self)


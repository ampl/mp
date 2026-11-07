import unittest
from pathlib import Path
import subprocess
import sys
import time
from unittest.mock import Mock, patch

from AMPLRunner import AMPLRunner


class ScriptTimeoutTest(unittest.TestCase):
    def make_runner(self, timeout):
        runner = AMPLRunner.__new__(AMPLRunner)
        runner._timeout = timeout
        runner._ampl = Mock()
        runner.readModel = Mock(return_value="read result")
        return runner

    def make_model(self, script):
        model = Mock()
        model.getFilePath.return_value = "/tmp/cases/example.run"
        model.isScript.return_value = script
        model.doCd.return_value = False
        return model

    @patch("AMPLRunner.TimeoutFunc.run_with_timeout_cleanup_multiprocessing_watcher")
    def test_script_runs_inside_timeout(self, watch):
        runner, model = self.make_runner(30), self.make_model(True)
        watch.side_effect = lambda target, cleanup, first, second: target()
        self.assertEqual("read result", runner.doReadModel(model))
        runner.readModel.assert_called_once_with(model)
        runner._ampl.cd.assert_called_once_with(str(Path("/tmp/cases").resolve()))
        self.assertEqual((30, 5), watch.call_args.args[2:])

    @patch("AMPLRunner.TimeoutFunc.run_with_timeout_cleanup_multiprocessing_watcher")
    def test_unlimited_script_does_not_start_watcher(self, watch):
        runner, model = self.make_runner(None), self.make_model(True)
        self.assertEqual("read result", runner.doReadModel(model))
        watch.assert_not_called()

    @patch("AMPLRunner.TimeoutFunc.run_with_timeout_cleanup_multiprocessing_watcher")
    def test_model_read_is_not_solve_execution(self, watch):
        runner, model = self.make_runner(30), self.make_model(False)
        self.assertEqual("read result", runner.doReadModel(model))
        watch.assert_not_called()
        runner._ampl.cd.assert_not_called()

    def test_script_timeout_stops_a_solver_descendant(self):
        runner, model = self.make_runner(1), self.make_model(True)
        # Match AMPL -> solver: the sleeper is not an immediate child of
        # this test process, so the normal descendant cleanup must stop it.
        wrapper = (
            "import subprocess,sys; subprocess.run([sys.executable, '-c', "
            "'import time; time.sleep(8)'])")
        runner.readModel.side_effect = lambda _: subprocess.run(
            [sys.executable, "-c", wrapper], check=True)
        start = time.monotonic()
        runner.doReadModel(model)
        self.assertLess(time.monotonic() - start, 6)


if __name__ == "__main__":
    unittest.main()

"""Regression tests for the AMPL model layout and NL comparison tools."""

import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import nlcheck
from amplfmt import format_model


class FormatterTest(unittest.TestCase):
    def test_comment_does_not_consume_next_statement(self):
        result = format_model("var x; # comment with ; and (\nmaximize obj: x;\n")
        self.assertIn("# comment with ; and (\nmaximize obj: x;", result)

    def test_literal_contents_and_delimiters_are_preserved(self):
        text = "param p symbolic := 'one  two; # it''s text';\n"
        self.assertEqual(text, format_model(text))


class NLCheckTest(unittest.TestCase):
    def check_model(self, runner, text="var  x;\n"):
        with tempfile.TemporaryDirectory() as directory:
            model = Path(directory) / "model.mod"
            model.write_text(text)
            with patch.object(nlcheck.subprocess, "run", side_effect=runner):
                passed = nlcheck.check_dir("ampl", directory, apply=True)
            return passed, model.read_text()

    def test_ampl_failure_does_not_apply(self):
        def failed(argv, **kwargs):
            Path(argv[1]).with_suffix(".nl").write_bytes(b"header\npartial\n")
            return subprocess.CompletedProcess(argv, 1, "", "AMPL failed")

        self.assertEqual((False, "var  x;\n"), self.check_model(failed))

    def test_missing_nl_does_not_apply(self):
        def no_nl(argv, **kwargs):
            return subprocess.CompletedProcess(argv, 0, "", "")

        text = "param  constant := 1;\n"
        self.assertEqual((False, text), self.check_model(no_nl, text))

    def test_matching_nl_is_applied(self):
        def matched(argv, **kwargs):
            Path(argv[1]).with_suffix(".nl").write_bytes(b"header\nsame model\n")
            return subprocess.CompletedProcess(argv, 0, "", "")

        self.assertEqual((True, "var x;\n"), self.check_model(matched))

    def test_source_edit_after_original_nl_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            model = Path(directory) / "model.mod"
            model.write_text("var  x;\n")
            original_writer = nlcheck.write_nls

            def matched(argv, **kwargs):
                Path(argv[1]).with_suffix(".nl").write_bytes(b"header\nsame model\n")
                return subprocess.CompletedProcess(argv, 0, "", "")

            def edit_after_original(ampl, models, outdir, tag):
                runs = original_writer(ampl, models, outdir, tag)
                if tag == "o":
                    model.write_text("var y;\n")
                return runs

            with (patch.object(nlcheck.subprocess, "run", side_effect=matched),
                  patch.object(nlcheck, "write_nls", side_effect=edit_after_original)):
                self.assertFalse(nlcheck.check_dir("ampl", directory, apply=True))
            self.assertEqual("var y;\n", model.read_text())


if __name__ == "__main__":
    unittest.main()

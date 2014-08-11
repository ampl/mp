import os, sys, tempfile, util, shutil
from nose2.tests._common import TestCase

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'support'))
from download import Downloader

def readfile(path):
  with open(path, 'r') as f:
    return f.read()

class TestException(Exception):
  pass

class DownloadTest(TestCase):
  def test_download(self):
    d = Downloader()
    with util.CaptureStdout():
      with d.download('file://' + __file__) as f:
        self.assertEqual(readfile(__file__), readfile(f))

  def test_file_removed_on_exception(self):
    try:
      d = Downloader()
      with util.CaptureStdout():
        with d.download('file://' + __file__) as f:
          filename = f
          raise TestException()
    except TestException:
      pass
    self.assertFalse(os.path.exists(filename))

  def test_download_to_temp_dir(self):
    d = Downloader()
    with util.CaptureStdout():
      with d.download('file://' + __file__) as f:
        filename = f
    self.assertEqual(tempfile.gettempdir(), os.path.dirname(filename))

  def test_download_to_dir(self):
    dir = tempfile.mkdtemp()
    try:
      d = Downloader(dir)
      with util.CaptureStdout():
        with d.download('file://' + __file__) as f:
          filename = f
    finally:
      shutil.rmtree(dir)
    self.assertEqual(dir, os.path.dirname(filename))

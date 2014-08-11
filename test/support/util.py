import sys
from cStringIO import StringIO

# Captures output to stdout in a block.
# Usage:
#   stdout = CaptureStdout()
#   with stdout:
#     print('something')
#   output = stdout.str()
class CaptureStdout:
  def __enter__(self):
    self.old_stdout = sys.stdout
    sys.stdout = self.out = StringIO()

  def __exit__(self, type, value, traceback):
    sys.stdout = self.old_stdout

  def str(self):
    return self.out.getvalue()

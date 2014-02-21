# File utils.

import errno, shutil

# Delete an entire directory tree if it exists.
def rmtree_if_exists(path):
  try:
    shutil.rmtree(path)
  except OSError as e:
    if e.errno != errno.ENOENT:
      raise

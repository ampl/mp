# File utils.

import errno, os, shutil

# Delete an entire directory tree if it exists.
def rmtree_if_exists(path):
  try:
    shutil.rmtree(path)
  except OSError as e:
    if e.errno != errno.ENOENT:
      raise

# Delete a file if it exists.
def remove_if_exists(path):
  try:
    os.remove(path)
  except OSError as e:
    if e.errno != errno.ENOENT:
      raise

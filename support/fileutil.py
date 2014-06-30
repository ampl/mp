# File utils.

from __future__ import print_function
import errno, os, shutil, zipfile

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

UNIX = 3

# Archives a directory with zip preserving symlinks. The created archive may not
# be supported by some zip utilities.
def make_archive(archive_name, dirname):
  with zipfile.ZipFile(archive_name, 'w', compression=zipfile.ZIP_DEFLATED) as zip:
    for root, dirs, files in os.walk(dirname):
      for filename in files:
        path = os.path.join(root, filename)
        if os.path.islink(path):
          zipinfo = zipfile.ZipInfo(path)
          zipinfo.create_system = UNIX
          zipinfo.external_attr = 2716663808L
          zip.writestr(zipinfo, os.readlink(path))
        else:
          zip.write(path, path)

def move(filename, target_dir):
  print('Moving', filename, 'to', target_dir)
  os.rename(filename, os.path.join(target_dir, filename))

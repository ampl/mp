# Common bootstrap functionality.

from __future__ import print_function
import os, platform, re, shutil, sys, tarfile, urllib2, urlparse, zipfile

class TempFile:
  def __init__(self, filename):
    self.filename = filename
  def __enter__(self):
    return self.filename
  def __exit__(self, type, value, traceback):
    os.remove(self.filename)

# Downloads into a temporary file.
def download(url, cookie = None):
  filename = os.path.basename(urlparse.urlsplit(url)[2])
  print('Downloading', url, 'to', filename)
  sys.stdout.flush()
  opener = urllib2.build_opener()
  if cookie:
    opener.addheaders.append(('Cookie', cookie))
  with open(filename, 'wb') as f:
    shutil.copyfileobj(opener.open(url), f)
  return TempFile(filename)

windows = platform.system() == 'Windows'
opt_dir = r'\Program Files' if windows else '/opt'

# Return the path to an executable which would be run if the given cmd was called.
def which(cmd):
  if windows:
    cmd += '.exe'
  for path in os.environ['PATH'].split(os.pathsep):
    filename = os.path.join(path.strip('"'), cmd)
    if os.path.isfile(filename) and os.access(filename, os.X_OK):
      return filename
  return None

# Adds filename to search paths.
def add_to_path(filename):
  if windows:
    return
  path = '/usr/local/bin'
  # Create /usr/local/bin directory if it doesn't exist (can happen on OS X).
  if not os.path.exists(path):
    os.makedirs(path)
  linkname = os.path.join(path, os.path.basename(filename))
  print('Creating a symlink from', linkname, 'to', filename)
  os.symlink(filename, linkname)

# Downloads and installs CMake.
# filename: The name of a CMake archive, e.g. 'cmake-2.8.12.2-Linux-i386.tar.gz'.
def install_cmake(filename):
  path = which('cmake')
  if path:
    print('CMake already installed in', path)
    return
  dir, version, minor = re.match(
    r'(cmake-(\d+\.\d+)\.(\d+)\.[^\.]+)\..*', filename).groups()
  archiver = zipfile.ZipFile if filename.endswith('zip') else tarfile
  # extractall overwrites existing files, so no need to prepare the destination.
  with download('http://www.cmake.org/files/v{}/{}'.format(version, filename)) as f:
    with archiver.open(f, "r:gz") as archive:
      archive.extractall(opt_dir)
  if platform.system() == 'Darwin':
    dir = os.path.join(dir, 'CMake {}-{}.app'.format(version, minor), 'Contents')
  add_to_path(os.path.join(opt_dir, dir, 'bin', 'cmake'))

# Returns true iff module exists.
def module_exists(module):
  try:
    importlib.import_module(module)
    return True
  except ImportError:
    return False

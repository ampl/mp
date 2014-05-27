# Common bootstrap functionality.

from __future__ import print_function
import os, platform, re, shutil, sys, tarfile, urllib2, urlparse, zipfile
from subprocess import check_call

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

# If we are in a VM managed by Vagrant, then do everything in the shared
# /vagrant directory to avoid growth of the VM drive.
def bootstrap_init():
  vagrant_dir = '/vagrant'
  if os.path.exists(vagrant_dir):
    os.chdir(vagrant_dir)
    return True
  return False

# Returns true if executable is installed on the path.
def installed(name):
  filename = name + '.exe' if windows else name
  for path in os.environ['PATH'].split(os.pathsep):
    path = os.path.join(path.strip('"'), filename)
    if os.path.isfile(path) and os.access(path, os.X_OK):
      print(name, 'is installed in', path)
      return True
  return False

# Adds filename to search paths.
def add_to_path(filename, linkname = None):
  if windows:
    path = os.path.basename(filename)
    print('Adding', path, 'to PATH')
    paths = os.getenv('PATH').split(os.pathsep)
    paths.append(path)
    check_call(['setx', 'PATH', os.pathsep.join(paths)])
    return
  path = '/usr/local/bin'
  # Create /usr/local/bin directory if it doesn't exist which can happen
  # on OS X.
  if not os.path.exists(path):
    os.makedirs(path)
  linkname = os.path.join(path, linkname or os.path.basename(filename))
  print('Creating a symlink from', linkname, 'to', filename)
  os.symlink(filename, linkname)

# Downloads and installs CMake.
# filename: The name of a CMake archive,
# e.g. 'cmake-2.8.12.2-Linux-i386.tar.gz'.
def install_cmake(filename):
  if installed('cmake'):
    return
  dir, version, minor = re.match(
    r'(cmake-(\d+\.\d+)\.(\d+)\.[^\.]+)\..*', filename).groups()
  # extractall overwrites existing files, so no need to prepare the destination.
  url = 'http://www.cmake.org/files/v{}/{}'.format(version, filename)
  with download(url) as f:
    iszip = filename.endswith('zip')
    with zipfile.ZipFile(f) if iszip else tarfile.open(f, 'r:gz') as archive:
      archive.extractall(opt_dir)
  if platform.system() == 'Darwin':
    dir = os.path.join(
      dir, 'CMake {}-{}.app'.format(version, minor), 'Contents')
  add_to_path(os.path.join(opt_dir, dir, 'bin', 'cmake'))

# Installs symlinks for ccache.
def install_ccache_links():
  for name in ['gcc', 'cc', 'g++', 'c++']:
    add_to_path('/usr/bin/ccache', name)

# Install f90cache.
def install_f90cache():
  if not installed('f90cache'):
    f90cache = 'f90cache-0.95'
    with download(
        'http://people.irisa.fr/Edouard.Canot/f90cache/' + f90cache + '.tar.bz2') as f:
      with tarfile.open(f, "r:bz2") as archive:
        archive.extractall('.')
    check_call(['./configure'], cwd=f90cache)
    check_call(['make', 'all', 'install'], cwd=f90cache)
    shutil.rmtree(f90cache)
    add_to_path('/usr/local/bin/f90cache', 'gfortran-4.9')

# Returns true iff module exists.
def module_exists(module):
  try:
    importlib.import_module(module)
    return True
  except ImportError:
    return False

# Installs pip.
def install_pip():
  if not module_exists('pip'):
    with download('https://bootstrap.pypa.io/get-pip.py') as f:
      check_call(['python', f])

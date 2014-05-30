# Common bootstrap functionality.

from __future__ import print_function
import glob, os, platform, re, shutil, sys
import tarfile, tempfile, urllib2, urlparse, zipfile
from contextlib import closing
from subprocess import check_call

# Downloads into a temporary file.
def download(url, cookie = None):
  suffix = os.path.splitext(urlparse.urlsplit(url)[2])[1]
  file = tempfile.NamedTemporaryFile(suffix=suffix, dir='')
  print('Downloading', url, 'to', file.name)
  sys.stdout.flush()
  opener = urllib2.build_opener()
  if cookie:
    opener.addheaders.append(('Cookie', cookie))
  shutil.copyfileobj(opener.open(url), file)
  file.flush()
  file.seek(0)
  return file

windows = platform.system() == 'Windows'
opt_dir = r'\Program Files' if windows else '/opt'

# If we are in a VM managed by Vagrant, then do everything in the shared
# /vagrant directory to avoid growth of the VM drive.
vagrant_dir = '/vagrant'
def bootstrap_init():
  if os.path.exists(vagrant_dir):
    os.chdir(vagrant_dir)
    return True
  return False

# Returns executable location.
def which(name):
  filename = name + '.exe' if windows else name
  for path in os.environ['PATH'].split(os.pathsep):
    path = os.path.join(path.strip('"'), filename)
    if os.path.isfile(path) and os.access(path, os.X_OK):
      return path
  return None

# Returns true if executable is installed on the path.
def installed(name):
  path = which(name)
  if path:
    print(name, 'is installed in', path)
  return path != None

# Adds filename to search paths.
def add_to_path(filename, linkname = None):
  paths = os.environ['PATH'].split(os.pathsep)
  path = os.path.dirname(filename)
  print('Adding', path, 'to PATH')
  paths.append(path)
  os.environ['PATH'] = os.pathsep.join(paths)
  if windows:
    check_call(['setx', 'PATH', os.environ['PATH']])
    return
  path = '/usr/local/bin'
  # Create /usr/local/bin directory if it doesn't exist which can happen
  # on OS X.
  if not os.path.exists(path):
    os.makedirs(path)
  linkname = os.path.join(path, linkname or os.path.basename(filename))
  if not os.path.exists(linkname):
    print('Creating a symlink from', linkname, 'to', filename)
    os.symlink(filename, linkname)
  else:
    print('File already exists:', linkname)

# Downloads and installs CMake.
# filename: The name of a CMake archive,
# e.g. 'cmake-2.8.12.2-Linux-i386.tar.gz'.
def install_cmake(filename):
  if installed('cmake'):
    return
  dir, version, minor = re.match(
    r'(cmake-(\d+\.\d+)\.(\d+)\.[^\.]+)\..*', filename).groups()
  # extractall overwrites existing files, so no need to prepare the
  # destination.
  url = 'http://www.cmake.org/files/v{0}/{1}'.format(version, filename)
  with download(url) as f:
    iszip = filename.endswith('zip')
    with zipfile.ZipFile(f) if iszip \
         else closing(tarfile.open(None, 'r:gz', fileobj=f)) as archive:
      archive.extractall(opt_dir)
  if platform.system() == 'Darwin':
    dir = os.path.join(
      dir, 'CMake {0}-{1}.app'.format(version, minor), 'Contents')
  add_to_path(os.path.join(opt_dir, dir, 'bin', 'cmake'))

# Install f90cache.
def install_f90cache():
  if not installed('f90cache'):
    f90cache = 'f90cache-0.95'
    with download(
        'http://people.irisa.fr/Edouard.Canot/f90cache/' +
        f90cache + '.tar.bz2') as f:
      with closing(tarfile.open(None, "r:bz2", fileobj=f)) as archive:
        archive.extractall('.')
    check_call(['sh', 'configure'], cwd=f90cache)
    check_call(['make', 'all', 'install'], cwd=f90cache)
    shutil.rmtree(f90cache)
    add_to_path('/usr/local/bin/f90cache', 'gfortran-4.9')

# Returns true iff module exists.
def module_exists(module):
  try:
    __import__(module)
    return True
  except ImportError:
    return False

# Install package using pip if it hasn't been installed already.
def pip_install(package, test_module=None):
  if not test_module:
    test_module = package
  if module_exists(test_module):
    return
  # If pip doesn't exist install it first.
  if not module_exists('pip'):
    with download('https://bootstrap.pypa.io/get-pip.py') as f:
      check_call(['python', f.name])
  # Install the package.
  print('Installing', package)
  from pip.index import PackageFinder
  from pip.req import InstallRequirement, RequirementSet
  from pip.locations import build_prefix, src_prefix
  requirement_set = RequirementSet(
      build_dir=build_prefix,
      src_dir=src_prefix,
      download_dir=None)
  requirement_set.add_requirement(InstallRequirement.from_line(package, None))
  finder = PackageFinder(
    find_links=[], index_urls=['http://pypi.python.org/simple/'])
  requirement_set.prepare_files(finder, force_root_egg_info=False, bundle=False)
  requirement_set.install([], [])

# Installs buildbot slave.
def install_buildbot_slave(name, path=None, script_dir='', shell=False):
  path = path or os.path.expanduser('~vagrant/slave')
  if not os.path.exists(vagrant_dir) or os.path.exists(path):
    return
  pip_install('buildbot-slave', 'buildbot')
  # The password is insecure but it doesn't matter as the buildslaves are
  # not publicly accessible.
  command = [os.path.join(script_dir, 'buildslave'),
             'create-slave', path, '10.0.2.2', name, 'pass']
  if not windows:
    command = ['sudo', '-u', 'vagrant'] + command
  check_call(command, shell=shell)
  if windows:
    return
  pip_install('python-crontab', 'crontab')
  from crontab import CronTab
  cron = CronTab('vagrant')
  cron.new('PATH={0}:/usr/local/bin buildslave start {1}'.format(
    os.environ['PATH'], path))
  cron.write()
  check_call(['sudo', '-H', '-u', 'vagrant', 'buildslave', 'start', path])

# Copies optional dependencies from opt/<platform> to /opt.
def copy_optional_dependencies(platform):
  if os.path.exists('opt'):
    for dir in glob.glob('opt/' + platform + '/*'):
      dest = '/opt/' + dir
      if not os.path.exists(dest):
        shutil.copytree(dir, dest)

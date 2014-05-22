# Set up build environment on Windows.

from __future__ import print_function
import importlib, os, sys, shutil, urllib2, urlparse
from glob import glob
from subprocess import check_call
from zipfile import ZipFile

class TempFile:
  def __init__(self, filename):
    self.filename = filename
  def __enter__(self):
    return self.filename
  def __exit__(self, type, value, traceback):
    os.remove(self.filename)

# Returns true iff module exists.
def module_exists(module):
  try:
    importlib.import_module(module)
    return True
  except ImportError:
    return False

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

def unzip(filename, path):
  with ZipFile(filename) as zip:
    zip.extractall(path)

# Install CMake.
cmake = 'cmake-2.8.12.2-win32-x86'
cmake_dir = os.path.join(r'C:\Program Files', cmake)
if not os.path.exists(cmake_dir):
  with download('http://www.cmake.org/files/v2.8/' + cmake + '.zip') as f:
    unzip(f, r'C:\Program Files')

# Add Python and CMake to PATH.
python_dir = r'C:\Python27'
check_call(['setx', 'PATH',
  os.getenv('PATH') + ';' + python_dir + ';' + os.path.join(cmake_dir, 'bin')])

# Install .NET Framework 4 for msbuild.
if not os.path.exists(r'\Windows\Microsoft.NET\Framework64\v4.0.30319'):
  with download(
      'http://download.microsoft.com/download/9/5/A/' +
      '95A9616B-7A37-4AF6-BC36-D6EA96C8DAAE/dotNetFx40_Full_x86_x64.exe') as f:
    check_call([f, '/q', '/norestart'])

# Install 7zip.
sevenzip = r'C:\Program Files (x86)\7-Zip\7z.exe'
if not os.path.exists(sevenzip):
  with download('http://downloads.sourceforge.net/sevenzip/7z920.exe') as f:
    check_call([f, '/S'])

# Install Windows SDK.
if not os.path.exists(r'\Program Files\Microsoft SDKs\Windows\v7.1'):
  # Extract ISO.
  with download(
       'http://download.microsoft.com/download/F/1/0/'
       'F10113F5-B750-4969-A255-274341AC6BCE/GRMSDKX_EN_DVD.iso') as f:
    check_call([sevenzip, 'x', '-tudf', '-owinsdk', f])
  # Install SDK.
  check_call([r'winsdk\setup.exe', '-q'])
  shutil.rmtree('winsdk')

# Install MinGW.
def install_mingw(arch):
  bits = '64' if arch.endswith('64') else '32'
  if os.path.exists(r'\mingw' + bits):
    return
  with download(
      'http://sourceforge.net/projects/mingw-w64/files/' +
      'Toolchains%20targetting%20Win' + bits + '/Personal%20Builds/' +
      'mingw-builds/4.8.2/threads-win32/sjlj/' + arch +
      '-4.8.2-release-win32-sjlj-rt_v3-rev4.7z/download') as f:
    check_call([sevenzip, 'x', '-oC:\\', f])

install_mingw('i686')
install_mingw('x86_64')

# Install pywin32 - buildbot dependency.
if not module_exists('win32api'):
  shutil.rmtree('pywin32', True)
  with download(
      'http://sourceforge.net/projects/pywin32/files/pywin32/Build%20219/' +
      'pywin32-219.win-amd64-py2.7.exe/download') as f:
    check_call([sevenzip, 'x', '-opywin32', f])
  site_packages_dir = os.path.join(python_dir, r'lib\site-packages')
  for path in glob('pywin32/PLATLIB/*') + glob('pywin32/SCRIPTS/*'):
    shutil.move(path, site_packages_dir)
  shutil.rmtree('pywin32')
  import pywin32_postinstall
  pywin32_postinstall.install()
  os.remove(site_packages_dir + r'\pywin32_postinstall.py')

# Install pip.
if not module_exists('pip'):
  with download('https://bootstrap.pypa.io/get-pip.py') as f:
    check_call(['python', f])

from pip.index import PackageFinder
from pip.req import InstallRequirement, RequirementSet
from pip.locations import build_prefix, src_prefix

# Install package using pip if it hasn't been installed already.
def pip_install(package, test_module=None):
  if not test_module:
    test_module = package
  if module_exists(test_module):
    return
  print('Installing', package)
  requirement_set = RequirementSet(
      build_dir=build_prefix,
      src_dir=src_prefix,
      download_dir=None)
  requirement_set.add_requirement(InstallRequirement.from_line(package, None))
  finder = PackageFinder(
    find_links=[], index_urls=['http://pypi.python.org/simple/'])
  requirement_set.prepare_files(finder, force_root_egg_info=False, bundle=False)
  requirement_set.install([], [])

# Install buildbot dependencies.
pip_install('twisted')
pip_install('zope.interface')

# Install buildbot slave.
pip_install('buildbot-slave', 'buildbot')

buildslave_dir = r'\buildslave'
if not os.path.exists(buildslave_dir):
  # Create buildbot slave.
  check_call([os.path.join(python_dir, r'scripts\buildslave.bat'),
              'create-slave', buildslave_dir, '10.0.2.2', 'win2008', 'pass'])

  # Grant the user the right to "log on as a service".
  import win32api, win32security
  username = win32api.GetUserNameEx(win32api.NameSamCompatible)
  domain, username = username.split('\\')
  policy = win32security.LsaOpenPolicy(domain, win32security.POLICY_ALL_ACCESS)
  sid_obj, domain, tmp = win32security.LookupAccountName(domain, username)
  win32security.LsaAddAccountRights(policy, sid_obj, ('SeServiceLogonRight',))
  win32security.LsaClose(policy)

  # Set buildslave parameters.
  import _winreg as reg
  with reg.CreateKey(reg.HKEY_LOCAL_MACHINE,
      r'System\CurrentControlSet\services\BuildBot\Parameters') as key:
    reg.SetValueEx(key, 'directories', 0, reg.REG_SZ, buildslave_dir)

  # Install buildbot service.
  check_call([
    os.path.join(python_dir, 'python'),
    os.path.join(python_dir, r'Scripts\buildbot_service.py'),
    '--user', '.\\' + username, '--password', 'vagrant',
    '--startup', 'auto', 'install'])
  import win32serviceutil
  win32serviceutil.StartService('BuildBot')

# Install 32-bit JDK.
cookie = 'oraclelicense=accept-securebackup-cookie'
if not os.path.exists(r'\Program Files (x86)\Java\jdk1.7.0_55'):
  with download(
      'http://download.oracle.com/otn-pub/java/jdk/7u55-b13/' +
      'jdk-7u55-windows-i586.exe', cookie) as f:
    check_call([f, '/s'])

# Install 64-bit JDK.
if not os.path.exists(r'\Program Files\Java\jdk1.7.0_55'):
  with download(
      'http://download.oracle.com/otn-pub/java/jdk/7u55-b13/' +
      'jdk-7u55-windows-x64.exe', cookie) as f:
    check_call([f, '/s'])

# Copy optional dependencies.
opt_dir = r'opt\win64'
if os.path.exists(opt_dir):
  for entry in os.listdir(opt_dir):
    subdir = os.path.join(opt_dir, entry)
    for subentry in os.listdir(subdir):
      dest = os.path.join('C:\\', entry, subentry)
      if not os.path.exists(dest):
        shutil.copytree(os.path.join(subdir, subentry), dest)

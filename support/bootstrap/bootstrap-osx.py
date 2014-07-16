#!/usr/bin/env python
# Set up build environment on OS X Moutain Lion.

from __future__ import print_function
from bootstrap import *
import glob, os, sys, tempfile
from subprocess import check_call

vagrant = bootstrap_init()

install_cmake('cmake-2.8.12.2-Darwin64-universal.tar.gz')

# Installs an OS X package.
def install_pkg(filename):
  print('Installing', filename)
  check_call(['installer', '-pkg', filename, '-target', '/'])

# Installs a package from a .dmg file.
def install_dmg(filename):
  dir = tempfile.mkdtemp()
  check_call(['hdiutil', 'attach', filename, '-mountpoint', dir])
  install_pkg(glob.glob(dir + '/*pkg')[0])
  check_call(['hdiutil', 'detach', dir])
  os.rmdir(dir)

# Install command-line tools for Xcode.
if not installed('clang'):
  with download(
      'http://devimages.apple.com/downloads/xcode/' +
      'command_line_tools_for_xcode_os_x_mountain_lion_april_2013.dmg') as f:
    install_dmg(f)

# Install MacPorts.
if not installed('port'):
  with download(
      'https://distfiles.macports.org/MacPorts/' +
      'MacPorts-2.2.0-10.8-MountainLion.pkg') as f:
    install_pkg(f)
  # Get rid of "No Xcode installation was found" error.
  macports_conf = '/opt/local/etc/macports/macports.conf'
  macports_conf_new = macports_conf + '.new'
  if os.path.exists(macports_conf):
    with open(macports_conf, 'r') as f:
      conf = f.read()
    with open(macports_conf_new, 'w') as f:
      f.write(re.sub(r'#\s*(developer_dir\s*).*', r'\1/', conf))
    os.remove(macports_conf)
  # The new config may also exists if the script was interrupted just after
  # removing the old config and then restarted.
  if os.path.exists(macports_conf_new):
    os.rename(macports_conf_new, macports_conf)
  add_to_path('/opt/local/bin/port')

# Install ccache.
if not installed('ccache'):
  check_call(['port', 'install', 'ccache'])
  home = os.path.expanduser('~vagrant' if vagrant else '~')
  with open(home + '/.profile', 'a') as f:
    f.write('export PATH=/opt/local/libexec/ccache:$PATH')
  add_to_path('/opt/local/libexec/ccache', None, True)

# Install gfortran.
if not installed('gfortran-4.9'):
  check_call(['port', 'install', 'gcc49', '+gfortran'])
  add_to_path('/opt/local/bin/gfortran-mp-4.9', 'gfortran-4.9')

install_f90cache()
create_symlink('/usr/local/bin/f90cache',
               '/opt/local/libexec/ccache/gfortran-4.9')

if not installed('localsolver'):
  with download('http://www.localsolver.com/downloads/' +
      'LocalSolver_4_0_20140127_MacOS64.pkg') as f:
    install_dmg(f)

copy_optional_dependencies('osx')
install_buildbot_slave('osx-ml')

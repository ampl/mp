#!/usr/bin/env python
# Set up build environment on OS X Moutain Lion.

from __future__ import print_function
from bootstrap import *
import glob, mmap, os, sys, tempfile
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
  with open('/opt/local/etc/macports/macports.conf', 'r+b') as f:
    m = mmap.mmap(f.fileno(), 0)
    pos = m.find('# developer_dir')
    if pos != -1:
      m.seek(pos)
      m.write('developer_dir /\n#')
  add_to_path('/opt/local/bin/port')

# Install ccache.
if not installed('ccache'):
  check_call(['port', 'install', 'ccache'])
  home = os.path.expanduser('~vagrant' if vagrant else '~')
  with open(home + '/.profile', 'a') as f:
    f.write('export PATH=/opt/local/libexec/ccache:$PATH')

# Install gfortran.
if not installed('gfortran'):
  check_call(['port', 'install', 'gcc49', '+gfortran'])

install_f90cache()

copy_optional_dependencies('osx')

if vagrant:
  buildslave_dir = '/home/vagrant/slave'
  if not os.path.exists(buildslave_dir):
    install_buildbot_slave('osx-ml', buildslave_dir)

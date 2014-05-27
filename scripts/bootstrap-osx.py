#!/usr/bin/env python
# Set up build environment on OS X Moutain Lion.

import os, sys, tempfile
from glob import glob
from subprocess import check_call

sys.path.append('/vagrant/scripts')
from bootstrap import *
vagrant = bootstrap_init()

install_cmake('cmake-2.8.12.2-Darwin64-universal.tar.gz')

# Installs an OS X package.
def install_pkg(filename):
  check_call(['installer', '-pkg', filename, '-target', '/'])

# Installs a package from a .dmg file.
def install_dmg(filename):
  dir = tempfile.mkdtemp()
  check_call(['hdiutil', 'attach', filename, '-mountpoint', dir])
  install_pkg(glob(dir + '/*pkg')[0])
  check_call(['hdiutil', 'detach', dir])
  os.rmdir(dir)

# Install command-line tools for Xcode.
if not installed('clang'):
  with download(
      'http://devimages.apple.com/downloads/xcode/' +
      'command_line_tools_for_xcode_os_x_mountain_lion_april_2013.dmg') as f:
    install_dmg(f)

# Install MacPorts.
with download(
    'https://distfiles.macports.org/MacPorts/' +
    'MacPorts-2.2.0-10.8-MountainLion.pkg') as f:
  install_pkg(f)
  if vagrant:
    with open('/Users/vagrant', 'a') as f:
      f.write('export PATH=/opt/local/bin:/opt/local/sbin:$PATH\n')

# Install gfortran
if not installed('gfortran'):
  check_call(['port', 'install', 'gcc49', '+gfortran'])

install_f90cache()

# TODO: install buildbot, ccache, fortran cache

#!/usr/bin/env python
# Set up build environment on OS X.

import os, tempfile
from subprocess import check_call

# If we are in a VM managed by Vagrant, then do everything in the shared
# /vagrant directory to avoid growth of the VM drive.
vagrant_dir = '/vagrant'
if os.path.exists(vagrant_dir):
  os.chdir(vagrant_dir)
  import sys
  sys.path.append(vagrant_dir + '/scripts')

import bootstrap

bootstrap.install_cmake('cmake-2.8.12.2-Darwin64-universal.tar.gz')

# Install command-line tools for Xcode.
if not os.path.exists('/usr/bin/clang'):
  with bootstrap.download(
    'http://devimages.apple.com/downloads/xcode/' +
    'command_line_tools_for_xcode_os_x_mountain_lion_april_2013.dmg') as f:
    dir = tempfile.mkdtemp()
    check_call(['hdiutil', 'attach', f, '-mountpoint', dir])
    check_call(['installer', '-pkg',
                dir + '/Command Line Tools (Mountain Lion).mpkg', '-target', '/'])
    check_call(['hdiutil', 'detach', dir])
    os.rmdir(dir)

with download('http://coudert.name/software/gfortran-4.8.2-MountainLion.dmg') as f:
  pass # TODO

# TODO: install buildbot, gfortran, ccache

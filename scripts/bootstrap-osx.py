#!/usr/bin/env python
# Set up build environment on OS X.

import os, tempfile

vagrant_dir = '/vagrant'
if os.path.exists(vagrant_dir):
  import sys
  sys.path.append(vagrant_dir + '/scripts')

from download import download

# Install command-line tools for Xcode.
with download(
  'http://devimages.apple.com/downloads/xcode/' +
  'command_line_tools_for_xcode_os_x_mountain_lion_april_2013.dmg') as f:
  dir = tempfile.mkdtemp()
  check_call(['hdiutil', 'attach', f, '-mountpoint', dir])
  check_call(['installer', '-pkg',
              'Command Line Tools (Mountain Lion).mpkg', '-target', '/'])
  check_call(['hdiutil', 'detach', dir])
  os.rmdir(dir)

# TODO: install buildbot

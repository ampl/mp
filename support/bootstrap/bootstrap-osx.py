#!/usr/bin/env python
# Set up build environment on OS X Moutain Lion.

from __future__ import print_function
from bootstrap import *
import glob, os, sys, tempfile
from subprocess import call, check_call

vagrant = bootstrap_init()

class flushfile:
  def __init__(self, f):
    self.f = f

  def __getattr__(self, name):
    return object.__getattribute__(self.f, name)

  def write(self, x):
    self.f.write(x)
    self.f.flush()

# Flush the buffer on every output.
sys.stdout = flushfile(sys.stdout)

def run(cmd, **kwargs):
  print('Running', cmd)
  check_call(cmd, **kwargs)

sudo = ['sudo', '-H', '-u', 'vagrant'] if vagrant else []

if vagrant:
  # Disable the screensaver because it breaks GUI tests and consumes resources.
  run(sudo + ['defaults', 'write', 'com.apple.screensaver', 'idleTime', '0'])
  # Kill the screensaver if it has started already.
  call(['killall', 'ScreenSaverEngine'])

install_cmake('cmake-3.3.0-Darwin-x86_64.tar.gz')
install_maven()

# Install command-line tools for Xcode.
if not installed('clang'):
  with download(
      'http://devimages.apple.com/downloads/xcode/' +
      'command_line_tools_for_xcode_os_x_mountain_lion_april_2013.dmg') as f:
    install_dmg(f, True)

# Install Homebrew.
if not installed('brew'):
  with download(
      'https://raw.githubusercontent.com/Homebrew/install/' +
      'b43a27f13a2f6750608cd01a6e724fb5f956b089/install') as f:
    run(sudo + ['ruby', f])

# Install ccache.
if not installed('ccache'):
  run(sudo + ['brew', 'install', 'ccache'])
  home = os.path.expanduser('~vagrant' if vagrant else '~')
  with open(home + '/.profile', 'a') as f:
    f.write('export PATH=/usr/local/opt/ccache/libexec:$PATH')
  add_to_path('/usr/local/opt/ccache/libexec', None, True)

# Install gfortran which is a part of the gcc package in Homebrew.
gfortran = 'gfortran-5.2'
if not installed('gfortran'):
  # Install newer version of curl because curl 7.24.0 fails with error
  # "Unknown SSL protocol error in connection to gmplib.org:443".
  run(sudo + ['brew', 'install', 'curl'])
  env = os.environ.copy()
  env['PATH'] = glob.glob('/usr/local/Cellar/curl/*/bin')[0] + ':' + env['PATH']
  run(sudo + ['brew', 'install', 'gcc'], env=env)
  # Create a symlink in the form gfortran-<major>.<minor> for f90cache.
  create_symlink('/usr/local/bin/gfortran', '/usr/local/bin/' + gfortran)

install_f90cache()
create_symlink('/usr/local/bin/f90cache',
               '/usr/local/opt/ccache/libexec/' + gfortran)

# Install JDK.
if len(glob.glob('/Library/Java/JavaVirtualMachines/jdk1.7.0_79.jdk')) == 0:
  with download(
      'http://download.oracle.com/otn-pub/java/jdk/7u79-b15/' +
      'jdk-7u79-macosx-x64.dmg', jdk_cookie) as f:
    install_dmg(f)

# Install LocalSolver.
if not installed('localsolver'):
  with download(LOCALSOLVER_DOWNLOADS_URL +
      'LocalSolver_{0}_MacOS64.pkg'.format(LOCALSOLVER_VERSION)) as f:
    install_pkg(f)

copy_optional_dependencies('osx', '/mnt')
for dir in ['Xcode.app', 'MATLAB_R2014a.app']:
  if os.path.exists('/opt/' + dir):
    create_symlink('/opt/' + dir, '/Applications/' + dir)

run([vagrant_dir + '/support/bootstrap/accept-xcode-license'])

buildbot_path = install_buildbot_slave('osx-ml')
if buildbot_path:
  # Add buildslave launch agent.
  plist_name = 'buildslave.plist'
  with open(vagrant_dir + '/support/buildbot/' + plist_name, 'r') as f:
    plist_content = f.read()
  plist_content = plist_content.replace('$PATH', os.environ['PATH'])
  launchagents_dir = '/Users/vagrant/Library/LaunchAgents'
  try:
    os.makedirs(launchagents_dir)
  except OSError as e:
    if not (e.errno == errno.EEXIST and os.path.isdir(launchagents_dir)):
      raise
  plist_path = os.path.join(launchagents_dir, plist_name)
  with open(plist_path, 'w') as f:
    f.write(plist_content)
  run(['launchctl', 'load', plist_path])

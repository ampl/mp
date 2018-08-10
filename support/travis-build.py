#!/usr/bin/env python
# Build the project on Travis CI.

from __future__ import print_function
import os, re, shutil, tarfile, tempfile
from bootstrap import bootstrap
from contextlib import closing
from download import Downloader
from subprocess import call, check_call, check_output, Popen, PIPE, STDOUT

build = os.environ['BUILD']
if build == 'doc':
  returncode = 1
  travis = 'TRAVIS' in os.environ
  workdir = tempfile.mkdtemp()
  try:
    doxygen_url = 'http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.10.linux.bin.tar.gz'
    dir = os.path.dirname(os.path.realpath(__file__))
    with Downloader().download(doxygen_url) as f:
      with closing(tarfile.open(f, 'r:gz')) as archive:
        archive.extractall(dir)
    doxygen = os.path.join(dir, 'doxygen-1.8.10/bin/doxygen')
    returncode, repo_dir = __import__('build-docs').build_docs(workdir, doxygen)
    if returncode == 0 and os.environ['TRAVIS_BRANCH'] == 'master':
      # Push docs to GitHub pages if this is a master branch.
      if travis:
        check_call(['git', 'config', '--global', 'user.name', 'amplbot'])
        check_call(['git', 'config', '--global', 'user.email', 'bot@ampl.com'])
      check_call(['git', 'add', '--all'], cwd=repo_dir)
      if call(['git', 'diff-index', '--quiet', 'HEAD'], cwd=repo_dir):
        check_call(['git', 'commit', '-m', 'Update documentation'], cwd=repo_dir)
        cmd = 'git push https://$KEY@github.com/ampl/ampl.github.io.git master'
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, cwd=repo_dir)
        # Print the output without the key.
        print(p.communicate()[0].replace(os.environ['KEY'], '$KEY'))
        returncode = p.returncode
  finally:
    # Don't remove workdir on Travis because the VM is discarded anyway.
    if not travis:
      shutil.rmtree(workdir)
  exit(returncode)

cmake_flags = ['-DBUILD=all']

# This is not used right now because container-based infrastructure on Travis
# doesn't support sudo.
if build == 'cross':
  cmake_flags = [
    '-DCMAKE_SYSTEM_NAME=Windows',
    '-DCMAKE_SYSTEM_PROCESSOR=x86_64',
    '-DBUILD_SHARED_LIBS:BOOL=ON',
    '-DCMAKE_C_COMPILER=/opt/mingw64/bin/x86_64-w64-mingw32-gcc',
    '-DCMAKE_CXX_COMPILER=/opt/mingw64/bin/x86_64-w64-mingw32-g++',
    '-DCMAKE_RC_COMPILER=/opt/mingw64/bin/x86_64-w64-mingw32-windres',
    '-DCMAKE_FIND_ROOT_PATH=/opt/mingw64',
    '-DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=ONLY',
    '-DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=ONLY',
    '-DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER']
  ubuntu_packages = ['mingw64-x-gcc']
  for ppa in ['tobydox/mingw-x-precise', 'ubuntu-wine/ppa']:
    check_call(['sudo', 'add-apt-repository', 'ppa:' + ppa, '-y'])

# Install dependencies.
os_name = os.environ['TRAVIS_OS_NAME']
if os_name == 'linux':
  cmake_system = 'Linux-x86_64'
else:
  cmake_system = 'Darwin-x86_64'
  # Install Java as a workaround for bug
  # http://bugs.java.com/bugdatabase/view_bug.do?bug_id=7131356.
  # java_url = 'https://support.apple.com/downloads/DL1572/en_US/javaforosx.dmg'
  # with Downloader().download(java_url) as f:
  #   bootstrap.install_dmg(f)

# Install newer version of CMake.
cmake_path = bootstrap.install_cmake(
  'cmake-3.4.0-{}.tar.gz'.format(cmake_system), check_installed=False,
  download_dir=None, install_dir='.')

env = os.environ.copy()
env['PATH'] = '/usr/lib/ccache:' + env['PATH']
check_call([cmake_path] + cmake_flags + ['.'], env=env)
check_call(['make', '-j3'], env=env)

# Install test dependencies.
if build == 'cross':
  check_call(['sudo', 'apt-get', 'install', 'wine1.7'])
  if check_output(['objdump', '-p', 'bin/libasl.dll']).find('write_sol_ASL') < 0:
    print('ASL symbols not exported')
    exit(1)

# Run tests.
check_call(['ctest', '-V'])

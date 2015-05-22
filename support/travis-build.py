#!/usr/bin/env python
# Build the project on Travis CI.

from __future__ import print_function
import os, re, shutil, tempfile
from bootstrap import bootstrap
from download import Downloader
from subprocess import call, check_call, check_output, Popen, PIPE, STDOUT

build = os.environ['BUILD']
if build == 'doc':
  returncode = 1
  travis = 'TRAVIS' in os.environ
  workdir = tempfile.mkdtemp()
  try:
    # Install dependencies.
    if travis:
      check_call(['sudo', 'apt-get', 'install', 'python-virtualenv', 'doxygen'])
    returncode, repo_dir = __import__('build-docs').build_docs(workdir)
    if returncode == 0:
      if travis:
        check_call(['git', 'config', '--global', 'user.name', 'amplbot'])
        check_call(['git', 'config', '--global', 'user.email', 'viz@ampl.com'])
      # Push docs to GitHub pages.
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
ubuntu_packages = ['gfortran', 'unixodbc-dev']

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
  # Install newer version of CMake.
  check_call(['sudo', 'apt-get', 'update'])
  check_call(['sudo', 'apt-get', 'install'] + ubuntu_packages)
  cmake_path = bootstrap.install_cmake(
    'cmake-3.1.1-Linux-x86_64.tar.gz', check_installed=False,
    download_dir=None, install_dir='.')
else:
  # Install Java as a workaround for bug
  # http://bugs.java.com/bugdatabase/view_bug.do?bug_id=7131356.
  java_url = 'http://support.apple.com/downloads/DL1572/en_US/JavaForOSX2014-001.dmg'
  with Downloader().download(java_url) as f:
    bootstrap.install_dmg(f)
  cmake_path = 'cmake'

check_call([cmake_path] + cmake_flags + ['.'])
check_call(['make', '-j3'])

# Install test dependencies.
if build == 'cross':
  check_call(['sudo', 'apt-get', 'install', 'wine1.7'])
  if check_output(['objdump', '-p', 'bin/libasl.dll']).find('write_sol_ASL') < 0:
    print('ASL symbols not exported')
    exit(1)

# Run tests.
check_call(['ctest', '-V'])

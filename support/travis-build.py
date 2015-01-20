#!/usr/bin/env python
# Set up a build environment on Travis.

from __future__ import print_function
import os, tarfile
from bootstrap import bootstrap
from contextlib import closing
from subprocess import check_call, check_output

cmake_flags = ['-DBUILD=all']
ubuntu_packages = ['gfortran', 'unixodbc-dev']

build = os.environ['BUILD']
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
  for ppa in ['tobydox/mingw', 'ubuntu-wine/ppa']:
    check_call(['sudo', 'add-apt-repository', 'ppa:' + ppa, '-y'])

# Install dependencies.
os_name = os.environ['TRAVIS_OS_NAME']
if os_name == 'linux':
  check_call(['sudo', 'apt-get', 'update'])
  check_call(['sudo', 'apt-get', 'install', 'libc6:i386'] + ubuntu_packages)
  cmake_package = 'cmake-2.8.12.2-Linux-i386.tar.gz'
else:
  cmake_package = 'cmake-2.8.12.2-Darwin-universal.tar.gz'

# Install newer version of CMake.
cmake_path = bootstrap.install_cmake(
  cmake_package, check_installed=False, download_dir=None, install_dir='.')

check_call([cmake_path] + cmake_flags + ['.'])
check_call(['make', '-j3'])

# Install test dependencies.
if build == 'cross':
  check_call(['sudo', 'apt-get', 'install', 'wine1.7'])
  if check_output(['objdump', '-p', 'bin/libasl.dll']).find('write_sol_ASL') >= 0:
    print('ASL symbols not exported')
    exit(1)

# Run tests.
check_call(['ctest', '-V'])

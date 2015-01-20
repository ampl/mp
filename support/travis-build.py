#!/usr/bin/env python
# Set up a build environment on Travis.

import os, tarfile
from contextlib import closing
from download import Downloader
from subprocess import check_call

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
  '''
      TEST_SETUP='if [ -z "$(objdump -p bin/libasl.dll | grep write_sol_ASL)" ]; then
                    echo ASL symbols not exported;
                    exit 1;
                  fi;
                  sudo apt-get install wine1.7'
  '''

os_name = os.environ['TRAVIS_OS_NAME']
if os_name == 'linux':
  check_call(['sudo', 'apt-get', 'update'])
  check_call(['sudo', 'apt-get', 'install', 'libc6:i386'] + packages)

# Install newer version of CMake.
with Downloader().download('http://www.cmake.org/files/v2.8/cmake-2.8.12.2-Linux-i386.tar.gz') as f:
  with closing(tarfile.open(f, 'r:gz')) as archive:
    archive.extractall('.')

# TODO: download and extract cmake
check_call(['cmake'] + cmake_flags + ['.'])
check_call(['make', '-j3'])
  
'''      
      TEST_SETUP=""

install:
  - eval $TEST_SETUP
  - ctest -V
'''

#!/bin/sh
# Set up build environment on Ubuntu or other Debian-based
# Linux distribution.

import platform
from bootstrap import *
from subprocess import check_call

vagrant = bootstrap_init()

# Install build tools.
check_call(['apt-get', 'update'])
packages = [
  'git-core', 'gcc', 'g++', 'gfortran' 'ccache', 'make',
  'python-pip', 'python-dev', 'default-jdk', 'unixodbc-dev'
]
if platform.machine() == 'x86_64':
  packages.append('libc6-i386')
check_call(['apt-get', 'install', '-y', packages])

install_cmake('cmake-2.8.12.2-Linux-i386.tar.gz')

# Installs symlinks for ccache.
for name in ['gcc', 'cc', 'g++', 'c++']:
  add_to_path(which(ccache), name)

install_f90cache()

if vagrant:
  buildslave_name = 'lucid64' if platform.machine() == 'x86_64' else 'lucid32'
  buildslave_dir = '/home/vagrant/slave'
  if not os.path.exists(buildslave_dir):
    install_buildbot_slave(buildslave_name, buildslave_dir)

# Copy optional dependencies.
if os.path.exists('opt'):
  for dir in glob('opt/linux' + platform.machine() + '/*'):
    dest = '/opt/' + dir
    if not os.path.exists(dest):
      shutil.copytree(dir, dest)

#!/usr/bin/env python
# Set up a build environment on Travis.

from subprocess import check_call

cmake_flags = ''
check_call(['cmake', cmake_flags, '.'])
check_call(['make', '-j3'])

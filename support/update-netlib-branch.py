#!/usr/bin/env python
# This script updates the "netlib" branch from the AMPL Solver Library
# repository at rsync ampl.com::ampl/

import shutil, tempfile
from subprocess import call, check_call

workdir = tempfile.mkdtemp()
try:
  check_call(['git', 'clone', '--depth=1', '-b', 'netlib',
              'git@github.com:ampl/ampl.git', workdir])

  # Get changes via rsync.
  check_call(['rsync', '-avzq', '--delete', '--exclude=.git',
              '--exclude=student', '--exclude=*.tgz', '--exclude=*.gz',
              'ampl.com::ampl/', workdir])

  # Fix permissions.
  check_call(['chmod', '+x', workdir + '/solvers/configurehere'])

  check_call(['git', 'add', '--all', '.'], cwd=workdir)

  # Commit changes if there are any.
  if call(['git', 'diff', '--quiet', '--staged', '--exit-code'], cwd=workdir):
    check_call(['git', 'commit', '-m', 'Upstream changes'], cwd=workdir)
finally:
  shutil.rmtree(workdir)

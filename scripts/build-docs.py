#!/usr/bin/env python
# This script builds the HTML documentation and commits it to the
# ampl.github.io repository.

import os
from subprocess import check_call
from shutil import rmtree
from glob import glob

REPO = 'ampl.github.io'

rmtree('html')
rmtree(REPO)

check_call(['cmake', '.'])
check_call(['make', 'doc'])

check_call(['git', 'clone', 'git@github.com:ampl/' + REPO + '.git'])
for entry in glob(REPO + '/*'):
  if entry == REPO + '/demo' or entry == REPO + '/ampl-book.pdf':
    continue
  if os.path.isfile(entry):
    os.remove(entry)
  else:
    shutil.rmtree(entry)

for entry in glob('html/*'):
  os.copytree(entry, REPO)
os.copytree('models', REPO)

check_call(['git', 'add', '--all'], cwd=REPO)
check_call(['git', 'commit', '-m', 'Update documentation'], cwd=REPO)

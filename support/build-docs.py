#!/usr/bin/env python
# This script builds the HTML documentation and commits it to the
# ampl.github.io repository.

import os, shutil
from subprocess import call, check_call
from glob import glob
from fileutil import rmtree_if_exists

REPO = 'ampl.github.io'

rmtree_if_exists('html')
rmtree_if_exists(REPO)

check_call(['cmake', '.'])
check_call(['make', 'doc'])

# Clone the repo and remove generated content.
check_call(['git', 'clone', 'git@github.com:ampl/' + REPO + '.git'])
for entry in glob(REPO + '/*'):
  basename = os.path.basename(entry)
  if basename == 'demo' or basename == 'ampl-book.pdf':
    continue
  if os.path.isdir(entry):
    shutil.rmtree(entry)
  else:
    os.remove(entry)

# Recursively copy an entire directory tree rooted at src.
# A directory with the same basename as src is created in dst
# which must name an existing directory.
def copytodir(src, dst):
  shutil.copytree(src, os.path.join(dst, os.path.basename(src)))

for entry in glob('html/*'):
  if os.path.isdir(entry):
    copytodir(entry, REPO)
  else:
    shutil.copy(entry, REPO)
copytodir('models', REPO)

check_call(['git', 'add', '--all'], cwd=REPO)
if call(["git", "diff-index", "--quiet", "HEAD"], cwd=REPO):
  check_call(['git', 'commit', '-m', 'Update documentation'], cwd=REPO)

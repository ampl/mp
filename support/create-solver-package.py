#!/usr/bin/env python
# Create a solver source package for AMPL use.

from __future__ import print_function
import os, zipfile

mp_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
solver = 'ilogcp'

def add_to_archive(archive, source_dir, target_dir=None, **kwargs):
  """
    Adds all files from source_dir to target_dir in archive.
    source_dir is relative to mp_dir.
  """
  if not target_dir:
    target_dir = source_dir
  exclude = kwargs.get('exclude', {})
  dirpath = os.path.join(mp_dir, source_dir)
  for entry in os.listdir(dirpath):
    path = os.path.join(dirpath, entry)
    if not (os.path.isdir(path) or entry in exclude):
      archive.write(path, os.path.join(solver, target_dir, entry))

if __name__ == '__main__':
  with zipfile.ZipFile(solver + '.zip', 'w') as archive:
    add_to_archive(archive, 'include/mp')
    add_to_archive(archive, 'src')
    add_to_archive(archive, 'src/asl', exclude={'CMakeLists.txt'})
    add_to_archive(archive, os.path.join('solvers', solver), 'src')

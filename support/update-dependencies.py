#!/usr/bin/env python
# Updates dependencies integrated into the source tree such as C++ Format.

import download, os, re, shutil, zipfile, fileutil

project_dir = os.path.join(os.path.dirname(__file__), '..')

include_files = ['format.h', 'posix.h']
src_files = ['format.cc', 'posix.cc']

def copyfile(src, dst):
  for line in src:
    # Fix includes.
    dst.write(re.sub(r'#include "(.*)"', r'#include "mp/\1"', line))

def extract(archive, filenames, dest, archive_dir, **kwargs):
  dest = os.path.join(project_dir, dest)
  if kwargs.get('clean'):
    fileutil.rmtree_if_exists(dest)
    os.mkdir(dest)
  for filename in filenames:
    with archive.open(archive_dir + filename) as src:
      with open(os.path.join(dest, filename), 'w') as dst:
        copyfile(src, dst)

d = download.Downloader()
with d.download('https://github.com/cppformat/cppformat/archive/master.zip') as f:
  with zipfile.ZipFile(f, 'r') as zf:
    root = 'cppformat-master/'
    extract(zf, include_files, 'include/mp', root)
    extract(zf, src_files, 'src', root)
    test_files = []
    archive_test_dir = root + 'test/'
    for name in zf.namelist():
      if name.startswith(archive_test_dir) and name != archive_test_dir:
        test_files.append(name[len(archive_test_dir):])
    extract(zf, test_files, 'test/format', archive_test_dir, clean=True)

#!/usr/bin/env python
# This script creates demo packages of ampl.

from __future__ import print_function
import gzip, os, shutil, stat, tarfile, urllib, zipfile
import fileutil
from glob import glob
from sets import Set
from StringIO import StringIO
from urlparse import urlparse

# URL for downloading student versions of AMPL binaries.
student_url = 'http://ampl.com/netlib/ampl/student/'

# URL for downloading a command-line version of AMPL.
amplcml_url = 'http://www.ampl.com/NEW/TABLES/amplcml.zip'

# URL for downloading AMPL table handler.
googlecode_url = 'https://ampl.googlecode.com/files/'

# Files to download.
download_files = [
  'README', 'ampl.gz', 'cplex.gz', 'gjh.gz',
  'gurobi.tgz', 'minos.gz', 'snopt.gz']

# Files that need executable permission.
executables = Set(['ampl', 'cplex', 'gjh', 'gurobix', 'minos', 'snopt'])

# Paths to files and directories to copy from amplcml.zip to the demo package.
extra_paths = ['MODELS/', 'kestrel', 'modinc']

# Writes a file object f to the file with specified name.
def writefile(f, filename):
  with open(filename, 'wb') as out:
    out.write(f.read())

# Retrieve the url or use cached version of the file if available.
cache_dir = 'cache'
def retrieve_cached(url, system = None):
  filename = os.path.basename(urlparse(url).path)
  cached_path = cache_dir
  if system is not None:
    cached_path = os.path.join(cache_dir, system)
    if not os.path.exists(cached_path):
      os.mkdir(cached_path)
  cached_path = os.path.join(cached_path, filename)
  if os.path.exists(cached_path):
    print('Using cached version of', filename)
  else:
    print('Downloading', filename)
    urllib.urlretrieve(url, cached_path)
  return cached_path

amplcml = zipfile.ZipFile(retrieve_cached(amplcml_url))

ampl_demo_dir = 'ampl-demo'

# Extract files from amplcml.zip.
def extract_amplcml(extra_paths = None):
  for name in amplcml.namelist():
    if extra_paths is not None:
      found = False
      for path in extra_paths:
        if name.startswith('amplcml/' + path):
          found = True
          break
      if not found:
        continue
    outname = name.replace('amplcml/', ampl_demo_dir + '/')
    if name.endswith('/'):
      os.makedirs(outname)
    else:
      writefile(amplcml.open(name), outname)

# Prepare a demo package for UNIX systems.
def prepare_unix_package(system):  
  os.mkdir(ampl_demo_dir)
  extract_amplcml(extra_paths)

  # Download ampl and solvers.
  for filename in download_files:
    sysdir = system
    if system == 'macosx':
      sysdir = system + '/x86_32'
    elif system == 'linux32':
      sysdir = 'linux'
    retrieved_file = retrieve_cached('{}/{}/{}'.format(student_url, sysdir, filename), system)
    # Unpack if necessary.
    outfilename = filename
    if filename.endswith('.gz'):
      outfilename = filename.replace('.gz', '')
      with gzip.GzipFile(retrieved_file) as f:
        writefile(f, os.path.join(ampl_demo_dir, outfilename))
    elif filename.endswith('.tgz'):
      with tarfile.open(retrieved_file) as tar:
        tar.extractall(ampl_demo_dir)
    else:
      shutil.copy(retrieved_file, os.path.join(ampl_demo_dir, filename))
    # Add executable permissions.
    if outfilename in executables:
      path = os.path.join(ampl_demo_dir, outfilename)
      st = os.stat(path)
      os.chmod(path, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

  # Replace libgurobi*.so link with the library ligurobi.so* because some
  # programs don't support symlinks in zip archives.
  libs = glob(os.path.join(ampl_demo_dir, 'libgurobi.so*'))
  if len(libs) > 0:
    libgurobi = libs[0]
    libgurobi_link = glob(os.path.join(ampl_demo_dir, 'libgurobi*.so'))[0]
    os.remove(libgurobi_link)
    shutil.move(libgurobi, libgurobi_link)

  # Download ampltabl.dll.
  ampltabl_url = googlecode_url + 'ampltabl-20131212-{}.zip'.format(system)
  with zipfile.ZipFile(retrieve_cached(ampltabl_url)) as zip:
    writefile(zip.open('ampltabl.dll'), os.path.join(ampl_demo_dir, 'ampltabl.dll'))

# Prepare a demo package for Windows.
def prepare_windows_package():
  fileutil.rmtree_if_exists(ampl_demo_dir)
  extract_amplcml()

# Map from system name to IDE package suffix.
sys2ide = {
  'linux32':  'linux32.tgz',
  'linux64':  'linux32.tgz',
  'macosx': 'mac64.tgz',
  'mswin':  'win32.zip'
}

for system in ['linux32', 'linux64', 'macosx', 'mswin']:
  # Prepare the command-line demo package.
  fileutil.rmtree_if_exists(ampl_demo_dir)
  if system != 'mswin':
    archive_format = 'gztar'
    prepare_unix_package(system)
  else:
    archive_format = 'zip'
    prepare_windows_package()
  basename = 'ampl-demo-' + system
  shutil.make_archive(basename, archive_format, '.', ampl_demo_dir)

  # Prepare the IDE demo package.
  amplide_demo_dir = 'amplide-demo'
  fileutil.rmtree_if_exists(amplide_demo_dir)
  amplide_url = 'http://www.ampl.com/dl/IDE/amplide.' + sys2ide[system]
  amplide = retrieve_cached(amplide_url)
  archive_open = zipfile.ZipFile if amplide_url.endswith('zip') else tarfile.open
  with archive_open(amplide) as archive:
    archive.extractall()
  shutil.move('amplide', amplide_demo_dir)
  shutil.move(ampl_demo_dir, os.path.join(amplide_demo_dir, 'ampl'))
  basename = 'amplide-demo-' + system
  shutil.make_archive(basename, archive_format, '.', amplide_demo_dir)

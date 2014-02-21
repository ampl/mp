#!/usr/bin/env python
# This script creates demo packages of ampl.

from __future__ import print_function
import gzip, os, shutil, stat, tarfile, urllib, zipfile
import fileutil
from glob import glob
from sets import Set
from StringIO import StringIO

# URL for downloading student versions of AMPL binaries.
student_url = 'http://ampl.com/netlib/ampl/student/'

amplcml_filename = 'amplcml.zip'

# URL for downloading a command-line version of AMPL.
amplcml_url = 'http://www.ampl.com/NEW/TABLES/' + amplcml_filename

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

# Download amplcml.zip or use cached version if available.
if os.path.exists(amplcml_filename):
  amplcml = zipfile.ZipFile(amplcml_filename)
else:
  print('Downloading', amplcml_filename)
  response = urllib.urlopen(amplcml_url).read()
  writefile(StringIO(response), amplcml_filename)
  amplcml = zipfile.ZipFile(StringIO(response))

for system in ['linux', 'macosx']:
  dirname = 'ampl-demo'
  fileutil.rmtree_if_exists(dirname)
  os.mkdir(dirname)

  # Extract files from amplcml.zip.
  for name in amplcml.namelist():
    found = False
    for path in extra_paths:
      if name.startswith('amplcml/' + path):
        found = True
        break
    if not found:
      continue
    outname = name.replace('amplcml/', dirname + '/')
    if name.endswith('/'):
      os.makedirs(outname)
    else:
      writefile(amplcml.open(name), outname)

  # Download ampl and solvers.
  for filename in download_files:
    print('Downloading', filename)
    fullsys = system if system != 'macosx' else system + '/x86_32'
    response = urllib.urlopen('{}/{}/{}'.format(student_url, fullsys, filename))
    # Unpack if necessary.
    outfilename = filename
    if filename.endswith('.gz'):
      outfilename = filename.replace('.gz', '')
      with gzip.GzipFile(fileobj=StringIO(response.read())) as f:
        writefile(f, os.path.join(dirname, outfilename))
    elif filename.endswith('.tgz'):
      with tarfile.open(fileobj=StringIO(response.read())) as tar:
        tar.extractall(dirname)
    else:
      writefile(response, os.path.join(dirname, filename))
    # Add executable permissions.
    if outfilename in executables:
      path = os.path.join(dirname, outfilename)
      st = os.stat(path)
      os.chmod(path, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

  # Replace libgurobi*.so link with the library ligurobi.so* because some
  # programs don't support symlinks in zip archives.
  libgurobi_link = glob(os.path.join(dirname, 'libgurobi*.so'))[0]
  libgurobi = glob(os.path.join(dirname, 'libgurobi.so*'))[0]
  os.remove(libgurobi_link)
  shutil.move(libgurobi, libgurobi_link)

  # Download ampltabl.dll.
  suffix = system + '32' if system != 'macosx' else system
  ampltabl_url = googlecode_url + 'ampltabl-20131212-{}.zip'.format(suffix)
  print('Downloading', ampltabl_url)
  response = urllib.urlopen(ampltabl_url)
  with zipfile.ZipFile(StringIO(response.read())) as zip:
    writefile(zip.open('ampltabl.dll'), os.path.join(dirname, 'ampltabl.dll'))

  # Create an archive.
  shutil.make_archive('ampl-demo-' + system, 'zip', '.', dirname)

# TODO: windows
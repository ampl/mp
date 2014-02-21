#!/usr/bin/env python
# This script creates demo packages of ampl.

from __future__ import print_function
import gzip, os, stat, tarfile, urllib, zipfile
from fileutil import rmtree_if_exists
from sets import Set
from StringIO import StringIO

# URL for downloading student versions of AMPL binaries.
student_url = 'http://ampl.com/netlib/ampl/student/linux'

amplcml_filename = 'amplcml.zip'

# URL for downloading a command-line version of AMPL.
amplcml_url = 'http://www.ampl.com/NEW/TABLES/' + amplcml_filename

# URL for downloading AMPL table handler.
ampltabl_url = 'https://ampl.googlecode.com/files/ampltabl-20131212-linux32.zip'

# Files to download.
files = [
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

dirname = 'ampl-demo'
rmtree_if_exists(dirname)
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
for filename in files:
  print('Downloading', filename)
  response = urllib.urlopen(student_url + '/' + filename)
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

# Download ampltabl.dll.
print('Downloading', ampltabl_url)
response = urllib.urlopen(ampltabl_url)
with zipfile.ZipFile(StringIO(response.read())) as zip:
  writefile(zip.open('ampltabl.dll'), os.path.join(dirname, 'ampltabl.dll'))

# Create an archive.
with zipfile.ZipFile('ampl-demo-linux.zip', 'w') as zip:
  for root, dirs, files in os.walk(dirname):
    for file in files:
      zip.write(os.path.join(root, file))

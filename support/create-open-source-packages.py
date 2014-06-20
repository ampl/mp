#!/usr/bin/env python

"""Create packages of open-source AMPL solvers and libraries.

Usage:
  create-open-source-packages.py [upload | simulate]
  
When run in "simulate" mode, this script doesn't uploading anything.
"""

from __future__ import print_function
import datetime, os, re, shutil, sys, tempfile, zipfile
from docopt import docopt
from subprocess import check_call

project = "ampl"

class Package:
  def __init__(self, name, files, **args):
    self.name = name
    self._files = files
    self.license = args.get('project', name) + '-license.txt'
    self._winfiles = args.get('winfiles', [])

  def getfiles(self, system, versions):
    """Get the list of files for the given system."""
    version = versions[self.name].split('-')[0]
    files = self._files + (self._winfiles if system.startswith('win') else [])
    for i in range(len(files)):
      files[i] = files[i].format(version=version)
      if system.startswith('win') and os.path.splitext(files[i])[1] == '':
        files[i] += '.exe'
    return files

packages = [
  Package('amplgsl', ['amplgsl.dll', 'gsl.ampl'], project='gsl'),
  Package('bonmin',  ['bonmin'], project='coin', winfiles=['libipoptfort.dll']),
  Package('cbc',     ['cbc'], project='coin'),
  Package('gecode',  ['gecode', 'gecode.ampl']),
  Package('ipopt',   ['ipopt'], project='coin', winfiles=['libipoptfort.dll']),
  Package('jacop',   ['jacop', 'ampljacop.jar', 'jacop-{version}.jar'])
]

def create_packages(system, workdir, simulate):
  """Create packages and upload them to the server."""

  # Download build artifacts.
  artifact_dir = os.path.join(workdir, 'artifacts')
  print("Downloading files for {}:".format(system))
  cmd = "scp -r ampl.com:/var/lib/buildbot/upload/{} {}".format(system, artifact_dir)
  check_call(cmd, shell=True)

  # Read versions.
  versions = {}
  for filename in ['versions.txt', 'coin-versions.txt']:
    with open(os.path.join(artifact_dir, filename)) as f:
      for line in f:
        items = line.rstrip().split(' ')
        if len(items) < 2:
          continue
        version = items[1].rstrip(',')
        m = re.match(r'.*(driver|library)\(([0-9]+)\)', line)
        if m:
          version += '-' + m.group(2)
        versions[items[0].lower()] = version

  # Create individual packages.
  for package in packages:
    archive_name = '{}-{}.zip'.format(package.name, system)
    with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED) as zip:
      for f in package.getfiles(system, versions):
        if not f:
          continue
        zip.write(os.path.join(artifact_dir, f), f)
      zip.write(os.path.join('licenses', package.license), package.license)

  #date = datetime.datetime.today()
  #date = "{}{:02}{:02}".format(date.year, date.month, date.day)
  # TODO: create full package
  shutil.rmtree(artifact_dir)

  return
  # Upload all in one archive.
  archive_name = os.path.join(tempdir, "ampl-open-{}-{}.zip".format(date, platform))
  with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED) as zip:
    for path in paths:
      zip.write(path, path[dirlen:])
    zip.write("LICENSE", "LICENSE")
  upload(archive_name, "Open-source AMPL solvers and libraries", simulate)

if __name__ == '__main__':
  args = docopt(__doc__)
  simulate = args['simulate']
  if not args['upload'] and not simulate:
    print(__doc__)
    sys.exit()
  workdir = tempfile.mkdtemp()
  try:
    for system in ['linux32', 'linux64', 'osx', 'win32', 'win64']:
      create_packages(system, workdir, simulate)
  finally:
    shutil.rmtree(workdir)

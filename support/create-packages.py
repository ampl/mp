#!/usr/bin/env python

"""Create packages of AMPL solvers and libraries.

Usage:
  create-packages.py [update]
"""

from __future__ import print_function
import docopt, fileutil, os, re, shutil, subprocess, tempfile, zipfile

project = "ampl"

class Package:
  def __init__(self, name, files, **args):
    self.name = name
    self._files = files
    self.is_open = args.get('is_open', True)
    self.version_name = args.get('version_name', name)
    self.license = args.get('project', name) + '-license.txt'
    self._winfiles = args.get('winfiles', [])

  def getfiles(self, system, versions):
    """Get the list of files for the given system."""
    version = versions[self.version_name].split('-')[0]
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
  Package('couenne', ['couenne'], project='coin', winfiles=['libipoptfort.dll']),
  Package('gecode',  ['gecode', 'gecode.ampl']),
  Package('ipopt',   ['ipopt'], project='coin', winfiles=['libipoptfort.dll']),
  Package('jacop',   ['jacop', 'ampljacop.jar', 'jacop-{version}.jar']),
  Package('path',    ['path'], is_open=False, version_name='ampl/path')
]

def get_archive_name(system, package=None):
  if package:
    return '{}-{}.zip'.format(package.name, system)
  return 'ampl-open-{}.zip'.format(system)

def create_packages(system, workdir):
  """Create packages and upload them to the server."""

  # Download build artifacts.
  artifact_dir = os.path.join(workdir, 'artifacts')
  print("Downloading files for {}:".format(system))
  cmd = "scp -r ampl.com:/var/lib/buildbot/upload/{} {}".format(system, artifact_dir)
  subprocess.check_call(cmd, shell=True)

  # Read versions.
  versions = {}
  for filename in ['versions.txt', 'coin-versions.txt']:
    with open(os.path.join(artifact_dir, filename)) as f:
      for line in f:
        m = re.match(r'([^ ]+) ([^\s]+)(.*)', line)
        if not m:
          continue
        name, version, rest = m.groups()
        name = name.lower()
        if name in versions:
          continue
        m = re.match(r'.*(driver|library)\(([0-9]+)\)', rest)
        if m:
          version += '-' + m.group(2)
        versions[name] = version
      print(versions)

  # Create individual packages.
  paths = set()
  for package in packages:
    archive_name = get_archive_name(system, package)
    with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED) as zip:
      for f in package.getfiles(system, versions):
        path = os.path.join(artifact_dir, f)
        zip.write(path, f)
        if package.is_open:
          paths.add(path)
      path = os.path.join('licenses', package.license)
      zip.write(path, package.license)
      if package.is_open:
        paths.add(path)

  # Careate a full package.
  with zipfile.ZipFile(get_archive_name(system), 'w', zipfile.ZIP_DEFLATED) as zip:
    for path in paths:
      zip.write(path, os.path.join('ampl-open', os.path.basename(path)))
  shutil.rmtree(artifact_dir)

systems = ['linux32', 'linux64', 'osx', 'win32', 'win64']

if __name__ == '__main__':
  args = docopt.docopt(__doc__)
  workdir = tempfile.mkdtemp()
  try:
    for system in systems:
      create_packages(system, workdir)
    if args['update']:
      target_dir = '/var/www/dl/open'
      for system in systems:
        for package in packages:
          dest = os.path.join(target_dir, package.name)
          if not os.path.exists(dest):
            os.mkdir(dest)
          fileutil.move(get_archive_name(system, package), dest)
        fileutil.move(get_archive_name(system), target_dir)
  finally:
    shutil.rmtree(workdir)

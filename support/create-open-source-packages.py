#!/usr/bin/env python

"""Create packages of open-source AMPL solvers and libraries.

Usage:
  create-open-source-packages.py [upload | simulate]
  
When run in "simulate" mode, this script only simulates upload without
actually uploading anything.
"""

import datetime, os, re, shutil, sys, tempfile, zipfile
from docopt import docopt
from subprocess import call, check_call

server = "ampl.com"
project = "ampl"

summaries = {
  "amplgsl" : "AMPL bindings for the GNU Scientific Library",
  "ampltabl": "ODBC table handler",
  "cbc"     : "COIN-OR Cbc solver",
  "gecode"  : "Gecode solver",
  "ipopt"   : "COIN-OR Ipopt solver",
  "jacop"   : "JaCoP solver"
}

# A mapping from the package name to the list of files it contains.
package_files = {
  'amplgsl': ['amplgsl.dll', 'gsl.ampl'],
  'cbc'    : ['bonmin'],
  'cbc'    : ['cbc'],
  'gecode' : ['gecode', 'gecode.ampl'],
  'ipopt'  : ['ipopt', 'win:libipoptfort.dll'],
  'jacop'  : ['jacop', 'ampljacop.jar', 'jacop-${version}.jar'],
}

def rmtree_if_exists(path):
  "Delete an entire directory tree if it exists."
  if os.path.exists(path):
    shutil.rmtree(path)

def upload(filename, summary, simulate):
  "Upload a file to the server."
  print("Uploading {}".format(filename))
  return
  # TODO: implement
  from netrc import netrc
  authenticators = netrc().authenticators("code.google.com")
  username = authenticators[0]
  password = authenticators[2]
  if simulate:
    return

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

  date = datetime.datetime.today()
  date = "{}{:02}{:02}".format(date.year, date.month, date.day)

  # Create individual packages.
  for package, files in package_files.iteritems():
    version = versions[package].split('-')[0]
    print package, version
    # TODO
  
  return
  paths = []
  for base, dirs, files in os.walk(artifact_dir):
    for file in files:
      path = os.path.join(base, file)
      name = path[dirlen:]
      if name == "versions" or name in file_package or 'ampl-lic' in name:
        continue
      basename = os.path.splitext(name)[0]
      suffix = versions.get(basename, date)
      archive_name = os.path.join(tempdir, "{}-{}-{}.zip".format(basename, suffix, platform))
      with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED) as zip:
        zip.write(path, name)
        files = extra_files.get(basename)
        if files:
          for f in files:
            items = f.split(':')
            filename = items[0]
            if len(items) > 1 and not platform.startswith(items[1]):
              continue
            extra_path = os.path.join(base, filename)
            zip.write(extra_path, filename)
            paths.append(extra_path)
      upload(archive_name, summaries[basename], simulate)
      paths.append(path)
  # Upload all in one archive.
  archive_name = os.path.join(tempdir, "ampl-open-{}-{}.zip".format(date, platform))
  with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED) as zip:
    for path in paths:
      zip.write(path, path[dirlen:])
    zip.write("LICENSE", "LICENSE")
  upload(archive_name, "Open-source AMPL solvers and libraries", simulate)
  shutil.rmtree(filedir)

def update_redir_page(repo_path, name, versions):
  filename = os.path.join(repo_path, name + ".html")
  with open(filename, 'r') as f:
    content = f.read()
  version = re.sub(r".*-()", r"\1", versions.get(name))
  print 'changing version: {}, {}'.format(versions.get(name), version)
  content = re.sub(r"q={}\+(\d+)".format(name), "q={}+{}".format(name, version), content)
  with open(filename, 'w') as f:
    f.write(content)

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

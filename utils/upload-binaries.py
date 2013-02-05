#!/usr/bin/env python

import datetime, os, shutil, sys, zipfile
from googlecode_upload import upload_find_auth

project = "ampl"

summaries = {
  "amplgsl.dll": "AMPL bindings for the GNU Scientific Library",
  "ampltabl.dll": "ODBC table handler",
  "gecode": "Gecode solver",
  "gecode.exe": "Gecode solver",
  "cbc": "COIN-OR CBC solver",
  "cbc.exe": "COIN-OR CBC solver"
}

def upload(filename, summary):
  print("Uploading {}".format(filename))
  status, reason, url = upload_find_auth(
    filename, project, summary, None, None, None)
  os.remove(filename)
  if not url:
    print('Google Code upload error: {} ({})'.format(reason, status))
    sys.exit(1)

dir = "ampl-open"
for platform in reversed(["linux32", "linux64", "macosx", "win32", "win64"]):
  if os.path.exists(dir):
    shutil.rmtree(dir)
  print("Downloading binaries for {}:".format(platform))
  os.system("scp -r callisto.local:/var/lib/buildbot/upload/{} {}".format(platform, dir))
  date = datetime.datetime.today()
  date = "{}{:02}{:02}".format(date.year, date.month, date.day)
  dirlen = len(dir) + 1
  # Upload individual files.
  paths = []
  for base, dirs, files in os.walk(dir):
    for file in files:
      path = os.path.join(base, file)
      name = path[dirlen:]
      archive_name = "{}-{}-{}.zip".format(name, date, platform)
      with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED) as zip:
        zip.write(path, name)
      upload(archive_name, summaries[name])
      paths.append(path)
  # Upload all in one archive.
  archive_name = "ampl-open-{}-{}.zip".format(date, platform)
  with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED) as zip:
    for path in paths:
      path = os.path.join(base, file)
      zip.write(path, path[dirlen:])
    zip.write("LICENSE", "LICENSE")
  upload(archive_name, "Open-source AMPL solvers and libraries")

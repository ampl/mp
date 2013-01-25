#!/usr/bin/env python

import datetime, os, shutil, sys, zipfile
from googlecode_upload import upload_find_auth

project = "ampl"
summary = "Open-source AMPL solvers and libraries"

dir = "ampl-open"
for platform in reversed(["linux32", "linux64", "macosx", "win32", "win64"]):
  if os.path.exists(dir):
    shutil.rmtree(dir)
  print("Downloading binaries for {}:".format(platform))
  os.system("scp -r callisto.local:/var/lib/buildbot/upload/{} {}".format(platform, dir))
  date = datetime.datetime.today()
  date = "{}{:02}{:02}".format(date.year, date.month, date.day)
  file_path = "ampl-open-{}-{}.zip".format(date, platform)
  with zipfile.ZipFile(file_path, 'w', zipfile.ZIP_DEFLATED) as zip:
    dirlen = len(dir) + 1
    for base, dirs, files in os.walk(dir):
      for file in files:
        fn = os.path.join(base, file)
        zip.write(fn, fn[dirlen:])
  print("Uploading {}".format(file_path))
  status, reason, url = upload_find_auth(
    file_path, project, summary, None, None, None)
  os.remove(file_path)
  if not url:
    print('Google Code upload error: {} ({})'.format(reason, status))
    sys.exit(1)

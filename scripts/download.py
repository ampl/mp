# Download with automatic cleanup.

from __future__ import print_function
import os, shutil, sys, urllib2, urlparse

class TempFile:
  def __init__(self, filename):
    self.filename = filename
  def __enter__(self):
    return self.filename
  def __exit__(self, type, value, traceback):
    os.remove(self.filename)

# Downloads into a temporary file.
def download(url, cookie = None):
  filename = os.path.basename(urlparse.urlsplit(url)[2])
  print('Downloading', url, 'to', filename)
  sys.stdout.flush()
  opener = urllib2.build_opener()
  if cookie:
    opener.addheaders.append(('Cookie', cookie))
  with open(filename, 'wb') as f:
    shutil.copyfileobj(opener.open(url), f)
  return TempFile(filename)

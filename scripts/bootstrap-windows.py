# Set up build environment on Windows.

from __future__ import print_function
import os, sys, shutil, urllib
from zipfile import ZipFile
from subprocess import check_call

def download(url, filename):
  print('Downloading', url, 'to', filename)
  sys.stdout.flush()
  urllib.urlretrieve(url, filename)

def unzip(filename, path):
  with ZipFile(filename) as zip:
    zip.extractall(path)

# Install CMake.
cmake = 'cmake-2.8.12.2-win32-x86'
cmake_install_dir = 'C:\\Program Files\\' + cmake
if not os.path.exists(cmake_install_dir):
  filename = "cmake.zip"
  download('http://www.cmake.org/files/v2.8/' + cmake + '.zip', filename)
  unzip(filename, 'C:\\Program Files')
  os.remove(filename)

# Add Python and CMake to PATH.
check_call(['setx', 'PATH',
  os.getenv('PATH') + ';C:\\Python27;' + cmake_install_dir + '\\bin', '/m'])

# Install .NET Framework 4 for msbuild.
# This requires vagrant-windows plugin version 1.7.0.pre.2 or later.
# See https://github.com/WinRb/vagrant-windows/pull/189
if not os.path.exists('\\Windows\\Microsoft.NET\\Framework64\\v4.0.30319'):
  filename = "dotnetfx40.exe"
  download(
    'http://download.microsoft.com/download/9/5/A/' +
    '95A9616B-7A37-4AF6-BC36-D6EA96C8DAAE/dotNetFx40_Full_x86_x64.exe',
    filename)
  check_call([filename, '/q', '/norestart'])
  os.remove(filename)

# Install 7zip.
sevenzip = 'C:\\Program Files (x86)\\7-Zip\\7z.exe'
if not os.path.exists(sevenzip):
  filename = '7z.exe'
  download('http://downloads.sourceforge.net/sevenzip/7z920.exe', filename)
  check_call([filename, '/S'])
  os.remove(filename)

# Install Windows SDK.
if not os.path.exists('\\Program Files\\Microsoft SDKs\\Windows\\v7.1'):
  # Extract ISO.
  filename = 'winsdk.iso'
  download(
     'http://download.microsoft.com/download/F/1/0/'
     'F10113F5-B750-4969-A255-274341AC6BCE/GRMSDKX_EN_DVD.iso',
     filename)
  check_call([sevenzip, 'x', '-tudf', '-owinsdk', filename])
  os.remove(filename)
  # Install SDK.
  check_call(['winsdk\\setup.exe', '-q'])
  shutil.rmtree('winsdk')

# Install MinGW.
def install_mingw(arch):
  bits = '64' if arch.endswith('64') else '32'
  if os.path.exists('\\mingw' + bits):
    return
  filename = 'mingw.7z'
  download(
    'http://sourceforge.net/projects/mingw-w64/files/' +
    'Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/4.8.2/' +
    'threads-win32/sjlj/' + arch +
    '-4.8.2-release-win32-sjlj-rt_v3-rev4.7z/download', filename)
  check_call([sevenzip, '-o', 'C:\\', 'x', filename])

install_mingw('i686')
install_mingw('x86_64')

# TODO: install buildbot

# rem Install 64-bit JDK.
# if not exist "\Program Files\Java\jdk1.7.0_55" (
#   if exist opt\win64\jdk-7u55-windows-x64.exe (
#     opt\win64\jdk-7u55-windows-x64.exe /s
#   )
# )
#
# rem Install 32-bit JDK.
# if not exist "\Program Files (x86)\Java\jdk1.7.0_55" (
#   if exist opt\win32\jdk-7u55-windows-i586.exe (
#     opt\win32\jdk-7u55-windows-i586.exe /s
#   )
# )
#
# rem Install other dependencies using xcopy.
# if not exist "\Program Files\CMake" (
#   if exist opt\win64\root (
#     xcopy opt\win64\root C:\ /s /e
#   )
# )

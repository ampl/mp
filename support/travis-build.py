#!/usr/bin/env python
# Build the project on Travis CI.

from __future__ import print_function
import mmap, os, re, shutil, tempfile, zipfile
from bootstrap import bootstrap
from download import Downloader
from subprocess import call, check_call, check_output, Popen, PIPE, STDOUT

def extract_docs(output_dir):
  "Extract the AMPLGSL documentation from the code."
  output = None
  dir = os.path.dirname(__file__)
  output_dir = os.path.join(output_dir, 'amplgsl')
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  with open(os.path.join(dir, '../src/gsl/amplgsl.c'), 'r+b') as input:
    map = mmap.mmap(input.fileno(), 0)
    for i in re.finditer(r'/\*\*(.*?)\*/', map, re.DOTALL):
      s = re.sub(r'\n +\* ?', r'\n', i.group(1))
      s = re.sub(r'\$(.+?)\$', r':math:`\1`', s, flags=re.DOTALL)
      m = re.search(r'@file (.*)', s)
      if m:
        filename = m.group(1)
        if output:
          output.close()
        output = open(os.path.join(output_dir, filename + '.rst'), 'w')
        s = s[:m.start()] + s[m.end():]
      output.write(s.rstrip(' '))
    map.close()

def get_mp_version():
  filename = os.path.join(os.path.dirname(__file__), '..', 'CMakeLists.txt')
  with open(filename, 'r') as f:
    for line in f:
      m = re.search(r'MP_VERSION (\d+\.\d+\.\d+)', line)
      if m:
        return m.group(1)

def pip_install(package, commit):
  "Install package from GitHub using pip."
  check_call(['pip', 'install',
              'git+git://github.com/{0}.git@{1}'.format(package, commit)])

def build_docs(workdir, travis):
  # Install dependencies.
  if travis:
    check_call(['sudo', 'apt-get', 'install', 'python-virtualenv', 'doxygen'])
  # Create virtualenv.
  check_call(['virtualenv', 'venv'])
  activate_this_file = 'venv/bin/activate_this.py'
  execfile(activate_this_file, dict(__file__=activate_this_file))
  # Install Sphinx and Breathe.
  pip_install('sphinx-doc/sphinx', 'a1a80ab509fbf01aa459e0ec5a5c9b66f011ee47')
  pip_install('michaeljones/breathe', '173c3d7f523f5d69af2b0b32543c90dd943340df')
  # Clone the ampl.github.io repo.
  repo = 'ampl.github.io'
  check_call(['git', 'clone', 'https://github.com/ampl/{}.git'.format(repo)],
             cwd=workdir)
  repo_dir = os.path.join(workdir, repo)
  # Copy API docs and the database connection guides to the build directory.
  # The guides are not stored in the mp repo to avoid polluting history with
  # image blobs.
  build_dir = os.path.join(workdir, 'build')
  shutil.copytree(os.path.join(repo_dir, 'src'), build_dir)
  for entry in os.listdir('doc'):
    src = os.path.join('doc', entry)
    dst = os.path.join(build_dir, entry)
    if os.path.isdir(src):
      shutil.copytree(src, dst)
    else:
      shutil.copyfile(src, dst)
  # Remove generated content.
  keep = {
    '.git', 'demo', 'models', 'src', 'ampl-book.pdf', 'nlwrite.pdf',
    '.gitignore', '.nojekyll'
  }
  for entry in os.listdir(repo_dir):
    if entry in keep:
      continue
    path = os.path.join(repo_dir, entry)
    if os.path.isdir(path):
      shutil.rmtree(path)
    else:
      os.remove(path)
  # Build docs.
  extract_docs(build_dir)
  p = Popen(['doxygen', '-'], stdin=PIPE, cwd=build_dir)
  p.communicate(input=r'''
      PROJECT_NAME      = MP
      INPUT             = {0}/include/mp/common.h \
                          {0}/include/mp/nl.h
      EXCLUDE_SYMBOLS   = mp::internal::*
      GENERATE_LATEX    = NO
      GENERATE_MAN      = NO
      GENERATE_RTF      = NO
      GENERATE_HTML     = NO
      GENERATE_XML      = YES
      XML_OUTPUT        = doxyxml
      QUIET             = YES
      JAVADOC_AUTOBRIEF = YES
      AUTOLINK_SUPPORT  = NO
      ALIASES           = "rst=\verbatim embed:rst"
      ALIASES          += "endrst=\endverbatim"
    '''.format(os.path.abspath('.')))
  if p.returncode != 0:
    return p.returncode
  check_call(['sphinx-build', '-D', 'version=' + get_mp_version(),
              '-b', 'html', build_dir, repo_dir])
  # Push docs to GitHub pages.
  if travis:
    check_call(['git', 'config', '--global', 'user.name', 'amplbot'])
    check_call(['git', 'config', '--global', 'user.email', 'viz@ampl.com'])
  check_call(['git', 'add', '--all'], cwd=repo_dir)
  returncode = 0
  if call(['git', 'diff-index', '--quiet', 'HEAD'], cwd=repo_dir):
    check_call(['git', 'commit', '-m', 'Update documentation'], cwd=repo_dir)
    cmd = 'git push https://$KEY@github.com/ampl/{}.git master'.format(repo)
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, cwd=repo_dir)
    # Remove URL from output because it may contain a token.
    print(re.sub(r'https:.*\.git', '<url>', p.communicate()[0]))
    returncode = p.returncode
  return returncode

build = os.environ['BUILD']
if build == 'doc':
  returncode = 1
  travis = 'TRAVIS' in os.environ
  workdir = tempfile.mkdtemp()
  try:
    returncode = build_docs(workdir, travis=travis)
  finally:
    # Don't remove workdir on Travis because the VM is discarded anyway.
    if not travis:
      shutil.rmtree(workdir)
  exit(returncode)

cmake_flags = ['-DBUILD=all']
ubuntu_packages = ['gfortran', 'unixodbc-dev']

if build == 'cross':
  cmake_flags = [
    '-DCMAKE_SYSTEM_NAME=Windows',
    '-DCMAKE_SYSTEM_PROCESSOR=x86_64',
    '-DBUILD_SHARED_LIBS:BOOL=ON',
    '-DCMAKE_C_COMPILER=/opt/mingw64/bin/x86_64-w64-mingw32-gcc',
    '-DCMAKE_CXX_COMPILER=/opt/mingw64/bin/x86_64-w64-mingw32-g++',
    '-DCMAKE_RC_COMPILER=/opt/mingw64/bin/x86_64-w64-mingw32-windres',
    '-DCMAKE_FIND_ROOT_PATH=/opt/mingw64',
    '-DCMAKE_FIND_ROOT_PATH_MODE_LIBRARY=ONLY',
    '-DCMAKE_FIND_ROOT_PATH_MODE_INCLUDE=ONLY',
    '-DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER']
  ubuntu_packages = ['mingw64-x-gcc']
  for ppa in ['tobydox/mingw-x-precise', 'ubuntu-wine/ppa']:
    check_call(['sudo', 'add-apt-repository', 'ppa:' + ppa, '-y'])

# Install dependencies.
os_name = os.environ['TRAVIS_OS_NAME']
if os_name == 'linux':
  # Install newer version of CMake.
  check_call(['sudo', 'apt-get', 'update'])
  check_call(['sudo', 'apt-get', 'install'] + ubuntu_packages)
  cmake_path = bootstrap.install_cmake(
    'cmake-3.1.1-Linux-x86_64.tar.gz', check_installed=False,
    download_dir=None, install_dir='.')
else:
  # Install Java as a workaround for bug
  # http://bugs.java.com/bugdatabase/view_bug.do?bug_id=7131356.
  java_url = 'http://support.apple.com/downloads/DL1572/en_US/JavaForOSX2014-001.dmg'
  with Downloader().download(java_url) as f:
    bootstrap.install_dmg(f)
  cmake_path = 'cmake'

check_call([cmake_path] + cmake_flags + ['.'])
check_call(['make', '-j3'])

# Install test dependencies.
if build == 'cross':
  check_call(['sudo', 'apt-get', 'install', 'wine1.7'])
  if check_output(['objdump', '-p', 'bin/libasl.dll']).find('write_sol_ASL') < 0:
    print('ASL symbols not exported')
    exit(1)

# Run tests.
check_call(['ctest', '-V'])

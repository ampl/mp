#!/usr/bin/env python
"""Build documentation.

Usage: build-docs.py [extract-docs <file>]
"""

from __future__ import print_function
import fileutil, mmap, os, re, shutil
from subprocess import call, check_call, check_output, Popen, PIPE
from docopt import docopt

mp_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

def run(*args, **kwargs):
  check_call(args, **kwargs)

def create_virtualenv(venv_dir):
  "Create and activate virtualenv in the given directory."
  # File "check" is used to make sure that we don't have
  # a partial environment in case virtualenv was interrupted.
  check_path = os.path.join(venv_dir, 'check')
  if not os.path.exists(check_path):
    run('virtualenv', venv_dir)
    os.mknod(check_path)
  # Activate virtualenv.
  import sysconfig
  scripts_dir = os.path.basename(sysconfig.get_path('scripts'))
  activate_this_file = os.path.join(venv_dir, scripts_dir, 'activate_this.py')
  with open(activate_this_file) as f:
    exec(f.read(), dict(__file__=activate_this_file))

def extract_docs(filename, output_dir):
  "Extract the AMPLGSL documentation from the code."
  output = None
  output_dir = os.path.join(output_dir, 'amplgsl')
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  with open(filename, 'r+b') as input:
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

def pip_install(package, **kwargs):
  "Install package using pip."
  commit = kwargs.get('commit')
  if commit:
    package = 'git+git://github.com/{0}.git@{1}'.format(package, commit)
  run('pip', 'install', '-q', package)

def copy_content(src_dir, dst_dir):
  "Copy content of the src_dir to dst_dir recursively."
  for entry in os.listdir(src_dir):
    src = os.path.join(src_dir, entry)
    dst = os.path.join(dst_dir, entry)
    if os.path.isdir(src):
      fileutil.rmtree_if_exists(dst)
      shutil.copytree(src, dst)
    else:
      shutil.copyfile(src, dst)

def build_docs(workdir, doxygen='doxygen'):
  create_virtualenv(os.path.join(workdir, 'build', 'virtualenv'))
  # Install Sphinx and Breathe.
  pip_install('sphinx==1.3.1')
  pip_install('michaeljones/breathe',
              commit='07b6e501fbe71917ec7919982ee9e5b71d318e38')

  # Clone the ampl.github.io repo.
  repo = 'ampl.github.io'
  repo_dir = os.path.join(workdir, repo)
  if not os.path.exists(repo_dir):
    run('git', 'clone', 'https://github.com/ampl/{}.git'.format(repo), cwd=workdir)

  # Copy API docs and the database connection guides to the build directory.
  # The guides are not stored in the mp repo to avoid polluting history with
  # image blobs.
  build_dir = os.path.join(workdir, 'build')
  copy_content(os.path.join(repo_dir, 'src'), build_dir)
  doc_dir = os.path.join(mp_dir, 'doc')
  copy_content(doc_dir, build_dir)
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
  dir = os.path.dirname(__file__)
  extract_docs(os.path.join(dir, '../src/gsl/amplgsl.cc'), build_dir)
  p = Popen([doxygen, '-'], stdin=PIPE, cwd=build_dir)
  p.communicate(input=r'''
      PROJECT_NAME      = MP
      INPUT             = {0}/include/mp/common.h \
                          {0}/include/mp/nl-reader.h \
                          {0}/include/mp/problem.h
      EXCLUDE_SYMBOLS   = mp::internal::*
      GENERATE_LATEX    = NO
      GENERATE_MAN      = NO
      GENERATE_RTF      = NO
      GENERATE_HTML     = NO
      GENERATE_XML      = YES
      XML_OUTPUT        = doxyxml
      QUIET             = YES
      JAVADOC_AUTOBRIEF = YES
      ALIASES           = "rst=\verbatim embed:rst"
      ALIASES          += "endrst=\endverbatim"
    '''.format(mp_dir))
  returncode = p.returncode
  if returncode == 0:
    # Pass the MP version via environment variables rather than command-line
    # -D version=value option because the latter doesn't work when version
    # is used in conf.py.
    env = os.environ.copy()
    env['MP_VERSION'] = get_mp_version()
    run('sphinx-build', '-b', 'html', build_dir, repo_dir, env=env)
  return (returncode, repo_dir)

if __name__ == '__main__':
  args = docopt(__doc__)
  if args['extract-docs']:
    extract_docs(args['<file>'], '.')
  else:
    build_docs('.')

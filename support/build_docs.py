# Build documentation

from __future__ import print_function
import mmap, os, re, shutil
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

def pip_install(package, **kwargs):
  "Install package using pip."
  commit = kwargs.get('commit')
  if commit:
    package = 'git+git://github.com/{0}.git@{1}'.format(package, commit)
  check_call(['pip', 'install', '-q', package])

def build_docs(workdir, travis):
  # Install dependencies.
  if travis:
    check_call(['sudo', 'apt-get', 'install', 'python-virtualenv', 'doxygen'])
  # Create virtualenv.
  virtualenv_dir = os.path.join(workdir, 'virtualenv')
  check_call(['virtualenv', virtualenv_dir])
  activate_this_file = os.path.join(virtualenv_dir, 'bin', 'activate_this.py')
  execfile(activate_this_file, dict(__file__=activate_this_file))
  # Install Sphinx and Breathe.
  pip_install('sphinx==1.3.1')
  pip_install('michaeljones/breathe', commit='18bd461b4e29dde0adf5df4b3da7e5473e2c2983')
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
      EXCLUDE_SYMBOLS   = mp::internal::*
      ALIASES           = "rst=\verbatim embed:rst"
      ALIASES          += "endrst=\endverbatim"
    '''.format(os.path.abspath('.')))
  if p.returncode != 0:
    return p.returncode
  # Pass the MP version via environment variables rather than command-line
  # -D version=value option because the latter doesn't work when version
  # is used in conf.py.
  env = os.environ.copy()
  env['MP_VERSION'] = get_mp_version()
  check_call(['sphinx-build', '-b', 'html', build_dir, repo_dir], env=env)
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

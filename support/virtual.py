# Virtualenv support

import os, subprocess, sysconfig
from distutils.version import LooseVersion

def run(*args):
  subprocess.check_call(args)

def create(virtualenv_dir):
  "Create and activate virtualenv in the given directory."
  # File "check" is used to make sure that we don't have
  # a partial environment in case virtualenv was interrupted.
  check_path = os.path.join(virtualenv_dir, 'check')
  if not os.path.exists(check_path):
    run('virtualenv', virtualenv_dir)
    os.mknod(check_path)
  # Activate virtualenv.
  scripts_dir = os.path.basename(sysconfig.get_path('scripts'))
  activate_this_file = os.path.join(virtualenv_dir, scripts_dir,
                                    'activate_this.py')
  with open(activate_this_file) as f:
    exec(f.read(), dict(__file__=activate_this_file))
  # Upgrade pip because installation of sphinx with pip 1.1 available on Travis
  # is broken and it doesn't support the show command.
  #from pkg_resources import get_distribution, DistributionNotFound
  #pip_version = get_distribution('pip').version
  #print('setuptools: ' + get_distribution('setuptools').version)
  #if LooseVersion(pip_version) < LooseVersion('1.5.4'):
  #  print("Upgrading pip")
  #  run('pip', 'install', '--upgrade', 'setuptools')
  #  run('pip', 'install', '--upgrade', 'pip')

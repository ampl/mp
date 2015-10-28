# Virtualenv support

import os, subprocess, sysconfig
from distutils.version import LooseVersion

def run(*args):
  subprocess.check_call(args)

def create(virtualenv_dir):
  "Create and activate virtualenv in the given directory."
  if not os.path.exists(virtualenv_dir):
    # Use a temporary directory and then rename it atomically to prevent leaving
    # a partial environment in case virtualenv is interrupted.
    tmp_dir = virtualenv_dir + '.tmp'
    run('virtualenv', tmp_dir)
    run('virtualenv', '--relocatable', tmp_dir)
    os.rename(tmp_dir, virtualenv_dir)
  # Activate virtualenv.
  scripts_dir = os.path.basename(sysconfig.get_path('scripts'))
  activate_this_file = os.path.join(virtualenv_dir, scripts_dir,
                                    'activate_this.py')
  with open(activate_this_file) as f:
    exec(f.read(), dict(__file__=activate_this_file))
  run('which', 'pip')
  run('pip', '--version')
  # Upgrade pip because installation of sphinx with pip 1.1 available on Travis
  # is broken and it doesn't support the show command.
  from pkg_resources import get_distribution, DistributionNotFound
  pip_version = get_distribution('pip').version
  if LooseVersion(pip_version) < LooseVersion('1.5.4'):
    print("Upgrading pip")
    run('pip', 'install', '--upgrade', 'pip')

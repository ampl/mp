#!/usr/bin/env python
"""
Set up build environment on Ubuntu or other Debian-based
Linux distribution.

Usage:
  bootstrap-linux.py <version> <name> <mirror> [buildbot | docker]

buildbot: install buildbot slave
"""

import platform, re, os, shutil
from bootstrap import *
from subprocess import check_call, Popen, PIPE

if __name__ == '__main__':
  vagrant = bootstrap_init()

  import docopt
  args = docopt.docopt(__doc__)

  version = args['<version>']
  image_name = args['<name>']
  mirror = args['<mirror>']
  x86_64 = platform.architecture()[0] == '64bit'
  docker = args['docker']

  env = os.environ.copy()
  env['DEBIAN_FRONTEND'] = 'noninteractive'
  
  check_call(['cp', '-r', '/support/bootstrap/cache', '/var/'])

  # Install AMPL.
  if x86_64:
    install_ampl('amplide.linux64.tgz')
  else:
    install_ampl('amplide.linux32.tgz')

  # Install build tools.
  if not installed('cmake'):
    # Install python-software-properties or software-properties-common for apt-add-repository.
    try:
      check_call(['apt-get', 'install', '-qy', 'python-software-properties'], env=env)
    except:
      check_call(['apt-get', 'install', '-qy', 'software-properties-common'], env=env)
    # Add git-core PPA for newer version of Git because version 1.7 available
    # in Lucid cannot access private repos on GitHub via a token.
    check_call(['add-apt-repository', 'ppa:git-core/ppa'], env=env)
    # Add webupd8team java PPA for Java 7.
    check_call(['add-apt-repository', 'ppa:webupd8team/java'], env=env)
    # Suppress license dialog.
    cmd = 'echo debconf shared/accepted-oracle-license-v1-1 {0} true | ' + \
          'debconf-set-selections'
    check_call(cmd.format('select'), shell=True, env=env)
    check_call(cmd.format('seen'), shell=True, env=env)
    # Install packages.
    check_call(['apt-get', 'update', '-q'], env=env)
    packages = [
      'doxygen', 'git-core', 'gcc', 'g++', 'gfortran', 'ccache', 'make',
      'oracle-java7-installer', 'oracle-java7-set-default',
      'libgtk2.0-0', 'libxrender1', 'libxtst6', # Java/Eclipse requirements
      'python-dev', 'unixodbc-dev'
    ]
    if x86_64:
      packages.append('libc6-i386')
    check_call(['apt-get', 'install', '-qy'] + packages, env=env)
    shutil.rmtree('/var/cache/oracle-jdk7-installer')

    install_cmake('cmake-3.3.0-Linux-i386.tar.gz')
    install_maven()

  # Upgrade gcc to 4.8 if necessary
  output = Popen(['gcc', '--version'], stdout=PIPE).communicate()[0]
  gcc_version = re.match(r'.* (\d+\.\d+)\.\d+', output).group(1)
  if gcc_version.split('.')[:2] < ['4', '8']:
    # Add ubuntu-toolchain-r PPA for gcc-4.8 or superior
    check_call(['add-apt-repository', 'ppa:ubuntu-toolchain-r/test'], env=env)
    check_call(['apt-get', 'update', '-q'], env=env)
    check_call(['apt-get', 'install', '-qy', 'gcc-4.8', 'g++-4.8'], env=env)
    check_call(['sudo', 'update-alternatives', '--install', 
                '/usr/bin/gcc', 'gcc', '/usr/bin/gcc-4.8', '50'], env=env)
    check_call(['sudo', 'update-alternatives', '--install', 
                '/usr/bin/g++', 'g++', '/usr/bin/g++-4.8', '50'], env=env)

  # Installs symlinks for ccache.
  for name in ['gcc', 'cc', 'g++', 'c++']:
    add_to_path(which('ccache'), name)

  install_f90cache()
  output = Popen(['gfortran', '--version'], stdout=PIPE).communicate()[0]
  gf_version = re.match(r'.* (\d+\.\d+)\.\d+', output).group(1)
  if gf_version.split('.')[:2] <= ['5', '2']:
    # f90cache 0.96 is not compatible with gfortran > 5.2
    add_to_path('/usr/local/bin/f90cache', 'gfortran-' + gf_version)

  if docker:
    # Install xvfb and miwm (a window manager) for GUI tests and x11vnc to
    # be able to remotely connect to the X server for debugging.
    check_call(['add-apt-repository',
                'deb {0} {1} universe'.format(mirror, version)], env=env)
    check_call(['apt-get', 'update', '-q'], env=env)
    check_call(['apt-get', 'install', '-qy', 'xvfb', 'x11vnc', 'miwm'], env=env)

  check_call(['apt-get', 'clean'], env=env)

  # Install LocalSolver.
  if not installed('localsolver'):
    ls_filename = 'LocalSolver_{0}_Linux{1}.run'.format(
      LOCALSOLVER_VERSION, 64 if x86_64 else 32)
    with download(LOCALSOLVER_DOWNLOADS_URL + ls_filename) as f:
      check_call(['sh', f])

  copy_optional_dependencies('linux-' + platform.machine())

  if args['buildbot'] or docker:
    ip = '172.17.42.1' if docker else None
    path = install_buildbot_slave(image_name, ip=ip)
    if not docker and path:
      pip_install('python-crontab', 'crontab')
      from crontab import CronTab
      username = 'buildbot'
      cron = CronTab(username)
      cron.new('PATH={0}:/usr/local/bin buildslave start {1}'.format(
        os.environ['PATH'], path)).every_reboot()
      cron.write()
      # Ignore errors from buildslave as the buildbot may not be accessible.
      call(['sudo', '-H', '-u', username, 'buildslave', 'start', path])

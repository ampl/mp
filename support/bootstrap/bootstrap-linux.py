#!/usr/bin/env python
"""
Set up build environment on Ubuntu or other Debian-based
Linux distribution.

Usage:
  bootstrap-linux.py [buildbot | docker]

buildbot: install buildbot slave
"""

import platform, re, shutil
from bootstrap import *
from subprocess import check_call, Popen, PIPE

if __name__ == '__main__':
  vagrant = bootstrap_init()

  import docopt
  args = docopt.docopt(__doc__)

  x86_64 = platform.architecture()[0] == '64bit'

  # Install build tools.
  if not installed('cmake'):
    # Install python-software-properties for apt-add-repository.
    check_call(['apt-get', 'install', '-qy', 'python-software-properties'])
    # Add git-core PPA for newer version of Git because version 1.7 available
    # in Lucid cannot access private repos on GitHub via a token.
    check_call(['add-apt-repository', 'ppa:git-core/ppa'])
    # Add webupd8team java PPA for Java 7.
    check_call(['add-apt-repository', 'ppa:webupd8team/java'])
    # Suppress license dialog.
    cmd = 'echo debconf shared/accepted-oracle-license-v1-1 {0} true | ' + \
          'debconf-set-selections'
    check_call(cmd.format('select'), shell=True)
    check_call(cmd.format('seen'), shell=True)
    # Install packages.
    check_call(['apt-get', 'update', '-q'])
    packages = [
      'doxygen', 'git-core', 'gcc', 'g++', 'gfortran', 'ccache', 'make',
      'oracle-java7-installer', 'oracle-java7-set-default',
      'libgtk2.0-0', 'libxrender1', 'libxtst6', # Java/Eclipse requirements
      'python-dev', 'unixodbc-dev'
    ]
    if x86_64:
      packages.append('libc6-i386')
    check_call(['apt-get', 'install', '-qy'] + packages)
    shutil.rmtree('/var/cache/oracle-jdk7-installer')

    install_cmake('cmake-3.3.0-Linux-i386.tar.gz')
    install_maven()

  # Installs symlinks for ccache.
  for name in ['gcc', 'cc', 'g++', 'c++']:
    add_to_path(which('ccache'), name)

  install_f90cache()
  output = Popen(['gfortran', '--version'], stdout=PIPE).communicate()[0]
  version = re.match(r'.* (\d+\.\d+)\.\d+', output).group(1)
  add_to_path('/usr/local/bin/f90cache', 'gfortran-' + version)

  docker = args['docker']
  if docker:
    # Install xvfb and miwm (a window manager) for GUI tests and x11vnc to
    # be able to remotely connect to the X server for debugging.
    check_call(['add-apt-repository',
                'deb http://archive.ubuntu.com/ubuntu lucid universe'])
    check_call(['apt-get', 'update', '-q'])
    check_call(['apt-get', 'install', '-qy', 'xvfb', 'x11vnc', 'miwm'])

  check_call(['apt-get', 'clean'])

  # Install LocalSolver.
  if not installed('localsolver'):
    ls_filename = 'LocalSolver_{0}_Linux{1}.run'.format(
      LOCALSOLVER_VERSION, 64 if x86_64 else 32)
    with download(LOCALSOLVER_DOWNLOADS_URL + ls_filename) as f:
      check_call(['sh', f])

  copy_optional_dependencies('linux-' + platform.machine())

  if args['buildbot'] or docker:
    ip = '172.17.42.1' if docker else None
    path = install_buildbot_slave(
      'lucid64' if x86_64 else 'lucid32', ip=ip)
    if not docker and path:
      pip_install('python-crontab', 'crontab')
      from crontab import CronTab
      cron = CronTab(username)
      cron.new('PATH={0}:/usr/local/bin buildslave start {1}'.format(
        os.environ['PATH'], path)).every_reboot()
      cron.write()
      # Ignore errors from buildslave as the buildbot may not be accessible.
      call(['sudo', '-H', '-u', username, 'buildslave', 'start', path])

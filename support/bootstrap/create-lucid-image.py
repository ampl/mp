#!/usr/bin/env python
"""
Create a base docker image for Ubuntu Lucid.

Usage:
  create-lucid-image.py [32]
"""

import os, sys
from subprocess import check_call, Popen, PIPE

if __name__ == '__main__':
  parent_dir = os.path.join(os.path.dirname(__file__), '..')
  sys.path.append(parent_dir)
  import docopt
  args = docopt.docopt(__doc__)

  image = 'lucid64'
  cmd = ['sudo', 'debootstrap']
  if args['32']:
    image = 'lucid32'
    cmd.append('--arch=i386')
  cmd += ['lucid', image]
  try:
    check_call(cmd)
    p1 = Popen(['sudo', 'tar', '-C', image, '-c', '.'], stdout=PIPE)
    p2 = Popen(['docker', 'import', '-', image], stdin=p1.stdout)
    p2.communicate()
  finally:
    check_call(['sudo', 'rm', '-rf', image])
  check_call(['docker', 'run', '-v', os.path.abspath(parent_dir) + ':/support', '--name',
              image, image, 'support/bootstrap/bootstrap-linux.py', 'docker'])
  check_call(['docker', 'commit', image, 'vitaut/ampl:' + image])

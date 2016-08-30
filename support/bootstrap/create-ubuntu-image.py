#!/usr/bin/env python
"""
Create a base docker image for Ubuntu Lucid.

Usage:
  create-ubuntu-image.py [32|64] <version> <name> <mirror>
"""

import os, sys
from subprocess import check_call, call, Popen, PIPE

if __name__ == '__main__':
  parent_dir = os.path.join(os.path.dirname(__file__), '..')
  sys.path.append(parent_dir)
  import docopt
  args = docopt.docopt(__doc__)

  version = args['<version>']
  image_name = args['<name>']
  mirror = args['<mirror>']
  nbits = 32 if args['32'] else 64   

  try:
    check_call(['docker', 'inspect', 'debootstrap:' + image_name])
  except:
    try:
      cmd = ['sudo', 'debootstrap']
      if nbits == 32:    
        cmd.append('--arch=i386')
        cmd += [version, image_name, mirror] 
      check_call(cmd)
      call(['docker', 'rmi', 'debootstrap:' + image_name])  
      p1 = Popen(['sudo', 'tar', '-C', image_name, '-c', '.'], stdout=PIPE)
      p2 = Popen(['docker', 'import', '-', 'debootstrap:' + image_name], stdin=p1.stdout)
      p2.communicate()
    finally:
      check_call(['sudo', 'rm', '-rf', image_name])

  call(['docker', 'rm', image_name + '_tmp'])  
  check_call(['docker', 'run', '-v', 
              os.path.abspath(parent_dir) + ':/support',
              '--name', image_name + '_tmp', 'debootstrap:' + image_name, '/bin/sh', '-c',
              'export ftp_proxy=http://172.17.42.1:3128;'
              'export http_proxy=http://172.17.42.1:3128;' 
              'export https_proxy=http://172.17.42.1:3128;'
              'apt-get install -y python;'
              'support/bootstrap/bootstrap-linux.py'
              ' {0} {1} {2} docker'.format(version, image_name, mirror)])
  call(['docker', 'rmi', 'buildbot:' + image_name])
  check_call(['docker', 'commit', image_name + '_tmp', 
              'buildbot:' + image_name])
  call(['docker', 'rm', image_name + '_tmp'])


#!/usr/bin/env python
# This scripts generates test .nl files from AMPL model and data files.

from subprocess import Popen, PIPE

DATA_DIR = 'data/'
MODELS_DIR = '../models/logic/'

def generate_nl_file(stub, *files, **kwargs):
  code = ''
  dir = MODELS_DIR
  for key, value in kwargs.iteritems():
    if key == 'code':
      code = value
    elif key == 'dir':
      dir = value
    else:
      raise Exception('Invalid parameter: ' + key)
  input = ''
  for f in files:
    f = f.replace('$', stub)
    if f.endswith('.dat'):
      input += 'data ' + dir + f + ';\n';
    else:
      input += 'model ' + dir + f + ';\n';
  if 'code' in kwargs:
    input += kwargs['code']
  input += 'write g' + DATA_DIR + stub + ';\n';
  ampl = Popen(['ampl'], stdin=PIPE)
  ampl.communicate(input)

generate_nl_file('assign0', '$.mod', '$.dat')
generate_nl_file('assign1', '$.mod', '$.dat')
generate_nl_file('assign1a', '$.mod', '$.dat')
generate_nl_file('balassign0', '$.mod', '$.dat')
generate_nl_file('balassign1', '$.mod', '$.dat')
generate_nl_file('magic', '$.mod', code='let n := 8;')
generate_nl_file('nqueens', '$.mod', code='let n := 8;')
generate_nl_file('nqueens0', '$.mod', code='let n := 8;')
generate_nl_file('party1', '$.mod', '$.dat')
generate_nl_file('party2', '$.mod', '$.dat')
generate_nl_file('objconst', '$.ampl', dir='data/')
generate_nl_file('numberof', '$.ampl', dir='data/')

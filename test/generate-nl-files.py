#!/usr/bin/env python
# This scripts generates test .nl files from AMPL model and data files.

from subprocess import Popen, PIPE

def generate_nl_file(stub, *files, **kwargs):
  code = ''
  for key, value in kwargs.iteritems():
    if key == 'code':
      code = value
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
  input += 'write g' + outdir + stub + ';\n';
  ampl = Popen(['ampl'], stdin=PIPE)
  ampl.communicate(input)

outdir = 'data/'

dir = 'data/'
generate_nl_file('test', '$.ampl')
generate_nl_file('numberof', '$.ampl')
generate_nl_file('feasible', '$.ampl')
generate_nl_file('infeasible', '$.ampl')
generate_nl_file('unbounded', '$.ampl')
generate_nl_file('noobj', '$.ampl')
generate_nl_file('simple', '$.ampl')
generate_nl_file('ssd', '$.ampl')
generate_nl_file('suffix', '$.ampl')
generate_nl_file('element', '$.ampl')

dir = outdir = 'data/smps/'
generate_nl_file('inconsistent-probabilities', '$.ampl')
generate_nl_file('int-var', '$.ampl')
generate_nl_file('nonlinear', '$.ampl')
generate_nl_file('random-bound', '$.ampl')
generate_nl_file('random-con-matrix', '$.ampl')
generate_nl_file('random-con-matrix2', '$.ampl')
generate_nl_file('random-rhs', '$.ampl')
generate_nl_file('range-con', '$.ampl')
generate_nl_file('single-scenario', '$.ampl')
generate_nl_file('single-stage', '$.ampl')
generate_nl_file('three-stage', '$.ampl')
generate_nl_file('vars-not-in-stage-order', '$.ampl')
generate_nl_file('zero-core-coefs', '$.ampl')
generate_nl_file('zero-core-con', '$.ampl')

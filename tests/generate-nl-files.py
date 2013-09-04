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

dir = '../models/logic/'
outdir = 'data/'
generate_nl_file('assign0', '$.mod', 'assign.dat')
generate_nl_file('assign1', '$.mod', 'assign.dat')
generate_nl_file('assign2', '$.mod', 'assign.dat')
generate_nl_file('assign3', '$.mod', 'assign.dat')
generate_nl_file('balassign0', '$.mod', '$.dat')
generate_nl_file('balassign1', '$.mod', '$.dat')
generate_nl_file('flowshp0', '$.mod', 'flowshp.dat')
generate_nl_file('flowshp1', '$.mod', 'flowshp.dat')
generate_nl_file('flowshp2', '$.mod', 'flowshp.dat')
generate_nl_file('grpassign0', '$.mod', '$.dat')
generate_nl_file('magic', '$.mod', code='let n := 8;')
generate_nl_file('mapcoloring', '$.mod')
generate_nl_file('nqueens', '$.mod', code='let n := 8;')
generate_nl_file('nqueens0', '$.mod', code='let n := 8;')
generate_nl_file('openshop', '$.mod', '$.dat')
generate_nl_file('party1', '$.mod', '$.dat')
generate_nl_file('sched0', '$.mod', 'sched.dat')
generate_nl_file('sched1', '$.mod', 'sched.dat')
generate_nl_file('sched2', '$.mod', 'sched.dat')
generate_nl_file('send-more-money', '$.mod')
generate_nl_file('send-most-money', '$.mod')
generate_nl_file('seq0', '$.mod')
generate_nl_file('seq0a', '$.mod')
generate_nl_file('seq1', '$.mod')
generate_nl_file('seq1a', '$.mod')
generate_nl_file('seq2', '$.mod')
generate_nl_file('sudokuHard', 'sudoku.mod', '$.dat')
generate_nl_file('sudokuVeryEasy', 'sudoku.mod', '$.dat')

dir = 'data/'
generate_nl_file('test', '$.ampl')
generate_nl_file('objconst', '$.ampl')
generate_nl_file('objconstint', '$.ampl')
generate_nl_file('numberof', '$.ampl')
generate_nl_file('feasible', '$.ampl')
generate_nl_file('infeasible', '$.ampl')
generate_nl_file('unbounded', '$.ampl')
generate_nl_file('noobj', '$.ampl')
generate_nl_file('simple', '$.ampl')
generate_nl_file('ssd', '$.ampl')
generate_nl_file('suffix', '$.ampl')

dir = outdir = 'data/miplib/'
generate_nl_file('assign1', '$.mod', '$.dat')

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

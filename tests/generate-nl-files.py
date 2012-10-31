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
generate_nl_file('grpassign1', '$.mod', '$.dat')
generate_nl_file('grpassign1a', '$.mod', '$.dat')
generate_nl_file('magic', '$.mod', code='let n := 8;')
generate_nl_file('mapcoloring', '$.mod')
generate_nl_file('nqueens', '$.mod', code='let n := 8;')
generate_nl_file('nqueens0', '$.mod', code='let n := 8;')
generate_nl_file('party1', '$.mod', '$.dat')
generate_nl_file('party2', '$.mod', '$.dat')
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
generate_nl_file('objconst', '$.ampl', dir='data/')
generate_nl_file('numberof', '$.ampl', dir='data/')
generate_nl_file('feasible', '$.ampl', dir='data/')
generate_nl_file('infeasible', '$.ampl', dir='data/')
generate_nl_file('unbounded', '$.ampl', dir='data/')

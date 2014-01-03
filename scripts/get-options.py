#!/usr/bin/env python
# This scripts extract options from solvers and formats them in HTML.

import os, errno
from ctypes import cdll, Structure, c_char_p, c_void_p, byref
from docutils.core import publish_string

class SolverOption:
  def __init__(self, name, description):
    self.name = name
    self.description = description

class ASL_SolverOption(Structure):
  _fields_ = [("name", c_char_p),
              ("description", c_char_p)]

class Solver:
  def __init__(self, libpath):
    self.lib = cdll.LoadLibrary(libpath)
    error = c_void_p()
    self.solver = self.lib.ASL_CreateSolver(byref(error))
    self.lib.ASL_GetErrorMessage.restype = c_char_p
    if not self.solver:
      raise Exception(self.lib.ASL_GetErrorMessage(error))

  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    self.lib.ASL_DestroySolver(self.solver)

  def _check(self, success):
    if not success:
      error = self.lib.ASL_GetLastError(self.solver)
      raise Exception(self.lib.ASL_GetErrorMessage(error))

  def get_options(self):
    num_options = self.lib.ASL_GetSolverOptions(self.solver, None, 0)
    self._check(num_options != -1)
    asl_options = (ASL_SolverOption * num_options)()
    num_options = self.lib.ASL_GetSolverOptions(self.solver, asl_options, num_options)
    self._check(num_options != -1)
    options = []
    for opt in asl_options:
      options.append(SolverOption(opt.name, opt.description))
    return options

solver_names = ['gecode', 'ilogcp', 'jacop', 'ssdsolver', 'sulum']
for solver_name in solver_names:
  html_filename = solver_name + '-options.html'
  try:
    os.remove(html_filename)
  except OSError as e:
    if e.errno != errno.ENOENT:
      raise
  libname = '../solvers/{0}/libampl{0}.so'.format(solver_name)
  if not os.path.exists(libname):
    continue
  with Solver(libname) as solver:
    rst = ''
    for opt in solver.get_options():
      rst += '``{}``\n'.format(opt.name)
      for line in opt.description.split('\n'):
        rst += '  ' + line + '\n'
    with open(html_filename, 'w') as f:
      f.write(publish_string(rst, writer_name='html'))

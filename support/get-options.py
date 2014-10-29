#!/usr/bin/env python
# This scripts extract options from solvers and formats them in HTML.

import errno, os, sys
from ctypes import cdll, Structure, c_int, c_char_p, c_void_p, byref
from docutils.core import publish_parts
from docutils.parsers.rst import directives
from docutils.parsers import rst
from docutils import nodes

class SolverOption:
  def __init__(self, name, description, values):
    self.name = name
    self.description = description
    self.values = values

class MP_SolverOptionInfo(Structure):
  _fields_ = [("name", c_char_p),
              ("description", c_char_p),
              ("flags", c_int),
              ("option", c_void_p)]

class OptionValueInfo:
  def __init__(self, value, description):
    self.value = value
    self.description = description

class MP_OptionValueInfo(Structure):
  _fields_ = [("value", c_char_p),
              ("description", c_char_p)]
  
class Solver:
  def __init__(self, libpath):
    self.lib = cdll.LoadLibrary(libpath)
    error = c_void_p()
    self.solver = self.lib.MP_CreateSolver(byref(error))
    self.lib.MP_GetErrorMessage.restype = c_char_p
    if not self.solver:
      raise Exception(self.lib.MP_GetErrorMessage(error))

  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    self.lib.MP_DestroySolver(self.solver)

  def _check(self, success):
    if not success:
      error = self.lib.MP_GetLastError(self.solver)
      raise Exception(self.lib.MP_GetErrorMessage(error))
    
  def option_header(self):
    self.lib.MP_GetOptionHeader.restype = c_char_p
    header = self.lib.MP_GetOptionHeader(self.solver)
    self._check(header is not None)
    return header

  def get_options(self):
    num_options = self.lib.MP_GetSolverOptions(self.solver, None, 0)
    self._check(num_options != -1)
    info = (MP_SolverOptionInfo * num_options)()
    num_options = self.lib.MP_GetSolverOptions(self.solver, info, num_options)
    self._check(num_options != -1)
    options = []
    for i in info:
      num_values = 0
      if (i.flags & 1) != 0:
        num_values = self.lib.MP_GetOptionValues(self.solver, i.option, None, 0)
        self._check(num_values != -1)
      values = []
      if num_values != 0:
        asl_values = (MP_OptionValueInfo * num_values)()
        num_values = self.lib.MP_GetOptionValues(self.solver, i.option, asl_values, num_values)
        self._check(num_values != -1)
        for v in asl_values:
          values.append(OptionValueInfo(v.value, v.description))
      options.append(SolverOption(i.name, i.description, values))
    return options

class ValueTableDirective(rst.Directive):
  values = []
  
  def run(self):
    if ValueTableDirective.values[0].description is None:
      list = nodes.bullet_list()
      for v in ValueTableDirective.values:
        item = nodes.list_item()
        item += nodes.literal(v.value, v.value)
        list += item
      return [list]
    table = nodes.table()
    tgroup = nodes.tgroup()
    tbody = nodes.tbody()
    for v in ValueTableDirective.values:
      row = nodes.row()
      entry = nodes.entry()
      entry += nodes.literal(v.value, v.value)
      row += entry
      entry = nodes.entry()
      entry += nodes.paragraph(text=v.description)
      row += entry
      tbody += row
    tgroup += nodes.colspec(colwidth=10)
    tgroup += nodes.colspec(colwidth=90)
    tgroup += tbody
    table += tgroup
    return [table]

directives.register_directive('value-table', ValueTableDirective)

solver_names = ['Gecode', 'IlogCP', 'Jacop', 'LocalSolver', 'SSDSolver', 'Sulum']
for solver_name in solver_names:
  solver_name_lower = solver_name.lower()
  html_filename = solver_name_lower + '-options.html'
  try:
    os.remove(html_filename)
  except OSError as e:
    if e.errno != errno.ENOENT:
      raise
  libname = '{}/libampl{}.so'.format(sys.argv[1], solver_name_lower)
  if not os.path.exists(libname):
    continue
  with Solver(libname) as solver:
    with open(html_filename, 'w') as f:
      f.write(publish_parts(
        solver.option_header(), writer_name='html')['html_body'].encode('utf-8'))
      f.write('<h2>Options</h2>\n')
      for opt in solver.get_options():
        rst = '``{}``\n'.format(opt.name)
        ValueTableDirective.values = opt.values
        for line in opt.description.split('\n'):
          rst += '  ' + line + '\n'
        f.write(publish_parts(rst, writer_name='html')['body'].encode('utf-8'))

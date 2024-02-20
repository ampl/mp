'''
 NL Writer Python API test.

 Copyright (C) 2024 AMPL Optimization Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
'''

import sys
import numpy as np
from scipy.sparse import csr_matrix

import nlwpy as m

assert m.__version__ == "0.0.1"

nlwo = m.NLW2_MakeNLOptionsBasic_Default()
assert 0 == nlwo.n_text_mode_
assert 0 == nlwo.want_nl_comments_
assert 1 == nlwo.flags_

b = 500
if b>400:
  b=b+400

print(b)

## ---------------------------------------------------------------
class ModelBuilder:
  def GetModel(self):
    nlme = m.NLW2_NLModel(self.prob_name_)

    nlme.SetCols(self.var_lb_, self.var_ub_, self.var_type_)
    nlme.SetColNames(self.var_names_)

    if self.A_ is not None:
      self.A_ = csr_matrix(self.A_)
      nlme.SetRows(self.row_lb_, self.row_ub_,
                   self.A_format_,
                   self.A_.indptr, self.A_.indices, self.A_.data)
    nlme.SetRowNames(self.row_names_)

    nlme.SetLinearObjective(self.obj_sense_, self.obj_c0_,
                            self.obj_c_)

    if self.Q_ is not None:
      self.Q_ = csr_matrix(self.Q_)
      nlme.SetHessian(self.Q_format_,
                      self.Q_.indptr, self.Q_.indices, self.Q_.data)
    nlme.SetObjName(self.obj_name_)

    return nlme

  def Check(self, sol):
    if not self.ApproxEqual(sol.obj_val_, self.obj_val_ref_):
      print("MIQP 1: wrong obj val ({:.17} !~ {:.17})".format(
             sol.obj_val_, self.obj_val_ref_))
      return False

    for i in range(len(sol.x_)):
      if not self.ApproxEqual(self.x_ref_[i], sol.x_[i]):
        print("MIQP 1: wrong x[{}] ({:.17} !~ {:.17})".format(
               i+1, sol.x_[i], self.x_ref_[i]))
        return False

    print("MIQP 1: solution check ok, obj={:.17}.".format(sol.obj_val_))
    return True

  def ApproxEqual(self, n, m):
    return abs(n-m) \
        <= 1e-5 * min(1.0, abs(n)+abs(m))

  def __init__(self):
    self.prob_name_ = "nlwpy_prob"
    self.var_lb_ = [0, -3, 0, -1, -1, -2]
    self.var_ub_ = [0, 20, 1, 1e20, -1, 10]
    self.var_type_ = [0, 1, 1, 1, 0, 0]
    self.var_names_ = \
      ["x1_4", "x2_6", "x3_5", "x4_3", "x5_1", "x6_2"]
    self.A_format_ = m.NLW2_MatrixFormat.Rowwise
    self.A_ = np.array([
      [0,1,1,1,0,1],
      [0,1,-1,-1,0,1]])
    self.row_lb_ = [15, 10]
    self.row_ub_ = [15, np.inf]
    self.row_names_ = ["C1", "C2"]
    self.obj_sense_ = m.NLW2_ObjSense.Minimize
    self.obj_c0_ = 3.24
    self.obj_c_ = [0,1,0,0,0,0]
    self.Q_format_ = m.NLW2_HessianFormat.Square
    self.Q_ = np.zeros([6, 6])
    self.Q_[3, 3] = 10
    self.Q_[3, 5] = 12
    self.Q_[4, 4] = 14
    self.obj_name_ = "obj[1]"

    ### Solution
    self.x_ref_ = [0, 5, 1, -1, -1, 10]
    self.obj_val_ref_ = -39.76

def SolveAndCheck(solver, sopts, binary, stub):
  mb = ModelBuilder()
  nlme = mb.GetModel()
  nlse = m.NLW2_NLSolver()
  nlopts = m.NLW2_MakeNLOptionsBasic_Default()
  nlopts.n_text_mode_ = not binary
  nlopts.want_nl_comments_ = 1
  nlse.SetNLOptions(nlopts)
  nlse.SetFileStub(stub)
  sol = nlse.Solve(nlme, solver, sopts)
  if sol.solve_result_ > -2:      ## Some method for this?
    if (not mb.Check(sol)):
      print("Solution check failed.")
      return False
  else:
    print(nlse.GetErrorMessage())
    return False

  return True

argc=len(sys.argv)
argv=sys.argv

if argc<2:
  print("AMPL NL Writer Python API example.\n"
  "Usage:\n"
  "  python <this_script> <solver> [\"<solver_options>\" [binary/text [<stub>]]],\n\n"
  "where <solver> is ipopt, gurobi, minos, ...;\n"
  "binary/text is the NL format (default: binary.)\n"
  "Examples:\n"
  "  python <this_script> ipopt \"\" text /tmp/stub\n"
  "  python <this_script> gurobi \"nonconvex=2 funcpieces=-2 funcpieceratio=1e-4\"")
  sys.exit(0)

solver = argv[1] if (argc>1) else "minos"
sopts =  argv[2] if argc>2 else ""
binary = ((argc<=3) or "text" == argv[3])
stub = argv[4] if argc>4 else ""

if not SolveAndCheck(solver, sopts, binary, stub):
  print("SolveAndCheck() failed.")
  sys.exit(1)

## ---------------------------------------------------------------
print("Test finished.")

/*
 Utilities for writing AMPL solver drivers.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef SOLVERS_UTIL_DRIVER_H_
#define SOLVERS_UTIL_DRIVER_H_

#include "solvers/util/expr.h"

namespace ampl {

// An AMPL solver driver.
class Driver {
 private:
  ASL_fg *asl_;

  // Do not implement.
  Driver(const Driver&);
  Driver& operator=(const Driver&);

 public:
  Driver();
  virtual ~Driver();

  // Gets the options.
  int GetOptions(char **argv, Option_Info *oi);

  // Reads the problem form a .nl file.
  bool Read(char **&argv, Option_Info *oi);

  // Returns the number of variables.
  int num_vars() const { return asl_->i.n_var_; }

  // Returns the number of objectives.
  int num_objs() const { return asl_->i.n_obj_; }

  // Returns the number of constraints.
  int num_cons() const { return asl_->i.n_con_; }

  // Returns the number of integer variables including binary.
  int num_integer_vars() const {
    return asl_->i.nbv_ + asl_->i.niv_ + asl_->i.nlvbi_ +
        asl_->i.nlvci_ + asl_->i.nlvoi_;
  }

  // Returns the number of continuous variables.
  int num_continuous_vars() const {
    return num_vars() - num_integer_vars();
  }

  // Returns the number of nonlinear objectives.
  int num_nonlinear_objs() const { return asl_->i.nlo_; }

  // Returns the number of nonlinear constraints.
  int num_nonlinear_cons() const { return asl_->i.nlc_; }

  // Returns the number of logical constraints.
  int num_logical_cons() const { return asl_->i.n_lcon_; }

  // Returns the variable lower bound.
  double GetVarLB(int var_index) const {
    assert(var_index >= 0 && var_index < num_vars());
    return asl_->i.LUv_[var_index];
  }

  // Returns the variable lower bound.
  double GetVarUB(int var_index) const {
    assert(var_index >= 0 && var_index < num_vars());
    return asl_->i.Uvx_[var_index];
  }

  // Returns the constraint lower bound.
  double GetConLB(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return asl_->i.LUrhs_[con_index];
  }

  // Returns the constraint lower bound.
  double GetConUB(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return asl_->i.Urhsx_[con_index];
  }

  enum ObjType { MIN = 0, MAX = 1 };

  // Returns the objective type.
  ObjType GetObjType(int obj_index) const {
    assert(obj_index >= 0 && obj_index < num_objs());
    return static_cast<ObjType>(asl_->i.objtype_[obj_index]);
  }

  // Returns the linear part of the objective gradient.
  ograd *GetObjGradient(int obj_index) const {
    assert(obj_index >= 0 && obj_index < num_objs());
    return asl_->i.Ograd_[obj_index];
  }

  // Returns the linear part of the constraint gradient.
  cgrad *GetConGradient(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return asl_->i.Cgrad_[con_index];
  }

  // Returns the nonlinear part of an objective expression.
  NumericExpr GetNonlinearObjExpr(int obj_index) const {
    assert(obj_index >= 0 && obj_index < num_objs());
    return Expr::Create<NumericExpr>(asl_->I.obj_de_[obj_index].e);
  }

  // Returns the nonlinear part of a constraint expression.
  NumericExpr GetNonlinearConExpr(int con_index) const {
    assert(con_index >= 0 && con_index < num_cons());
    return Expr::Create<NumericExpr>(asl_->I.con_de_[con_index].e);
  }

  // Returns a logical constraint expression.
  LogicalExpr GetLogicalConExpr(int lcon_index) const {
    assert(lcon_index >= 0 && lcon_index < num_logical_cons());
    return Expr::Create<LogicalExpr>(asl_->I.lcon_de_[lcon_index].e);
  }

  // Returns the solve code.
  int solve_code() const { return asl_->p.solve_code_; }

  // Sets the solve code.
  void SetSolveCode(int value) {
    asl_->p.solve_code_ = value;
  }

  // Writes the solution.
  void WriteSolution(char *msg, double *x, double *y, Option_Info* oi) {
    write_sol_ASL(reinterpret_cast<ASL*>(asl_), msg, x, y, oi);
  }
};
}

#endif  // SOLVERS_UTIL_DRIVER_H_


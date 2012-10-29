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

 public:
  Driver() : asl_(reinterpret_cast<ASL_fg*>(ASL_alloc(ASL_read_fg))) {}
  virtual ~Driver() { ASL_free(reinterpret_cast<ASL**>(&asl_)); }

  ASL_fg *asl() const { return asl_; }

  // Returns the nonlinear part of an objective expression.
  NumericExpr GetNonlinearObjExpr(int obj_index) const {
    assert(obj_index >= 0 && obj_index < asl_->i.n_obj_);
    return Expr::Create<NumericExpr>(asl_->I.obj_de_[obj_index].e);
  }

  // Returns the nonlinear part of a constraint expression.
  NumericExpr GetNonlinearConExpr(int con_index) const {
    assert(con_index >= 0 && con_index < asl_->i.n_con_);
    return Expr::Create<NumericExpr>(asl_->I.con_de_[con_index].e);
  }

  // Returns a logical constraint expression.
  LogicalExpr GetLogicalConExpr(int lcon_index) const {
    assert(lcon_index >= 0 && lcon_index < asl_->i.n_lcon_);
    return Expr::Create<LogicalExpr>(asl_->I.lcon_de_[lcon_index].e);
  }
};
}

#endif  // SOLVERS_UTIL_DRIVER_H_


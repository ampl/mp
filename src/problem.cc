/*
 Optimization problem

 Copyright (C) 2014 AMPL Optimization Inc

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

 Author: Victor Zverovich
 */

#include "mp/problem.h"

using mp::Problem;

int Problem::GetSuffixSize(int suffix_type) {
  switch (suffix_type) {
  default:
    MP_ASSERT(false, "invalid suffix type");
    // Fall through.
  case suf::VAR:
    return vars_.capacity();
  case suf::CON:
    return algebraic_cons_.capacity();
  case suf::OBJ:
    return linear_objs_.capacity();
  case suf::PROBLEM:
    return 1;
  }
}

Problem::LinearObjBuilder Problem::AddObj(
    mp::obj::Type type, mp::NumericExpr expr, int num_linear_terms) {
  MP_ASSERT(linear_objs_.size() < MP_MAX_INDEX, "too many objectives");
  obj_types_.push_back(type != obj::MIN);
  linear_objs_.push_back(LinearExpr());
  LinearExpr &linear_expr = linear_objs_.back();
  linear_expr.Reserve(num_linear_terms);
  if (expr) {
    if (nonlinear_objs_.empty()) {
      nonlinear_objs_.reserve(linear_objs_.capacity());
      nonlinear_objs_.resize(linear_objs_.size() - 1);
    }
    nonlinear_objs_.push_back(expr);
  }
  return LinearObjBuilder(&linear_expr);
}

Problem::LinearConBuilder Problem::AddCon(
    mp::NumericExpr expr, double lb, double ub, int num_linear_terms) {
  MP_ASSERT(algebraic_cons_.size() < MP_MAX_INDEX,
            "too many algebraic constraints");
  algebraic_cons_.push_back(AlgebraicCon(lb, ub));
  AlgebraicCon &con = algebraic_cons_.back();
  con.linear_expr.Reserve(num_linear_terms);
  if (expr) {
    if (nonlinear_cons_.empty()) {
      nonlinear_cons_.reserve(algebraic_cons_.capacity());
      nonlinear_cons_.resize(algebraic_cons_.size() - 1);
    }
    nonlinear_cons_.push_back(expr);
  }
  return LinearConBuilder(&con.linear_expr);
}

void Problem::SetComplement(int con_index, int var_index, int flags) {
  MP_ASSERT(0 <= con_index && con_index <= num_algebraic_cons(),
            "invalid index");
  if (compl_vars_.size() <= con_index) {
    compl_vars_.reserve(algebraic_cons_.capacity());
    compl_vars_.resize(algebraic_cons_.size());
  }
  compl_vars_[con_index] = var_index + 1u;
  double inf = std::numeric_limits<double>::infinity();
  AlgebraicCon &con = algebraic_cons_[con_index];
  con.lb = (flags & comp::INF_LB) != 0 ? -inf : 0;
  con.ub = (flags & comp::INF_UB) != 0 ?  inf : 0;
}

void Problem::SetInfo(const mp::ProblemInfo &info) {
  vars_.reserve(info.num_vars);
  var_int_.reserve(info.num_vars);
  obj_types_.reserve(info.num_objs);
  linear_objs_.reserve(info.num_objs);
  if (info.num_nl_objs != 0)
    nonlinear_objs_.reserve(info.num_objs);
  algebraic_cons_.reserve(info.num_algebraic_cons);
  if (info.num_compl_conds != 0)
    compl_vars_.reserve(info.num_algebraic_cons);
  if (info.num_nl_cons != 0)
    nonlinear_cons_.reserve(info.num_algebraic_cons);
  logical_cons_.reserve(info.num_logical_cons);
  int num_common_exprs = info.num_common_exprs();
  linear_exprs_.reserve(num_common_exprs);
  nonlinear_exprs_.reserve(num_common_exprs);
  funcs_.reserve(info.num_funcs);
}

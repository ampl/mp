/*
 AMPL solver for problems with second-order stochastic dominance (SSD)
 constraints.

 Copyright (C) 2013 AMPL Optimization Inc

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

#ifndef AMPL_SOLVERS_SSDSOLVER_H
#define AMPL_SOLVERS_SSDSOLVER_H

#include <vector>

#include "solvers/util/solver.h"

#define SSDSOLVER_VERSION 20130226

namespace ampl {

// Expression visitor that extracts SSD constraints.
class SSDExtractor : public ExprVisitor<SSDExtractor, void, void> {
 private:
  // A matrix of variable coefficients in the SSD constraint with one
  // row per scenario and one column per variable. The matrix is stored in
  // a row-major order.
  std::vector<double> coefs_;

  unsigned num_vars_;
  unsigned con_index_;
  double sign_;

  // A vector of right-hand sides (reference distribution) in the SSD
  // constraint with one element per scenario.
  std::vector<double> rhs_;

  friend class ExprVisitor<SSDExtractor, void, void>;

  void VisitMult(BinaryExpr e) {
    NumericConstant coef = Cast<NumericConstant>(e.lhs());
    Variable var = Cast<Variable>(e.rhs());
    if (!coef || !var)
      throw UnsupportedExprError::CreateFromExprString("nonlinear *");
    coefs_[con_index_ * num_vars_ + var.index()] = sign_ * coef.value();
  }

  void VisitSum(SumExpr e) {
    for (SumExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
      Visit(*i);
  }

 public:
  SSDExtractor(unsigned num_scenarios, unsigned num_vars)
  : coefs_(num_scenarios * num_vars), num_vars_(num_vars),
    con_index_(0), sign_(0), rhs_(num_scenarios) {}

  void Extract(CallExpr call) {
    for (CallExpr::arg_expr_iterator
        i = call.arg_expr_begin(), e = call.arg_expr_end(); i != e; ++i) {
      sign_ = call.arg_index(i) == 0 ? 1 : -1;
      Visit(*i);
    }
    rhs_[con_index_] = call.arg_constant(1) - call.arg_constant(0);
    ++con_index_;
  }

  const double *coefs() const { return &coefs_[0]; }
  const std::vector<double> &rhs() const { return rhs_; }
};

class SSDSolver : public Solver {
 private:
  bool output_;
  bool scaled_;
  double abs_tolerance_;
  std::string solver_name_;

  int GetBoolOption(const char *, bool *opt) const { return *opt; }
  void SetBoolOption(const char *name, int value, bool *opt) {
    if (value != 0 && value != 1)
      throw InvalidOptionValue(name, value);
    *opt = value != 0;
  }

  double GetAbsTolerance(const char *) const { return abs_tolerance_; }
  void SetAbsTolerance(const char *name, double value) {
    if (value < 0)
      throw InvalidOptionValue(name, value);
    abs_tolerance_ = value;
  }

  std::string GetSolverName(const char *) const { return solver_name_; }
  void SetSolverName(const char *, const char *value) { solver_name_ = value; }

 protected:
  void DoSolve(Problem &p);

 public:
  SSDSolver();
};
}

#endif // AMPL_SOLVERS_SSDSOLVER_H

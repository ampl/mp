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

#ifndef MP_SOLVERS_SSDSOLVER_H_
#define MP_SOLVERS_SSDSOLVER_H_

#include <vector>

#include "asl/expr.h"
#include "asl/aslsolver.h"

#define SSDSOLVER_VERSION 20130226

namespace mp {

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

  void VisitNumericConstant(NumericConstant n) {
    rhs_[con_index_] -= sign_ * n.value();
  }

 public:
  SSDExtractor(unsigned num_scenarios, unsigned num_vars)
  : coefs_(num_scenarios * num_vars), num_vars_(num_vars),
    con_index_(0), sign_(0), rhs_(num_scenarios) {}

  void Extract(CallExpr call) {
    assert(call.num_args() == 2);
    sign_ = 1;
    Visit(Cast<NumericExpr>(call[0]));
    sign_ = -1;
    Visit(Cast<NumericExpr>(call[1]));
    ++con_index_;
  }

  const double *coefs() const { return &coefs_[0]; }
  const std::vector<double> &rhs() const { return rhs_; }
};

class SSDSolver : public ASLSolver {
 private:
  bool output_;
  bool scaled_;
  double abs_tolerance_;
  std::string solver_name_;

  int GetBoolOption(const SolverOption &, bool *ptr) const { return *ptr; }
  void SetBoolOption(const SolverOption &opt, int value, bool *ptr) {
    if (value != 0 && value != 1)
      throw InvalidOptionValue(opt, value);
    *ptr = value != 0;
  }

  double GetAbsTolerance(const SolverOption &) const { return abs_tolerance_; }
  void SetAbsTolerance(const SolverOption &opt, double value) {
    if (value < 0)
      throw InvalidOptionValue(opt, value);
    abs_tolerance_ = value;
  }

  std::string GetSolverName(const SolverOption &) const { return solver_name_; }
  void SetSolverName(const SolverOption &, const char *value) {
    solver_name_ = value;
  }

 protected:
  int DoSolve(Problem &p, SolutionHandler &sh);

 public:
  SSDSolver();
};
}

#endif  // MP_SOLVERS_SSDSOLVER_H_

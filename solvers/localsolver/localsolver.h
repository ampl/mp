/*
 AMPL solver interface to LocalSolver.

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

#ifndef AMPL_SOLVERS_LOCALSOLVER_H
#define AMPL_SOLVERS_LOCALSOLVER_H

#include <localsolver.h>
#include "solvers/util/solver.h"

namespace ampl {

namespace ls = localsolver;

// Converter of optimization problems from NL to LocalSolver format.
class NLToLocalSolverConverter :
  public ExprVisitor<NLToLocalSolverConverter,
                     ls::LSExpression*, ls::LSExpression*> {
 private:
  ls::LSModel &model_;
  std::vector<ls::LSExpression*> vars_;

  template<typename Term>
  ls::LSExpression *ConvertExpr(LinearExpr<Term> linear, NumericExpr nonlinear);

 public:
  NLToLocalSolverConverter(ls::LSModel &model) : model_(model) {}

  void Convert(const Problem &p);

  ls::LSExpression *const *vars() const { return &vars_[0]; }
};

class LocalSolver : public Solver {
 private:
  ls::LocalSolver solver_;
  int timelimit_;

  int GetTimeLimit(const SolverOption &) const {
    return timelimit_;
  }

  void SetTimeLimit(const SolverOption &opt, int value) {
    if (value <= 0)
      throw InvalidOptionValue(opt, value);
    timelimit_ = value;
  }

 protected:
  void DoSolve(Problem &p);

 public:
  LocalSolver();
};
}

#endif // AMPL_SOLVERS_LOCALSOLVER_H

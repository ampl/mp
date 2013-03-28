/*
 Solver test suite.

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

#ifndef TESTS_SOLVER_TEST_H_
#define TESTS_SOLVER_TEST_H_

#include "solvers/util/solver.h"
#include "tests/expr_builder.h"
#include "gtest/gtest.h"

typedef std::auto_ptr<ampl::BasicSolver> (*SolverFactory)();

// Abstract solver test.
class SolverTest
    : private ampl::Noncopyable,
      public ::testing::TestWithParam<SolverFactory>, public ampl::ExprBuilder {
 protected:
  std::auto_ptr<ampl::BasicSolver> solver_;
  ampl::Variable x;
  ampl::Variable y;
  ampl::Variable z;

  // TODO: get min an max for a solver
  int min() { return -100; }
  int max() { return 100; }

  double Eval(ampl::NumericExpr e, int var1 = 0, int var2 = 0, int var3 = 0);

  double Eval(ampl::LogicalExpr e, int var1 = 0, int var2 = 0, int var3 = 0) {
    // TODO: handle AllDiff specially
    return Eval(AddIf(e, AddNum(1), AddNum(0)), var1, var2, var3);
  }

 public:
  SolverTest()
  : solver_(GetParam()()), x(AddVar(1)), y(AddVar(2)), z(AddVar(3)) {}
};

#endif  // TESTS_SOLVER_TEST_H_

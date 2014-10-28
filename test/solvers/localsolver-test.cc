/*
 LocalSolver tests

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

#include "localsolver/localsolver.h"
#include "feature.h"

typedef mp::LocalSolver Solver;
enum { FEATURES = ~feature::PLTERM };

#define MP_THREAD_INTERRUPT 1

// Demo version of LocalSolver cannot handle a TSP problem with n > 9.
#define MP_TSP_SIZE 9

#include "nl-solver-test.h"

using localsolver::LSParam;

// Creates and solves a test problem.
void SolveTestProblem(mp::LocalSolver &s) {
  mp::LSProblemBuilder pb(s.GetProblemBuilder(""));
  MakeAllDiffProblem(pb);
  TestSolutionHandler sh;
  s.Solve(pb, sh);
}

// Creates a test LocalSolver model.
void CreateTestModel(localsolver::LocalSolver &s) {
  s.getParam().setVerbosity(0);
  localsolver::LSModel model = s.getModel();
  model.addObjective(model.createConstant(localsolver::lsint(0)),
                     localsolver::OD_Minimize);
  model.close();
}

struct IntOption {
  const char *name;
  int default_value;
  int lb;
  int ub;
  int (LSParam::*get)() const;
  void (LSParam::*set)(int);
};

class OptionTest : public testing::TestWithParam<IntOption> {
};

// Test that the default value agrees with LocalSolver.
TEST_P(OptionTest, DefaultValue) {
  auto opt = GetParam();
  localsolver::LocalSolver ls;
  CreateTestModel(ls);
  auto param = ls.getParam();
  EXPECT_EQ((param.*opt.get)(), mp::LocalSolver().GetIntOption(opt.name));
}

// Test that the range agrees with LocalSolver.
TEST_P(OptionTest, Range) {
  auto opt = GetParam();
  localsolver::LocalSolver ls;
  CreateTestModel(ls);
  auto param = ls.getParam();
  mp::LocalSolver solver;
  (param.*opt.set)(opt.lb);
  solver.SetIntOption(opt.name, opt.lb);
  (param.*opt.set)(opt.ub);
  solver.SetIntOption(opt.name, opt.ub);
  EXPECT_THROW((param.*opt.set)(opt.lb - 1), localsolver::LSException);
  EXPECT_THROW(solver.SetIntOption(opt.name, opt.lb - 1),
               mp::InvalidOptionValue);
  if (opt.ub != INT_MAX) {
    EXPECT_THROW((param.*opt.set)(opt.ub + 1), localsolver::LSException);
    EXPECT_THROW(solver.SetIntOption(opt.name, opt.ub + 1),
                 mp::InvalidOptionValue);
  }
}

// Test that the option value is passed to LocalSolver.
TEST_P(OptionTest, PassValue) {
  auto opt = GetParam();
  enum {TEST_VALUE = 7};
  struct TestLocalSolver : mp::LocalSolver {
    IntOption opt;
    void DoSolve(localsolver::LocalSolver &s) {
      EXPECT_EQ(TEST_VALUE, (s.getParam().*opt.get)());
    }
    TestLocalSolver(IntOption opt) : opt(opt) {}
  } solver(opt);
  solver.SetIntOption(opt.name, TEST_VALUE);
  EXPECT_EQ(TEST_VALUE, solver.GetIntOption(opt.name));
  SolveTestProblem(solver);
  EXPECT_THROW(solver.SetStrOption(opt.name, "oops"), mp::OptionError);
}

IntOption options[] = {
  {"seed", 0, 0, INT_MAX, &LSParam::getSeed, &LSParam::setSeed},
  {"threads", 2, 1, 1024, &LSParam::getNbThreads, &LSParam::setNbThreads},
  {"annealing_level", 1, 0, 9,
   &LSParam::getAnnealingLevel, &LSParam::setAnnealingLevel},
};
INSTANTIATE_TEST_CASE_P(, OptionTest, testing::ValuesIn(options));

// TODO: test solver options verbosity,
//       time_between_displays, logfile, timelimit, iterlimit

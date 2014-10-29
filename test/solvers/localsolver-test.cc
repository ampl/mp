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

template <typename T = int>
struct TestLocalSolver : mp::LocalSolver {
  T (LSParam::*get)() const;
  T value;
  void DoSolve(localsolver::LocalSolver &s) {
    value = (s.getParam().*get)();
  }
  TestLocalSolver() : value() {}
};

// Creates and solves a test problem.
template <typename T>
T GetOption(TestLocalSolver<T> &solver, T (LSParam::*get)() const) {
  mp::LSProblemBuilder pb(solver.GetProblemBuilder(""));
  MakeAllDiffProblem(pb);
  TestSolutionHandler sh;
  solver.get = get;
  solver.Solve(pb, sh);
  return solver.value;
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

class OptionTest : public testing::TestWithParam<IntOption> {};

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
  TestLocalSolver<> solver;
  enum {TEST_VALUE = 7};
  solver.SetIntOption(opt.name, TEST_VALUE);
  EXPECT_EQ(TEST_VALUE, solver.GetIntOption(opt.name));
  EXPECT_EQ(TEST_VALUE, GetOption(solver, opt.get));
  EXPECT_THROW(solver.SetStrOption(opt.name, "oops"), mp::OptionError);
}

IntOption options[] = {
  {"seed", 0, 0, INT_MAX, &LSParam::getSeed, &LSParam::setSeed},
  {"threads", 2, 1, 1024, &LSParam::getNbThreads, &LSParam::setNbThreads},
  {"annealing_level", 1, 0, 9,
   &LSParam::getAnnealingLevel, &LSParam::setAnnealingLevel},
  {"time_between_displays", 1, 1, 65535,
   &LSParam::getTimeBetweenDisplays, &LSParam::setTimeBetweenDisplays},
};
INSTANTIATE_TEST_CASE_P(, OptionTest, testing::ValuesIn(options));

TEST(LocalSolverTest, VerbosityOption) {
  TestLocalSolver<> solver;

  // Test default value.
  EXPECT_EQ("quiet", solver.GetStrOption("verbosity"));
  mp::LSProblemBuilder builder(solver.GetProblemBuilder(""));
  EXPECT_EQ(0, builder.solver().getParam().getVerbosity());

  // Test option values.
  auto values = solver.FindOption("verbosity")->values();
  EXPECT_EQ(3, values.size());
  struct ValueInfo {
    const char *value;
    int data;
  };
  ValueInfo value_info[] = {{"quiet", 0}, {"normal", 1}, {"detailed", 10}};
  int index = 0;
  for (auto i = values.begin(), end = values.end(); i != end; ++i) {
    const auto &info = value_info[index++];
    EXPECT_STREQ(info.value, i->value);
    EXPECT_EQ(info.data, i->data);
    solver.SetStrOption("verbosity", info.value);
    EXPECT_EQ(info.value, solver.GetStrOption("verbosity"));
    solver.SetStrOption("verbosity", fmt::format("{}", info.data));
    EXPECT_EQ(info.value, solver.GetStrOption("verbosity"));
  }

  // Test that value is passed to LocalSolver.
  solver.SetStrOption("verbosity", "normal");
  EXPECT_EQ(1, GetOption(solver, &LSParam::getVerbosity));

  EXPECT_THROW(solver.SetDblOption("verbosity", 1.2), mp::OptionError);
  EXPECT_THROW(solver.SetStrOption("verbosity", "oops"),
               mp::InvalidOptionValue);
}

TEST(LocalSolverTest, LogFileOption) {
  TestLocalSolver<std::string> solver;
  EXPECT_EQ("", solver.GetStrOption("logfile"));
  solver.SetStrOption("logfile", "testlog");
  EXPECT_EQ("testlog", GetOption(solver, &LSParam::getLogFile));
  EXPECT_THROW(solver.SetIntOption("logfile", 1), mp::OptionError);
}

// TODO: test solver options timelimit, iterlimit

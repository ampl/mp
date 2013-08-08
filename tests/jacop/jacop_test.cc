/*
 JaCoP solver tests.

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

#include <algorithm>
#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "solvers/jacop/jacop.h"

#include "tests/solver_test.h"
#include "tests/util.h"

using std::string;
using ampl::InvalidOptionValue;

namespace {

// ----------------------------------------------------------------------------
// Solver tests

std::auto_ptr<ampl::BasicSolver> CreateSolver() {
  return std::auto_ptr<ampl::BasicSolver>(new ampl::JaCoPSolver());
}

INSTANTIATE_TEST_CASE_P(JaCoP, SolverTest,
    ::testing::Values(SolverTestParam(CreateSolver, feature::POW)));

TEST_P(SolverTest, SolveFlowshp2) {
  EXPECT_EQ(22, Solve("flowshp2").obj);
}

class JaCoPSolverTest : public ::testing::Test {
 protected:
  ampl::JaCoPSolver solver_;

  SolveResult Solve(const char *stub, const char *opt1 = nullptr,
      const char *opt2 = nullptr, const char *opt3 = nullptr) {
    return SolverTest::Solve(solver_, stub, opt1, opt2, opt3);
  }
};

// ----------------------------------------------------------------------------
// Option tests

TEST_F(JaCoPSolverTest, BacktrackLimitOption) {
  Solve("miplib/assign1", "backtracklimit=42");
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_EQ(42, solver_.GetIntOption("backtracklimit"));
  EXPECT_THROW(solver_.SetIntOption("backtracklimit", -1), InvalidOptionValue);
}

TEST_F(JaCoPSolverTest, DecisionLimitOption) {
  Solve("miplib/assign1", "decisionlimit=42");
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_EQ(42, solver_.GetIntOption("decisionlimit"));
  EXPECT_THROW(solver_.SetIntOption("decisionlimit", -1), InvalidOptionValue);
}

TEST_F(JaCoPSolverTest, FailLimitOption) {
  string message = Solve("miplib/assign1", "faillimit=42").message;
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find(" 43 fails") != string::npos);
  EXPECT_EQ(42, solver_.GetIntOption("faillimit"));
  EXPECT_THROW(solver_.SetIntOption("faillimit", -1), InvalidOptionValue);
}

TEST_F(JaCoPSolverTest, NodeLimitOption) {
  string message = Solve("miplib/assign1", "nodelimit=42").message;
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find("43 nodes") != string::npos);
  EXPECT_EQ(42, solver_.GetIntOption("nodelimit"));
  EXPECT_THROW(solver_.SetIntOption("nodelimit", -1), InvalidOptionValue);
}

TEST_F(JaCoPSolverTest, TimeLimitOption) {
  Solve("miplib/assign1", "timelimit=1");
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_EQ(1, solver_.GetIntOption("timelimit"));
  EXPECT_THROW(solver_.SetIntOption("timelimit", -1), InvalidOptionValue);
}

const char *const VAL_SELECT[] = {
  "IndomainMax",
  "IndomainMedian",
  "IndomainMiddle",
  "IndomainMin",
  "IndomainRandom",
  "IndomainSimpleRandom",
  0
};

TEST_F(JaCoPSolverTest, ValSelectOption) {
  EXPECT_EQ("indomainmin", solver_.GetStrOption("val_select"));
  unsigned count = 0;
  for (const char *const *s = VAL_SELECT; *s; ++s, ++count) {
    std::string value = *s;
    std::transform(value.begin(), value.end(), value.begin(), ::tolower);
    solver_.SetStrOption("val_select", value.c_str());
    EXPECT_EQ(value, solver_.GetStrOption("val_select"));
  }
  EXPECT_EQ(6u, count);
}

const char *const VAR_SELECT[] = {
  "LargestDomain",
  "LargestMax",
  "LargestMin",
  "MaxRegret",
  "MinDomainOverDegree",
  "MostConstrainedDynamic",
  "MostConstrainedStatic",
  "SmallestDomain",
  "SmallestMax",
  "SmallestMin",
  "WeightedDegree",
  0
};

TEST_F(JaCoPSolverTest, VarSelectOption) {
  EXPECT_EQ("smallestdomain", solver_.GetStrOption("var_select"));
  unsigned count = 0;
  for (const char *const *s = VAR_SELECT; *s; ++s, ++count) {
    std::string value = *s;
    std::transform(value.begin(), value.end(), value.begin(), ::tolower);
    solver_.SetStrOption("var_select", value.c_str());
    EXPECT_EQ(value, solver_.GetStrOption("var_select"));
  }
  EXPECT_EQ(11u, count);
}

/*TEST_F(JaCoPSolverTest, OutLevOption) {
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve("objconstint");
    fclose(f);
    _exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("", ReadFile("out"));

  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve("objconstint", "outlev=1");
    fclose(f);
    _exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("outlev=1\n"
      " Max Depth      Nodes      Fails      Best Obj\n"
      "                                            42\n", ReadFile("out"));

  solver_.SetIntOption("outlev", 0);
  EXPECT_EQ(0, solver_.GetIntOption("outlev"));
  solver_.SetIntOption("outlev", 1);
  EXPECT_EQ(1, solver_.GetIntOption("outlev"));
  EXPECT_THROW(solver_.SetIntOption("outlev", -1), InvalidOptionValue);
  EXPECT_THROW(solver_.SetIntOption("outlev", 2), InvalidOptionValue);
}

TEST_F(JaCoPSolverTest, OutFreqOption) {
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve("party1", "outlev=1", "outfreq=0.4", "timelimit=1");
    fclose(f);
    _exit(0);
  }, ::testing::ExitedWithCode(0), "");
  string out = ReadFile("out");
  EXPECT_EQ(6, std::count(out.begin(), out.end(), '\n'));

  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve("party1", "outlev=1", "outfreq=0.8", "timelimit=1");
    fclose(f);
    _exit(0);
  }, ::testing::ExitedWithCode(0), "");
  out = ReadFile("out");
  EXPECT_EQ(5, std::count(out.begin(), out.end(), '\n'));

  solver_.SetDblOption("outfreq", 1.23);
  EXPECT_EQ(1.23, solver_.GetDblOption("outfreq"));
  EXPECT_THROW(solver_.SetDblOption("outfreq", -1), InvalidOptionValue);
  EXPECT_THROW(solver_.SetDblOption("outfreq", 0), InvalidOptionValue);
}*/
}

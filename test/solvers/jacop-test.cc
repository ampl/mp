/*
 JaCoP solver tests.

 Copyright (C) 2012 AMPL Optimization Inc

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

#include <algorithm>
#include <memory>
#include <string>

#include "jacop/jacop.h"
#include "feature.h"
#include "../util.h"

typedef mp::JaCoPSolver Solver;
enum {FEATURES = feature::POW};

#include "solver-impl-test.h"

using std::string;
using mp::InvalidOptionValue;
using mp::Problem;

// ----------------------------------------------------------------------------
// Option tests

TEST_F(SolverImplTest, BacktrackLimitOption) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeTSP(pb);
  solver_.SetIntOption("backtracklimit", 42);
  EXPECT_EQ(42, solver_.GetIntOption("backtracklimit"));
  EXPECT_EQ(400, Solve(pb).solve_code());
  EXPECT_THROW(solver_.SetIntOption("backtracklimit", -1), InvalidOptionValue);
}

TEST_F(SolverImplTest, DecisionLimitOption) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeTSP(pb);
  solver_.SetIntOption("decisionlimit", 42);
  EXPECT_EQ(400, Solve(pb).solve_code());
  EXPECT_EQ(42, solver_.GetIntOption("decisionlimit"));
  EXPECT_THROW(solver_.SetIntOption("decisionlimit", -1), InvalidOptionValue);
}

TEST_F(SolverImplTest, FailLimitOption) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeTSP(pb);
  solver_.SetIntOption("faillimit", 42);
  TestSolutionHandler sh;
  solver_.Solve(pb, sh);
  EXPECT_EQ(400, sh.status());
  EXPECT_TRUE(sh.message().find(" 43 fails") != string::npos);
  EXPECT_EQ(42, solver_.GetIntOption("faillimit"));
  EXPECT_THROW(solver_.SetIntOption("faillimit", -1), InvalidOptionValue);
}

TEST_F(SolverImplTest, NodeLimitOption) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeTSP(pb);
  solver_.SetIntOption("nodelimit", 42);
  TestSolutionHandler sh;
  solver_.Solve(pb, sh);
  EXPECT_EQ(400, sh.status());
  EXPECT_TRUE(sh.message().find("43 nodes") != string::npos);
  EXPECT_EQ(42, solver_.GetIntOption("nodelimit"));
  EXPECT_THROW(solver_.SetIntOption("nodelimit", -1), InvalidOptionValue);
}

TEST_F(SolverImplTest, TimeLimitOption) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeTSP(pb);
  solver_.SetIntOption("timelimit", 1);
  EXPECT_EQ(400, Solve(pb).solve_code());
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

TEST_F(SolverImplTest, ValSelectOption) {
  EXPECT_EQ("indomainmin", solver_.GetStrOption("val_select"));
  unsigned count = 0;
  for (const char *const *s = VAL_SELECT; *s; ++s, ++count) {
    string value = *s;
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

TEST_F(SolverImplTest, VarSelectOption) {
  EXPECT_EQ("smallestdomain", solver_.GetStrOption("var_select"));
  unsigned count = 0;
  for (const char *const *s = VAR_SELECT; *s; ++s, ++count) {
    string value = *s;
    std::transform(value.begin(), value.end(), value.begin(), ::tolower);
    solver_.SetStrOption("var_select", value.c_str());
    EXPECT_EQ(value, solver_.GetStrOption("var_select"));
  }
  EXPECT_EQ(11u, count);
}

struct TestOutputHandler : public mp::OutputHandler {
  string output;
  void HandleOutput(fmt::StringRef output) { this->output += output; }
};

TEST_F(SolverImplTest, OutLevOption) {
  TestOutputHandler h;
  solver_.set_output_handler(&h);
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/objconstint.nl");
  mp::BasicSolutionHandler sh;
  solver_.Solve(p, sh);
  EXPECT_EQ("", h.output);

  h.output.clear();
  solver_.SetIntOption("outlev", 1);
  solver_.Solve(p, sh);
  EXPECT_EQ(
      " Max Depth      Nodes      Fails      Best Obj\n"
      "                                            42\n", h.output);

  solver_.SetIntOption("outlev", 0);
  EXPECT_EQ(0, solver_.GetIntOption("outlev"));
  solver_.SetIntOption("outlev", 1);
  EXPECT_EQ(1, solver_.GetIntOption("outlev"));
  EXPECT_THROW(solver_.SetIntOption("outlev", -1), InvalidOptionValue);
  EXPECT_THROW(solver_.SetIntOption("outlev", 2), InvalidOptionValue);
}

TEST_F(SolverImplTest, OutFreqOption) {
  TestOutputHandler h;
  solver_.set_output_handler(&h);
  Problem p;
  p.Read(MP_TEST_DATA_DIR "/party1.nl");
  solver_.SetIntOption("outlev", 1);
  solver_.SetIntOption("timelimit", 1);

  solver_.SetDblOption("outfreq", 0.4);
  mp::BasicSolutionHandler sh;
  solver_.Solve(p, sh);
  string out = h.output;
  EXPECT_EQ(3, std::count(out.begin(), out.end(), '\n'));

  h.output.clear();
  solver_.SetDblOption("outfreq", 0.8);
  solver_.Solve(p, sh);
  out = h.output;
  EXPECT_EQ(2, std::count(out.begin(), out.end(), '\n'));

  solver_.SetDblOption("outfreq", 1.23);
  EXPECT_EQ(1.23, solver_.GetDblOption("outfreq"));
  EXPECT_THROW(solver_.SetDblOption("outfreq", -1), InvalidOptionValue);
  EXPECT_THROW(solver_.SetDblOption("outfreq", 0), InvalidOptionValue);
}

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

#include <csignal>
#include <cstdlib>

#include "gtest/gtest.h"

#include "solvers/jacop/jacop.h"
#include "solvers/util/expr.h"

extern "C" {
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
}

#include "tests/args.h"
#include "tests/expr_builder.h"
#include "tests/solution_handler.h"
#include "tests/solver_test.h"
#include "tests/util.h"
#include "tests/config.h"

#ifdef HAVE_THREADS
# include <thread>
#endif

using std::string;

namespace {

// ----------------------------------------------------------------------------
// Solver tests

std::auto_ptr<ampl::BasicSolver> CreateSolver() {
  return std::auto_ptr<ampl::BasicSolver>(new ampl::JaCoPSolver());
}

INSTANTIATE_TEST_CASE_P(JaCoP, SolverTest,
    ::testing::Values(SolverTestParam(CreateSolver, feature::POW)));

// TODO: fix stackoverflow
/*TEST_P(SolverTest, SolveBalassign0) {
  EXPECT_EQ(14, Solve("balassign0").obj);
}

TEST_P(SolverTest, SolveBalassign1) {
  EXPECT_EQ(14, Solve("balassign1").obj);
}*/

TEST_P(SolverTest, SolveFlowshp2) {
  EXPECT_EQ(22, Solve("flowshp2").obj);
}

/*class JaCoPSolverTest : public ::testing::Test {
 protected:
  ampl::JaCoPSolver solver_;

  class ParseResult {
   private:
    bool result_;
    std::string error_;

    void True() const {}
    typedef void (ParseResult::*SafeBool)() const;

   public:
    ParseResult(bool result, std::string error)
    : result_(result), error_(error) {}

    operator SafeBool() const { return result_ ? &ParseResult::True : 0; }

    std::string error() const {
      EXPECT_FALSE(result_);
      return error_;
    }
  };

  struct TestErrorHandler : ampl::ErrorHandler {
    std::string error;
    void HandleError(fmt::StringRef message) {
      error += message.c_str();
    }
  };

  ParseResult ParseOptions(const char *opt1, const char *opt2 = nullptr) {
    TestErrorHandler eh;
    solver_.set_error_handler(&eh);
    bool result = solver_.ParseOptions(
        Args(opt1, opt2), solver_, ampl::BasicSolver::NO_OPTION_ECHO);
    if (result)
      EXPECT_EQ("", eh.error);
    return ParseResult(result, eh.error);
  }
};

// ----------------------------------------------------------------------------
// Solve code tests

TEST_F(JaCoPSolverTest, OptimalSolveCode) {
  Solve("objconstint");
  EXPECT_EQ(0, solver_.problem().solve_code());
}

TEST_F(JaCoPSolverTest, FeasibleSolveCode) {
  Solve("feasible");
  EXPECT_EQ(100, solver_.problem().solve_code());
}

TEST_F(JaCoPSolverTest, InfeasibleSolveCode) {
  Solve("infeasible");
  EXPECT_EQ(200, solver_.problem().solve_code());
}

// ----------------------------------------------------------------------------
// Interrupt tests

#ifdef HAVE_THREADS
void Interrupt() {
  // Wait until started.
  while (ampl::SignalHandler::stop())
    std::this_thread::yield();
  std::raise(SIGINT);
}

TEST_F(JaCoPSolverTest, InterruptSolution) {
  std::thread t(Interrupt);
  std::string message = Solve("miplib/assign1").message;
  t.join();
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}
#endif

// ----------------------------------------------------------------------------
// Option tests

TEST_F(JaCoPSolverTest, ADOption) {
  EXPECT_EQ(Gecode::Search::Options().a_d, solver_.options().a_d);
  EXPECT_TRUE(ParseOptions("a_d=42"));
  EXPECT_EQ(42u, solver_.options().a_d);
  EXPECT_EQ("Invalid value -1 for option a_d", ParseOptions("a_d=-1").error());
}

TEST_F(JaCoPSolverTest, CDOption) {
  EXPECT_EQ(Gecode::Search::Options().c_d, solver_.options().c_d);
  EXPECT_TRUE(ParseOptions("c_d=42"));
  EXPECT_EQ(42u, solver_.options().c_d);
  EXPECT_EQ("Invalid value -1 for option c_d", ParseOptions("c_d=-1").error());
}

TEST_F(JaCoPSolverTest, FailLimitOption) {
  std::string message =
      Solve("miplib/assign1", "faillimit=10").message;
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find(" 11 fails") != string::npos);
  EXPECT_EQ("Invalid value -1 for option faillimit",
      ParseOptions("faillimit=-1").error());
}

TEST_F(JaCoPSolverTest, MemoryLimitOption) {
  Solve("miplib/assign1", "memorylimit=1000000");
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_EQ("Invalid value -1 for option memorylimit",
      ParseOptions("memorylimit=-1").error());
}

TEST_F(JaCoPSolverTest, NodeLimitOption) {
  std::string message =
      Solve("miplib/assign1", "nodelimit=10").message;
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find("11 nodes") != string::npos);
  EXPECT_EQ("Invalid value -1 for option nodelimit",
      ParseOptions("nodelimit=-1").error());
}

TEST_F(JaCoPSolverTest, TimeLimitOption) {
  Solve("miplib/assign1", "timelimit=0.1");
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_EQ("Invalid value -1 for option timelimit",
      ParseOptions("timelimit=-1").error());
}

TEST_F(JaCoPSolverTest, ThreadsOption) {
  EXPECT_EQ(Gecode::Search::Options().threads, solver_.options().threads);
  EXPECT_TRUE(ParseOptions("threads=0.5"));
  EXPECT_EQ(0.5, solver_.options().threads);
  EXPECT_TRUE(ParseOptions("threads=-10"));
  EXPECT_EQ(-10.0, solver_.options().threads);
}

template <typename T>
struct OptionValue {
  const char *name;
  T value;
};

const OptionValue<Gecode::IntValBranch> VAL_BRANCHINGS[] = {
    {"min",        Gecode::INT_VAL_MIN},
    {"med",        Gecode::INT_VAL_MED},
    {"max",        Gecode::INT_VAL_MAX},
    {"rnd",        Gecode::INT_VAL_RND},
    {"split_min",  Gecode::INT_VAL_SPLIT_MIN},
    {"split_max",  Gecode::INT_VAL_SPLIT_MAX},
    {"range_min",  Gecode::INT_VAL_RANGE_MIN},
    {"range_max",  Gecode::INT_VAL_RANGE_MAX},
    {"values_min", Gecode::INT_VALUES_MIN},
    {"values_max", Gecode::INT_VALUES_MAX},
    {}
};

TEST_F(JaCoPSolverTest, ValBranchingOption) {
  EXPECT_EQ(Gecode::INT_VAL_MIN, solver_.val_branching());
  unsigned count = 0;
  for (const OptionValue<Gecode::IntValBranch>
      *p = VAL_BRANCHINGS; p->name; ++p, ++count) {
    EXPECT_TRUE(ParseOptions(
        c_str(fmt::Format("val_branching={}") << p->name)));
    EXPECT_EQ(p->value, solver_.val_branching());
  }
  EXPECT_EQ(10u, count);
}

const OptionValue<Gecode::IntVarBranch> VAR_BRANCHINGS[] = {
    {"none",            Gecode::INT_VAR_NONE},
    {"rnd",             Gecode::INT_VAR_RND},
    {"degree_min",      Gecode::INT_VAR_DEGREE_MIN},
    {"degree_max",      Gecode::INT_VAR_DEGREE_MAX},
    {"afc_min",         Gecode::INT_VAR_AFC_MIN},
    {"afc_max",         Gecode::INT_VAR_AFC_MAX},
    {"min_min",         Gecode::INT_VAR_MIN_MIN},
    {"min_max",         Gecode::INT_VAR_MIN_MAX},
    {"max_min",         Gecode::INT_VAR_MAX_MIN},
    {"max_max",         Gecode::INT_VAR_MAX_MAX},
    {"size_min",        Gecode::INT_VAR_SIZE_MIN},
    {"size_max",        Gecode::INT_VAR_SIZE_MAX},
    {"size_degree_min", Gecode::INT_VAR_SIZE_DEGREE_MIN},
    {"size_degree_max", Gecode::INT_VAR_SIZE_DEGREE_MAX},
    {"size_afc_min",    Gecode::INT_VAR_SIZE_AFC_MIN},
    {"size_afc_max",    Gecode::INT_VAR_SIZE_AFC_MAX},
    {"regret_min_min",  Gecode::INT_VAR_REGRET_MIN_MIN},
    {"regret_min_max",  Gecode::INT_VAR_REGRET_MIN_MAX},
    {"regret_max_min",  Gecode::INT_VAR_REGRET_MAX_MIN},
    {"regret_max_max",  Gecode::INT_VAR_REGRET_MAX_MAX},
    {}
};

TEST_F(JaCoPSolverTest, VarBranchingOption) {
  EXPECT_EQ(Gecode::INT_VAR_SIZE_MIN, solver_.var_branching());
  unsigned count = 0;
  for (const OptionValue<Gecode::IntVarBranch>
      *p = VAR_BRANCHINGS; p->name; ++p, ++count) {
    EXPECT_TRUE(ParseOptions(
        c_str(fmt::Format("var_branching={}") << p->name)));
    EXPECT_EQ(p->value, solver_.var_branching());
  }
  EXPECT_EQ(20u, count);
}

TEST_F(JaCoPSolverTest, OutLevOption) {
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve("objconstint");
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("", ReadFile("out"));

  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve("objconstint", "outlev=1");
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("outlev=1\n"
      " Max Depth      Nodes      Fails      Best Obj\n"
      "                                            42\n", ReadFile("out"));

  EXPECT_TRUE(ParseOptions("outlev=0"));
  EXPECT_TRUE(ParseOptions("outlev=1"));
  EXPECT_EQ("Invalid value -1 for option outlev",
      ParseOptions("outlev=-1").error());
  EXPECT_EQ("Invalid value 2 for option outlev",
      ParseOptions("outlev=2").error());
}

TEST_F(JaCoPSolverTest, OutFreqOption) {
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve("party1", "outlev=1", "outfreq=1", "timelimit=2.5");
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  std::string out = ReadFile("out");
  EXPECT_EQ(6, std::count(out.begin(), out.end(), '\n'));

  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve("party1", "outlev=1", "outfreq=2", "timelimit=2.5");
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  out = ReadFile("out");
  EXPECT_EQ(5, std::count(out.begin(), out.end(), '\n'));

  EXPECT_EQ("Invalid value -1 for option outfreq",
      ParseOptions("outfreq=-1").error());
  EXPECT_EQ("Invalid value 0 for option outfreq",
      ParseOptions("outfreq=0").error());
}*/
}

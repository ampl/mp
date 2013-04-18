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
#include <fstream>
#include <memory>
#include <string>
#include <vector>

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

#if defined(_MSC_VER)
# define isnan _isnan
#else
# define isnan std::isnan
#endif

using std::size_t;
using std::string;
using std::vector;

using ampl::Class;
using ampl::ExprBuilder;
using ampl::NLToJaCoPConverter;

#define DATA_DIR "../data/"

namespace {

std::auto_ptr<ampl::BasicSolver> CreateSolver() {
  return std::auto_ptr<ampl::BasicSolver>(new ampl::JaCoPSolver());
}

INSTANTIATE_TEST_CASE_P(JaCoP, SolverTest, ::testing::Values(CreateSolver));

class JaCoPConverterTest : public ::testing::Test, public ExprBuilder {
 protected:
  ampl::Env env_;
  int min_int_;
  int max_int_;

  JaCoPConverterTest() : env_(ampl::JVM::env()) {
    jclass domain_class = env_.FindClass("JaCoP/core/IntDomain");
    min_int_ = env_.GetStaticIntField(
        domain_class, env_.GetStaticFieldID(domain_class, "MinInt", "I"));
    max_int_ = env_.GetStaticIntField(
        domain_class, env_.GetStaticFieldID(domain_class, "MaxInt", "I"));
  }

  double min() const { return min_int_; }
  double max() const { return max_int_; }
};

// ----------------------------------------------------------------------------
// Solver tests

// TODO
/*class JaCoPSolverTest : public ::testing::Test, public ExprBuilder {
 protected:
  ampl::GecodeSolver solver_;

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

  int RunSolver(const char *stub = nullptr, const char *opt1 = nullptr,
      const char *opt2 = nullptr, const char *opt3 = nullptr) {
    return solver_.Run(Args("gecode", "-s", stub, opt1, opt2, opt3));
  }

  SolveResult Solve(const char *stub, const char *opt1 = nullptr,
      const char *opt2 = nullptr, const char *opt3 = nullptr);
};

SolveResult JaCoPSolverTest::Solve(const char *stub,
    const char *opt1, const char *opt2, const char *opt3) {
  TestSolutionHandler sh;
  solver_.set_solution_handler(&sh);
  RunSolver(stub, opt1, opt2, opt3);
  const string &message = sh.message();
  int solve_code = sh.solve_code();
  EXPECT_GE(solve_code, 0);
  bool solved = true;
  if (solve_code < 100)
    EXPECT_TRUE(message.find("optimal solution") != string::npos);
  else if (solve_code < 200)
    EXPECT_TRUE(message.find("feasible solution") != string::npos);
  else
    solved = false;
  return SolveResult(solved, sh.obj_value(), message);
}

TEST_F(JaCoPSolverTest, ContinuousVarsNotSupported) {
  EXPECT_THROW(RunSolver(DATA_DIR "objconst"), std::runtime_error);
}

TEST_F(JaCoPSolverTest, SolveAssign0) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign0").obj);
}

TEST_F(JaCoPSolverTest, SolveAssign1) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign1").obj);
}

// Disabled because it takes too long to solve.
TEST_F(JaCoPSolverTest, DISABLED_SolveBalassign0) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign0").obj);
}

// Disabled because it takes too long to solve.
TEST_F(JaCoPSolverTest, DISABLED_SolveBalassign1) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign1").obj);
}

TEST_F(JaCoPSolverTest, SolveFlowshp0) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp0").obj);
}

TEST_F(JaCoPSolverTest, SolveFlowshp1) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp1").obj);
}

TEST_F(JaCoPSolverTest, SolveFlowshp2) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp2").obj);
}

TEST_F(JaCoPSolverTest, SolveGrpassign0) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign0").obj);
}

// Disabled because variables in subscripts are not yet allowed.
TEST_F(JaCoPSolverTest, DISABLED_SolveGrpassign1) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1").obj);
}

// Disabled because object-valued variables are not yet allowed.
TEST_F(JaCoPSolverTest, DISABLED_SolveGrpassign1a) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1a").obj);
}

TEST_F(JaCoPSolverTest, SolveMagic) {
  EXPECT_TRUE(Solve(DATA_DIR "magic").solved);
}

TEST_F(JaCoPSolverTest, SolveMapcoloring) {
  EXPECT_TRUE(Solve(DATA_DIR "mapcoloring").solved);
}

TEST_F(JaCoPSolverTest, SolveNQueens) {
  EXPECT_TRUE(Solve(DATA_DIR "nqueens").solved);
}

TEST_F(JaCoPSolverTest, SolveNQueens0) {
  EXPECT_TRUE(Solve(DATA_DIR "nqueens0").solved);
}

// Disabled because it takes somewhat long (compared to other tests).
TEST_F(JaCoPSolverTest, DISABLED_SolveOpenShop) {
  EXPECT_EQ(1955, Solve(DATA_DIR "openshop").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(JaCoPSolverTest, DISABLED_SolveParty1) {
  EXPECT_EQ(61, Solve(DATA_DIR "party1").obj);
}

// Disabled because Gecode doesn't support 'alldiff' as a subexpression.
TEST_F(JaCoPSolverTest, DISABLED_SolveParty2) {
  EXPECT_EQ(3, Solve(DATA_DIR "party2").obj);
}

// Disabled because it takes somewhat long (compared to other tests).
TEST_F(JaCoPSolverTest, DISABLED_SolvePhoto9) {
  EXPECT_EQ(10, Solve(DATA_DIR "photo9").obj);
}

// Disabled because it takes somewhat long (compared to other tests).
TEST_F(JaCoPSolverTest, DISABLED_SolvePhoto11) {
  EXPECT_EQ(12, Solve(DATA_DIR "photo11").obj);
}

TEST_F(JaCoPSolverTest, SolveSched0) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched0").obj);
}

TEST_F(JaCoPSolverTest, SolveSched1) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched1").obj);
}

TEST_F(JaCoPSolverTest, SolveSched2) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched2").obj);
}

TEST_F(JaCoPSolverTest, SolveSendMoreMoney) {
  EXPECT_TRUE(Solve(DATA_DIR "send-more-money").solved);
}

TEST_F(JaCoPSolverTest, SolveSendMostMoney) {
  EXPECT_NEAR(10876, Solve(DATA_DIR "send-most-money").obj, 1e-5);
}

TEST_F(JaCoPSolverTest, SolveSeq0) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0").obj, 1e-5);
}

TEST_F(JaCoPSolverTest, SolveSeq0a) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0a").obj, 1e-5);
}

TEST_F(JaCoPSolverTest, SolveSudokuHard) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuHard").solved);
}

TEST_F(JaCoPSolverTest, SolveSudokuVeryEasy) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuVeryEasy").solved);
}

// ----------------------------------------------------------------------------
// Solve code tests

TEST_F(JaCoPSolverTest, OptimalSolveCode) {
  Solve(DATA_DIR "objconstint");
  EXPECT_EQ(0, solver_.problem().solve_code());
}

TEST_F(JaCoPSolverTest, FeasibleSolveCode) {
  Solve(DATA_DIR "feasible");
  EXPECT_EQ(100, solver_.problem().solve_code());
}

TEST_F(JaCoPSolverTest, InfeasibleSolveCode) {
  Solve(DATA_DIR "infeasible");
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
  std::string message = Solve(DATA_DIR "miplib/assign1").message;
  t.join();
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}
#endif

// ----------------------------------------------------------------------------
// Option tests

TEST_F(JaCoPSolverTest, VersionOption) {
  EXPECT_FALSE((solver_.flags() & ASL_OI_show_version) != 0);
  EXPECT_TRUE(ParseOptions("version"));
  EXPECT_TRUE((solver_.flags() & ASL_OI_show_version) != 0);
}

TEST_F(JaCoPSolverTest, WantsolOption) {
  EXPECT_EQ(0, solver_.wantsol());
  EXPECT_TRUE(ParseOptions("wantsol=1"));
  EXPECT_EQ(1, solver_.wantsol());
  EXPECT_TRUE(ParseOptions("wantsol=5"));
  EXPECT_EQ(5, solver_.wantsol());
}

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
      Solve(DATA_DIR "miplib/assign1", "faillimit=10").message;
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find(" 11 fails") != string::npos);
  EXPECT_EQ("Invalid value -1 for option faillimit",
      ParseOptions("faillimit=-1").error());
}

TEST_F(JaCoPSolverTest, MemoryLimitOption) {
  Solve(DATA_DIR "miplib/assign1", "memorylimit=1000000");
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_EQ("Invalid value -1 for option memorylimit",
      ParseOptions("memorylimit=-1").error());
}

TEST_F(JaCoPSolverTest, NodeLimitOption) {
  std::string message =
      Solve(DATA_DIR "miplib/assign1", "nodelimit=10").message;
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find("11 nodes") != string::npos);
  EXPECT_EQ("Invalid value -1 for option nodelimit",
      ParseOptions("nodelimit=-1").error());
}

TEST_F(JaCoPSolverTest, TimeLimitOption) {
  Solve(DATA_DIR "miplib/assign1", "timelimit=0.1");
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
    Solve(DATA_DIR "objconstint");
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("", ReadFile("out"));

  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve(DATA_DIR "objconstint", "outlev=1");
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
    Solve(DATA_DIR "party1", "outlev=1", "outfreq=1", "timelimit=2.5");
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  std::string out = ReadFile("out");
  EXPECT_EQ(6, std::count(out.begin(), out.end(), '\n'));

  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve(DATA_DIR "party1", "outlev=1", "outfreq=2", "timelimit=2.5");
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

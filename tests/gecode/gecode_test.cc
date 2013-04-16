/*
 Gecode solver tests.

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

#include "solvers/gecode/gecode.h"
#include "solvers/util/expr.h"

extern "C" {
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
}

#include "tests/args.h"
#include "tests/expr_builder.h"
#include "tests/solver_test.h"
#include "tests/solution_handler.h"
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

using std::ifstream;
using std::size_t;
using std::string;
using std::vector;

using ampl::Variable;
using ampl::ExprBuilder;
using ampl::GecodeProblem;
using ampl::LogicalExpr;
using ampl::NumericExpr;
using ampl::UnsupportedExprError;

#define DATA_DIR "../data/"

namespace {

// ----------------------------------------------------------------------------
// Solver tests

std::auto_ptr<ampl::BasicSolver> CreateSolver() {
  return std::auto_ptr<ampl::BasicSolver>(new ampl::GecodeSolver());
}

INSTANTIATE_TEST_CASE_P(Gecode, SolverTest, ::testing::Values(CreateSolver));

TEST_P(SolverTest, FloorSqrt) {
  EXPECT_EQ(6, Eval(AddUnary(FLOOR, AddUnary(OP_sqrt, x)), 42));
}

class GecodeSolverTest : public ::testing::Test, public ExprBuilder {
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

SolveResult GecodeSolverTest::Solve(const char *stub,
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

TEST_F(GecodeSolverTest, ObjConst) {
  EXPECT_EQ(42, Solve(DATA_DIR "objconstint").obj);
}

TEST_F(GecodeSolverTest, ContinuousVarsNotSupported) {
  EXPECT_THROW(RunSolver(DATA_DIR "objconst"), std::runtime_error);
}

TEST_F(GecodeSolverTest, SolveAssign0) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign0").obj);
}

TEST_F(GecodeSolverTest, SolveAssign1) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign1").obj);
}

// Disabled because it takes too long to solve.
TEST_F(GecodeSolverTest, DISABLED_SolveBalassign0) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign0").obj);
}

// Disabled because it takes too long to solve.
TEST_F(GecodeSolverTest, DISABLED_SolveBalassign1) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign1").obj);
}

TEST_F(GecodeSolverTest, SolveFlowshp0) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp0").obj);
}

TEST_F(GecodeSolverTest, SolveFlowshp1) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp1").obj);
}

TEST_F(GecodeSolverTest, SolveFlowshp2) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp2").obj);
}

TEST_F(GecodeSolverTest, SolveGrpassign0) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign0").obj);
}

// Disabled because variables in subscripts are not yet allowed.
TEST_F(GecodeSolverTest, DISABLED_SolveGrpassign1) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1").obj);
}

// Disabled because object-valued variables are not yet allowed.
TEST_F(GecodeSolverTest, DISABLED_SolveGrpassign1a) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1a").obj);
}

TEST_F(GecodeSolverTest, SolveMagic) {
  EXPECT_TRUE(Solve(DATA_DIR "magic").solved);
}

TEST_F(GecodeSolverTest, SolveMapcoloring) {
  EXPECT_TRUE(Solve(DATA_DIR "mapcoloring").solved);
}

TEST_F(GecodeSolverTest, SolveNQueens) {
  EXPECT_TRUE(Solve(DATA_DIR "nqueens").solved);
}

TEST_F(GecodeSolverTest, SolveNQueens0) {
  EXPECT_TRUE(Solve(DATA_DIR "nqueens0").solved);
}

// Disabled because it takes somewhat long (compared to other tests).
TEST_F(GecodeSolverTest, DISABLED_SolveOpenShop) {
  EXPECT_EQ(1955, Solve(DATA_DIR "openshop").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(GecodeSolverTest, DISABLED_SolveParty1) {
  EXPECT_EQ(61, Solve(DATA_DIR "party1").obj);
}

// Disabled because Gecode doesn't support 'alldiff' as a subexpression.
TEST_F(GecodeSolverTest, DISABLED_SolveParty2) {
  EXPECT_EQ(3, Solve(DATA_DIR "party2").obj);
}

// Disabled because it takes somewhat long (compared to other tests).
TEST_F(GecodeSolverTest, DISABLED_SolvePhoto9) {
  EXPECT_EQ(10, Solve(DATA_DIR "photo9").obj);
}

// Disabled because it takes somewhat long (compared to other tests).
TEST_F(GecodeSolverTest, DISABLED_SolvePhoto11) {
  EXPECT_EQ(12, Solve(DATA_DIR "photo11").obj);
}

TEST_F(GecodeSolverTest, SolveSched0) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched0").obj);
}

TEST_F(GecodeSolverTest, SolveSched1) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched1").obj);
}

TEST_F(GecodeSolverTest, SolveSched2) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched2").obj);
}

TEST_F(GecodeSolverTest, SolveSendMoreMoney) {
  EXPECT_TRUE(Solve(DATA_DIR "send-more-money").solved);
}

TEST_F(GecodeSolverTest, SolveSendMostMoney) {
  EXPECT_NEAR(10876, Solve(DATA_DIR "send-most-money").obj, 1e-5);
}

TEST_F(GecodeSolverTest, SolveSeq0) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0").obj, 1e-5);
}

TEST_F(GecodeSolverTest, SolveSeq0a) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0a").obj, 1e-5);
}

TEST_F(GecodeSolverTest, SolveSudokuHard) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuHard").solved);
}

TEST_F(GecodeSolverTest, SolveSudokuVeryEasy) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuVeryEasy").solved);
}

// ----------------------------------------------------------------------------
// Solve code tests

TEST_F(GecodeSolverTest, OptimalSolveCode) {
  Solve(DATA_DIR "objconstint");
  EXPECT_EQ(0, solver_.problem().solve_code());
}

TEST_F(GecodeSolverTest, FeasibleSolveCode) {
  Solve(DATA_DIR "feasible");
  EXPECT_EQ(100, solver_.problem().solve_code());
}

TEST_F(GecodeSolverTest, InfeasibleSolveCode) {
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

TEST_F(GecodeSolverTest, InterruptSolution) {
  std::thread t(Interrupt);
  std::string message = Solve(DATA_DIR "miplib/assign1").message;
  t.join();
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}
#endif

// ----------------------------------------------------------------------------
// Option tests

TEST_F(GecodeSolverTest, VersionOption) {
  EXPECT_FALSE((solver_.flags() & ASL_OI_show_version) != 0);
  EXPECT_TRUE(ParseOptions("version"));
  EXPECT_TRUE((solver_.flags() & ASL_OI_show_version) != 0);
}

TEST_F(GecodeSolverTest, WantsolOption) {
  EXPECT_EQ(0, solver_.wantsol());
  EXPECT_TRUE(ParseOptions("wantsol=1"));
  EXPECT_EQ(1, solver_.wantsol());
  EXPECT_TRUE(ParseOptions("wantsol=5"));
  EXPECT_EQ(5, solver_.wantsol());
}

TEST_F(GecodeSolverTest, ADOption) {
  EXPECT_EQ(Gecode::Search::Options().a_d, solver_.options().a_d);
  EXPECT_TRUE(ParseOptions("a_d=42"));
  EXPECT_EQ(42u, solver_.options().a_d);
  EXPECT_EQ("Invalid value -1 for option a_d", ParseOptions("a_d=-1").error());
}

TEST_F(GecodeSolverTest, CDOption) {
  EXPECT_EQ(Gecode::Search::Options().c_d, solver_.options().c_d);
  EXPECT_TRUE(ParseOptions("c_d=42"));
  EXPECT_EQ(42u, solver_.options().c_d);
  EXPECT_EQ("Invalid value -1 for option c_d", ParseOptions("c_d=-1").error());
}

TEST_F(GecodeSolverTest, FailLimitOption) {
  std::string message =
      Solve(DATA_DIR "miplib/assign1", "faillimit=10").message;
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find(" 11 fails") != string::npos);
  EXPECT_EQ("Invalid value -1 for option faillimit",
      ParseOptions("faillimit=-1").error());
}

TEST_F(GecodeSolverTest, MemoryLimitOption) {
  Solve(DATA_DIR "miplib/assign1", "memorylimit=100000");
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_EQ("Invalid value -1 for option memorylimit",
      ParseOptions("memorylimit=-1").error());
}

TEST_F(GecodeSolverTest, NodeLimitOption) {
  std::string message =
      Solve(DATA_DIR "miplib/assign1", "nodelimit=10").message;
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find("11 nodes") != string::npos);
  EXPECT_EQ("Invalid value -1 for option nodelimit",
      ParseOptions("nodelimit=-1").error());
}

TEST_F(GecodeSolverTest, TimeLimitOption) {
  Solve(DATA_DIR "miplib/assign1", "timelimit=0.1");
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_EQ("Invalid value -1 for option timelimit",
      ParseOptions("timelimit=-1").error());
}

TEST_F(GecodeSolverTest, ThreadsOption) {
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

const OptionValue<Gecode::IntConLevel> INT_CON_LEVELS[] = {
    {"val", Gecode::ICL_VAL},
    {"bnd", Gecode::ICL_BND},
    {"dom", Gecode::ICL_DOM},
    {"def", Gecode::ICL_DEF},
    {}
};

TEST_F(GecodeSolverTest, IntConLevelOption) {
  EXPECT_EQ(Gecode::ICL_DEF, solver_.icl());
  unsigned count = 0;
  for (const OptionValue<Gecode::IntConLevel>
      *p = INT_CON_LEVELS; p->name; ++p, ++count) {
    EXPECT_TRUE(ParseOptions(c_str(fmt::Format("icl={}") << p->name)));
    EXPECT_EQ(p->value, solver_.icl());
  }
  EXPECT_EQ(4u, count);
}

const OptionValue<Gecode::IntValBranch> VAL_BRANCHINGS[] = {
    {"min",        Gecode::INT_VAL_MIN()},
    {"med",        Gecode::INT_VAL_MED()},
    {"max",        Gecode::INT_VAL_MAX()},
    {"rnd",        Gecode::INT_VAL_RND(Gecode::Rnd(0))},
    {"split_min",  Gecode::INT_VAL_SPLIT_MIN()},
    {"split_max",  Gecode::INT_VAL_SPLIT_MAX()},
    {"range_min",  Gecode::INT_VAL_RANGE_MIN()},
    {"range_max",  Gecode::INT_VAL_RANGE_MAX()},
    {"values_min", Gecode::INT_VALUES_MIN()},
    {"values_max", Gecode::INT_VALUES_MAX()},
    {}
};

TEST_F(GecodeSolverTest, ValBranchingOption) {
  EXPECT_EQ(Gecode::INT_VAL_MIN().select(), solver_.val_branching().select());
  unsigned count = 0;
  for (const OptionValue<Gecode::IntValBranch>
      *p = VAL_BRANCHINGS; p->name; ++p, ++count) {
    EXPECT_TRUE(ParseOptions(
        c_str(fmt::Format("val_branching={}") << p->name)));
    EXPECT_EQ(p->value.select(), solver_.val_branching().select());
  }
  EXPECT_EQ(10u, count);
}

const OptionValue<Gecode::IntVarBranch> VAR_BRANCHINGS[] = {
    {"none",            Gecode::INT_VAR_NONE()},
    {"rnd",             Gecode::INT_VAR_RND(Gecode::Rnd(0))},
    {"degree_min",      Gecode::INT_VAR_DEGREE_MIN()},
    {"degree_max",      Gecode::INT_VAR_DEGREE_MAX()},
    {"afc_min",         Gecode::INT_VAR_AFC_MIN()},
    {"afc_max",         Gecode::INT_VAR_AFC_MAX()},
    {"min_min",         Gecode::INT_VAR_MIN_MIN()},
    {"min_max",         Gecode::INT_VAR_MIN_MAX()},
    {"max_min",         Gecode::INT_VAR_MAX_MIN()},
    {"max_max",         Gecode::INT_VAR_MAX_MAX()},
    {"size_min",        Gecode::INT_VAR_SIZE_MIN()},
    {"size_max",        Gecode::INT_VAR_SIZE_MAX()},
    {"degree_size_min", Gecode::INT_VAR_DEGREE_SIZE_MIN()},
    {"degree_size_max", Gecode::INT_VAR_DEGREE_SIZE_MAX()},
    {"afc_size_min",    Gecode::INT_VAR_AFC_SIZE_MIN()},
    {"afc_size_max",    Gecode::INT_VAR_AFC_SIZE_MAX()},
    {"regret_min_min",  Gecode::INT_VAR_REGRET_MIN_MIN()},
    {"regret_min_max",  Gecode::INT_VAR_REGRET_MIN_MAX()},
    {"regret_max_min",  Gecode::INT_VAR_REGRET_MAX_MIN()},
    {"regret_max_max",  Gecode::INT_VAR_REGRET_MAX_MAX()},
    {}
};

TEST_F(GecodeSolverTest, VarBranchingOption) {
  EXPECT_EQ(Gecode::INT_VAR_SIZE_MIN().select(),
      solver_.var_branching().select());
  unsigned count = 0;
  for (const OptionValue<Gecode::IntVarBranch>
      *p = VAR_BRANCHINGS; p->name; ++p, ++count) {
    EXPECT_TRUE(ParseOptions(
        c_str(fmt::Format("var_branching={}") << p->name)));
    EXPECT_EQ(p->value.select(), solver_.var_branching().select());
  }
  EXPECT_EQ(20u, count);
}

TEST_F(GecodeSolverTest, OutLevOption) {
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

TEST_F(GecodeSolverTest, OutFreqOption) {
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
}
}

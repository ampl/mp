/*
 IlogCP solver tests.

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

#include <ilconcert/ilodiffi.h>
#include <ilconcert/ilopathi.h>
#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

#include "gtest/gtest.h"

#include "solvers/ilogcp/ilogcp.h"
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
#include "tests/config.h"

#ifdef HAVE_THREADS
# include <thread>
#endif

using std::ifstream;
using std::size_t;
using std::string;
using std::vector;

using ampl::CPLEXOptimizer;
using ampl::CPOptimizer;
using ampl::ExprBuilder;
using ampl::IlogCPSolver;
using ampl::LogicalExpr;
using ampl::NumericExpr;
using ampl::UnsupportedExprError;

#define DATA_DIR "../data/"

namespace {

std::auto_ptr<ampl::BasicSolver> CreateSolver() {
  return std::auto_ptr<ampl::BasicSolver>(new ampl::IlogCPSolver());
}

INSTANTIATE_TEST_CASE_P(IlogCP, SolverTest, ::testing::Values(CreateSolver));

bool AreBothSpaces(char lhs, char rhs) { return lhs == ' ' && rhs == ' '; }

// Replace all occurrences of string old_s in s with new_s.
void Replace(string &s, const string &old_s, const string &new_s) {
  size_t pos = 0;
  while ((pos = s.find(old_s, pos)) != string::npos) {
    s.replace(pos, old_s.length(), new_s);
    pos += new_s.length();
  }
}

// Returns a string representation of the argument.
template <typename T>
string str(T t) {
  std::ostringstream ss;
  ss << t;
  string s = ss.str();

  // Replace adjacent duplicate spaces and possible trailing space.
  string::iterator end = std::unique(s.begin(), s.end(), AreBothSpaces);
  if (*(end - 1) == ' ') --end;
  s.erase(end, s.end());

  // Normalize representation of infinity.
  Replace(s, "1.#INF", "inf");
  return s;
}

// Evaluates a Concert expression.
double eval(IloExpr e) {
  return e.getImpl()->eval(IloAlgorithm());
}

struct EnumValue {
  const char *name;
  IloCP::ParameterValues value;
};

class IlogCPTest : public ::testing::Test, public ExprBuilder {
 protected:
  IlogCPSolver s;
  ampl::NLToConcertConverter converter_;

  IloNumVarArray CreateVars() {
    IloEnv env = s.env();
    IloNumVarArray vars = IloNumVarArray(env, 3);
    vars[0] = IloIntVar(env, 0, 1, "x");
    vars[1] = IloNumVar(env, 0, 1, "y");
    vars[2] = IloNumVar(env, 0, 1, "theta");
    return vars;
  }

  IlogCPTest() :
    converter_(s.env(), CreateVars(), false, false) {}

  int RunSolver(const char *stub = nullptr, const char *opt = nullptr) {
    try {
      return s.Run(Args("ilogcp", "-s", stub, opt));
    } catch (const IloException &e) {  // NOLINT(whitespace/parens)
      throw std::runtime_error(e.getMessage());
    }
    return 0;
  }

  bool ParseOptions(const char *opt1, const char *opt2 = nullptr) {
    try {
      return s.ParseOptions(Args(opt1, opt2));
    } catch (const IloException &e) {  // NOLINT(whitespace/parens)
      throw std::runtime_error(e.getMessage());
    }
    return false;
  }

  SolveResult Solve(const char *stub, const char *opt = nullptr);

  template <typename T>
  static string Option(const char *name, T value) {
    std::ostringstream os;
    os << name << "=" << value;
    return os.str();
  }

  void CheckIntCPOption(const char *option, IloCP::IntParam param,
      int start, int end, int offset = 0, bool accepts_auto = false,
      const EnumValue *values = 0);

  template <typename ParamT>
  void CheckIntCPLEXOption(const char *option,
      ParamT param, int start, int end) {
    for (int i = std::max(start, -9); i <= std::min(end, 9); ++i) {
      EXPECT_TRUE(ParseOptions("optimizer=cplex", Option(option, i).c_str()));
      CPLEXOptimizer *opt = dynamic_cast<CPLEXOptimizer*>(s.optimizer());
      ASSERT_TRUE(opt != nullptr);
      EXPECT_EQ(i, opt->cplex().getParam(param))
        << "Failed option: " << option;
    }
    if (end != INT_MAX) {
      EXPECT_FALSE(ParseOptions("optimizer=cplex",
          Option(option, end + 1).c_str()));
    }
    if (start != INT_MIN) {
      int small = start - 1;
      EXPECT_FALSE(ParseOptions("optimizer=cplex",
          Option(option, small).c_str()));
      EXPECT_FALSE(ParseOptions("optimizer=cp", Option(option, start).c_str()));
    }
  }

  void CheckDblCPOption(const char *option, IloCP::NumParam param,
      double good, double bad);
};

SolveResult IlogCPTest::Solve(const char *stub, const char *opt) {
  TestSolutionHandler sh;
  s.set_solution_handler(&sh);
  RunSolver(stub, opt);
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

void IlogCPTest::CheckIntCPOption(const char *option,
    IloCP::IntParam param, int start, int end, int offset, bool accepts_auto,
    const EnumValue *values) {
  for (int i = start; i <= std::min(end, 9); ++i) {
    EXPECT_TRUE(ParseOptions("optimizer=cp", Option(option, i).c_str()));
    CPOptimizer *opt = dynamic_cast<CPOptimizer*>(s.optimizer());
    ASSERT_TRUE(opt != nullptr);
    EXPECT_EQ(offset + i, opt->solver().getParameter(param))
      << "Failed option: " << option;
  }
  if (end != INT_MAX)
    EXPECT_FALSE(ParseOptions("optimizer=cp", Option(option, end + 1).c_str()));
  if (accepts_auto) {
    EXPECT_TRUE(ParseOptions("optimizer=cp", Option(option, -1).c_str()));
    CPOptimizer *opt = dynamic_cast<CPOptimizer*>(s.optimizer());
    ASSERT_TRUE(opt != nullptr);
    EXPECT_EQ(IloCP::Auto, opt->solver().getParameter(param));

    EXPECT_TRUE(ParseOptions("optimizer=cp", Option(option, "auto").c_str()));
    opt = dynamic_cast<CPOptimizer*>(s.optimizer());
    ASSERT_TRUE(opt != nullptr);
    EXPECT_EQ(IloCP::Auto, opt->solver().getParameter(param));
  }
  int small = start - 1;
  if (accepts_auto && small == -1)
    --small;
  EXPECT_FALSE(ParseOptions("optimizer=cp", Option(option, small).c_str()));
  EXPECT_FALSE(ParseOptions("optimizer=cplex", Option(option, start).c_str()));
  if (values) {
    int count = 0;
    for (const EnumValue *v = values; v->name; ++v, ++count) {
      EXPECT_TRUE(ParseOptions(
          "optimizer=cp", Option(option, v->name).c_str()));
      CPOptimizer *opt = dynamic_cast<CPOptimizer*>(s.optimizer());
      ASSERT_TRUE(opt != nullptr);
      EXPECT_EQ(v->value, opt->solver().getParameter(param))
        << "Failed option: " << option;
    }
    EXPECT_EQ(end - start + 1, count);
  }
}

void IlogCPTest::CheckDblCPOption(const char *option,
    IloCP::NumParam param, double good, double bad) {
  EXPECT_TRUE(ParseOptions("optimizer=cp", Option(option, good).c_str()));
  CPOptimizer *opt = dynamic_cast<CPOptimizer*>(s.optimizer());
  ASSERT_TRUE(opt != nullptr);
  EXPECT_EQ(good, opt->solver().getParameter(param))
    << "Failed option: " << option;

  EXPECT_FALSE(ParseOptions("optimizer=cp", Option(option, bad).c_str()));
  EXPECT_FALSE(ParseOptions("optimizer=cplex", Option(option, good).c_str()));
}

TEST_F(IlogCPTest, IloArrayCopyingIsCheap) {
  IloIntArray array(s.env());
  array.add(42);
  EXPECT_TRUE(array.getImpl() != nullptr);
  EXPECT_EQ(array.getImpl(), IloIntArray(array).getImpl());
}

TEST_F(IlogCPTest, ConvertSingleNumberOfToIloDistribute) {
  s.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  EXPECT_EQ("IloIntVar(4)" + bounds, str(converter_.Visit(
      AddNumberOf(AddNum(42), AddVar(0), AddVar(1)))));
  converter_.FinishBuildingNumberOf();
  IloModel::Iterator iter(converter_.model());
  ASSERT_NE(0, iter.ok());
  EXPECT_EQ("IloIntVar(6)" + bounds + " == x", str(*iter));
  ++iter;
  ASSERT_NE(0, iter.ok());
  EXPECT_EQ("IloIntVar(9)" + bounds + " == y", str(*iter));
  ++iter;
  ASSERT_NE(0, iter.ok());
  IloDistributeI *dist = dynamic_cast<IloDistributeI*>((*iter).getImpl());
  ASSERT_TRUE(dist != nullptr);
  EXPECT_EQ("[IloIntVar(4)" + bounds + " ]", str(dist->getCardVarArray()));
  EXPECT_EQ("[IloIntVar(6)" + bounds + " , IloIntVar(9)" + bounds + " ]",
      str(dist->getVarArray()));
  EXPECT_EQ("[42]", str(dist->getValueArray()));
  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithSameValuesToIloDistribute) {
  s.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  NumericExpr expr(AddNumberOf(AddNum(42), AddVar(0), AddVar(1)));
  EXPECT_EQ("IloIntVar(4)" + bounds, str(converter_.Visit(expr)));
  EXPECT_EQ("IloIntVar(4)" + bounds, str(converter_.Visit(
      AddNumberOf(AddNum(42), AddVar(0), AddVar(1)))));
  converter_.FinishBuildingNumberOf();
  IloModel::Iterator iter(converter_.model());
  ASSERT_NE(0, iter.ok());
  EXPECT_EQ("IloIntVar(7)" + bounds + " == x", str(*iter));
  ++iter;
  ASSERT_NE(0, iter.ok());
  EXPECT_EQ("IloIntVar(10)" + bounds + " == y", str(*iter));
  ++iter;
  ASSERT_NE(0, iter.ok());
  IloDistributeI *dist = dynamic_cast<IloDistributeI*>((*iter).getImpl());
  ASSERT_TRUE(dist != nullptr);
  EXPECT_EQ("[IloIntVar(4)" + bounds + " ]", str(dist->getCardVarArray()));
  EXPECT_EQ("[IloIntVar(7)" + bounds + " , IloIntVar(10)" + bounds + " ]",
      str(dist->getVarArray()));
  EXPECT_EQ("[42]", str(dist->getValueArray()));
  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithDiffValuesToIloDistribute) {
  s.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  EXPECT_EQ("IloIntVar(4)" + bounds,
      str(converter_.Visit(AddNumberOf(AddNum(42), AddVar(0), AddVar(1)))));
  EXPECT_EQ("IloIntVar(6)" + bounds,
      str(converter_.Visit(AddNumberOf(AddNum(43), AddVar(0), AddVar(1)))));
  converter_.FinishBuildingNumberOf();
  IloModel::Iterator iter(converter_.model());
  ASSERT_NE(0, iter.ok());
  EXPECT_EQ("IloIntVar(8)" + bounds + " == x", str(*iter));
  ++iter;
  ASSERT_NE(0, iter.ok());
  EXPECT_EQ("IloIntVar(11)" + bounds + " == y", str(*iter));
  ++iter;
  ASSERT_NE(0, iter.ok());
  IloDistributeI *dist = dynamic_cast<IloDistributeI*>((*iter).getImpl());
  ASSERT_TRUE(dist != nullptr);
  EXPECT_EQ("[IloIntVar(4)" + bounds + " , IloIntVar(6)" + bounds + " ]",
      str(dist->getCardVarArray()));
  EXPECT_EQ("[IloIntVar(8)" + bounds + " , IloIntVar(11)" + bounds + " ]",
      str(dist->getVarArray()));
  EXPECT_EQ("[42, 43]", str(dist->getValueArray()));
  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithDiffExprs) {
  s.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  EXPECT_EQ("IloIntVar(4)" + bounds,
      str(converter_.Visit(AddNumberOf(AddNum(42), AddVar(0), AddVar(1)))));
  EXPECT_EQ("IloIntVar(6)" + bounds,
      str(converter_.Visit(AddNumberOf(AddNum(42), AddVar(2)))));
  converter_.FinishBuildingNumberOf();
  IloModel::Iterator iter(converter_.model());
  ASSERT_NE(0, iter.ok());
  EXPECT_EQ("IloIntVar(8)" + bounds + " == x", str(*iter));
  ++iter;
  ASSERT_NE(0, iter.ok());
  EXPECT_EQ("IloIntVar(11)" + bounds + " == y", str(*iter));
  ++iter;
  ASSERT_NE(0, iter.ok());
  IloDistributeI *dist = dynamic_cast<IloDistributeI*>((*iter).getImpl());
  ASSERT_TRUE(dist != nullptr);
  EXPECT_EQ("[IloIntVar(4)" + bounds + " ]", str(dist->getCardVarArray()));
  EXPECT_EQ("[IloIntVar(8)" + bounds + " , IloIntVar(11)" + bounds + " ]",
      str(dist->getVarArray()));
  EXPECT_EQ("[42]", str(dist->getValueArray()));
  ++iter;
  ASSERT_NE(0, iter.ok());
  EXPECT_EQ("IloIntVar(15)" + bounds + " == theta", str(*iter));
  ++iter;
  ASSERT_NE(0, iter.ok());
  dist = dynamic_cast<IloDistributeI*>((*iter).getImpl());
  ASSERT_TRUE(dist != nullptr);
  EXPECT_EQ("[IloIntVar(6)" + bounds + " ]", str(dist->getCardVarArray()));
  EXPECT_EQ("[IloIntVar(15)" + bounds + " ]",
      str(dist->getVarArray()));
  EXPECT_EQ("[42]", str(dist->getValueArray()));
  ++iter;
  EXPECT_FALSE(iter.ok());
}

// ----------------------------------------------------------------------------
// Solver tests

TEST_F(IlogCPTest, CPOptimizerDoesntSupportContinuousVars) {
  EXPECT_THROW(RunSolver(DATA_DIR "objconst", "optimizer=cp"), ampl::Error);
}

TEST_F(IlogCPTest, SolveNumberOfCplex) {
  s.use_numberof(false);
  RunSolver(DATA_DIR "numberof", "optimizer=cplex");
}

TEST_F(IlogCPTest, SolveAssign0) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign0").obj);
}

TEST_F(IlogCPTest, SolveAssign1) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign1").obj);
}

TEST_F(IlogCPTest, SolveBalassign0) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign0").obj);
}

TEST_F(IlogCPTest, SolveBalassign1) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign1").obj);
}

TEST_F(IlogCPTest, SolveFlowshp0) {
  EXPECT_NEAR(22, Solve(DATA_DIR "flowshp0").obj, 1e-5);
}

TEST_F(IlogCPTest, SolveFlowshp1) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp1").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(IlogCPTest, DISABLED_SolveFlowshp2) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp2").obj);
}

TEST_F(IlogCPTest, SolveGrpassign0) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign0").obj);
}

// Disabled because variables in subscripts are not yet allowed.
TEST_F(IlogCPTest, DISABLED_SolveGrpassign1) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1").obj);
}

// Disabled because object-valued variables are not yet allowed.
TEST_F(IlogCPTest, DISABLED_SolveGrpassign1a) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1a").obj);
}

TEST_F(IlogCPTest, SolveMagic) {
  EXPECT_TRUE(Solve(DATA_DIR "magic").solved);
}

TEST_F(IlogCPTest, SolveMapcoloring) {
  EXPECT_TRUE(Solve(DATA_DIR "mapcoloring").solved);
}

TEST_F(IlogCPTest, SolveNQueens) {
  EXPECT_TRUE(Solve(DATA_DIR "nqueens").solved);
}

TEST_F(IlogCPTest, SolveNQueens0) {
  EXPECT_TRUE(Solve(DATA_DIR "nqueens0").solved);
}

TEST_F(IlogCPTest, SolveOpenShop) {
  EXPECT_NEAR(1955, Solve(DATA_DIR "openshop", "optimizer=cplex").obj, 1e-5);
}

// Disabled because it's too difficult to solve.
TEST_F(IlogCPTest, DISABLED_SolveParty1) {
  EXPECT_EQ(61, Solve(DATA_DIR "party1").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(IlogCPTest, DISABLED_SolveParty2) {
  EXPECT_EQ(3, Solve(DATA_DIR "party2").obj);
}

TEST_F(IlogCPTest, SolvePhoto9) {
  EXPECT_EQ(10, Solve(DATA_DIR "photo9").obj);
}

// Disabled because it takes somewhat long (compared to other tests).
TEST_F(IlogCPTest, DISABLED_SolvePhoto11) {
  EXPECT_EQ(12, Solve(DATA_DIR "photo11").obj);
}

TEST_F(IlogCPTest, SolveSched0) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched0").obj);
}

TEST_F(IlogCPTest, SolveSched1) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched1").obj);
}

TEST_F(IlogCPTest, SolveSched2) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched2").obj);
}

TEST_F(IlogCPTest, SolveSendMoreMoney) {
  EXPECT_TRUE(Solve(DATA_DIR "send-more-money").solved);
}

TEST_F(IlogCPTest, SolveSendMostMoney) {
  EXPECT_NEAR(10876, Solve(DATA_DIR "send-most-money",
      "relativeoptimalitytolerance=1e-5").obj, 1e-5);
}

TEST_F(IlogCPTest, SolveSeq0) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0").obj, 1e-5);
}

TEST_F(IlogCPTest, SolveSeq0a) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0a").obj, 1e-5);
}

TEST_F(IlogCPTest, SolveSudokuHard) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuHard").solved);
}

TEST_F(IlogCPTest, SolveSudokuVeryEasy) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuVeryEasy").solved);
}

// ----------------------------------------------------------------------------
// Option tests

TEST_F(IlogCPTest, VersionOption) {
  EXPECT_FALSE((s.flags() & ASL_OI_show_version) != 0);
  EXPECT_TRUE(ParseOptions("version"));
  EXPECT_TRUE((s.flags() & ASL_OI_show_version) != 0);
}

TEST_F(IlogCPTest, WantsolOption) {
  EXPECT_EQ(0, s.wantsol());
  EXPECT_TRUE(ParseOptions("wantsol=1"));
  EXPECT_EQ(1, s.wantsol());
  EXPECT_TRUE(ParseOptions("wantsol=5"));
  EXPECT_EQ(5, s.wantsol());
}

TEST_F(IlogCPTest, DebugExprOption) {
  EXPECT_TRUE(ParseOptions("debugexpr=0"));
  EXPECT_EQ(0, s.GetOption(IlogCPSolver::DEBUGEXPR));
  EXPECT_TRUE(ParseOptions("debugexpr=1"));
  EXPECT_EQ(1, s.GetOption(IlogCPSolver::DEBUGEXPR));
  EXPECT_FALSE(ParseOptions("debugexpr=42"));
  EXPECT_FALSE(ParseOptions("debugexpr=oops"));
}

TEST_F(IlogCPTest, OptimizerOption) {
  EXPECT_EQ(IlogCPSolver::AUTO, s.GetOption(IlogCPSolver::OPTIMIZER));

  EXPECT_TRUE(ParseOptions("optimizer=cplex"));
  EXPECT_EQ(IlogCPSolver::CPLEX, s.GetOption(IlogCPSolver::OPTIMIZER));
  EXPECT_TRUE(dynamic_cast<IloCplexI*>(s.alg().getImpl()) != nullptr);

  EXPECT_TRUE(ParseOptions("optimizer=cp"));
  EXPECT_EQ(IlogCPSolver::CP, s.GetOption(IlogCPSolver::OPTIMIZER));
  EXPECT_TRUE(dynamic_cast<IloCplexI*>(s.alg().getImpl()) == nullptr);
}

TEST_F(IlogCPTest, TimingOption) {
  EXPECT_TRUE(ParseOptions("timing=0"));
  EXPECT_EQ(0, s.GetOption(IlogCPSolver::TIMING));
  EXPECT_TRUE(ParseOptions("timing=1"));
  EXPECT_EQ(1, s.GetOption(IlogCPSolver::TIMING));
  EXPECT_FALSE(ParseOptions("timing=42"));
  EXPECT_FALSE(ParseOptions("timing=oops"));
}

TEST_F(IlogCPTest, UseNumberOfOption) {
  EXPECT_TRUE(ParseOptions("usenumberof=0"));
  EXPECT_EQ(0, s.GetOption(IlogCPSolver::USENUMBEROF));
  EXPECT_TRUE(ParseOptions("usenumberof=1"));
  EXPECT_EQ(1, s.GetOption(IlogCPSolver::USENUMBEROF));
  EXPECT_FALSE(ParseOptions("usenumberof=42"));
  EXPECT_FALSE(ParseOptions("timing=oops"));
}

TEST_F(IlogCPTest, CPFlagOptions) {
  const EnumValue flags[] = {
      {"off", IloCP::Off},
      {"on",  IloCP::On},
      {nullptr}
  };
  CheckIntCPOption("constraintaggregation", IloCP::ConstraintAggregation,
      0, 1, IloCP::Off, false, flags);
  CheckIntCPOption("dynamicprobing", IloCP::DynamicProbing,
      0, 1, IloCP::Off, true, flags);
  CheckIntCPOption("temporalrelaxation", IloCP::TemporalRelaxation,
      0, 1, IloCP::Off, false, flags);
}

TEST_F(IlogCPTest, CPInferenceLevelOptions) {
  const EnumValue inf_levels[] = {
      {"default",  IloCP::Default},
      {"low",      IloCP::Low},
      {"basic",    IloCP::Basic},
      {"medium",   IloCP::Medium},
      {"extended", IloCP::Extended},
      {nullptr}
  };
  CheckIntCPOption("alldiffinferencelevel", IloCP::AllDiffInferenceLevel,
      0, 4, IloCP::Default, false, inf_levels);
  CheckIntCPOption("defaultinferencelevel", IloCP::DefaultInferenceLevel,
      1, 4, IloCP::Default, false, inf_levels + 1);
  CheckIntCPOption("distributeinferencelevel", IloCP::DistributeInferenceLevel,
      0, 4, IloCP::Default, false, inf_levels);
}

TEST_F(IlogCPTest, CPDefaultVerbosityQuiet) {
  EXPECT_TRUE(ParseOptions("optimizer=cp"));
  CPOptimizer *opt = dynamic_cast<CPOptimizer*>(s.optimizer());
  ASSERT_TRUE(opt != nullptr);
  EXPECT_EQ(IloCP::Quiet, opt->solver().getParameter(IloCP::LogVerbosity));
}

TEST_F(IlogCPTest, CPVerbosityOptions) {
  const EnumValue verbosities[] = {
      {"quiet",   IloCP::Quiet},
      {"terse",   IloCP::Terse},
      {"normal",  IloCP::Normal},
      {"verbose", IloCP::Verbose},
      {nullptr}
  };
  CheckIntCPOption("logverbosity", IloCP::LogVerbosity,
      0, 3, IloCP::Quiet, false, verbosities);
  CheckIntCPOption("outlev", IloCP::LogVerbosity,
      0, 3, IloCP::Quiet, false, verbosities);
  CheckIntCPOption("propagationlog", IloCP::PropagationLog,
      0, 3, IloCP::Quiet, false, verbosities);
}

TEST_F(IlogCPTest, CPSearchTypeOption) {
  const EnumValue types[] = {
      {"depthfirst", IloCP::DepthFirst},
      {"restart",    IloCP::Restart},
      {"multipoint", IloCP::MultiPoint},
      {nullptr}
  };
  CheckIntCPOption("searchtype", IloCP::SearchType,
      0, 2, IloCP::DepthFirst, CPX_VERSION > 1220, types);
}

TEST_F(IlogCPTest, CPTimeModeOption) {
  const EnumValue modes[] = {
      {"cputime",     IloCP::CPUTime},
      {"elapsedtime", IloCP::ElapsedTime},
      {nullptr}
  };
  CheckIntCPOption("timemode", IloCP::TimeMode,
      0, 1, IloCP::CPUTime, false, modes);
}

TEST_F(IlogCPTest, CPOptions) {
  CheckIntCPOption("branchlimit", IloCP::BranchLimit, 0, INT_MAX);
  CheckIntCPOption("choicepointlimit", IloCP::ChoicePointLimit, 0, INT_MAX);
  CheckDblCPOption("dynamicprobingstrength",
      IloCP::DynamicProbingStrength, 42, -1);
  CheckIntCPOption("faillimit", IloCP::FailLimit, 0, INT_MAX);
  CheckIntCPOption("logperiod", IloCP::LogPeriod, 1, INT_MAX);
  CheckIntCPOption("multipointnumberofsearchpoints",
      IloCP::MultiPointNumberOfSearchPoints, 2, INT_MAX);
  CheckDblCPOption("optimalitytolerance", IloCP::OptimalityTolerance, 42, -1);
  CheckIntCPOption("randomseed", IloCP::RandomSeed, 0, INT_MAX);
  CheckDblCPOption("relativeoptimalitytolerance",
      IloCP::RelativeOptimalityTolerance, 42, -1);
  CheckDblCPOption("restartgrowthfactor", IloCP::RestartGrowthFactor, 42, -1);
  CheckIntCPOption("restartfaillimit", IloCP::RestartFailLimit, 1, INT_MAX);
  CheckIntCPOption("solutionlimit", IloCP::SolutionLimit, 0, INT_MAX);
  CheckDblCPOption("timelimit", IloCP::TimeLimit, 42, -1);
  if (CPX_VERSION > 1220)
    CheckIntCPOption("workers", IloCP::Workers, 0, INT_MAX, 0, true);
  else
    CheckIntCPOption("workers", IloCP::Workers, 1, 4, 0, false);
}

TEST_F(IlogCPTest, CPLEXDefaultMIPDisplayZero) {
  EXPECT_TRUE(ParseOptions("optimizer=cplex"));
  CPLEXOptimizer *opt = dynamic_cast<CPLEXOptimizer*>(s.optimizer());
  ASSERT_TRUE(opt != nullptr);
  EXPECT_EQ(0, opt->cplex().getParam(IloCplex::MIPDisplay));
}

TEST_F(IlogCPTest, CPLEXOptions) {
  CheckIntCPLEXOption("mipdisplay", IloCplex::MIPDisplay, 0, 5);
  CheckIntCPLEXOption("mipinterval",
      IloCplex::MIPInterval, INT_MIN, INT_MAX);
}

// ----------------------------------------------------------------------------
// Solve code tests

TEST_F(IlogCPTest, OptimalSolveCode) {
  Solve(DATA_DIR "objconst");
  EXPECT_EQ(0, s.problem().solve_code());
}

TEST_F(IlogCPTest, FeasibleSolveCode) {
  Solve(DATA_DIR "feasible");
  EXPECT_EQ(100, s.problem().solve_code());
}

TEST_F(IlogCPTest, InfeasibleSolveCode) {
  Solve(DATA_DIR "infeasible");
  EXPECT_EQ(200, s.problem().solve_code());
}

TEST_F(IlogCPTest, InfeasibleOrUnboundedSolveCode) {
  Solve(DATA_DIR "unbounded");
  EXPECT_EQ(201, s.problem().solve_code());
}

// ----------------------------------------------------------------------------

#ifdef HAVE_THREADS
void Interrupt() {
  // Wait until started.
  while (ampl::SignalHandler::stop())
    std::this_thread::yield();
  std::raise(SIGINT);
}

TEST_F(IlogCPTest, InterruptCPLEX) {
  std::thread t(Interrupt);
  std::string message =
      Solve(DATA_DIR "miplib/assign1", "optimizer=cplex").message;
  t.join();
  EXPECT_EQ(600, s.problem().solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}

TEST_F(IlogCPTest, InterruptCP) {
  std::thread t(Interrupt);
  std::string message =
      Solve(DATA_DIR "miplib/assign1", "optimizer=cp").message;
  t.join();
  EXPECT_EQ(600, s.problem().solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}
#endif
}

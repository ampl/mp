/*
 IlogCP solver tests.

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

#if __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunused-parameter"
# pragma clang diagnostic ignored "-Wunused-private-field"
#endif

#include <ilconcert/ilodiffi.h>
#include <ilconcert/ilopathi.h>
#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#if __clang__
# pragma clang diagnostic pop
#endif

#include <algorithm>
#include <memory>
#include <string>
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

using std::size_t;
using std::string;

using ampl::CPLEXOptimizer;
using ampl::CPOptimizer;
using ampl::IlogCPSolver;
using ampl::Problem;
using ampl::UnsupportedExprError;

namespace {

SolverPtr CreateSolver() { return SolverPtr(new ampl::IlogCPSolver()); }

INSTANTIATE_TEST_CASE_P(IlogCP, SolverTest,
    ::testing::Values(SolverTestParam(CreateSolver, feature::ALL)));

TEST_P(SolverTest, CPOptimizerDoesntSupportContinuousVars) {
  EXPECT_THROW(Solve("objconst", "optimizer=cp"), ampl::Error);
}

TEST_P(SolverTest, SolveBalassign0) {
  EXPECT_EQ(14, Solve("balassign0").obj);
}

TEST_P(SolverTest, SolveBalassign1) {
  EXPECT_EQ(14, Solve("balassign1").obj);
}

TEST_P(SolverTest, SolveFlowshp2) {
  EXPECT_EQ(22, Solve("flowshp2", "optimizer=cplex").obj);
}

TEST_P(SolverTest, SolveOpenShop) {
  EXPECT_NEAR(1955, Solve("openshop", "optimizer=cplex").obj, 1e-5);
}

struct EnumValue {
  const char *name;
  IloCP::ParameterValues value;
};

class IlogCPTest : public ::testing::Test, public ampl::ExprBuilder {
 protected:
  IlogCPSolver s;

  int CountIloDistribute();

  bool ParseOptions(const char *opt1,
      const char *opt2 = nullptr, const Problem *p = nullptr) {
    try {
      return s.ParseOptions(Args(opt1, opt2), 0, p);
    } catch (const IloException &e) {  // NOLINT(whitespace/parens)
      throw std::runtime_error(e.getMessage());
    }
    return false;
  }

  SolveResult Solve(Problem &p, const char *stub, const char *opt = nullptr) {
    return SolverTest::Solve(s, p, stub, opt);
  }

  template <typename T>
  static string Option(const char *name, T value) {
    return str(fmt::Format("{}={}") << name << value);
  }

  void CheckIntCPOption(const char *option, IloCP::IntParam param,
      int start, int end, int offset = 0, bool accepts_auto = false,
      const EnumValue *values = 0);

  template <typename ParamT>
  void CheckIntCPLEXOption(
      const char *option, ParamT param, int start, int end);

  void CheckDblCPOption(const char *option,
      IloCP::NumParam param, double good, double bad);
};

int IlogCPTest::CountIloDistribute() {
  int count = 0;
  for (IloModel::Iterator i(s.alg().getModel()); i.ok(); ++i) {
    if (dynamic_cast<IloDistributeI*>((*i).getImpl()))
      ++count;
  }
  return count;
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
    if (!values) {
      if (accepts_auto)
        EXPECT_EQ(str(fmt::Format("{}") << i), s.GetStrOption(option));
      else
        EXPECT_EQ(i, s.GetIntOption(option));
    }
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
    EXPECT_EQ("auto", s.GetStrOption(option));
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
      EXPECT_EQ(v->name, s.GetStrOption(option)) << "Failed option: " << option;
    }
    EXPECT_EQ(end - start + 1, count);
  }
}

template <typename ParamT>
void IlogCPTest::CheckIntCPLEXOption(const char *option,
    ParamT param, int start, int end) {
  for (int i = std::max(start, -9); i <= std::min(end, 9); ++i) {
    EXPECT_TRUE(ParseOptions("optimizer=cplex", Option(option, i).c_str()));
    CPLEXOptimizer *opt = dynamic_cast<CPLEXOptimizer*>(s.optimizer());
    ASSERT_TRUE(opt != nullptr);
    EXPECT_EQ(i, opt->cplex().getParam(param))
      << "Failed option: " << option;
    EXPECT_EQ(i, s.GetIntOption(option));
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

void IlogCPTest::CheckDblCPOption(const char *option,
    IloCP::NumParam param, double good, double bad) {
  EXPECT_TRUE(ParseOptions("optimizer=cp", Option(option, good).c_str()));
  CPOptimizer *opt = dynamic_cast<CPOptimizer*>(s.optimizer());
  ASSERT_TRUE(opt != nullptr);
  EXPECT_EQ(good, opt->solver().getParameter(param))
    << "Failed option: " << option;
  EXPECT_EQ(good, s.GetDblOption(option));

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
  ampl::Problem p;
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddCon(AddRelational(EQ, AddNum(0),
      AddNumberOf(AddNum(42), AddVar(0), AddVar(1))));
  s.Solve(p);
  ASSERT_EQ(1, CountIloDistribute());
}

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithSameValuesToIloDistribute) {
  s.use_numberof();
  ampl::Problem p;
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddCon(AddRelational(EQ, AddNum(0),
      AddNumberOf(AddNum(42), AddVar(0), AddVar(1))));
  p.AddCon(AddRelational(EQ, AddNum(0),
      AddNumberOf(AddNum(42), AddVar(0), AddVar(1))));
  s.Solve(p);
  ASSERT_EQ(1, CountIloDistribute());
}

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithDiffValuesToIloDistribute) {
  s.use_numberof();
  ampl::Problem p;
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddCon(AddRelational(EQ, AddNum(0),
      AddNumberOf(AddNum(42), AddVar(0), AddVar(1))));
  p.AddCon(AddRelational(EQ, AddNum(0),
      AddNumberOf(AddNum(43), AddVar(0), AddVar(1))));
  s.Solve(p);
  ASSERT_EQ(1, CountIloDistribute());
}

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithDiffExprs) {
  s.use_numberof();
  ampl::Problem p;
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddCon(AddRelational(EQ, AddNum(0),
      AddNumberOf(AddNum(42), AddVar(0), AddVar(1))));
  p.AddCon(AddRelational(EQ, AddNum(0),
      AddNumberOf(AddNum(42), AddVar(1))));
  s.Solve(p);
  ASSERT_EQ(2, CountIloDistribute());
}

struct TestSolutionHandler : ampl::SolutionHandler {
  int num_solutions;
  TestSolutionHandler() : num_solutions(0) {}
  void HandleSolution(Problem &, fmt::StringRef,
        const double *values, const double *, double) {
    ++num_solutions;
    for (int i = 0; i < 3; ++i) {
      EXPECT_GE(values[i], 1);
      EXPECT_LE(values[i], 3);
    }
    EXPECT_NE(values[0], values[1]);
    EXPECT_NE(values[0], values[2]);
    EXPECT_NE(values[1], values[2]);
  }
};

TEST_F(IlogCPTest, SolutionLimit) {
  ampl::Problem p;
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddCon(AddAllDiff(AddVar(0), AddVar(1), AddVar(2)));
  TestSolutionHandler sh;
  s.set_solution_handler(&sh);
  ParseOptions("optimizer=cp", nullptr, &p);
  s.Solve(p);
  EXPECT_EQ(1, sh.num_solutions);
}

// TODO: this test is currently disabled because CP reports identical solutions
TEST_F(IlogCPTest, DISABLED_SolutionLimit) {
  ampl::Problem p;
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddCon(AddAllDiff(AddVar(0), AddVar(1), AddVar(2)));
  TestSolutionHandler sh;
  ParseOptions("optimizer=cp", "solutionlimit=100");
  s.set_solution_handler(&sh);
  s.Solve(p);
  EXPECT_EQ(6, sh.num_solutions);
}

TEST_F(IlogCPTest, SolveNumberOfCplex) {
  s.use_numberof(false);
  Problem p;
  Solve(p, "numberof", "optimizer=cplex");
}

TEST_F(IlogCPTest, InfeasibleOrUnboundedSolveCode) {
  Problem p;
  Solve(p, "unbounded");
  EXPECT_EQ(201, p.solve_code());
}

// ----------------------------------------------------------------------------
// Option tests

TEST_F(IlogCPTest, DebugExprOption) {
  EXPECT_TRUE(ParseOptions("debugexpr=0"));
  EXPECT_EQ(0, s.GetOption(IlogCPSolver::DEBUGEXPR));
  EXPECT_EQ(0, s.GetIntOption("debugexpr"));
  EXPECT_TRUE(ParseOptions("debugexpr=1"));
  EXPECT_EQ(1, s.GetOption(IlogCPSolver::DEBUGEXPR));
  EXPECT_EQ(1, s.GetIntOption("debugexpr"));
  EXPECT_FALSE(ParseOptions("debugexpr=42"));
  EXPECT_FALSE(ParseOptions("debugexpr=oops"));
}

TEST_F(IlogCPTest, OptimizerOption) {
  EXPECT_EQ(IlogCPSolver::AUTO, s.GetOption(IlogCPSolver::OPTIMIZER));
  EXPECT_EQ("auto", s.GetStrOption("optimizer"));

  EXPECT_TRUE(ParseOptions("optimizer=cplex"));
  EXPECT_EQ("cplex", s.GetStrOption("optimizer"));
  EXPECT_EQ(IlogCPSolver::CPLEX, s.GetOption(IlogCPSolver::OPTIMIZER));
  EXPECT_TRUE(dynamic_cast<IloCplexI*>(s.alg().getImpl()) != nullptr);

  EXPECT_TRUE(ParseOptions("optimizer=cp"));
  EXPECT_EQ("cp", s.GetStrOption("optimizer"));
  EXPECT_EQ(IlogCPSolver::CP, s.GetOption(IlogCPSolver::OPTIMIZER));
  EXPECT_TRUE(dynamic_cast<IloCplexI*>(s.alg().getImpl()) == nullptr);
}

TEST_F(IlogCPTest, UseNumberOfOption) {
  EXPECT_TRUE(ParseOptions("usenumberof=0"));
  EXPECT_EQ(0, s.GetOption(IlogCPSolver::USENUMBEROF));
  EXPECT_EQ(0, s.GetIntOption("usenumberof"));
  EXPECT_TRUE(ParseOptions("usenumberof=1"));
  EXPECT_EQ(1, s.GetOption(IlogCPSolver::USENUMBEROF));
  EXPECT_EQ(1, s.GetIntOption("usenumberof"));
  EXPECT_FALSE(ParseOptions("usenumberof=42"));
  EXPECT_FALSE(ParseOptions("usenumberof=oops"));
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
  EXPECT_EQ("quiet", s.GetStrOption("logverbosity"));
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
  EXPECT_EQ(0, s.GetIntOption("mipdisplay"));
}

TEST_F(IlogCPTest, CPLEXOptions) {
  CheckIntCPLEXOption("mipdisplay", IloCplex::MIPDisplay, 0, 5);
  CheckIntCPLEXOption("mipinterval", IloCplex::MIPInterval, INT_MIN, INT_MAX);
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
  Problem p;
  std::string message = Solve(p, "miplib/assign1", "optimizer=cplex").message;
  t.join();
  EXPECT_EQ(600, p.solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}

TEST_F(IlogCPTest, InterruptCP) {
  std::thread t(Interrupt);
  Problem p;
  std::string message = Solve(p, "miplib/assign1", "optimizer=cp").message;
  t.join();
  EXPECT_EQ(600, p.solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}
#endif
}

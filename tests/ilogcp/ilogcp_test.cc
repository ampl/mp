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

#include <ilconcert/ilopathi.h>
#include <ilcp/cp.h>
#include <ilcplex/cpxconst.h>

#if __clang__
# pragma clang diagnostic pop
#endif

#include <algorithm>
#include <stdexcept>
#include <string>

#include "gtest/gtest.h"

#include "solvers/ilogcp/ilogcp.h"

extern "C" {
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
}

#include "tests/args.h"
#include "tests/expr_builder.h"
#include "tests/solver_test.h"

using ampl::IlogCPSolver;
using ampl::InvalidOptionValue;
using ampl::OptionError;
using ampl::Problem;

using std::string;

namespace {

SolverPtr CreateSolver() { return SolverPtr(new IlogCPSolver()); }

INSTANTIATE_TEST_CASE_P(IlogCP, SolverTest,
    ::testing::Values(SolverTestParam(CreateSolver, feature::ALL)));

TEST_P(SolverTest, SolveBalassign0) {
  EXPECT_EQ(14, Solve("balassign0").obj);
}

TEST_P(SolverTest, SolveBalassign1) {
  EXPECT_EQ(14, Solve("balassign1").obj);
}

// ----------------------------------------------------------------------------
// element constraint tests

TEST_P(SolverTest, ElementConstraint) {
  EXPECT_EQ(22, Eval(AddCall("element", 11, 22, 33, x), 1));
}

TEST_P(SolverTest, TooFewArgsToElementConstraint) {
  EXPECT_THROW_MSG(Eval(AddCall("element", x), 0),
      ampl::Error, "element: too few arguments");
}

TEST_P(SolverTest, ElementConstantIndexOutOfBounds) {
  EXPECT_THROW_MSG(Eval(AddCall("element", 11, 22, 2)),
        ampl::Error, "element: index 2 is out of bounds");
}

TEST_P(SolverTest, ElementConstantAtConstantIndex) {
  EXPECT_EQ(22, Eval(AddCall("element", 11, 22, 1)));
}

TEST_P(SolverTest, ElementExprAtConstantIndex) {
  EXPECT_EQ(42, Eval(AddCall("element", x, 22, 0), 42));
}

TEST_P(SolverTest, ElementExprPlusConstantAtConstantIndex) {
  EXPECT_EQ(44, Eval(AddCall("element", 11, ampl::CallArg(2, x), 1), 42));
}

TEST_P(SolverTest, ElementVariableIndexOutOfBounds) {
  EXPECT_EQ(ampl::INFEASIBLE,
      Eval(AddCall("element", 11, 22, x), 2).solve_code());
}

TEST_P(SolverTest, ElementConstantAtVariableIndex) {
  EXPECT_EQ(22, Eval(AddCall("element", 11, 22, x), 1));
}

TEST_P(SolverTest, ElementExprAtVariableIndex) {
  EXPECT_EQ(42, Eval(AddCall("element", x, 22, y), 42, 0));
}

TEST_P(SolverTest, ElementExprPlusConstantAtVariableIndex) {
  EXPECT_EQ(44, Eval(AddCall("element", 11, ampl::CallArg(2, x), y), 42, 1));
}

// ----------------------------------------------------------------------------
// in_relation constraint tests

TEST_P(SolverTest, InRelationConstraint) {
  Problem p;
  p.AddVar(0, 100, ampl::INTEGER);
  p.AddObj(ampl::MIN, AddVar(0));
  p.AddCon(AddRelational(NE, AddCall("in_relation", AddVar(0), 42), AddNum(0)));
  EXPECT_EQ(42, Solve(p).obj_value());
}

TEST_P(SolverTest, NestedInRelationNotSupported) {
  EXPECT_THROW_MSG(Eval(AddBinary(OPPLUS,
      AddCall("in_relation", AddVar(0), 42), AddNum(1)));,
      ampl::UnsupportedExprError,
      "unsupported expression: nested 'in_relation'");
}

TEST_P(SolverTest, TooFewArgsToInRelationConstraint) {
  Problem p;
  p.AddVar(0, 100, ampl::INTEGER);
  p.AddObj(ampl::MIN, AddVar(0));
  p.AddCon(AddRelational(NE, AddCall("in_relation"), AddNum(0)));
  EXPECT_THROW_MSG(Solve(p), ampl::Error, "in_relation: too few arguments");
}

TEST_P(SolverTest, InRelationSizeIsNotMultipleOfArity) {
  Problem p;
  p.AddVar(0, 100, ampl::INTEGER);
  p.AddVar(0, 100, ampl::INTEGER);
  p.AddObj(ampl::MIN, AddVar(0));
  p.AddCon(AddRelational(NE,
      AddCall("in_relation", AddVar(0), AddVar(1), 1, 2, 3), AddNum(0)));
  EXPECT_THROW_MSG(Solve(p), ampl::Error,
      "in_relation: the number of arguments 5 is not a multiple of arity 2");
}

TEST_P(SolverTest, InRelationTuple) {
  Problem p;
  p.AddVar(0, 100, ampl::INTEGER);
  p.AddVar(0, 100, ampl::INTEGER);
  p.AddObj(ampl::MIN, AddBinary(OPPLUS, AddVar(0), AddVar(1)));
  p.AddCon(AddRelational(NE,
      AddCall("in_relation", AddVar(0), AddVar(1), 11, 22), AddNum(0)));
  EXPECT_EQ(33, Solve(p).obj_value());
}

TEST_P(SolverTest, InRelationEmptySet) {
  Problem p;
  p.AddVar(0, 100, ampl::INTEGER);
  p.AddCon(AddRelational(NE, AddCall("in_relation", AddVar(0)), AddNum(0)));
  EXPECT_EQ(ampl::INFEASIBLE, Solve(p).solve_code());
}

TEST_P(SolverTest, InRelationNonConstantSetElement) {
  Problem p;
  p.AddVar(0, 100, ampl::INTEGER);
  p.AddCon(AddRelational(NE,
      AddCall("in_relation", AddVar(0), 0, AddVar(0)), AddNum(0)));
  EXPECT_THROW_MSG(Solve(p), ampl::Error,
      "in_relation: argument 3 is not constant");
}

// ----------------------------------------------------------------------------
// Other test

struct EnumValue {
  const char *name;
  IloCP::ParameterValues value;
};

class IlogCPTest : public ::testing::Test, public ampl::ExprBuilder {
 protected:
  IlogCPSolver s;

  int CountIloDistribute();

  SolveResult Solve(Problem &p, const char *stub, const char *opt = nullptr) {
    return SolverTest::Solve(s, p, stub, opt);
  }

  template <typename T>
  static std::string Option(const char *name, T value) {
    return str(fmt::Format("{}={}") << name << value);
  }

  void CheckIntCPOption(const char *option, IloCP::IntParam param,
      int start, int end, int offset = 0, bool accepts_auto = false,
      const EnumValue *values = 0);

  void CheckDblCPOption(const char *option,
      IloCP::NumParam param, double good, double bad);
};

int IlogCPTest::CountIloDistribute() {
  int count = 0;
  for (IloModel::Iterator i(s.cp().getModel()); i.ok(); ++i) {
    if (dynamic_cast<IloDistributeI*>((*i).getImpl()))
      ++count;
  }
  return count;
}

void IlogCPTest::CheckIntCPOption(const char *option,
    IloCP::IntParam param, int start, int end, int offset, bool accepts_auto,
    const EnumValue *values) {
  IloCP cp = s.cp();
  for (int i = start; i <= std::min(end, 9); ++i) {
    if (accepts_auto || values)
      s.SetStrOption(option, c_str(fmt::Format("{}") << i));
    else
      s.SetIntOption(option, i);
    EXPECT_EQ(offset + i, cp.getParameter(param))
      << "Failed option: " << option;
    if (!values) {
      if (accepts_auto)
        EXPECT_EQ(str(fmt::Format("{}") << i), s.GetStrOption(option));
      else
        EXPECT_EQ(i, s.GetIntOption(option));
    }
  }
  if (end != INT_MAX) {
    if (accepts_auto || values) {
      EXPECT_THROW(s.SetStrOption(option, c_str(fmt::Format("{}") << end + 1)),
          InvalidOptionValue);
    } else {
      EXPECT_THROW(s.SetIntOption(option, end + 1), InvalidOptionValue);
    }
  }
  if (accepts_auto) {
    s.SetStrOption(option, "-1");
    EXPECT_EQ(IloCP::Auto, cp.getParameter(param));
    s.SetStrOption(option, "auto");
    EXPECT_EQ(IloCP::Auto, cp.getParameter(param));
    EXPECT_EQ("auto", s.GetStrOption(option));
  }
  int small = start - 1;
  if (accepts_auto && small == -1)
    --small;
  if (accepts_auto || values) {
    EXPECT_THROW(s.SetStrOption(option, c_str(fmt::Format("{}") << small)),
        InvalidOptionValue);
  } else {
    EXPECT_THROW(s.SetIntOption(option, small), InvalidOptionValue);
  }
  if (values) {
    int count = 0;
    for (const EnumValue *v = values; v->name; ++v, ++count) {
      s.SetStrOption(option, v->name);
      EXPECT_EQ(v->value, cp.getParameter(param))
        << "Failed option: " << option;
      EXPECT_EQ(v->name, s.GetStrOption(option)) << "Failed option: " << option;
    }
    EXPECT_EQ(end - start + 1, count);
  }
}

void IlogCPTest::CheckDblCPOption(const char *option,
    IloCP::NumParam param, double good, double bad) {
  s.SetDblOption(option, good);
  EXPECT_EQ(good, s.cp().getParameter(param))
    << "Failed option: " << option;
  EXPECT_EQ(good, s.GetDblOption(option));
  EXPECT_THROW(s.SetDblOption(option, bad), InvalidOptionValue);
}

TEST_F(IlogCPTest, IloArrayCopyingIsCheap) {
  IloIntArray array(s.env());
  array.add(42);
  EXPECT_TRUE(array.getImpl() != nullptr);
  EXPECT_EQ(array.getImpl(), IloIntArray(array).getImpl());
}

TEST_F(IlogCPTest, ConvertSingleNumberOfToIloDistribute) {
  s.use_numberof();
  Problem p;
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddCon(AddRelational(EQ, AddNum(0),
      AddNumberOf(AddNum(42), AddVar(0), AddVar(1))));
  s.Solve(p);
  ASSERT_EQ(1, CountIloDistribute());
}

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithSameValuesToIloDistribute) {
  s.use_numberof();
  Problem p;
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
  Problem p;
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
  Problem p;
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddCon(AddRelational(EQ, AddNum(0),
      AddNumberOf(AddNum(42), AddVar(0), AddVar(1))));
  p.AddCon(AddRelational(EQ, AddNum(0),
      AddNumberOf(AddNum(42), AddVar(1))));
  s.Solve(p);
  ASSERT_EQ(2, CountIloDistribute());
}

struct TestSolutionHandler : ampl::BasicSolutionHandler {
  int num_solutions;
  TestSolutionHandler() : num_solutions(0) {}
  virtual ~TestSolutionHandler() {}
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

TEST_F(IlogCPTest, DefaultSolutionLimit) {
  Problem p;
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddCon(AddAllDiff(AddVar(0), AddVar(1), AddVar(2)));
  TestSolutionHandler sh;
  s.set_solution_handler(&sh);
  s.SetIntOption("solutionlimit", 100);
  s.Solve(p);
  EXPECT_EQ(1, sh.num_solutions);
}

TEST_F(IlogCPTest, CPOptimizerDoesntSupportContinuousVars) {
  Problem p;
  p.AddVar(0, 1);
  p.AddObj(ampl::MIN, AddVar(0));
  EXPECT_THROW(s.Solve(p), ampl::Error);
}

// ----------------------------------------------------------------------------
// Option tests

TEST_F(IlogCPTest, OptimizerOption) {
  EXPECT_EQ("auto", s.GetStrOption("optimizer"));
  Problem p;
  p.AddVar(1, 2, ampl::INTEGER);
  p.AddVar(1, 2, ampl::INTEGER);
  p.AddCon(AddAllDiff(AddVar(0), AddVar(1)));
  s.SetStrOption("optimizer", "cp");
  s.Solve(p);
  EXPECT_EQ(0, p.solve_code());
  s.SetStrOption("optimizer", "cplex");
  EXPECT_THROW(s.Solve(p), ampl::UnsupportedExprError);
  s.SetStrOption("optimizer", "auto");
  s.Solve(p);
  EXPECT_EQ(0, p.solve_code());
}

// TODO: move to solver_test
TEST_P(SolverTest, ObjnoOption) {
  EXPECT_EQ(1, solver_->GetIntOption("objno"));
  Problem p;
  p.AddVar(11, 22, ampl::INTEGER);
  ampl::Variable var = AddVar(0);
  p.AddObj(ampl::MIN, var);
  p.AddObj(ampl::MAX, var);
  EXPECT_EQ(11, Solve(p).obj_value());
  solver_->SetIntOption("objno", 0);
  EXPECT_EQ(0, solver_->GetIntOption("objno"));
  EXPECT_EQ(0, Solve(p).obj_value());
  solver_->SetIntOption("objno", 2);
  EXPECT_EQ(2, solver_->GetIntOption("objno"));
  EXPECT_EQ(22, Solve(p).obj_value());
  solver_->SetStrOption("optimizer", "cplex");
  EXPECT_EQ(2, solver_->GetIntOption("objno"));
  EXPECT_EQ(22, Solve(p).obj_value());
  EXPECT_THROW(solver_->SetIntOption("objno", -1), InvalidOptionValue);
  solver_->SetIntOption("objno", 3);
  EXPECT_THROW(Solve(p), InvalidOptionValue);
}

TEST_F(IlogCPTest, UseCplexForLinearProblem) {
  Problem p;
  p.AddVar(1, 2);
  p.AddObj(ampl::MIN, AddNum(42));
  s.Solve(p);
  EXPECT_EQ(0, p.solve_code());
}

TEST_F(IlogCPTest, DebugExprOption) {
  s.SetIntOption("debugexpr", 0);
  EXPECT_EQ(0, s.GetOption(IlogCPSolver::DEBUGEXPR));
  EXPECT_EQ(0, s.GetIntOption("debugexpr"));
  s.SetIntOption("debugexpr", 1);
  EXPECT_EQ(1, s.GetOption(IlogCPSolver::DEBUGEXPR));
  EXPECT_EQ(1, s.GetIntOption("debugexpr"));
  EXPECT_THROW(s.SetIntOption("debugexpr", 42), InvalidOptionValue);
  EXPECT_THROW(s.SetStrOption("debugexpr", "oops"), OptionError);
}

TEST_F(IlogCPTest, UseNumberOfOption) {
  s.SetIntOption("usenumberof", 0);
  EXPECT_EQ(0, s.GetOption(IlogCPSolver::USENUMBEROF));
  EXPECT_EQ(0, s.GetIntOption("usenumberof"));
  s.SetIntOption("usenumberof", 1);
  EXPECT_EQ(1, s.GetOption(IlogCPSolver::USENUMBEROF));
  EXPECT_EQ(1, s.GetIntOption("usenumberof"));
  EXPECT_THROW(s.SetIntOption("usenumberof", 42), InvalidOptionValue);
  EXPECT_THROW(s.SetStrOption("usenumberof", "oops"), OptionError);
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
  EXPECT_EQ(IloCP::Quiet, s.cp().getParameter(IloCP::LogVerbosity));
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
  CheckDblCPOption("timelimit", IloCP::TimeLimit, 42, -1);
  if (CPX_VERSION > 1220)
    CheckIntCPOption("workers", IloCP::Workers, 0, INT_MAX, 0, true);
  else
    CheckIntCPOption("workers", IloCP::Workers, 1, 4, 0, false);
}

TEST_P(SolverTest, MultiObjOption) {
  Problem p;
  p.AddVar(0, 10, ampl::INTEGER);
  p.AddObj(ampl::MIN, AddBinary(OPREM, AddVar(0), AddNum(3)));
  p.AddObj(ampl::MIN,
      AddUnary(OP2POW, AddBinary(OPMINUS, AddVar(0), AddNum(5))));
  EXPECT_EQ(0, Solve(p));
  solver_->SetIntOption("multiobj", 1);
  EXPECT_EQ(6, Solve(p));
  solver_->SetIntOption("multiobj", 0);
  EXPECT_EQ(0, Solve(p));
}

TEST_F(IlogCPTest, SolutionLimitOption) {
  EXPECT_EQ(-1, s.GetOption(IlogCPSolver::SOLUTION_LIMIT));
  s.SetIntOption("solutionlimit", 0);
  EXPECT_EQ(0, s.GetOption(IlogCPSolver::SOLUTION_LIMIT));
  EXPECT_EQ(0, s.GetIntOption("solutionlimit"));
  s.SetIntOption("solutionlimit", 42);
  EXPECT_EQ(42, s.GetOption(IlogCPSolver::SOLUTION_LIMIT));
  EXPECT_EQ(42, s.GetIntOption("solutionlimit"));
  EXPECT_THROW(s.SetIntOption("solutionlimit", -1), InvalidOptionValue);
  EXPECT_THROW(s.SetStrOption("solutionlimit", "oops"), OptionError);
}

TEST_F(IlogCPTest, MIPDisplayOption) {
  EXPECT_EQ(0, s.GetIntOption("mipdisplay"));
  EXPECT_EQ(0, s.cplex().getParam(IloCplex::MIPDisplay));
  s.SetIntOption("mipdisplay", 4);
  EXPECT_EQ(4, s.GetIntOption("mipdisplay"));
  EXPECT_EQ(4, s.cplex().getParam(IloCplex::MIPDisplay));
  EXPECT_THROW(s.SetIntOption("mipdisplay", -1), InvalidOptionValue);
  EXPECT_THROW(s.SetStrOption("mipdisplay", "oops"), OptionError);
}

TEST_F(IlogCPTest, MIPIntervalOption) {
  EXPECT_EQ(0, s.GetIntOption("mipinterval"));
  EXPECT_EQ(0, s.cplex().getParam(IloCplex::MIPInterval));
  s.SetIntOption("mipinterval", 42);
  EXPECT_EQ(42, s.GetIntOption("mipinterval"));
  EXPECT_EQ(42, s.cplex().getParam(IloCplex::MIPInterval));
  EXPECT_THROW(s.SetStrOption("mipinterval", "oops"), OptionError);
}

// ----------------------------------------------------------------------------
// Interrupt tests

#ifdef HAVE_THREADS
TEST_P(SolverTest, CPInterruptSolution) {
  std::thread t(Interrupt);
  Problem p;
  solver_->SetStrOption("optimizer", "cp");
  string message = Solve(p, "miplib/assign1").message;
  t.join();
  EXPECT_EQ(600, p.solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}

TEST_P(SolverTest, CPLEXInterruptSolution) {
  std::thread t(Interrupt);
  Problem p;
  solver_->SetStrOption("optimizer", "cplex");
  string message = Solve(p, "miplib/assign1").message;
  t.join();
  EXPECT_EQ(600, p.solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}
#endif
}


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

#include "ilogcp/ilogcp.h"
#include "mp/nl.h"
#include "feature.h"

typedef mp::IlogCPSolver Solver;
enum {FEATURES = feature::ALL};

#include "nl-solver-test.h"

using mp::Expr;
using mp::NumericExpr;
using mp::IlogCPSolver;
using mp::InvalidOptionValue;
using mp::OptionError;
using mp::Problem;

namespace var = mp::var;
namespace obj = mp::obj;

using std::string;

class FunctionTest : public NLSolverTest {
 protected:
  template <typename IndexFactory>
  struct CallFactory : NumericExprFactory {
    int num_args_;
    IndexFactory factory_;

    CallFactory(int num_args, IndexFactory factory)
      : num_args_(num_args), factory_(factory) {}

    NumericExpr Create(ProblemBuilder &pb) const {
      auto args = pb.BeginCall(0, num_args_);
      for (int i = 1; i < num_args_; ++i)
        args.AddArg(pb.MakeNumericConstant(11 * i));
      args.AddArg(factory_(pb));
      return pb.EndCall(args);
    }
  };

  template <typename IndexFactory>
  CallFactory<IndexFactory> MakeCallFactory(
      int num_args, IndexFactory factory) {
    return CallFactory<IndexFactory>(num_args, factory);
  }

  virtual void SetInfo(ProblemBuilder &pb, mp::ProblemInfo &info) {
    info.num_funcs = 2;
    NLSolverTest::SetInfo(pb, info);
    // Create functions permitting less arguments than necessary.
    // This is done to be able to test calls with invalid arguments.
    pb.SetFunction(0, "element", -2, mp::func::NUMERIC);
    pb.SetFunction(1, "in_relation", -1, mp::func::NUMERIC);
  }
};

// ----------------------------------------------------------------------------
// element constraint tests

TEST_F(FunctionTest, ElementConstraint) {
  EXPECT_EQ(22, Eval(MakeCallFactory(4, VariableFactory(1)), 1));
}

TEST_F(FunctionTest, TooFewArgsToElementConstraint) {
  EXPECT_THROW_MSG(Eval(MakeCallFactory(1, VariableFactory(1)), 0),
      mp::Error, "element: too few arguments");
}

TEST_F(FunctionTest, ElementConstantIndexOutOfBounds) {
  EXPECT_THROW_MSG(Eval(MakeCallFactory(3, NumericConstantFactory(2))),
        mp::Error, "element: index 2 is out of bounds");
}

TEST_F(FunctionTest, ElementAtConstantIndex) {
  EXPECT_EQ(22, Eval(MakeCallFactory(3, NumericConstantFactory(1))));
}

TEST_F(FunctionTest, ElementExprAtConstantIndex) {
  struct Factory : NumericExprFactory {
    NumericExpr Create(ProblemBuilder &pb) const {
      auto args = pb.BeginCall(0, 3);
      args.AddArg(x);
      args.AddArg(pb.MakeNumericConstant(22));
      args.AddArg(pb.MakeNumericConstant(0));
      return pb.EndCall(args);
    }
  } factory;
  EXPECT_EQ(42, Eval(factory, 42));
}

TEST_F(FunctionTest, ElementExprPlusConstantAtConstantIndex) {
  struct Factory : NumericExprFactory {
    NumericExpr Create(ProblemBuilder &pb) const {
      auto args = pb.BeginCall(0, 3);
      args.AddArg(pb.MakeNumericConstant(11));
      args.AddArg(pb.MakeBinary(mp::expr::ADD, x, pb.MakeNumericConstant(2)));
      args.AddArg(pb.MakeNumericConstant(1));
      return pb.EndCall(args);
    }
  } factory;
  EXPECT_EQ(42, Eval(factory, 40));
}

TEST_F(FunctionTest, ElementVariableIndexOutOfBounds) {
  EXPECT_EQ(mp::sol::INFEASIBLE,
            Eval(MakeCallFactory(3, VariableFactory(1)), 2).solve_code());
}

TEST_F(FunctionTest, ElementConstantAtVariableIndex) {
  EXPECT_EQ(22, Eval(MakeCallFactory(3, VariableFactory(1)), 1));
}

TEST_F(FunctionTest, ElementExprAtVariableIndex) {
  struct Factory : NumericExprFactory {
    NumericExpr Create(ProblemBuilder &pb) const {
      auto args = pb.BeginCall(0, 3);
      args.AddArg(x);
      args.AddArg(pb.MakeNumericConstant(22));
      args.AddArg(y);
      return pb.EndCall(args);
    }
  } factory;
  EXPECT_EQ(42, Eval(factory, 42, 0));
}

TEST_F(FunctionTest, ElementExprPlusConstantAtVariableIndex) {
  struct Factory : NumericExprFactory {
    NumericExpr Create(ProblemBuilder &pb) const {
      auto args = pb.BeginCall(0, 3);
      args.AddArg(pb.MakeNumericConstant(11));
      args.AddArg(pb.MakeBinary(mp::expr::ADD, x, pb.MakeNumericConstant(2)));
      args.AddArg(y);
      return pb.EndCall(args);
    }
  } factory;
  EXPECT_EQ(42, Eval(factory, 40, 1));
}

// ----------------------------------------------------------------------------
// in_relation constraint tests

// Makes a problem for testing in_relation constraint.
void MakeInRelationProblem(mp::internal::ASLBuilder &pb) {
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_nl_integer_vars_in_cons = 1;
  info.num_objs = info.num_logical_cons = info.num_funcs = 1;
  pb.SetInfo(info);
  pb.SetFunction(0, "in_relation", -1, mp::func::NUMERIC);
  pb.AddVar(0, 100, var::INTEGER);
  pb.AddObj(obj::MIN, NumericExpr(), 0).AddTerm(0, 1);
}

TEST_F(FunctionTest, InRelationConstraint) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeInRelationProblem(pb);
  auto args = pb.BeginCall(0, 2);
  args.AddArg(pb.MakeVariable(0));
  args.AddArg(pb.MakeNumericConstant(42));
  pb.AddCon(pb.MakeRelational(
              mp::expr::NE, pb.EndCall(args), pb.MakeNumericConstant(0)));
  EXPECT_EQ(42, Solve(pb).obj_value());
}

TEST_F(FunctionTest, NestedInRelationNotSupported) {
  struct Factory : NumericExprFactory {
    NumericExpr Create(ProblemBuilder &pb) const {
      auto args = pb.BeginCall(1, 2);
      args.AddArg(x);
      args.AddArg(pb.MakeNumericConstant(42));
      return pb.MakeBinary(mp::expr::ADD, pb.EndCall(args), one);
    }
  } factory;
  EXPECT_THROW_MSG(Eval(factory), mp::UnsupportedExprError,
      "unsupported expression: nested 'in_relation'");
}

TEST_F(FunctionTest, TooFewArgsToInRelationConstraint) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeInRelationProblem(pb);
  auto args = pb.BeginCall(0, 0);
  pb.AddCon(pb.MakeRelational(
              mp::expr::NE, pb.EndCall(args), pb.MakeNumericConstant(0)));
  EXPECT_THROW_MSG(Solve(pb), mp::Error, "in_relation: too few arguments");
}

void MakeInRelationProblem2(mp::internal::ASLBuilder &pb, int num_const_args) {
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_nl_integer_vars_in_cons = 2;
  info.num_objs = info.num_logical_cons = info.num_funcs = 1;
  pb.SetInfo(info);
  pb.SetFunction(0, "in_relation", -1, mp::func::NUMERIC);
  pb.AddVar(0, 100, var::INTEGER);
  pb.AddVar(0, 100, var::INTEGER);
  auto obj = pb.AddObj(obj::MIN, NumericExpr(), 2);
  obj.AddTerm(0, 1);
  obj.AddTerm(1, 1);
  auto args = pb.BeginCall(0, num_const_args + 2);
  args.AddArg(pb.MakeVariable(0));
  args.AddArg(pb.MakeVariable(1));
  for (int i = 1; i <= num_const_args; ++i)
    args.AddArg(pb.MakeNumericConstant(11 * i));
  pb.AddCon(pb.MakeRelational(
              mp::expr::NE, pb.EndCall(args), pb.MakeNumericConstant(0)));
}

TEST_F(FunctionTest, InRelationSizeIsNotMultipleOfArity) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeInRelationProblem2(pb, 3);
  EXPECT_THROW_MSG(Solve(pb), mp::Error,
      "in_relation: the number of arguments 5 is not a multiple of arity 2");
}

TEST_F(FunctionTest, InRelationTuple) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeInRelationProblem2(pb, 2);
  EXPECT_EQ(33, Solve(pb).obj_value());
}

TEST_F(FunctionTest, InRelationEmptySet) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeInRelationProblem(pb);
  auto args = pb.BeginCall(0, 1);
  args.AddArg(pb.MakeVariable(0));
  pb.AddCon(pb.MakeRelational(
              mp::expr::NE, pb.EndCall(args), pb.MakeNumericConstant(0)));
  EXPECT_EQ(mp::sol::INFEASIBLE, Solve(pb).solve_code());
}

TEST_F(FunctionTest, InRelationNonConstantSetElement) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeInRelationProblem(pb);
  auto args = pb.BeginCall(0, 3);
  args.AddArg(pb.MakeVariable(0));
  args.AddArg(pb.MakeNumericConstant(0));
  args.AddArg(pb.MakeVariable(0));
  pb.AddCon(pb.MakeRelational(
              mp::expr::NE, pb.EndCall(args), pb.MakeNumericConstant(0)));
  EXPECT_THROW_MSG(Solve(pb), mp::Error,
      "in_relation: argument 3 is not constant");
}

// ----------------------------------------------------------------------------
// Other test

struct EnumValue {
  const char *name;
  IloCP::ParameterValues value;
};

class IlogCPTest : public ::testing::Test, public mp::internal::ASLBuilder {
 protected:
  IlogCPSolver s;

  IlogCPTest() {
    set_flags(mp::internal::ASL_STANDARD_OPCODES);
    mp::ProblemInfo info = mp::ProblemInfo();
    info.num_vars = 3;
    info.num_objs = 1;
    SetInfo(info);
  }

  EvalResult Solve(mp::Problem &p) {
    TestSolutionHandler sh(1);
    s.Solve(p, sh);
    const double *sol = sh.primal();
    EvalResult result(sol ? sol[0] : 0, sh.obj_value());
    result.set_solve_code(sh.status());
    return result;
  }

  int CountIloDistribute();

  mp::NumericConstant MakeConst(double value) {
    return MakeNumericConstant(value);
  }

  template <typename T>
  static std::string Option(const char *name, T value) {
    return fmt::format("{}={}", name, value);
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
      s.SetStrOption(option, fmt::format("{}", i));
    else
      s.SetIntOption(option, i);
    EXPECT_EQ(offset + i, cp.getParameter(param))
      << "Failed option: " << option;
    if (!values) {
      if (accepts_auto)
        EXPECT_EQ(fmt::format("{}", i), s.GetStrOption(option));
      else
        EXPECT_EQ(i, s.GetIntOption(option));
    }
  }
  if (end != INT_MAX) {
    if (accepts_auto || values) {
      EXPECT_THROW(s.SetStrOption(option, fmt::format("{}", end + 1)),
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
    EXPECT_THROW(s.SetStrOption(option, fmt::format("{}", small)),
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
  EXPECT_TRUE(array.getImpl() != 0);
  EXPECT_EQ(array.getImpl(), IloIntArray(array).getImpl());
}

TEST_F(IlogCPTest, ConvertSingleNumberOfToIloDistribute) {
  s.use_numberof();
  Problem p;
  p.AddVar(0, 0, var::INTEGER);
  p.AddVar(0, 0, var::INTEGER);
  NumericExpr args[] = {MakeConst(42), MakeVariable(0), MakeVariable(1)};
  p.AddCon(MakeRelational(mp::expr::EQ, MakeConst(0), MakeNumberOf(args)));
  mp::BasicSolutionHandler sh;
  s.Solve(p, sh);
  ASSERT_EQ(1, CountIloDistribute());
}

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithSameValuesToIloDistribute) {
  s.use_numberof();
  Problem p;
  p.AddVar(0, 0, var::INTEGER);
  p.AddVar(0, 0, var::INTEGER);
  NumericExpr args[] = {MakeConst(42), MakeVariable(0), MakeVariable(1)};
  p.AddCon(MakeRelational(mp::expr::EQ, MakeConst(0), MakeNumberOf(args)));
  p.AddCon(MakeRelational(mp::expr::EQ, MakeConst(0), MakeNumberOf(args)));
  mp::BasicSolutionHandler sh;
  s.Solve(p, sh);
  ASSERT_EQ(1, CountIloDistribute());
}

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithDiffValuesToIloDistribute) {
  s.use_numberof();
  Problem p;
  p.AddVar(0, 0, var::INTEGER);
  p.AddVar(0, 0, var::INTEGER);
  NumericExpr args[] = {MakeConst(42), MakeVariable(0), MakeVariable(1)};
  p.AddCon(MakeRelational(mp::expr::EQ, MakeConst(0), MakeNumberOf(args)));
  args[0] = MakeConst(43);
  p.AddCon(MakeRelational(mp::expr::EQ, MakeConst(0), MakeNumberOf(args)));
  mp::BasicSolutionHandler sh;
  s.Solve(p, sh);
  ASSERT_EQ(1, CountIloDistribute());
}

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithDiffExprs) {
  s.use_numberof();
  Problem p;
  p.AddVar(0, 0, var::INTEGER);
  p.AddVar(0, 0, var::INTEGER);
  NumericExpr args[] = {MakeConst(42), MakeVariable(0), MakeVariable(1)};
  p.AddCon(MakeRelational(mp::expr::EQ, MakeConst(0), MakeNumberOf(args)));
  NumericExpr args2[] = {MakeConst(42), MakeVariable(1)};
  p.AddCon(MakeRelational(mp::expr::EQ, MakeConst(0), MakeNumberOf(args2)));
  mp::BasicSolutionHandler sh;
  s.Solve(p, sh);
  ASSERT_EQ(2, CountIloDistribute());
}

TEST_F(IlogCPTest, DefaultSolutionLimit) {
  Problem p;
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeVariable(2)};
  p.AddCon(MakeAllDiff(args));
  s.SetIntOption("solutionlimit", 100);

  struct TestSolutionHandler : mp::BasicSolutionHandler {
    int num_solutions;
    TestSolutionHandler() : num_solutions(0) {}
    void HandleSolution(int, fmt::StringRef,
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
  } sh;
  s.Solve(p, sh);
  EXPECT_EQ(1, sh.num_solutions);
}

TEST_F(IlogCPTest, CPOptimizerDoesntSupportContinuousVars) {
  Problem p;
  p.AddVar(0, 1);
  p.AddObj(mp::obj::MIN, MakeVariable(0));
  mp::BasicSolutionHandler sh;
  EXPECT_THROW(s.Solve(p, sh), mp::Error);
}

// ----------------------------------------------------------------------------
// Option tests

TEST_F(IlogCPTest, OptimizerOption) {
  EXPECT_EQ("auto", s.GetStrOption("optimizer"));
  Problem p;
  p.AddVar(1, 2, var::INTEGER);
  p.AddVar(1, 2, var::INTEGER);
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1)};
  p.AddCon(MakeAllDiff(args));
  s.SetStrOption("optimizer", "cp");
  mp::BasicSolutionHandler sh;
  s.Solve(p, sh);
  EXPECT_EQ(0, p.solve_code());
  s.SetStrOption("optimizer", "cplex");
  EXPECT_THROW(s.Solve(p, sh), mp::UnsupportedExprError);
  s.SetStrOption("optimizer", "auto");
  s.Solve(p, sh);
  EXPECT_EQ(0, p.solve_code());
}

TEST_F(IlogCPTest, UseCplexForLinearProblem) {
  Problem p;
  p.AddVar(1, 2);
  p.AddObj(obj::MIN, MakeConst(42));
  mp::BasicSolutionHandler sh;
  s.Solve(p, sh);
  EXPECT_EQ(0, p.solve_code());
}

TEST_F(NLSolverTest, ObjnoOptionCP) {
  solver_.SetStrOption("optimizer", "cp");
  solver_.SetIntOption("objno", 2);
  EXPECT_EQ(22, SolveObjNoProblem());
}

TEST_F(NLSolverTest, ObjnoOptionCPLEX) {
  solver_.SetStrOption("optimizer", "cplex");
  solver_.SetIntOption("objno", 2);
  EXPECT_EQ(22, SolveObjNoProblem());
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
      {0,     IloCP::Off}
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
      {0,          IloCP::Off}
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
      {0,         IloCP::Off}
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
      {0,            IloCP::Off}
  };
  CheckIntCPOption("searchtype", IloCP::SearchType,
      0, 2, IloCP::DepthFirst, CPX_VERSION > 1220, types);
}

TEST_F(IlogCPTest, CPTimeModeOption) {
  const EnumValue modes[] = {
      {"cputime",     IloCP::CPUTime},
      {"elapsedtime", IloCP::ElapsedTime},
      {0,             IloCP::Off}
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

TEST_F(IlogCPTest, MultiObjOption) {
  ASLBuilder pb(s.GetProblemBuilder(""));
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_nl_integer_vars_in_objs = 1;
  info.num_objs = info.num_nl_objs = 2;
  pb.SetInfo(info);
  pb.AddVar(0, 10, var::INTEGER);
  pb.AddObj(mp::obj::MIN, pb.MakeBinary(
              mp::expr::MOD, pb.MakeVariable(0), pb.MakeNumericConstant(3)), 0);
  pb.AddObj(mp::obj::MIN, pb.MakeUnary(
              mp::expr::POW2, pb.MakeBinary(mp::expr::SUB, pb.MakeVariable(0),
                                            pb.MakeNumericConstant(5))), 0);
  Problem p(pb.GetProblem());
  EXPECT_EQ(0, Solve(p));
  s.SetIntOption("multiobj", 1);
  EXPECT_EQ(6, Solve(p));
  s.SetIntOption("multiobj", 0);
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

TEST_F(NLSolverTest, InterruptCP) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeTSP(pb);
  pb.EndBuild();
  solver_.SetStrOption("optimizer", "cp");
  TestInterrupter interrupter(solver_);
  TestSolutionHandler sh;
  solver_.Solve(pb, sh);
  EXPECT_EQ(600, sh.status());
  EXPECT_TRUE(sh.message().find("interrupted") != std::string::npos);
}

TEST_F(NLSolverTest, InterruptCPLEX) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeTSP(pb);
  pb.EndBuild();
  solver_.SetStrOption("optimizer", "cplex");
  TestInterrupter interrupter(solver_);
  TestSolutionHandler sh;
  solver_.Solve(pb, sh);
  EXPECT_EQ(600, sh.status());
  EXPECT_TRUE(sh.message().find("interrupted") != std::string::npos);
}

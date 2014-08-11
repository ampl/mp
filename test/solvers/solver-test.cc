/*
 Solver test suite.

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

#include "solver-test.h"
#include "../util.h"

#include "mp/nl.h"

#include <cmath>

using std::string;

using mp::LogicalExpr;
using mp::NumericExpr;
using mp::Problem;
using mp::UnsupportedExprError;
namespace var = mp::var;
namespace obj = mp::obj;

#ifndef MP_TEST_DATA_DIR
# define MP_TEST_DATA_DIR "../data"
#endif

SolverTest::EvalResult SolverTest::Solve(Problem &p) {
  struct TestSolutionHandler : mp::BasicSolutionHandler {
    EvalResult result;
    virtual ~TestSolutionHandler() {}
    void HandleSolution(mp::Problem &p, fmt::StringRef,
          const double *values, const double *, double obj_value) {
      result = values ? EvalResult(values[0], obj_value, p.solve_code())
          : EvalResult(p.solve_code());
    }
  };
  TestSolutionHandler sh;
  solver_->set_solution_handler(&sh);
  solver_->Solve(p);
  return sh.result;
}

SolverTest::EvalResult SolverTest::Solve(
    LogicalExpr e, int var1, int var2, int var3, bool need_result) {
  Problem p;
  p.AddVar(need_result ? negInfinity : 0,
      need_result ? Infinity : 0, var::INTEGER);
  p.AddVar(var1, var1, var::INTEGER);
  p.AddVar(var2, var2, var::INTEGER);
  p.AddVar(var3, var3, var::INTEGER);
  p.AddCon(e);
  return Solve(p);
}

SolverTest::SolverTest()
: solver_(GetParam().create_solver()), features_(GetParam().features) {
  mp::NLHeader header = {};
  header.num_vars = 4;
  header.num_objs = 1;
  header.num_funcs = 2;
  BeginBuild("", header);
  x = MakeVariable(1);
  y = MakeVariable(2);
  z = MakeVariable(3);
}

SolveResult SolverTest::Solve(
    mp::Solver &s, Problem &p, const char *stub, const char *opt) {
  TestSolutionHandler sh;
  s.set_solution_handler(&sh);
  const std::string DATA_DIR = MP_TEST_DATA_DIR "/";
  if (s.ProcessArgs(Args(s.name(), "-s", (DATA_DIR + stub).c_str(), opt), p))
    s.Solve(p);
  const string &message = sh.message();
  int solve_code = sh.solve_code();
  EXPECT_GE(solve_code, 0);
  bool solved = true;
  if (solve_code < 100) {
    EXPECT_TRUE(message.find("optimal solution") != string::npos ||
        (p.num_objs() == 0 &&
            message.find("feasible solution") != string::npos));
  } else
    solved = false;
  return SolveResult(solved, sh.obj_value(), message);
}

TEST_P(SolverTest, Plus) {
  NumericExpr e = MakeBinary(OPPLUS, x, y);
  EXPECT_EQ(25, Eval(e, 10, 15));
  EXPECT_EQ(12, Eval(e, 19, -7));
}

TEST_P(SolverTest, Minus) {
  NumericExpr e = MakeBinary(OPMINUS, x, y);
  EXPECT_EQ(-5, Eval(e, 10, 15));
  EXPECT_EQ(26, Eval(e, 19, -7));
}

TEST_P(SolverTest, Mult) {
  NumericExpr e = MakeBinary(OPMULT, x, y);
  EXPECT_EQ(150, Eval(e, 10, 15));
  EXPECT_EQ(-133, Eval(e, 19, -7));
}

TEST_P(SolverTest, Div) {
  NumericExpr e = MakeBinary(OPDIV, x, y);
  if (!HasFeature(feature::DIV)) {
    EXPECT_THROW(Eval(e, 150, 15);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(10, Eval(e, 150, 15));
  EXPECT_EQ(-7, Eval(e, -133, 19));
}

TEST_P(SolverTest, Rem) {
  NumericExpr e = MakeBinary(OPREM, x, y);
  EXPECT_EQ(0, Eval(e, 9, 3));
  EXPECT_EQ(2, Eval(e, 8, 3));
  EXPECT_EQ(-2, Eval(e, -8, 3));
  EXPECT_EQ(2, Eval(e, 8, -3));
  EXPECT_EQ(-2, Eval(e, -8, -3));
}

TEST_P(SolverTest, Pow) {
  NumericExpr e = MakeBinary(OPPOW, x, y);
  if (!HasFeature(feature::POW)) {
    EXPECT_THROW(Eval(e, 2, 3);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(8, Eval(e, 2, 3));
  EXPECT_EQ(81, Eval(e, 3, 4));
}

TEST_P(SolverTest, NumericLess) {
  NumericExpr e = MakeBinary(OPLESS, x, y);
  EXPECT_EQ(0, Eval(e, 10, 15));
  EXPECT_EQ(26, Eval(e, 19, -7));
}

TEST_P(SolverTest, Min) {
  NumericExpr args[] = {x, y, z};
  NumericExpr e = MakeVarArg(MINLIST, args);
  EXPECT_EQ(-7, Eval(e, 3, -7, 5));
  EXPECT_EQ(10, Eval(e, 10, 20, 30));
}

TEST_P(SolverTest, Max) {
  NumericExpr args[] = {x, y, z};
  NumericExpr e = MakeVarArg(MAXLIST, args);
  EXPECT_EQ(5, Eval(e, 3, -7, 5));
  EXPECT_EQ(30, Eval(e, 30, 20, 10));
}

TEST_P(SolverTest, Floor) {
  NumericExpr e = MakeUnary(FLOOR, x);
  EXPECT_EQ(-42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
  if (!HasFeature(feature::FLOAT_CONST)) return;
  EXPECT_EQ(4, Eval(MakeUnary(FLOOR, MakeConst(4.9))));
  EXPECT_EQ(-5, Eval(MakeUnary(FLOOR, MakeConst(-4.1))));
}

TEST_P(SolverTest, Ceil) {
  NumericExpr e = MakeUnary(CEIL, x);
  EXPECT_EQ(-42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
  if (!HasFeature(feature::FLOAT_CONST)) return;
  EXPECT_EQ(5, Eval(MakeUnary(CEIL, MakeConst(4.1))));
  EXPECT_EQ(-4, Eval(MakeUnary(CEIL, MakeConst(-4.9))));
}

TEST_P(SolverTest, Abs) {
  NumericExpr e = MakeUnary(ABS, x);
  EXPECT_EQ(42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
}

TEST_P(SolverTest, UnaryMinus) {
  NumericExpr e = MakeUnary(OPUMINUS, x);
  EXPECT_EQ(42, Eval(e, -42));
  EXPECT_EQ(-42, Eval(e, 42));
}

TEST_P(SolverTest, If) {
  NumericExpr e = MakeIf(MakeRelational(EQ, x, MakeConst(1)), y, z);
  EXPECT_EQ(42, Eval(e, 1, 42, 10));
  EXPECT_EQ(10, Eval(e, 0, 42, 10));
  EXPECT_EQ(42, Eval(e, 1, 42, 42));
}

TEST_P(SolverTest, Tanh) {
  double arg = (std::log(1.5) - std::log(0.5)) / 2;
  NumericExpr e = MakeBinary(OPMULT,
      MakeConst(2), MakeUnary(OP_tanh, MakeConst(arg)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(1, Eval(e));
}

TEST_P(SolverTest, Tan) {
  EXPECT_THROW(Eval(MakeUnary(OP_tan, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Sqrt) {
  NumericExpr e = MakeUnary(OP_sqrt, x);
  if (!HasFeature(feature::SQRT)) {
    EXPECT_THROW(Eval(e, 64);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(8, Eval(e, 64));
}

TEST_P(SolverTest, Sinh) {
  NumericExpr e = MakeUnary(OP_sinh, MakeConst(std::log(2 + std::sqrt(5.0))));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(2, Eval(e));
}

TEST_P(SolverTest, Sin) {
  EXPECT_THROW(Eval(MakeUnary(OP_sin, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Log10) {
  NumericExpr e = MakeUnary(OP_log10, MakeConst(1000));
  if (!HasFeature(feature::LOG)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(3, Eval(e));
}

TEST_P(SolverTest, Log) {
  NumericExpr e = MakeUnary(OP_log, MakeConst(std::exp(5.0)));
  if (!HasFeature(feature::LOG)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverTest, Exp) {
  NumericExpr e = MakeUnary(OP_exp, MakeConst(std::log(5.0)));
  if (!HasFeature(feature::EXP)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverTest, Cosh) {
  double x = 5;
  NumericExpr e = MakeUnary(OP_cosh,
      MakeConst(std::log(x + std::sqrt(x + 1) * std::sqrt(x - 1))));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverTest, Cos) {
  EXPECT_THROW(Eval(MakeUnary(OP_cos, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Atanh) {
  mp::UnaryExpr x = MakeUnary(OP_atanh, MakeConst(std::tanh(5.0)));
  NumericExpr e = MakeUnary(FLOOR, MakeBinary(OPPLUS,
      MakeConst(0.5), MakeBinary(OPMULT, MakeConst(1000000), x)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5000000, Eval(e));
}

TEST_P(SolverTest, Atan2) {
  EXPECT_THROW(Eval(MakeBinary(OP_atan2, x, y)), UnsupportedExprError);
}

TEST_P(SolverTest, Atan) {
  EXPECT_THROW(Eval(MakeUnary(OP_atan, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Asinh) {
  NumericExpr e = MakeUnary(OP_asinh, MakeConst(std::sinh(5.0)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverTest, Asin) {
  EXPECT_THROW(Eval(MakeUnary(OP_asin, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Acosh) {
  NumericExpr e = MakeUnary(OP_acosh, MakeConst(std::cosh(5.0)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverTest, Acos) {
  EXPECT_THROW(Eval(MakeUnary(OP_acos, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Sum) {
  NumericExpr args[] = {x, y, z};
  using mp::MakeArrayRef;
  EXPECT_EQ(0, Eval(MakeSum(MakeArrayRef(args, 0))));
  EXPECT_EQ(42, Eval(MakeSum(MakeArrayRef(args, 1)), 42));
  EXPECT_EQ(123, Eval(MakeSum(args), 100, 20, 3));
}

TEST_P(SolverTest, IntDiv) {
  NumericExpr e = MakeBinary(OPintDIV, x, y);
  EXPECT_EQ(3, Eval(e, 9, 3));
  EXPECT_EQ(2, Eval(e, 8, 3));
  EXPECT_EQ(-2, Eval(e, -8, 3));
  EXPECT_EQ(-2, Eval(e, 8, -3));
  EXPECT_EQ(2, Eval(e, -8, -3));
}

TEST_P(SolverTest, Precision) {
  EXPECT_THROW(Eval(MakeBinary(OPprecision, x, y)), UnsupportedExprError);
}

TEST_P(SolverTest, Round) {
  EXPECT_EQ(42, Eval(MakeBinary(OPround, x, MakeConst(0)), 42));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(4, Eval(MakeBinary(OPround, MakeConst(4.4), MakeConst(0))));
    EXPECT_EQ(5, Eval(MakeBinary(OPround, MakeConst(4.6), MakeConst(0))));
    EXPECT_EQ(-4, Eval(MakeBinary(OPround, MakeConst(-4.4), MakeConst(0))));
    EXPECT_EQ(-5, Eval(MakeBinary(OPround, MakeConst(-4.6), MakeConst(0))));
  }
  EXPECT_THROW(Eval(MakeBinary(OPround, x, MakeConst(1))), UnsupportedExprError);
  EXPECT_THROW(Eval(MakeBinary(OPround, x, y)), UnsupportedExprError);
}

TEST_P(SolverTest, Trunc) {
  EXPECT_EQ(42, Eval(MakeBinary(OPtrunc, x, MakeConst(0)), 42));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(4, Eval(MakeBinary(OPtrunc, MakeConst(4.4), MakeConst(0))));
    EXPECT_EQ(4, Eval(MakeBinary(OPtrunc, MakeConst(4.6), MakeConst(0))));
    EXPECT_EQ(-4, Eval(MakeBinary(OPtrunc, MakeConst(-4.4), MakeConst(0))));
    EXPECT_EQ(-4, Eval(MakeBinary(OPtrunc, MakeConst(-4.6), MakeConst(0))));
  }
  EXPECT_THROW(Eval(MakeBinary(OPtrunc, x, MakeConst(1))), UnsupportedExprError);
  EXPECT_THROW(Eval(MakeBinary(OPtrunc, x, y)), UnsupportedExprError);
}

TEST_P(SolverTest, Count) {
  LogicalExpr args[] = {
    MakeRelational(NE, x, MakeConst(0)),
    MakeRelational(NE, y, MakeConst(0)),
    MakeRelational(NE, z, MakeConst(0))
  };
  EXPECT_EQ(0, Eval(MakeCount(args)));
  EXPECT_EQ(1, Eval(MakeCount(args), 1));
  EXPECT_EQ(2, Eval(MakeCount(args), 0, 1, 1));
  EXPECT_EQ(3, Eval(MakeCount(args), 1, 1, 1));
}

TEST_P(SolverTest, NumberOf) {
  NumericExpr args[] = {MakeConst(42), x};
  EXPECT_EQ(0, Eval(MakeNumberOf(args)));
  EXPECT_EQ(1, Eval(MakeNumberOf(args), 42));
  NumericExpr args3[] = {MakeConst(42), x, y};
  EXPECT_EQ(0, Eval(MakeNumberOf(args3)));
  EXPECT_EQ(1, Eval(MakeNumberOf(args3), 0, 42));
  EXPECT_EQ(2, Eval(MakeNumberOf(args3), 42, 42));
}

TEST_P(SolverTest, PiecewiseLinear) {
  // Test on the following piecewise-linear function:
  //
  //     y^
  //    \ |           /
  //     \|  3  6  9 /      x
  //  ----\-->-->-->/------->
  //     0|\       /
  //      | \     /
  //    -3|  \___/
  //      |
  //
  // Breakpoints are at x = 3 and x = 6. Slopes are -1, 0 and 1.
  //
  double breakpoints[] = {3, 6};
  double slopes[] = {-1, 0, 1};
  NumericExpr e = MakePiecewiseLinear(2, breakpoints, slopes, MakeVariable(1));
  if (!HasFeature(feature::PLTERM)) {
    EXPECT_THROW(Eval(e, 42);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(33, Eval(e, 42));
  EXPECT_EQ(-3, Eval(e, 4));
  EXPECT_EQ(1, Eval(e, -1));
}

TEST_P(SolverTest, UnsupportedFunctionCall) {
  mp::Function f = AddFunction("foo", TestFunc, 2);
  mp::Expr args[] = {MakeConst(1), MakeConst(2)};
  EXPECT_THROW(Eval(MakeCall(f, args), 3);, UnsupportedExprError);
}

TEST_P(SolverTest, PowConstExp) {
  EXPECT_EQ(16, Eval(MakeBinary(OP1POW, x, MakeConst(4)), 2));
}

TEST_P(SolverTest, Pow2) {
  EXPECT_EQ(49, Eval(MakeUnary(OP2POW, x), 7));
}

TEST_P(SolverTest, PowConstBase) {
  NumericExpr e = MakeBinary(OPCPOW, MakeConst(5), x);
  if (!HasFeature(feature::POW)) {
    EXPECT_THROW(Eval(e, 3);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(125, Eval(e, 3));
}

TEST_P(SolverTest, NumericConstant) {
  EXPECT_EQ(42, Eval(MakeConst(42)));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(42, Eval(MakeBinary(OPMULT, MakeConst(0.42), MakeConst(100))));
    return;
  }
  EXPECT_THROW_MSG(Eval(MakeConst(0.42));, UnsupportedExprError,
    "value 0.42 can't be represented as int");
  EXPECT_THROW(Solve("objconst"), mp::Error);
}

TEST_P(SolverTest, Var) {
  EXPECT_EQ(11, Eval(x, 11, 22));
  EXPECT_EQ(22, Eval(y, 11, 22));
  EXPECT_EQ(33, Eval(x, 33));
}

TEST_P(SolverTest, Or) {
  NumericExpr one = MakeConst(1);
  LogicalExpr e = MakeBinaryLogical(
      OPOR, MakeRelational(EQ, x, one), MakeRelational(EQ, y, one));
  EXPECT_EQ(0, Eval(e, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_P(SolverTest, And) {
  NumericExpr one = MakeConst(1);
  LogicalExpr e = MakeBinaryLogical(
      OPAND, MakeRelational(EQ, x, one), MakeRelational(EQ, y, one));
  EXPECT_EQ(0, Eval(e, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_P(SolverTest, Less) {
  LogicalExpr e = MakeRelational(LT, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_P(SolverTest, LessEqual) {
  LogicalExpr e = MakeRelational(LE, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_P(SolverTest, Equal) {
  LogicalExpr e = MakeRelational(EQ, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_P(SolverTest, GreaterEqual) {
  LogicalExpr e = MakeRelational(GE, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_P(SolverTest, Greater) {
  LogicalExpr e = MakeRelational(GT, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_P(SolverTest, NotEqual) {
  LogicalExpr e = MakeRelational(NE, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_P(SolverTest, Not) {
  LogicalExpr e = MakeNot(MakeRelational(EQ, x, MakeConst(1)));
  EXPECT_EQ(1, Eval(e, 0));
  EXPECT_EQ(0, Eval(e, 1));
}

TEST_P(SolverTest, AtLeast) {
  LogicalExpr args[] = {
    MakeRelational(NE, y, MakeConst(0)),
    MakeRelational(NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(OPATLEAST, x, MakeCount(args));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, AtMost) {
  LogicalExpr args[] = {
    MakeRelational(NE, y, MakeConst(0)),
    MakeRelational(NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(OPATMOST, x, MakeCount(args));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, Exactly) {
  LogicalExpr args[] = {
    MakeRelational(NE, y, MakeConst(0)),
    MakeRelational(NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(OPEXACTLY, x, MakeCount(args));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, NotAtLeast) {
  LogicalExpr args[] = {
    MakeRelational(NE, y, MakeConst(0)),
    MakeRelational(NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(OPNOTATLEAST, x, MakeCount(args));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, NotAtMost) {
  LogicalExpr args[] = {
    MakeRelational(NE, y, MakeConst(0)),
    MakeRelational(NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(OPNOTATMOST, x, MakeCount(args));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, NotExactly) {
  LogicalExpr args[] = {
    MakeRelational(NE, y, MakeConst(0)),
    MakeRelational(NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(OPNOTEXACTLY, x, MakeCount(args));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, ForAll) {
  LogicalExpr args[] = {
    MakeRelational(EQ, x, MakeConst(1)),
    MakeRelational(EQ, y, MakeConst(1)),
    MakeRelational(EQ, z, MakeConst(1))
  };
  LogicalExpr e = MakeIteratedLogical(ANDLIST, args);
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 0, 1));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 1));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
}

TEST_P(SolverTest, Exists) {
  LogicalExpr args[] = {
    MakeRelational(EQ, x, MakeConst(1)),
    MakeRelational(EQ, y, MakeConst(1)),
    MakeRelational(EQ, z, MakeConst(1))
  };
  LogicalExpr e = MakeIteratedLogical(ORLIST, args);
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 0, 1));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 1));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
}

TEST_P(SolverTest, Implication) {
  LogicalExpr e = MakeImplication(
      MakeRelational(EQ, x, MakeConst(1)),
      MakeRelational(EQ, y, MakeConst(1)),
      MakeRelational(EQ, z, MakeConst(1)));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 0, 1));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 1));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
}

TEST_P(SolverTest, Iff) {
  LogicalExpr e = MakeBinaryLogical(OP_IFF,
      MakeRelational(EQ, x, MakeConst(1)),
      MakeRelational(EQ, y, MakeConst(1)));
  EXPECT_EQ(1, Eval(e, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_P(SolverTest, AllDiff) {
  NumericExpr args[] = {MakeConst(1), x, y};
  LogicalExpr e = MakeAllDiff(args);
  EXPECT_TRUE(Solve(e, 2, 3).has_value());
  EXPECT_FALSE(Solve(e, 2, 1).has_value());
  EXPECT_FALSE(Solve(e, 1, 1).has_value());
}

TEST_P(SolverTest, NestedAllDiff) {
  NumericExpr args[] = {MakeConst(1), x, y};
  EXPECT_THROW(Eval(MakeNot(MakeAllDiff(args)), 1, 2), UnsupportedExprError);
}

TEST_P(SolverTest, LogicalConstant) {
  EXPECT_EQ(0, Eval(MakeLogicalConstant(false)));
  EXPECT_EQ(1, Eval(MakeLogicalConstant(true)));
}

TEST_P(SolverTest, NonlinearObj) {
  Problem p;
  p.AddVar(2, 2, var::INTEGER);
  mp::Variable x = MakeVariable(0);
  p.AddObj(obj::MIN, MakeBinary(OPMULT, x, x));
  EXPECT_EQ(4, Solve(p).obj_value());
}

TEST_P(SolverTest, ObjConst) {
  Problem p;
  p.AddVar(0, 0, var::INTEGER);
  p.AddObj(obj::MIN, MakeConst(42));
  EXPECT_EQ(42, Solve(p).obj_value());
}

TEST_P(SolverTest, Minimize) {
  Problem p;
  p.AddVar(42, 100, var::INTEGER);
  p.AddObj(obj::MIN, MakeVariable(0));
  EXPECT_EQ(42, Solve(p).obj_value());
}

TEST_P(SolverTest, Maximize) {
  Problem p;
  p.AddVar(0, 42, var::INTEGER);
  p.AddObj(obj::MAX, MakeVariable(0));
  EXPECT_EQ(42, Solve(p).obj_value());
}

TEST_P(SolverTest, TimingOption) {
  struct TestOutputHandler : mp::OutputHandler {
    std::string output;

    virtual ~TestOutputHandler() {}
    void HandleOutput(fmt::StringRef output) {
      this->output += output;
    }
  };
  TestOutputHandler oh;
  solver_->set_output_handler(&oh);

  Problem p;
  p.AddVar(42, 100, var::INTEGER);
  p.AddObj(obj::MIN, MakeVariable(0));

  solver_->SetIntOption("timing", 0);
  solver_->Solve(p);
  EXPECT_TRUE(oh.output.find("Setup time = ") == std::string::npos);
  EXPECT_TRUE(oh.output.find("Solution time = ") == std::string::npos);
  EXPECT_TRUE(oh.output.find("Output time = ") == std::string::npos);

  solver_->SetIntOption("timing", 1);
  solver_->Solve(p);
  EXPECT_TRUE(oh.output.find("Setup time = ") != std::string::npos);
  EXPECT_TRUE(oh.output.find("Solution time = ") != std::string::npos);
  EXPECT_TRUE(oh.output.find("Output time = ") != std::string::npos);
}

// ---------------------------------------------------------------------------
// Solve test problems

TEST_P(SolverTest, SolveAssign0) {
  EXPECT_EQ(6, Solve("assign0").obj);
}

TEST_P(SolverTest, SolveAssign1) {
  EXPECT_EQ(6, Solve("assign1").obj);
}

TEST_P(SolverTest, SolveFlowshp1) {
  EXPECT_EQ(22, Solve("flowshp1").obj);
}

TEST_P(SolverTest, SolveGrpassign0) {
  EXPECT_EQ(61, Solve("grpassign0").obj);
}

TEST_P(SolverTest, SolveMagic) {
  EXPECT_TRUE(Solve("magic").solved);
}

TEST_P(SolverTest, SolveMapcoloring) {
  EXPECT_TRUE(Solve("mapcoloring").solved);
}

TEST_P(SolverTest, SolveNQueens) {
  EXPECT_TRUE(Solve("nqueens").solved);
}

TEST_P(SolverTest, SolveNQueens0) {
  EXPECT_TRUE(Solve("nqueens0").solved);
}

TEST_P(SolverTest, SolveSched0) {
  EXPECT_EQ(5, Solve("sched0").obj);
}

TEST_P(SolverTest, SolveSched1) {
  EXPECT_EQ(5, Solve("sched1").obj);
}

TEST_P(SolverTest, SolveSched2) {
  EXPECT_EQ(5, Solve("sched2").obj);
}

TEST_P(SolverTest, SolveSendMoreMoney) {
  EXPECT_TRUE(Solve("send-more-money").solved);
}

TEST_P(SolverTest, SolveSendMostMoney) {
  EXPECT_NEAR(10876, Solve("send-most-money").obj, 1);
}

TEST_P(SolverTest, SolveSeq0) {
  EXPECT_NEAR(332, Solve("seq0").obj, 1e-5);
}

TEST_P(SolverTest, SolveSeq0a) {
  EXPECT_NEAR(332, Solve("seq0a").obj, 1e-5);
}

TEST_P(SolverTest, SolveSudokuHard) {
  EXPECT_TRUE(Solve("sudokuHard").solved);
}

TEST_P(SolverTest, SolveSudokuVeryEasy) {
  EXPECT_TRUE(Solve("sudokuVeryEasy").solved);
}

// ----------------------------------------------------------------------------
// Solve code tests

TEST_P(SolverTest, OptimalSolveCode) {
  Problem p;
  EXPECT_TRUE(Solve(p, "objconstint").solved);
  EXPECT_EQ(0, p.solve_code());
}

TEST_P(SolverTest, FeasibleSolveCode) {
  Problem p;
  EXPECT_TRUE(Solve(p, "feasible").solved);
  EXPECT_EQ(0, p.solve_code());
}

TEST_P(SolverTest, InfeasibleSolveCode) {
  Problem p;
  EXPECT_FALSE(Solve(p, "infeasible").solved);
  EXPECT_EQ(200, p.solve_code());
}

// ----------------------------------------------------------------------------
// Interrupt tests

#ifdef HAVE_THREADS
void Interrupt() {
  // Wait until started.
  while (mp::SignalHandler::stop())
    std::this_thread::yield();
  std::raise(SIGINT);
}

TEST_P(SolverTest, InterruptSolution) {
  std::thread t(Interrupt);
  Problem p;
  string message = Solve(p, "miplib/assign1").message;
  t.join();
  EXPECT_EQ(600, p.solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}
#endif

struct SolutionCounter : mp::BasicSolutionHandler {
  int num_solutions;
  SolutionCounter() : num_solutions(0) {}
  void HandleFeasibleSolution(Problem &, fmt::StringRef,
        const double *, const double *, double) {
    ++num_solutions;
  }
};

TEST_P(SolverTest, CountSolutions) {
  mp::Problem p;
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeVariable(2)};
  p.AddCon(MakeAllDiff(args));
  SolutionCounter sc;
  solver_->SetIntOption("solutionlimit", 10);
  solver_->SetIntOption("countsolutions", 1);
  solver_->set_solution_handler(&sc);
  solver_->Solve(p);
  EXPECT_EQ(6, sc.num_solutions);
}

TEST_P(SolverTest, SatisfactionSolutionLimit) {
  mp::Problem p;
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeVariable(2)};
  p.AddCon(MakeAllDiff(args));
  solver_->SetIntOption("solutionlimit", 5);
  solver_->Solve(p);
  EXPECT_EQ(0, p.solve_code());
}

TEST_P(SolverTest, OptimizationSolutionLimit) {
  mp::Problem p;
  p.Read(MP_TEST_DATA_DIR "/photo9");
  solver_->SetIntOption("solutionlimit", 2);
  solver_->Solve(p);
  EXPECT_GE(p.solve_code(), 400);
  EXPECT_LT(p.solve_code(), 499);
}

TEST_P(SolverTest, MultipleSolutions) {
  mp::Problem p;
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeVariable(2)};
  p.AddCon(MakeAllDiff(args));
  SolutionCounter sc;
  solver_->set_solution_handler(&sc);
  solver_->SetIntOption("solutionlimit", 3);
  solver_->SetStrOption("solutionstub", "test");
  solver_->Solve(p);
  EXPECT_EQ(0, p.solve_code());
  EXPECT_EQ(3, sc.num_solutions);
}

TEST_P(SolverTest, OptionValues) {
  for (mp::Solver::option_iterator
    i = solver_->option_begin(), e = solver_->option_end(); i != e; ++i) {
    for (mp::ValueArrayRef::iterator j = i->values().begin(),
        value_end = i->values().end(); j != value_end; ++j) {
      EXPECT_TRUE(j->value != 0);
    }
  }
}

TEST_P(SolverTest, CreateSolver) {
  EXPECT_STREQ(solver_->name(), mp::CreateSolver(0)->name());
}

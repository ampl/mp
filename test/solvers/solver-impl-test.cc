/*
 Solver implementation test

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

#include "solver-impl-test.h"
#include "../gtest-extra.h"

#include "mp/nl.h"

#include <cmath>

using std::string;

using mp::LogicalExpr;
using mp::NumericExpr;
using mp::Problem;
using mp::UnsupportedExprError;
namespace var = mp::var;
namespace obj = mp::obj;

SolverImplTest::EvalResult SolverImplTest::Solve(Problem &p) {
  struct TestSolutionHandler : mp::BasicSolutionHandler {
    EvalResult result;
    virtual ~TestSolutionHandler() {}
    void HandleSolution(int status, fmt::StringRef,
          const double *values, const double *, double obj_value) {
      result = values ? EvalResult(values[0], obj_value) : EvalResult();
      result.set_solve_code(status);
    }
  };
  TestSolutionHandler sh;
  solver_->Solve(p, sh);
  return sh.result;
}

SolverImplTest::EvalResult SolverImplTest::Solve(
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

SolverImplTest::SolverImplTest()
: solver_(GetParam().create_solver()), features_(GetParam().features) {
  mp::NLHeader header = mp::NLHeader();
  header.num_vars = 4;
  header.num_objs = 1;
  header.num_funcs = 2;
  set_flags(mp::internal::ASL_STANDARD_OPCODES);
  SetInfo(header);
  x = MakeVariable(1);
  y = MakeVariable(2);
  z = MakeVariable(3);
}

SolveResult SolverImplTest::Solve(
    mp::ASLSolver &s, Problem &p, const char *stub, const char *opt) {
  TestSolutionHandler sh;
  const std::string DATA_DIR = MP_TEST_DATA_DIR "/";
  p.Read(DATA_DIR + stub);
  char *args[] = {const_cast<char*>(opt), 0};
  if (s.ParseOptions(args, 0, &p))
    s.Solve(p, sh);
  const string &message = sh.message();
  int solve_result = sh.status();
  EXPECT_GE(solve_result, 0);
  bool solved = true;
  if (solve_result < 100) {
    EXPECT_TRUE(message.find("optimal solution") != string::npos ||
        (p.num_objs() == 0 &&
            message.find("feasible solution") != string::npos));
  } else
    solved = false;
  return SolveResult(solved, sh.obj_value(), message);
}

TEST_P(SolverImplTest, Add) {
  NumericExpr e = MakeBinary(mp::expr::ADD, x, y);
  EXPECT_EQ(25, Eval(e, 10, 15));
  EXPECT_EQ(12, Eval(e, 19, -7));
}

TEST_P(SolverImplTest, Sub) {
  NumericExpr e = MakeBinary(mp::expr::SUB, x, y);
  EXPECT_EQ(-5, Eval(e, 10, 15));
  EXPECT_EQ(26, Eval(e, 19, -7));
}

TEST_P(SolverImplTest, Mul) {
  NumericExpr e = MakeBinary(mp::expr::MUL, x, y);
  EXPECT_EQ(150, Eval(e, 10, 15));
  EXPECT_EQ(-133, Eval(e, 19, -7));
}

TEST_P(SolverImplTest, Div) {
  NumericExpr e = MakeBinary(mp::expr::DIV, x, y);
  if (!HasFeature(feature::DIV)) {
    EXPECT_THROW(Eval(e, 150, 15);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(10, Eval(e, 150, 15));
  EXPECT_EQ(-7, Eval(e, -133, 19));
}

TEST_P(SolverImplTest, Mod) {
  NumericExpr e = MakeBinary(mp::expr::MOD, x, y);
  EXPECT_EQ(0, Eval(e, 9, 3));
  EXPECT_EQ(2, Eval(e, 8, 3));
  EXPECT_EQ(-2, Eval(e, -8, 3));
  EXPECT_EQ(2, Eval(e, 8, -3));
  EXPECT_EQ(-2, Eval(e, -8, -3));
}

TEST_P(SolverImplTest, Pow) {
  NumericExpr e = MakeBinary(mp::expr::POW, x, y);
  if (!HasFeature(feature::POW)) {
    EXPECT_THROW(Eval(e, 2, 3);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(8, Eval(e, 2, 3));
  EXPECT_EQ(81, Eval(e, 3, 4));
}

TEST_P(SolverImplTest, Less) {
  NumericExpr e = MakeBinary(mp::expr::LESS, x, y);
  EXPECT_EQ(0, Eval(e, 10, 15));
  EXPECT_EQ(26, Eval(e, 19, -7));
}

TEST_P(SolverImplTest, Min) {
  NumericExpr args[] = {x, y, z};
  NumericExpr e = MakeVarArg(mp::expr::MIN, args);
  EXPECT_EQ(-7, Eval(e, 3, -7, 5));
  EXPECT_EQ(10, Eval(e, 10, 20, 30));
}

TEST_P(SolverImplTest, Max) {
  NumericExpr args[] = {x, y, z};
  NumericExpr e = MakeVarArg(mp::expr::MAX, args);
  EXPECT_EQ(5, Eval(e, 3, -7, 5));
  EXPECT_EQ(30, Eval(e, 30, 20, 10));
}

TEST_P(SolverImplTest, Floor) {
  NumericExpr e = MakeUnary(mp::expr::FLOOR, x);
  EXPECT_EQ(-42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
  if (!HasFeature(feature::FLOAT_CONST)) return;
  EXPECT_EQ(4, Eval(MakeUnary(mp::expr::FLOOR, MakeConst(4.9))));
  EXPECT_EQ(-5, Eval(MakeUnary(mp::expr::FLOOR, MakeConst(-4.1))));
}

TEST_P(SolverImplTest, Ceil) {
  NumericExpr e = MakeUnary(mp::expr::CEIL, x);
  EXPECT_EQ(-42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
  if (!HasFeature(feature::FLOAT_CONST)) return;
  EXPECT_EQ(5, Eval(MakeUnary(mp::expr::CEIL, MakeConst(4.1))));
  EXPECT_EQ(-4, Eval(MakeUnary(mp::expr::CEIL, MakeConst(-4.9))));
}

TEST_P(SolverImplTest, Abs) {
  NumericExpr e = MakeUnary(mp::expr::ABS, x);
  EXPECT_EQ(42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
}

TEST_P(SolverImplTest, Minus) {
  NumericExpr e = MakeUnary(mp::expr::MINUS, x);
  EXPECT_EQ(42, Eval(e, -42));
  EXPECT_EQ(-42, Eval(e, 42));
}

TEST_P(SolverImplTest, If) {
  NumericExpr e = MakeIf(MakeRelational(mp::expr::EQ, x, MakeConst(1)), y, z);
  EXPECT_EQ(42, Eval(e, 1, 42, 10));
  EXPECT_EQ(10, Eval(e, 0, 42, 10));
  EXPECT_EQ(42, Eval(e, 1, 42, 42));
}

TEST_P(SolverImplTest, Tanh) {
  double arg = (std::log(1.5) - std::log(0.5)) / 2;
  NumericExpr e = MakeBinary(mp::expr::MUL,
      MakeConst(2), MakeUnary(mp::expr::TANH, MakeConst(arg)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(1, Eval(e));
}

TEST_P(SolverImplTest, Tan) {
  EXPECT_THROW(Eval(MakeUnary(mp::expr::TAN, x)), UnsupportedExprError);
}

TEST_P(SolverImplTest, Sqrt) {
  NumericExpr e = MakeUnary(mp::expr::SQRT, x);
  if (!HasFeature(feature::SQRT)) {
    EXPECT_THROW(Eval(e, 64);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(8, Eval(e, 64));
}

TEST_P(SolverImplTest, Sinh) {
  NumericExpr e = MakeUnary(
        mp::expr::SINH, MakeConst(std::log(2 + std::sqrt(5.0))));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(2, Eval(e));
}

TEST_P(SolverImplTest, Sin) {
  EXPECT_THROW(Eval(MakeUnary(mp::expr::SIN, x)), UnsupportedExprError);
}

TEST_P(SolverImplTest, Log10) {
  NumericExpr e = MakeUnary(mp::expr::LOG10, MakeConst(1000));
  if (!HasFeature(feature::LOG)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(3, Eval(e));
}

TEST_P(SolverImplTest, Log) {
  NumericExpr e = MakeUnary(mp::expr::LOG, MakeConst(std::exp(5.0)));
  if (!HasFeature(feature::LOG)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverImplTest, Exp) {
  NumericExpr e = MakeUnary(mp::expr::EXP, MakeConst(std::log(5.0)));
  if (!HasFeature(feature::EXP)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverImplTest, Cosh) {
  double x = 5;
  NumericExpr e = MakeUnary(mp::expr::COSH,
      MakeConst(std::log(x + std::sqrt(x + 1) * std::sqrt(x - 1))));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverImplTest, Cos) {
  EXPECT_THROW(Eval(MakeUnary(mp::expr::COS, x)), UnsupportedExprError);
}

TEST_P(SolverImplTest, Atanh) {
  mp::UnaryExpr x = MakeUnary(mp::expr::ATANH, MakeConst(std::tanh(5.0)));
  NumericExpr e = MakeUnary(mp::expr::FLOOR, MakeBinary(mp::expr::ADD,
      MakeConst(0.5), MakeBinary(mp::expr::MUL, MakeConst(1000000), x)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5000000, Eval(e));
}

TEST_P(SolverImplTest, Atan2) {
  EXPECT_THROW(Eval(MakeBinary(mp::expr::ATAN2, x, y)), UnsupportedExprError);
}

TEST_P(SolverImplTest, Atan) {
  EXPECT_THROW(Eval(MakeUnary(mp::expr::ATAN, x)), UnsupportedExprError);
}

TEST_P(SolverImplTest, Asinh) {
  NumericExpr e = MakeUnary(mp::expr::ASINH, MakeConst(std::sinh(5.0)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverImplTest, Asin) {
  EXPECT_THROW(Eval(MakeUnary(mp::expr::ASIN, x)), UnsupportedExprError);
}

TEST_P(SolverImplTest, Acosh) {
  NumericExpr e = MakeUnary(mp::expr::ACOSH, MakeConst(std::cosh(5.0)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverImplTest, Acos) {
  EXPECT_THROW(Eval(MakeUnary(mp::expr::ACOS, x)), UnsupportedExprError);
}

TEST_P(SolverImplTest, Sum) {
  NumericExpr args[] = {x, y, z};
  using mp::MakeArrayRef;
  EXPECT_EQ(0, Eval(MakeSum(MakeArrayRef(args, 0))));
  EXPECT_EQ(42, Eval(MakeSum(MakeArrayRef(args, 1)), 42));
  EXPECT_EQ(123, Eval(MakeSum(args), 100, 20, 3));
}

TEST_P(SolverImplTest, IntDiv) {
  NumericExpr e = MakeBinary(mp::expr::INT_DIV, x, y);
  EXPECT_EQ(3, Eval(e, 9, 3));
  EXPECT_EQ(2, Eval(e, 8, 3));
  EXPECT_EQ(-2, Eval(e, -8, 3));
  EXPECT_EQ(-2, Eval(e, 8, -3));
  EXPECT_EQ(2, Eval(e, -8, -3));
}

TEST_P(SolverImplTest, Precision) {
  EXPECT_THROW(Eval(MakeBinary(mp::expr::PRECISION, x, y)),
               UnsupportedExprError);
}

TEST_P(SolverImplTest, Round) {
  mp::expr::Kind round = mp::expr::ROUND;
  EXPECT_EQ(42, Eval(MakeBinary(round, x, MakeConst(0)), 42));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(4, Eval(MakeBinary(round, MakeConst(4.4), MakeConst(0))));
    EXPECT_EQ(5, Eval(MakeBinary(round, MakeConst(4.6), MakeConst(0))));
    EXPECT_EQ(-4, Eval(MakeBinary(round, MakeConst(-4.4), MakeConst(0))));
    EXPECT_EQ(-5, Eval(MakeBinary(round, MakeConst(-4.6), MakeConst(0))));
  }
  EXPECT_THROW(Eval(MakeBinary(round, x, MakeConst(1))), UnsupportedExprError);
  EXPECT_THROW(Eval(MakeBinary(round, x, y)), UnsupportedExprError);
}

TEST_P(SolverImplTest, Trunc) {
  mp::expr::Kind trunc = mp::expr::TRUNC;
  EXPECT_EQ(42, Eval(MakeBinary(trunc, x, MakeConst(0)), 42));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(4, Eval(MakeBinary(trunc, MakeConst(4.4), MakeConst(0))));
    EXPECT_EQ(4, Eval(MakeBinary(trunc, MakeConst(4.6), MakeConst(0))));
    EXPECT_EQ(-4, Eval(MakeBinary(trunc, MakeConst(-4.4), MakeConst(0))));
    EXPECT_EQ(-4, Eval(MakeBinary(trunc, MakeConst(-4.6), MakeConst(0))));
  }
  EXPECT_THROW(Eval(MakeBinary(trunc, x, MakeConst(1))), UnsupportedExprError);
  EXPECT_THROW(Eval(MakeBinary(trunc, x, y)), UnsupportedExprError);
}

TEST_P(SolverImplTest, Count) {
  LogicalExpr args[] = {
    MakeRelational(mp::expr::NE, x, MakeConst(0)),
    MakeRelational(mp::expr::NE, y, MakeConst(0)),
    MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  EXPECT_EQ(0, Eval(MakeCount(args)));
  EXPECT_EQ(1, Eval(MakeCount(args), 1));
  EXPECT_EQ(2, Eval(MakeCount(args), 0, 1, 1));
  EXPECT_EQ(3, Eval(MakeCount(args), 1, 1, 1));
}

TEST_P(SolverImplTest, NumberOf) {
  NumericExpr args[] = {MakeConst(42), x};
  EXPECT_EQ(0, Eval(MakeNumberOf(args)));
  EXPECT_EQ(1, Eval(MakeNumberOf(args), 42));
  NumericExpr args3[] = {MakeConst(42), x, y};
  EXPECT_EQ(0, Eval(MakeNumberOf(args3)));
  EXPECT_EQ(1, Eval(MakeNumberOf(args3), 0, 42));
  EXPECT_EQ(2, Eval(MakeNumberOf(args3), 42, 42));
}

TEST_P(SolverImplTest, PiecewiseLinear) {
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

TEST_P(SolverImplTest, UnsupportedFunctionCall) {
  mp::Function f = AddFunction("foo", TestFunc, 2);
  mp::Expr args[] = {MakeConst(1), MakeConst(2)};
  EXPECT_THROW(Eval(MakeCall(f, args), 3);, UnsupportedExprError);
}

TEST_P(SolverImplTest, PowConstExp) {
  EXPECT_EQ(16, Eval(MakeBinary(mp::expr::POW_CONST_EXP, x, MakeConst(4)), 2));
}

TEST_P(SolverImplTest, Pow2) {
  EXPECT_EQ(49, Eval(MakeUnary(mp::expr::POW2, x), 7));
}

TEST_P(SolverImplTest, PowConstBase) {
  NumericExpr e = MakeBinary(mp::expr::POW_CONST_BASE, MakeConst(5), x);
  if (!HasFeature(feature::POW)) {
    EXPECT_THROW(Eval(e, 3);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(125, Eval(e, 3));
}

TEST_P(SolverImplTest, NumericConstant) {
  EXPECT_EQ(42, Eval(MakeConst(42)));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(42, Eval(MakeBinary(mp::expr::MUL, MakeConst(0.42), MakeConst(100))));
    return;
  }
  EXPECT_THROW_MSG(Eval(MakeConst(0.42));, UnsupportedExprError,
    "value 0.42 can't be represented as int");
  EXPECT_THROW(Solve("objconst"), mp::Error);
}

TEST_P(SolverImplTest, Var) {
  EXPECT_EQ(11, Eval(x, 11, 22));
  EXPECT_EQ(22, Eval(y, 11, 22));
  EXPECT_EQ(33, Eval(x, 33));
}

TEST_P(SolverImplTest, Or) {
  NumericExpr one = MakeConst(1);
  LogicalExpr e = MakeBinaryLogical(
      mp::expr::OR, MakeRelational(mp::expr::EQ, x, one),
        MakeRelational(mp::expr::EQ, y, one));
  EXPECT_EQ(0, Eval(e, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_P(SolverImplTest, And) {
  NumericExpr one = MakeConst(1);
  LogicalExpr e = MakeBinaryLogical(
      mp::expr::AND, MakeRelational(mp::expr::EQ, x, one),
        MakeRelational(mp::expr::EQ, y, one));
  EXPECT_EQ(0, Eval(e, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_P(SolverImplTest, LT) {
  LogicalExpr e = MakeRelational(mp::expr::LT, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_P(SolverImplTest, LE) {
  LogicalExpr e = MakeRelational(mp::expr::LE, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_P(SolverImplTest, EQ) {
  LogicalExpr e = MakeRelational(mp::expr::EQ, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_P(SolverImplTest, GE) {
  LogicalExpr e = MakeRelational(mp::expr::GE, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_P(SolverImplTest, GT) {
  LogicalExpr e = MakeRelational(mp::expr::GT, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_P(SolverImplTest, NE) {
  LogicalExpr e = MakeRelational(mp::expr::NE, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_P(SolverImplTest, Not) {
  LogicalExpr e = MakeNot(MakeRelational(mp::expr::EQ, x, MakeConst(1)));
  EXPECT_EQ(1, Eval(e, 0));
  EXPECT_EQ(0, Eval(e, 1));
}

TEST_P(SolverImplTest, AtLeast) {
  LogicalExpr args[] = {
    MakeRelational(mp::expr::NE, y, MakeConst(0)),
    MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(mp::expr::ATLEAST, x, MakeCount(args));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_P(SolverImplTest, AtMost) {
  LogicalExpr args[] = {
    MakeRelational(mp::expr::NE, y, MakeConst(0)),
    MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(mp::expr::ATMOST, x, MakeCount(args));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_P(SolverImplTest, Exactly) {
  LogicalExpr args[] = {
    MakeRelational(mp::expr::NE, y, MakeConst(0)),
    MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(mp::expr::EXACTLY, x, MakeCount(args));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_P(SolverImplTest, NotAtLeast) {
  LogicalExpr args[] = {
    MakeRelational(mp::expr::NE, y, MakeConst(0)),
    MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(mp::expr::NOT_ATLEAST, x, MakeCount(args));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_P(SolverImplTest, NotAtMost) {
  LogicalExpr args[] = {
    MakeRelational(mp::expr::NE, y, MakeConst(0)),
    MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(mp::expr::NOT_ATMOST, x, MakeCount(args));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_P(SolverImplTest, NotExactly) {
  LogicalExpr args[] = {
    MakeRelational(mp::expr::NE, y, MakeConst(0)),
    MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = MakeLogicalCount(mp::expr::NOT_EXACTLY, x, MakeCount(args));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_P(SolverImplTest, ForAll) {
  LogicalExpr args[] = {
    MakeRelational(mp::expr::EQ, x, MakeConst(1)),
    MakeRelational(mp::expr::EQ, y, MakeConst(1)),
    MakeRelational(mp::expr::EQ, z, MakeConst(1))
  };
  LogicalExpr e = MakeIteratedLogical(mp::expr::FORALL, args);
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 0, 1));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 1));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
}

TEST_P(SolverImplTest, Exists) {
  LogicalExpr args[] = {
    MakeRelational(mp::expr::EQ, x, MakeConst(1)),
    MakeRelational(mp::expr::EQ, y, MakeConst(1)),
    MakeRelational(mp::expr::EQ, z, MakeConst(1))
  };
  LogicalExpr e = MakeIteratedLogical(mp::expr::EXISTS, args);
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 0, 1));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 1));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
}

TEST_P(SolverImplTest, Implication) {
  LogicalExpr e = MakeImplication(
      MakeRelational(mp::expr::EQ, x, MakeConst(1)),
      MakeRelational(mp::expr::EQ, y, MakeConst(1)),
      MakeRelational(mp::expr::EQ, z, MakeConst(1)));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 0, 1));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 1));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
}

TEST_P(SolverImplTest, Iff) {
  LogicalExpr e = MakeBinaryLogical(mp::expr::IFF,
      MakeRelational(mp::expr::EQ, x, MakeConst(1)),
      MakeRelational(mp::expr::EQ, y, MakeConst(1)));
  EXPECT_EQ(1, Eval(e, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_P(SolverImplTest, AllDiff) {
  NumericExpr args[] = {MakeConst(1), x, y};
  LogicalExpr e = MakeAllDiff(args);
  EXPECT_TRUE(Solve(e, 2, 3).has_value());
  EXPECT_FALSE(Solve(e, 2, 1).has_value());
  EXPECT_FALSE(Solve(e, 1, 1).has_value());
}

TEST_P(SolverImplTest, NestedAllDiff) {
  NumericExpr args[] = {MakeConst(1), x, y};
  EXPECT_THROW(Eval(MakeNot(MakeAllDiff(args)), 1, 2), UnsupportedExprError);
}

TEST_P(SolverImplTest, LogicalConstant) {
  EXPECT_EQ(0, Eval(MakeLogicalConstant(false)));
  EXPECT_EQ(1, Eval(MakeLogicalConstant(true)));
}

TEST_P(SolverImplTest, NonlinearObj) {
  Problem p;
  p.AddVar(2, 2, var::INTEGER);
  mp::Variable x = MakeVariable(0);
  p.AddObj(obj::MIN, MakeBinary(mp::expr::MUL, x, x));
  EXPECT_EQ(4, Solve(p).obj_value());
}

TEST_P(SolverImplTest, ObjConst) {
  Problem p;
  p.AddVar(0, 0, var::INTEGER);
  p.AddObj(obj::MIN, MakeConst(42));
  EXPECT_EQ(42, Solve(p).obj_value());
}

TEST_P(SolverImplTest, Minimize) {
  Problem p;
  p.AddVar(42, 100, var::INTEGER);
  p.AddObj(obj::MIN, MakeVariable(0));
  EXPECT_EQ(42, Solve(p).obj_value());
}

TEST_P(SolverImplTest, Maximize) {
  Problem p;
  p.AddVar(0, 42, var::INTEGER);
  p.AddObj(obj::MAX, MakeVariable(0));
  EXPECT_EQ(42, Solve(p).obj_value());
}

TEST_P(SolverImplTest, TimingOption) {
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
  mp::BasicSolutionHandler sol_handler;
  solver_->Solve(p, sol_handler);
  EXPECT_TRUE(oh.output.find("Setup time = ") == std::string::npos);
  EXPECT_TRUE(oh.output.find("Solution time = ") == std::string::npos);
  EXPECT_TRUE(oh.output.find("Output time = ") == std::string::npos);

  solver_->SetIntOption("timing", 1);
  solver_->Solve(p, sol_handler);
  EXPECT_TRUE(oh.output.find("Setup time = ") != std::string::npos);
  EXPECT_TRUE(oh.output.find("Solution time = ") != std::string::npos);
  EXPECT_TRUE(oh.output.find("Output time = ") != std::string::npos);
}

// ---------------------------------------------------------------------------
// Solve test problems

TEST_P(SolverImplTest, SolveAssign0) {
  EXPECT_EQ(6, Solve("assign0").obj);
}

TEST_P(SolverImplTest, SolveAssign1) {
  EXPECT_EQ(6, Solve("assign1").obj);
}

TEST_P(SolverImplTest, SolveFlowshp1) {
  EXPECT_EQ(22, Solve("flowshp1").obj);
}

TEST_P(SolverImplTest, SolveGrpassign0) {
  EXPECT_EQ(61, Solve("grpassign0").obj);
}

TEST_P(SolverImplTest, SolveMagic) {
  EXPECT_TRUE(Solve("magic").solved);
}

TEST_P(SolverImplTest, SolveMapcoloring) {
  EXPECT_TRUE(Solve("mapcoloring").solved);
}

TEST_P(SolverImplTest, SolveNQueens) {
  EXPECT_TRUE(Solve("nqueens").solved);
}

TEST_P(SolverImplTest, SolveNQueens0) {
  EXPECT_TRUE(Solve("nqueens0").solved);
}

TEST_P(SolverImplTest, SolveSched0) {
  EXPECT_EQ(5, Solve("sched0").obj);
}

TEST_P(SolverImplTest, SolveSched1) {
  EXPECT_EQ(5, Solve("sched1").obj);
}

TEST_P(SolverImplTest, SolveSched2) {
  EXPECT_EQ(5, Solve("sched2").obj);
}

TEST_P(SolverImplTest, SolveSendMoreMoney) {
  EXPECT_TRUE(Solve("send-more-money").solved);
}

TEST_P(SolverImplTest, SolveSendMostMoney) {
  EXPECT_NEAR(10876, Solve("send-most-money").obj, 1);
}

TEST_P(SolverImplTest, SolveSeq0) {
  EXPECT_NEAR(332, Solve("seq0").obj, 1e-5);
}

TEST_P(SolverImplTest, SolveSeq0a) {
  EXPECT_NEAR(332, Solve("seq0a").obj, 1e-5);
}

TEST_P(SolverImplTest, SolveSudokuHard) {
  EXPECT_TRUE(Solve("sudokuHard").solved);
}

TEST_P(SolverImplTest, SolveSudokuVeryEasy) {
  EXPECT_TRUE(Solve("sudokuVeryEasy").solved);
}

// ----------------------------------------------------------------------------
// Solve code tests

TEST_P(SolverImplTest, OptimalSolveCode) {
  Problem p;
  EXPECT_TRUE(Solve(p, "objconstint").solved);
  EXPECT_EQ(0, p.solve_code());
}

TEST_P(SolverImplTest, FeasibleSolveCode) {
  Problem p;
  EXPECT_TRUE(Solve(p, "feasible").solved);
  EXPECT_EQ(0, p.solve_code());
}

TEST_P(SolverImplTest, InfeasibleSolveCode) {
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

TEST_P(SolverImplTest, InterruptSolution) {
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
  void HandleFeasibleSolution(
      fmt::StringRef, const double *, const double *, double) {
    ++num_solutions;
  }
};

TEST_P(SolverImplTest, CountSolutions) {
  mp::Problem p;
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeVariable(2)};
  p.AddCon(MakeAllDiff(args));
  solver_->SetIntOption("solutionlimit", 10);
  solver_->SetIntOption("countsolutions", 1);
  SolutionCounter sc;
  solver_->Solve(p, sc);
  EXPECT_EQ(6, sc.num_solutions);
}

TEST_P(SolverImplTest, SatisfactionSolutionLimit) {
  mp::Problem p;
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeVariable(2)};
  p.AddCon(MakeAllDiff(args));
  solver_->SetIntOption("solutionlimit", 5);
  mp::BasicSolutionHandler sh;
  solver_->Solve(p, sh);
  EXPECT_EQ(0, p.solve_code());
}

TEST_P(SolverImplTest, OptimizationSolutionLimit) {
  mp::Problem p;
  p.Read(MP_TEST_DATA_DIR "/photo9");
  solver_->SetIntOption("solutionlimit", 2);
  mp::BasicSolutionHandler sh;
  solver_->Solve(p, sh);
  EXPECT_GE(p.solve_code(), 400);
  EXPECT_LT(p.solve_code(), 499);
}

TEST_P(SolverImplTest, MultipleSolutions) {
  mp::Problem p;
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  p.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeVariable(2)};
  p.AddCon(MakeAllDiff(args));
  solver_->SetIntOption("solutionlimit", 3);
  solver_->SetStrOption("solutionstub", "test");
  SolutionCounter sc;
  solver_->Solve(p, sc);
  EXPECT_EQ(0, p.solve_code());
  EXPECT_EQ(3, sc.num_solutions);
}

TEST_P(SolverImplTest, OptionValues) {
  for (mp::Solver::option_iterator
    i = solver_->option_begin(), e = solver_->option_end(); i != e; ++i) {
    for (mp::ValueArrayRef::iterator j = i->values().begin(),
        value_end = i->values().end(); j != value_end; ++j) {
      EXPECT_TRUE(j->value != 0);
    }
  }
}

TEST_P(SolverImplTest, CreateSolver) {
  EXPECT_STREQ(solver_->name(), mp::CreateSolver(0)->name());
}

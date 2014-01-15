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

#include "tests/solver-test.h"
#include "tests/util.h"

#include <cmath>

using std::string;

using ampl::LogicalExpr;
using ampl::NumericExpr;
using ampl::Problem;
using ampl::UnsupportedExprError;

SolverTest::EvalResult SolverTest::Solve(Problem &p) {
  struct TestSolutionHandler : ampl::BasicSolutionHandler {
    EvalResult result;
    virtual ~TestSolutionHandler() {}
    void HandleSolution(ampl::Problem &p, fmt::StringRef,
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
      need_result ? Infinity : 0, ampl::INTEGER);
  p.AddVar(var1, var1, ampl::INTEGER);
  p.AddVar(var2, var2, ampl::INTEGER);
  p.AddVar(var3, var3, ampl::INTEGER);
  p.AddCon(e);
  return Solve(p);
}

SolverTest::SolverTest()
: solver_(GetParam().create_solver()), features_(GetParam().features),
  x(AddVar(1)), y(AddVar(2)), z(AddVar(3)) {
}

SolveResult SolverTest::Solve(
    ampl::Solver &s, Problem &p, const char *stub, const char *opt) {
  TestSolutionHandler sh;
  s.set_solution_handler(&sh);
  const std::string DATA_DIR = "../data/";
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
  NumericExpr e = AddBinary(OPPLUS, x, y);
  EXPECT_EQ(25, Eval(e, 10, 15));
  EXPECT_EQ(12, Eval(e, 19, -7));
}

TEST_P(SolverTest, Minus) {
  NumericExpr e = AddBinary(OPMINUS, x, y);
  EXPECT_EQ(-5, Eval(e, 10, 15));
  EXPECT_EQ(26, Eval(e, 19, -7));
}

TEST_P(SolverTest, Mult) {
  NumericExpr e = AddBinary(OPMULT, x, y);
  EXPECT_EQ(150, Eval(e, 10, 15));
  EXPECT_EQ(-133, Eval(e, 19, -7));
}

TEST_P(SolverTest, Div) {
  NumericExpr e = AddBinary(OPDIV, x, y);
  if (!HasFeature(feature::DIV)) {
    EXPECT_THROW(Eval(e, 150, 15);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(10, Eval(e, 150, 15));
  EXPECT_EQ(-7, Eval(e, -133, 19));
}

TEST_P(SolverTest, Rem) {
  NumericExpr e = AddBinary(OPREM, x, y);
  EXPECT_EQ(0, Eval(e, 9, 3));
  EXPECT_EQ(2, Eval(e, 8, 3));
  EXPECT_EQ(-2, Eval(e, -8, 3));
  EXPECT_EQ(2, Eval(e, 8, -3));
  EXPECT_EQ(-2, Eval(e, -8, -3));
}

TEST_P(SolverTest, Pow) {
  NumericExpr e = AddBinary(OPPOW, x, y);
  if (!HasFeature(feature::POW)) {
    EXPECT_THROW(Eval(e, 2, 3);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(8, Eval(e, 2, 3));
  EXPECT_EQ(81, Eval(e, 3, 4));
}

TEST_P(SolverTest, NumericLess) {
  NumericExpr e = AddBinary(OPLESS, x, y);
  EXPECT_EQ(0, Eval(e, 10, 15));
  EXPECT_EQ(26, Eval(e, 19, -7));
}

TEST_P(SolverTest, Min) {
  NumericExpr e = AddVarArg(MINLIST, x, y, z);
  EXPECT_EQ(-7, Eval(e, 3, -7, 5));
  EXPECT_EQ(10, Eval(e, 10, 20, 30));
}

TEST_P(SolverTest, Max) {
  NumericExpr e = AddVarArg(MAXLIST, x, y, z);
  EXPECT_EQ(5, Eval(e, 3, -7, 5));
  EXPECT_EQ(30, Eval(e, 30, 20, 10));
}

TEST_P(SolverTest, Floor) {
  NumericExpr e = AddUnary(FLOOR, x);
  EXPECT_EQ(-42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
  if (!HasFeature(feature::FLOAT_CONST)) return;
  EXPECT_EQ(4, Eval(AddUnary(FLOOR, AddNum(4.9))));
  EXPECT_EQ(-5, Eval(AddUnary(FLOOR, AddNum(-4.1))));
}

TEST_P(SolverTest, Ceil) {
  NumericExpr e = AddUnary(CEIL, x);
  EXPECT_EQ(-42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
  if (!HasFeature(feature::FLOAT_CONST)) return;
  EXPECT_EQ(5, Eval(AddUnary(CEIL, AddNum(4.1))));
  EXPECT_EQ(-4, Eval(AddUnary(CEIL, AddNum(-4.9))));
}

TEST_P(SolverTest, Abs) {
  NumericExpr e = AddUnary(ABS, x);
  EXPECT_EQ(42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
}

TEST_P(SolverTest, UnaryMinus) {
  NumericExpr e = AddUnary(OPUMINUS, x);
  EXPECT_EQ(42, Eval(e, -42));
  EXPECT_EQ(-42, Eval(e, 42));
}

TEST_P(SolverTest, If) {
  NumericExpr e = AddIf(AddRelational(EQ, x, AddNum(1)), y, z);
  EXPECT_EQ(42, Eval(e, 1, 42, 10));
  EXPECT_EQ(10, Eval(e, 0, 42, 10));
  EXPECT_EQ(42, Eval(e, 1, 42, 42));
}

TEST_P(SolverTest, Tanh) {
  NumericExpr e = AddBinary(OPMULT, AddNum(2),
        AddUnary(OP_tanh, AddNum((std::log(1.5) - std::log(0.5)) / 2)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(1, Eval(e));
}

TEST_P(SolverTest, Tan) {
  EXPECT_THROW(Eval(AddUnary(OP_tan, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Sqrt) {
  NumericExpr e = AddUnary(OP_sqrt, x);
  if (!HasFeature(feature::SQRT)) {
    EXPECT_THROW(Eval(e, 64);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(8, Eval(e, 64));
}

TEST_P(SolverTest, Sinh) {
  NumericExpr e = AddUnary(OP_sinh, AddNum(std::log(2 + std::sqrt(5.0))));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(2, Eval(e));
}

TEST_P(SolverTest, Sin) {
  EXPECT_THROW(Eval(AddUnary(OP_sin, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Log10) {
  NumericExpr e = AddUnary(OP_log10, AddNum(1000));
  if (!HasFeature(feature::LOG)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(3, Eval(e));
}

TEST_P(SolverTest, Log) {
  NumericExpr e = AddUnary(OP_log, AddNum(std::exp(5.0)));
  if (!HasFeature(feature::LOG)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverTest, Exp) {
  NumericExpr e = AddUnary(OP_exp, AddNum(std::log(5.0)));
  if (!HasFeature(feature::EXP)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverTest, Cosh) {
  double x = 5;
  NumericExpr e = AddUnary(OP_cosh,
      AddNum(std::log(x + std::sqrt(x + 1) * std::sqrt(x - 1))));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverTest, Cos) {
  EXPECT_THROW(Eval(AddUnary(OP_cos, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Atanh) {
  ampl::UnaryExpr x = AddUnary(OP_atanh, AddNum(std::tanh(5.0)));
  NumericExpr e = AddUnary(FLOOR, AddBinary(OPPLUS,
      AddNum(0.5), AddBinary(OPMULT, AddNum(1000000), x)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5000000, Eval(e));
}

TEST_P(SolverTest, Atan2) {
  EXPECT_THROW(Eval(AddBinary(OP_atan2, x, y)), UnsupportedExprError);
}

TEST_P(SolverTest, Atan) {
  EXPECT_THROW(Eval(AddUnary(OP_atan, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Asinh) {
  NumericExpr e = AddUnary(OP_asinh, AddNum(std::sinh(5.0)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverTest, Asin) {
  EXPECT_THROW(Eval(AddUnary(OP_asin, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Acosh) {
  NumericExpr e = AddUnary(OP_acosh, AddNum(std::cosh(5.0)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW(Eval(e);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_P(SolverTest, Acos) {
  EXPECT_THROW(Eval(AddUnary(OP_acos, x)), UnsupportedExprError);
}

TEST_P(SolverTest, Sum) {
  EXPECT_EQ(0, Eval(AddSum()));
  EXPECT_EQ(42, Eval(AddSum(x), 42));
  EXPECT_EQ(123, Eval(AddSum(x, y, z), 100, 20, 3));
}

TEST_P(SolverTest, IntDiv) {
  NumericExpr e = AddBinary(OPintDIV, x, y);
  EXPECT_EQ(3, Eval(e, 9, 3));
  EXPECT_EQ(2, Eval(e, 8, 3));
  EXPECT_EQ(-2, Eval(e, -8, 3));
  EXPECT_EQ(-2, Eval(e, 8, -3));
  EXPECT_EQ(2, Eval(e, -8, -3));
}

TEST_P(SolverTest, Precision) {
  EXPECT_THROW(Eval(AddBinary(OPprecision, x, y)), UnsupportedExprError);
}

TEST_P(SolverTest, Round) {
  EXPECT_EQ(42, Eval(AddBinary(OPround, x, AddNum(0)), 42));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(4, Eval(AddBinary(OPround, AddNum(4.4), AddNum(0))));
    EXPECT_EQ(5, Eval(AddBinary(OPround, AddNum(4.6), AddNum(0))));
    EXPECT_EQ(-4, Eval(AddBinary(OPround, AddNum(-4.4), AddNum(0))));
    EXPECT_EQ(-5, Eval(AddBinary(OPround, AddNum(-4.6), AddNum(0))));
  }
  EXPECT_THROW(Eval(AddBinary(OPround, x, AddNum(1))), UnsupportedExprError);
  EXPECT_THROW(Eval(AddBinary(OPround, x, y)), UnsupportedExprError);
}

TEST_P(SolverTest, Trunc) {
  EXPECT_EQ(42, Eval(AddBinary(OPtrunc, x, AddNum(0)), 42));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(4, Eval(AddBinary(OPtrunc, AddNum(4.4), AddNum(0))));
    EXPECT_EQ(4, Eval(AddBinary(OPtrunc, AddNum(4.6), AddNum(0))));
    EXPECT_EQ(-4, Eval(AddBinary(OPtrunc, AddNum(-4.4), AddNum(0))));
    EXPECT_EQ(-4, Eval(AddBinary(OPtrunc, AddNum(-4.6), AddNum(0))));
  }
  EXPECT_THROW(Eval(AddBinary(OPtrunc, x, AddNum(1))), UnsupportedExprError);
  EXPECT_THROW(Eval(AddBinary(OPtrunc, x, y)), UnsupportedExprError);
}

TEST_P(SolverTest, Count) {
  LogicalExpr a(AddRelational(NE, x, AddNum(0)));
  LogicalExpr b(AddRelational(NE, y, AddNum(0)));
  LogicalExpr c(AddRelational(NE, z, AddNum(0)));
  EXPECT_EQ(0, Eval(AddCount(a, b, c)));
  EXPECT_EQ(1, Eval(AddCount(a, b, c), 1));
  EXPECT_EQ(2, Eval(AddCount(a, b, c), 0, 1, 1));
  EXPECT_EQ(3, Eval(AddCount(a, b, c), 1, 1, 1));
}

TEST_P(SolverTest, NumberOf) {
  ampl::NumericConstant n = AddNum(42);
  EXPECT_EQ(0, Eval(AddNumberOf(n, x)));
  EXPECT_EQ(1, Eval(AddNumberOf(n, x), 42));
  EXPECT_EQ(0, Eval(AddNumberOf(n, x, y)));
  EXPECT_EQ(1, Eval(AddNumberOf(n, x, y), 0, 42));
  EXPECT_EQ(2, Eval(AddNumberOf(n, x, y), 42, 42));
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
  double args[] = {-1, 3, 0, 6, 1};
  NumericExpr e = AddPL(5, args, 1);
  if (!HasFeature(feature::PLTERM)) {
    EXPECT_THROW(Eval(e, 42);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(33, Eval(e, 42));
  EXPECT_EQ(-3, Eval(e, 4));
  EXPECT_EQ(1, Eval(e, -1));
}

TEST_P(SolverTest, UnsupportedFunctionCall) {
  EXPECT_THROW(Eval(AddCall("foo", 1, 2), 3);, UnsupportedExprError);
}

TEST_P(SolverTest, PowConstExp) {
  EXPECT_EQ(16, Eval(AddBinary(OP1POW, x, AddNum(4)), 2));
}

TEST_P(SolverTest, Pow2) {
  EXPECT_EQ(49, Eval(AddUnary(OP2POW, x), 7));
}

TEST_P(SolverTest, PowConstBase) {
  NumericExpr e = AddBinary(OPCPOW, AddNum(5), x);
  if (!HasFeature(feature::POW)) {
    EXPECT_THROW(Eval(e, 3);, UnsupportedExprError);
    return;
  }
  EXPECT_EQ(125, Eval(e, 3));
}

TEST_P(SolverTest, NumericConstant) {
  EXPECT_EQ(42, Eval(AddNum(42)));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(42, Eval(AddBinary(OPMULT, AddNum(0.42), AddNum(100))));
    return;
  }
  EXPECT_THROW_MSG(Eval(AddNum(0.42));, UnsupportedExprError,
    "value 0.42 can't be represented as int");
  EXPECT_THROW(Solve("objconst"), ampl::Error);
}

TEST_P(SolverTest, Var) {
  EXPECT_EQ(11, Eval(x, 11, 22));
  EXPECT_EQ(22, Eval(y, 11, 22));
  EXPECT_EQ(33, Eval(x, 33));
}

TEST_P(SolverTest, Or) {
  NumericExpr one = AddNum(1);
  LogicalExpr e = AddBinaryLogical(
      OPOR, AddRelational(EQ, x, one), AddRelational(EQ, y, one));
  EXPECT_EQ(0, Eval(e, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_P(SolverTest, And) {
  NumericExpr one = AddNum(1);
  LogicalExpr e = AddBinaryLogical(
      OPAND, AddRelational(EQ, x, one), AddRelational(EQ, y, one));
  EXPECT_EQ(0, Eval(e, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_P(SolverTest, Less) {
  LogicalExpr e = AddRelational(LT, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_P(SolverTest, LessEqual) {
  LogicalExpr e = AddRelational(LE, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_P(SolverTest, Equal) {
  LogicalExpr e = AddRelational(EQ, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_P(SolverTest, GreaterEqual) {
  LogicalExpr e = AddRelational(GE, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_P(SolverTest, Greater) {
  LogicalExpr e = AddRelational(GT, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_P(SolverTest, NotEqual) {
  LogicalExpr e = AddRelational(NE, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_P(SolverTest, Not) {
  LogicalExpr e = AddNot(AddRelational(EQ, x, AddNum(1)));
  EXPECT_EQ(1, Eval(e, 0));
  EXPECT_EQ(0, Eval(e, 1));
}

TEST_P(SolverTest, AtLeast) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPATLEAST, AddVar(1), AddCount(a, b));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, AtMost) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPATMOST, AddVar(1), AddCount(a, b));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, Exactly) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPEXACTLY, AddVar(1), AddCount(a, b));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, NotAtLeast) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPNOTATLEAST, AddVar(1), AddCount(a, b));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, NotAtMost) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPNOTATMOST, AddVar(1), AddCount(a, b));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, NotExactly) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPNOTEXACTLY, AddVar(1), AddCount(a, b));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_P(SolverTest, ForAll) {
  LogicalExpr e = AddIteratedLogical(ANDLIST,
      AddRelational(EQ, x, AddNum(1)),
      AddRelational(EQ, y, AddNum(1)),
      AddRelational(EQ, z, AddNum(1)));
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
  LogicalExpr e = AddIteratedLogical(ORLIST,
      AddRelational(EQ, x, AddNum(1)),
      AddRelational(EQ, y, AddNum(1)),
      AddRelational(EQ, z, AddNum(1)));
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
  LogicalExpr e = AddImplication(
      AddRelational(EQ, x, AddNum(1)),
      AddRelational(EQ, y, AddNum(1)),
      AddRelational(EQ, z, AddNum(1)));
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
  LogicalExpr e = AddBinaryLogical(OP_IFF,
      AddRelational(EQ, x, AddNum(1)),
      AddRelational(EQ, y, AddNum(1)));
  EXPECT_EQ(1, Eval(e, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_P(SolverTest, AllDiff) {
  LogicalExpr e = AddAllDiff(AddNum(1), x, y);
  EXPECT_TRUE(Solve(e, 2, 3).has_value());
  EXPECT_FALSE(Solve(e, 2, 1).has_value());
  EXPECT_FALSE(Solve(e, 1, 1).has_value());
}

TEST_P(SolverTest, NestedAllDiff) {
  EXPECT_THROW(Eval(AddNot(AddAllDiff(AddNum(1), x, y)), 1, 2),
      UnsupportedExprError);
}

TEST_P(SolverTest, LogicalConstant) {
  EXPECT_EQ(0, Eval(AddBool(false)));
  EXPECT_EQ(1, Eval(AddBool(true)));
}

TEST_P(SolverTest, NonlinearObj) {
  Problem p;
  p.AddVar(2, 2, ampl::INTEGER);
  ampl::Variable x = AddVar(0);
  p.AddObj(ampl::MIN, AddBinary(OPMULT, x, x));
  EXPECT_EQ(4, Solve(p).obj_value());
}

TEST_P(SolverTest, ObjConst) {
  Problem p;
  p.AddVar(0, 0, ampl::INTEGER);
  p.AddObj(ampl::MIN, AddNum(42));
  EXPECT_EQ(42, Solve(p).obj_value());
}

TEST_P(SolverTest, Minimize) {
  Problem p;
  p.AddVar(42, 100, ampl::INTEGER);
  p.AddObj(ampl::MIN, AddVar(0));
  EXPECT_EQ(42, Solve(p).obj_value());
}

TEST_P(SolverTest, Maximize) {
  Problem p;
  p.AddVar(0, 42, ampl::INTEGER);
  p.AddObj(ampl::MAX, AddVar(0));
  EXPECT_EQ(42, Solve(p).obj_value());
}

TEST_P(SolverTest, TimingOption) {
  struct TestOutputHandler : ampl::OutputHandler {
    std::string output;

    virtual ~TestOutputHandler() {}
    void HandleOutput(fmt::StringRef output) {
      this->output += output;
    }
  };
  TestOutputHandler oh;
  solver_->set_output_handler(&oh);

  Problem p;
  p.AddVar(42, 100, ampl::INTEGER);
  p.AddObj(ampl::MIN, AddVar(0));

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
  while (ampl::SignalHandler::stop())
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

struct SolutionCounter : ampl::BasicSolutionHandler {
  int num_solutions;
  SolutionCounter() : num_solutions(0) {}
  void HandleFeasibleSolution(Problem &, fmt::StringRef,
        const double *, const double *, double) {
    ++num_solutions;
  }
};

TEST_P(SolverTest, CountSolutions) {
  ampl::Problem p;
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddCon(AddAllDiff(AddVar(0), AddVar(1), AddVar(2)));
  SolutionCounter sc;
  solver_->SetIntOption("solutionlimit", 10);
  solver_->SetIntOption("countsolutions", 1);
  solver_->set_solution_handler(&sc);
  solver_->Solve(p);
  EXPECT_EQ(6, sc.num_solutions);
}

TEST_P(SolverTest, SatisfactionSolutionLimit) {
  ampl::Problem p;
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddCon(AddAllDiff(AddVar(0), AddVar(1), AddVar(2)));
  solver_->SetIntOption("solutionlimit", 5);
  solver_->Solve(p);
  EXPECT_EQ(0, p.solve_code());
}

TEST_P(SolverTest, OptimizationSolutionLimit) {
  ampl::Problem p;
  p.Read("../data/photo9");
  solver_->SetIntOption("solutionlimit", 2);
  solver_->Solve(p);
  EXPECT_GE(p.solve_code(), 400);
  EXPECT_LT(p.solve_code(), 499);
}

TEST_P(SolverTest, MultipleSolutions) {
  ampl::Problem p;
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddVar(1, 3, ampl::INTEGER);
  p.AddCon(AddAllDiff(AddVar(0), AddVar(1), AddVar(2)));
  SolutionCounter sc;
  solver_->set_solution_handler(&sc);
  solver_->SetIntOption("solutionlimit", 3);
  solver_->SetStrOption("solutionstub", "test");
  solver_->Solve(p);
  EXPECT_EQ(0, p.solve_code());
  EXPECT_EQ(3, sc.num_solutions);
}

TEST_P(SolverTest, OptionValues) {
  for (ampl::Solver::option_iterator
    i = solver_->option_begin(), e = solver_->option_end(); i != e; ++i) {
    for (ampl::ValueArrayRef::iterator j = i->values().begin(),
        value_end = i->values().end(); j != value_end; ++j) {
      EXPECT_TRUE(j->value != 0);
    }
  }
}

TEST_P(SolverTest, CreateSolver) {
  EXPECT_STREQ(solver_->name(), ampl::CreateSolver()->name());
}

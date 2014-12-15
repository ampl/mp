/*
 Expression writer tests

 Copyright (C) 2014 AMPL Optimization Inc

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

#include <gtest/gtest.h>

#include "mp/arrayref.h"
#include "mp/expr.h"
#include "expr-writer.h"

namespace ex = mp::expr;

using mp::NumericExpr;
using mp::LogicalExpr;

class ExprWriterTest : public ::testing::Test, protected mp::ExprFactory {
 protected:
  mp::LogicalConstant l0, l1;

  ExprWriterTest() {
    l0 = MakeLogicalConstant(false);
    l1 = MakeLogicalConstant(true);
  }

  mp::NumericConstant MakeConst(double value) {
    return MakeNumericConstant(value);
  }

  mp::CallExpr MakeCall(mp::Function func, mp::ArrayRef<mp::Expr> args) {
    int num_args = args.size();
    auto builder = BeginCall(func, num_args);
    for (auto i = 0; i < num_args; ++i)
      builder.AddArg(args[i]);
    return EndCall(builder);
  }

  mp::IteratedExpr MakeIterated(ex::Kind kind, mp::ArrayRef<NumericExpr> args) {
    int num_args = args.size();
    auto builder = BeginIterated(kind, num_args);
    for (auto i = 0; i < num_args; ++i)
      builder.AddArg(args[i]);
    return EndIterated(builder);
  }

  mp::CountExpr MakeCount(mp::ArrayRef<LogicalExpr> args) {
    int num_args = args.size();
    auto builder = BeginCount(num_args);
    for (auto i = 0; i < num_args; ++i)
      builder.AddArg(args[i]);
    return EndCount(builder);
  }

  mp::PairwiseExpr MakePairwise(ex::Kind kind, mp::ArrayRef<NumericExpr> args) {
    int num_args = args.size();
    auto builder = BeginPairwise(kind, num_args);
    for (auto i = 0; i < num_args; ++i)
      builder.AddArg(args[i]);
    return EndPairwise(builder);
  }

  mp::IteratedLogicalExpr MakeIteratedLogical(
      ex::Kind kind, mp::ArrayRef<LogicalExpr> args) {
    int num_args = args.size();
    auto builder = BeginIteratedLogical(kind, num_args);
    for (auto i = 0; i < num_args; ++i)
      builder.AddArg(args[i]);
    return EndIteratedLogical(builder);
  }
};

// Checks if formatting expr produces expected output.
#define CHECK_WRITE(expected_output, expr) \
  EXPECT_EQ(expected_output, fmt::format("{}", static_cast<NumericExpr>(expr)))

TEST_F(ExprWriterTest, WriteNumericConstant) {
  CHECK_WRITE("0", MakeConst(0));
  CHECK_WRITE("42", MakeConst(42));
  CHECK_WRITE("12.34", MakeConst(12.34));
}

TEST_F(ExprWriterTest, WriteVariable) {
  CHECK_WRITE("x1", MakeVariable(0));
  CHECK_WRITE("x3", MakeVariable(2));
}

TEST_F(ExprWriterTest, WriteUnaryExpr) {
  auto x1 = MakeVariable(0);
  CHECK_WRITE("-x1", MakeUnary(ex::MINUS, x1));
  CHECK_WRITE("abs(x1)", MakeUnary(ex::ABS, x1));
  CHECK_WRITE("floor(x1)", MakeUnary(ex::FLOOR, x1));
  CHECK_WRITE("ceil(x1)", MakeUnary(ex::CEIL, x1));
  CHECK_WRITE("sqrt(x1)", MakeUnary(ex::SQRT, x1));
  CHECK_WRITE("x1 ^ 2", MakeUnary(ex::POW2, x1));
  CHECK_WRITE("exp(x1)", MakeUnary(ex::EXP, x1));
  CHECK_WRITE("log(x1)", MakeUnary(ex::LOG, x1));
  CHECK_WRITE("log10(x1)", MakeUnary(ex::LOG10, x1));
  CHECK_WRITE("sin(x1)", MakeUnary(ex::SIN, x1));
  CHECK_WRITE("sinh(x1)", MakeUnary(ex::SINH, x1));
  CHECK_WRITE("cos(x1)", MakeUnary(ex::COS, x1));
  CHECK_WRITE("cosh(x1)", MakeUnary(ex::COSH, x1));
  CHECK_WRITE("tan(x1)", MakeUnary(ex::TAN, x1));
  CHECK_WRITE("tanh(x1)", MakeUnary(ex::TANH, x1));
  CHECK_WRITE("asin(x1)", MakeUnary(ex::ASIN, x1));
  CHECK_WRITE("asinh(x1)", MakeUnary(ex::ASINH, x1));
  CHECK_WRITE("acos(x1)", MakeUnary(ex::ACOS, x1));
  CHECK_WRITE("acosh(x1)", MakeUnary(ex::ACOSH, x1));
  CHECK_WRITE("atan(x1)", MakeUnary(ex::ATAN, x1));
  CHECK_WRITE("atanh(x1)", MakeUnary(ex::ATANH, x1));
}

TEST_F(ExprWriterTest, WriteBinaryExpr) {
  auto x1 = MakeVariable(0);
  auto n42 = MakeConst(42);
  CHECK_WRITE("x1 + 42", MakeBinary(ex::ADD, x1, n42));
  CHECK_WRITE("x1 - 42", MakeBinary(ex::SUB, x1, n42));
  CHECK_WRITE("x1 * 42", MakeBinary(ex::MUL, x1, n42));
  CHECK_WRITE("x1 / 42", MakeBinary(ex::DIV, x1, n42));
  CHECK_WRITE("x1 mod 42", MakeBinary(ex::MOD, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(ex::POW, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(ex::POW, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(ex::POW_CONST_BASE, x1, n42));
  CHECK_WRITE("x1 ^ 42", MakeBinary(ex::POW_CONST_EXP, x1, n42));
  CHECK_WRITE("x1 less 42", MakeBinary(ex::LESS, x1, n42));
  CHECK_WRITE("x1 div 42", MakeBinary(ex::INT_DIV, x1, n42));
}

TEST_F(ExprWriterTest, WriteBinaryFunc) {
  auto x1 = MakeVariable(0);
  auto n42 = MakeConst(42);
  CHECK_WRITE("atan2(x1, 42)", MakeBinary(ex::ATAN2, x1, n42));
  CHECK_WRITE("precision(x1, 42)", MakeBinary(ex::PRECISION, x1, n42));
  CHECK_WRITE("round(x1, 42)", MakeBinary(ex::ROUND, x1, n42));
  CHECK_WRITE("trunc(x1, 42)", MakeBinary(ex::TRUNC, x1, n42));
}

TEST_F(ExprWriterTest, WriteIfExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if x1 = 0 then 1",
      MakeIf(MakeRelational(ex::EQ, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 = 0 then 0 else 1",
      MakeIf(MakeRelational(ex::EQ, MakeVariable(0), n0), n0, n1));
}

TEST_F(ExprWriterTest, WritePiecewiseLinearExpr) {
  auto plterm = BeginPLTerm(2);
  plterm.AddSlope(-1);
  plterm.AddBreakpoint(5);
  plterm.AddSlope(0);
  plterm.AddBreakpoint(10);
  plterm.AddSlope(1);
  CHECK_WRITE("<<5, 10; -1, 0, 1>> x43", EndPLTerm(plterm, MakeVariable(42)));
}

TEST_F(ExprWriterTest, WriteCallExpr) {
  auto f = AddFunction("foo", -1);
  mp::Expr args[] = {
      MakeConst(3),
      MakeBinary(ex::ADD, MakeVariable(0), MakeConst(5)),
      MakeConst(7),
      MakeVariable(1)
  };
  CHECK_WRITE("foo()", MakeCall(f, MakeArrayRef(args, 0)));
  CHECK_WRITE("foo(3)", MakeCall(f, MakeArrayRef(args, 1)));
  CHECK_WRITE("foo(3, x1 + 5, 7)", MakeCall(f, MakeArrayRef(args, 3)));
  CHECK_WRITE("foo(3, x1 + 5, 7, x2)", MakeCall(f, args));
}

TEST_F(ExprWriterTest, WriteVarArgExpr) {
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  CHECK_WRITE("min(x1, x2, 42)", MakeIterated(ex::MIN, args));
  CHECK_WRITE("max(x1, x2, 42)", MakeIterated(ex::MAX, args));
}

TEST_F(ExprWriterTest, WriteSumExpr) {
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  CHECK_WRITE("/* sum */ (x1 + x2 + 42)", MakeIterated(ex::SUM, args));
  NumericExpr args2[] = {
    MakeBinary(ex::ADD, MakeVariable(0), MakeVariable(1)), MakeConst(42)
  };
  CHECK_WRITE("/* sum */ ((x1 + x2) + 42)", MakeIterated(ex::SUM, args2));
}

TEST_F(ExprWriterTest, WriteCountExpr) {
  LogicalExpr args[] = {
    MakeRelational(ex::EQ, MakeVariable(0), MakeConst(0)), l1, l0
  };
  CHECK_WRITE("count(x1 = 0, 1, 0)", MakeCount(args));
}

TEST_F(ExprWriterTest, WriteNumberOfExpr) {
  NumericExpr args[] = {MakeConst(42), MakeConst(43), MakeConst(44)};
  CHECK_WRITE("numberof 42 in (43, 44)", MakeIterated(ex::NUMBEROF, args));
}

TEST_F(ExprWriterTest, WriteNotExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if !(x1 = 0) then 1",
      MakeIf(MakeNot(MakeRelational(ex::EQ, MakeVariable(0), n0)), n1, n0));
}

TEST_F(ExprWriterTest, WriteBinaryLogicalExpr) {
  auto e1 = MakeRelational(ex::GT, MakeVariable(0), MakeConst(0));
  auto e2 = MakeRelational(ex::LT, MakeVariable(0), MakeConst(10));
  CHECK_WRITE("if x1 > 0 || x1 < 10 then 1",
      MakeIf(MakeBinaryLogical(ex::OR, e1, e2), MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if x1 > 0 && x1 < 10 then 1",
      MakeIf(MakeBinaryLogical(ex::AND, e1, e2), MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if x1 > 0 <==> x1 < 10 then 1",
      MakeIf(MakeBinaryLogical(ex::IFF, e1, e2), MakeConst(1), MakeConst(0)));
}

TEST_F(ExprWriterTest, WriteRelationalExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if x1 < 0 then 1",
      MakeIf(MakeRelational(ex::LT, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 <= 0 then 1",
      MakeIf(MakeRelational(ex::LE, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 = 0 then 1",
      MakeIf(MakeRelational(ex::EQ, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 >= 0 then 1",
      MakeIf(MakeRelational(ex::GE, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 > 0 then 1",
      MakeIf(MakeRelational(ex::GT, MakeVariable(0), n0), n1, n0));
  CHECK_WRITE("if x1 != 0 then 1",
      MakeIf(MakeRelational(ex::NE, MakeVariable(0), n0), n1, n0));
}

TEST_F(ExprWriterTest, WriteLogicalCountExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1), value = MakeConst(42);
  LogicalExpr args[] = {
    MakeRelational(ex::EQ, MakeVariable(0), MakeConst(0)), l1, l0
  };
  auto count = MakeCount(args);
  CHECK_WRITE("if atleast 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::ATLEAST, value, count), n1, n0));
  CHECK_WRITE("if atmost 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::ATMOST, value, count), n1, n0));
  CHECK_WRITE("if exactly 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::EXACTLY, value, count), n1, n0));
  CHECK_WRITE("if !atleast 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_ATLEAST, value, count), n1, n0));
  CHECK_WRITE("if !atmost 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_ATMOST, value, count), n1, n0));
  CHECK_WRITE("if !exactly 42 (x1 = 0, 1, 0) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_EXACTLY, value, count), n1, n0));
}

TEST_F(ExprWriterTest, WriteImplicationExpr) {
  auto e1 = MakeRelational(ex::EQ, MakeVariable(0), MakeConst(0));
  CHECK_WRITE("if x1 = 0 ==> 1 then 1",
      MakeIf(MakeImplication(e1, l1, l0), MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if x1 = 0 ==> 0 else 1 then 1",
      MakeIf(MakeImplication(e1, l0, l1), MakeConst(1), MakeConst(0)));
}

TEST_F(ExprWriterTest, WriteIteratedLogicalExpr) {
  LogicalExpr args[] = {
    MakeRelational(ex::EQ, MakeVariable(0), MakeConst(0)), l1, l0
  };
  CHECK_WRITE("if /* forall */ (x1 = 0 && 1 && 0) then 1",
      MakeIf(MakeIteratedLogical(ex::FORALL, args),
             MakeConst(1), MakeConst(0)));
  CHECK_WRITE("if /* exists */ (x1 = 0 || 1 || 0) then 1",
      MakeIf(MakeIteratedLogical(ex::EXISTS, args),
             MakeConst(1), MakeConst(0)));
}

TEST_F(ExprWriterTest, WritePairwiseExpr) {
  NumericExpr args[] = {MakeConst(42), MakeConst(43), MakeConst(44)};
  CHECK_WRITE("if alldiff(42, 43, 44) then 1",
      MakeIf(MakePairwise(ex::ALLDIFF, args), MakeConst(1), MakeConst(0)));
}

TEST_F(ExprWriterTest, WriteStringLiteral) {
  mp::Expr args[] = {MakeStringLiteral("abc")};
  auto f = AddFunction("f", 1, mp::func::SYMBOLIC);
  CHECK_WRITE("f('abc')", MakeCall(f, args));
  args[0] = MakeStringLiteral("ab'c");
  CHECK_WRITE("f('ab''c')", MakeCall(f, args));
  args[0] = MakeStringLiteral("ab\nc");
  CHECK_WRITE("f('ab\\\nc')", MakeCall(f, args));
}

TEST_F(ExprWriterTest, UnaryExprPrecedence) {
  auto x1 = MakeVariable(0);
  CHECK_WRITE("--x1", MakeUnary(ex::MINUS, MakeUnary(ex::MINUS, x1)));
  CHECK_WRITE("-(x1 ^ x1)", MakeUnary(ex::MINUS, MakeBinary(ex::POW, x1, x1)));
}

TEST_F(ExprWriterTest, UnaryFuncPrecedence) {
  auto x1 = MakeVariable(0);
  for (int i = ex::FIRST_UNARY; i < ex::LAST_UNARY; ++i) {
    ex::Kind kind = static_cast<ex::Kind>(i);
    if (kind == ex::MINUS || kind == ex::POW2)
      continue;
    CHECK_WRITE(fmt::format("{0}({0}(x1))", str(kind)),
        MakeUnary(kind, MakeUnary(kind, x1)));
    CHECK_WRITE(fmt::format("{0}(x1 + x1)", str(kind)),
        MakeUnary(kind, MakeBinary(ex::ADD, x1, x1)));
  }
}

TEST_F(ExprWriterTest, Pow2Precedence) {
  auto x1 = MakeVariable(0);
  CHECK_WRITE("(x1 ^ 2) ^ 2", MakeUnary(ex::POW2, MakeUnary(ex::POW2, x1)));
  CHECK_WRITE("(x1 * x1) ^ 2",
              MakeUnary(ex::POW2, MakeBinary(ex::MUL, x1, x1)));
}

TEST_F(ExprWriterTest, AdditiveExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  CHECK_WRITE("x1 + x2 + x3",
      MakeBinary(ex::ADD, MakeBinary(ex::ADD, x1, x2), x3));
  CHECK_WRITE("x1 + x2 - x3",
      MakeBinary(ex::SUB, MakeBinary(ex::ADD, x1, x2), x3));
  CHECK_WRITE("x1 + x2 less x3",
      MakeBinary(ex::LESS, MakeBinary(ex::ADD, x1, x2), x3));
  CHECK_WRITE("x1 + (x2 + x3)",
      MakeBinary(ex::ADD, x1, MakeBinary(ex::ADD, x2, x3)));
  CHECK_WRITE("(x1 + x2) * x3",
      MakeBinary(ex::MUL, MakeBinary(ex::ADD, x1, x2), x3));
  CHECK_WRITE("if 1 then x1 else x2 + x3",
      MakeIf(l1, x1, MakeBinary(ex::ADD, x2, x3)));
}

TEST_F(ExprWriterTest, MultiplicativeExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  CHECK_WRITE("x1 * x2 * x3",
      MakeBinary(ex::MUL, MakeBinary(ex::MUL, x1, x2), x3));
  CHECK_WRITE("x1 * x2 / x3",
      MakeBinary(ex::DIV, MakeBinary(ex::MUL, x1, x2), x3));
  CHECK_WRITE("x1 * x2 div x3",
      MakeBinary(ex::INT_DIV, MakeBinary(ex::MUL, x1, x2), x3));
  CHECK_WRITE("x1 * x2 mod x3",
      MakeBinary(ex::MOD, MakeBinary(ex::MUL, x1, x2), x3));
  CHECK_WRITE("x1 * (x2 * x3)",
      MakeBinary(ex::MUL, x1, MakeBinary(ex::MUL, x2, x3)));
  CHECK_WRITE("(x1 * x2) ^ x3",
      MakeBinary(ex::POW, MakeBinary(ex::MUL, x1, x2), x3));
  CHECK_WRITE("(x1 + x2) * x3",
      MakeBinary(ex::MUL, MakeBinary(ex::ADD, x1, x2), x3));
}

TEST_F(ExprWriterTest, ExponentiationExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  CHECK_WRITE("x1 ^ x2 ^ x3",
      MakeBinary(ex::POW, x1, MakeBinary(ex::POW, x2, x3)));
  CHECK_WRITE("x1 ^ x2 ^ 3",
      MakeBinary(ex::POW, x1, MakeBinary(ex::POW_CONST_BASE, x2,
                                         MakeConst(3))));
  CHECK_WRITE("x1 ^ 3 ^ x2",
      MakeBinary(ex::POW, x1, MakeBinary(ex::POW_CONST_EXP, MakeConst(3), x2)));
  CHECK_WRITE("(x1 ^ 2) ^ 3",
      MakeBinary(ex::POW_CONST_BASE, MakeBinary(ex::POW_CONST_BASE, x1,
                                                MakeConst(2)), MakeConst(3)));
  CHECK_WRITE("-x1 ^ -x2",
      MakeBinary(ex::POW, MakeUnary(ex::MINUS, x1), MakeUnary(ex::MINUS, x2)));
  CHECK_WRITE("x1 ^ (x2 * x3)",
      MakeBinary(ex::POW, x1, MakeBinary(ex::MUL, x2, x3)));
}

TEST_F(ExprWriterTest, BinaryFuncPrecedence) {
  auto x1 = MakeVariable(0);
  auto e = MakeBinary(ex::ADD, x1, x1);
  CHECK_WRITE("atan2(atan2(x1, x1), x1 + x1)",
      MakeBinary(ex::ATAN2, MakeBinary(ex::ATAN2, x1, x1), e));
  CHECK_WRITE("precision(precision(x1, x1), x1 + x1)",
      MakeBinary(ex::PRECISION, MakeBinary(ex::PRECISION, x1, x1), e));
  CHECK_WRITE("round(round(x1, x1), x1 + x1)",
      MakeBinary(ex::ROUND, MakeBinary(ex::ROUND, x1, x1), e));
  CHECK_WRITE("trunc(trunc(x1, x1), x1 + x1)",
      MakeBinary(ex::TRUNC, MakeBinary(ex::TRUNC, x1, x1), e));
}

TEST_F(ExprWriterTest, IfExprPrecedence) {
  NumericExpr n0 = MakeConst(0), n1 = MakeConst(1), n2 = MakeConst(2);
  LogicalExpr e = MakeBinaryLogical(ex::OR, l0, l1);
  CHECK_WRITE("if 0 || 1 then if 0 || 1 then 1",
      MakeIf(e, MakeIf(e, n1, n0), n0));
  CHECK_WRITE("if 0 || 1 then if 0 || 1 then 1 else 2",
      MakeIf(e, MakeIf(e, n1, n2), n0));
  CHECK_WRITE("if 0 || 1 then (if 0 || 1 then 1) else 2",
      MakeIf(e, MakeIf(e, n1, n0), n2));
  CHECK_WRITE("if 0 || 1 then 0 else if 0 || 1 then 1 else 2",
      MakeIf(e, n0, MakeIf(e, n1, n2)));
  CHECK_WRITE("if !(0 || 1) then x1 + 1",
      MakeIf(MakeNot(e), MakeBinary(ex::ADD, MakeVariable(0), n1), n0));
}

TEST_F(ExprWriterTest, PiecewiseLinearExprPrecedence) {
  auto plterm = BeginPLTerm(2);
  plterm.AddSlope(-1);
  plterm.AddBreakpoint(5);
  plterm.AddSlope(0);
  plterm.AddBreakpoint(10);
  plterm.AddSlope(1);
  CHECK_WRITE("<<5, 10; -1, 0, 1>> x43 ^ 2",
      MakeBinary(ex::POW, EndPLTerm(plterm, MakeVariable(42)), MakeConst(2)));
}

TEST_F(ExprWriterTest, CallExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1);
  auto f = AddFunction("foo", -1);
  mp::Expr args[] = {
      MakeCall(f, mp::MakeArrayRef<mp::Expr>(0, 0)),
      MakeBinary(ex::ADD, x1, MakeConst(5)),
      MakeConst(7),
      MakeUnary(ex::FLOOR, x2)
  };
  CHECK_WRITE("foo(foo(), x1 + 5, 7, floor(x2))", MakeCall(f, args));
}

TEST_F(ExprWriterTest, VarArgExprPrecedence) {
  NumericExpr x1 = MakeVariable(0), x2 = MakeVariable(1);
  NumericExpr e = MakeBinary(ex::ADD, x1, x2);
  NumericExpr args[] = {e, e};
  CHECK_WRITE("min(x1 + x2, x1 + x2)", MakeIterated(ex::MIN, args));
  CHECK_WRITE("max(x1 + x2, x1 + x2)", MakeIterated(ex::MAX, args));
  NumericExpr args2[] = {x1}, args3[] = {x2};
  NumericExpr args4[] = {
    MakeIterated(ex::MIN, args2), MakeIterated(ex::MIN, args3)
  };
  CHECK_WRITE("min(min(x1), min(x2))", MakeIterated(ex::MIN, args4));
  NumericExpr args5[] = {
    MakeIterated(ex::MAX, args2), MakeIterated(ex::MAX, args3)
  };
  CHECK_WRITE("max(max(x1), max(x2))", MakeIterated(ex::MAX, args5));
}

TEST_F(ExprWriterTest, SumExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1), x3 = MakeVariable(2);
  NumericExpr args1[] = {x2, x3};
  NumericExpr args2[] = {x1, MakeIterated(ex::SUM, args1)};
  CHECK_WRITE("/* sum */ (x1 + /* sum */ (x2 + x3))",
              MakeIterated(ex::SUM, args2));
  NumericExpr args3[] = {x1, MakeBinary(ex::MUL, x2, x3)};
  CHECK_WRITE("/* sum */ (x1 + x2 * x3)", MakeIterated(ex::SUM, args3));
  NumericExpr args4[] = {MakeBinary(ex::ADD, x1, x2), x3};
  CHECK_WRITE("/* sum */ ((x1 + x2) + x3)", MakeIterated(ex::SUM, args4));
}

TEST_F(ExprWriterTest, CountExprPrecedence) {
  auto e = MakeBinaryLogical(ex::OR, l0, l1);
  LogicalExpr args[] = {e, e};
  CHECK_WRITE("count(0 || 1, 0 || 1)", MakeCount(args));
}

TEST_F(ExprWriterTest, NumberOfExprPrecedence) {
  auto x1 = MakeVariable(0), x2 = MakeVariable(1);
  NumericExpr args[] = {MakeConst(42), x1, x2};
  auto e = MakeIterated(ex::NUMBEROF, args);
  NumericExpr args2[] = {e, e, e};
  CHECK_WRITE("numberof numberof 42 in (x1, x2) in ("
      "numberof 42 in (x1, x2), numberof 42 in (x1, x2))",
              MakeIterated(ex::NUMBEROF, args2));
  auto e2 = MakeBinary(ex::ADD, x1, x2);
  NumericExpr args3[] = {e2, e2, e2};
  CHECK_WRITE("numberof x1 + x2 in (x1 + x2, x1 + x2)",
              MakeIterated(ex::NUMBEROF, args3));
}

TEST_F(ExprWriterTest, NotExprPrecedence) {
  NumericExpr n0 = MakeConst(0), n1 = MakeConst(1), x = MakeVariable(0);
  CHECK_WRITE("if !!(x1 = 0) then 1",
      MakeIf(MakeNot(MakeNot(MakeRelational(ex::EQ, x, n0))), n1, n0));
}

TEST_F(ExprWriterTest, LogicalOrExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 || 1 || 1 then 1",
      MakeIf(MakeBinaryLogical(ex::OR, MakeBinaryLogical(ex::OR, l0, l1), l1),
          n1, n0));
  CHECK_WRITE("if 0 || (1 || 1) then 1",
      MakeIf(MakeBinaryLogical(ex::OR, l0, MakeBinaryLogical(ex::OR, l1, l1)),
          n1, n0));
  CHECK_WRITE("if 0 || 1 && 1 then 1",
      MakeIf(MakeBinaryLogical(ex::OR, l0, MakeBinaryLogical(ex::AND, l1, l1)),
          n1, n0));
}

TEST_F(ExprWriterTest, LogicalAndExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 && 1 && 1 then 1",
      MakeIf(MakeBinaryLogical(ex::AND, MakeBinaryLogical(ex::AND, l0, l1), l1),
          n1, n0));
  CHECK_WRITE("if 0 && (1 && 1) then 1",
      MakeIf(MakeBinaryLogical(ex::AND, l0, MakeBinaryLogical(ex::AND, l1, l1)),
          n1, n0));
  CHECK_WRITE("if 0 <= 1 && 1 then 1",
      MakeIf(MakeBinaryLogical(ex::AND,
                               MakeRelational(ex::LE, n0, n1), l1), n1, n0));
}

TEST_F(ExprWriterTest, IffExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 <==> 1 <==> 1 then 1",
      MakeIf(MakeBinaryLogical(ex::IFF, MakeBinaryLogical(ex::IFF, l0, l1), l1),
          n1, n0));
  CHECK_WRITE("if 0 <==> (1 <==> 1) then 1",
      MakeIf(MakeBinaryLogical(ex::IFF, l0, MakeBinaryLogical(ex::IFF, l1, l1)),
          n1, n0));
  CHECK_WRITE("if (0 <==> 1) && 1 then 1",
      MakeIf(MakeBinaryLogical(ex::AND, MakeBinaryLogical(ex::IFF, l0, l1), l1),
          n1, n0));
}

TEST_F(ExprWriterTest, RelationalExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  auto e1 = MakeBinary(ex::ADD, MakeVariable(0), n1);
  auto e2 = MakeBinary(ex::ADD, MakeVariable(1), n1);
  CHECK_WRITE("if x1 + 1 < x2 + 1 then 1",
      MakeIf(MakeRelational(ex::LT, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 <= x2 + 1 then 1",
      MakeIf(MakeRelational(ex::LE, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 = x2 + 1 then 1",
      MakeIf(MakeRelational(ex::EQ, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 >= x2 + 1 then 1",
      MakeIf(MakeRelational(ex::GE, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 > x2 + 1 then 1",
      MakeIf(MakeRelational(ex::GT, e1, e2), n1, n0));
  CHECK_WRITE("if x1 + 1 != x2 + 1 then 1",
      MakeIf(MakeRelational(ex::NE, e1, e2), n1, n0));
}

TEST_F(ExprWriterTest, LogicalCountExprPrecedence) {
  NumericExpr n0 = MakeConst(0), n1 = MakeConst(1), lhs = MakeConst(42);
  LogicalExpr args[] = {
    MakeRelational(ex::EQ, MakeVariable(0), MakeConst(0)), l1
  };
  auto count1 = MakeLogicalCount(ex::ATLEAST, lhs, MakeCount(args));
  LogicalExpr args2[] = {count1, l1};
  auto count2 = MakeCount(args2);
  CHECK_WRITE("if atleast 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::ATLEAST, lhs, count2), n1, n0));
  CHECK_WRITE("if atmost 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::ATMOST, lhs, count2), n1, n0));
  CHECK_WRITE("if exactly 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::EXACTLY, lhs, count2), n1, n0));
  CHECK_WRITE("if !atleast 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_ATLEAST, lhs, count2), n1, n0));
  CHECK_WRITE("if !atmost 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_ATMOST, lhs, count2), n1, n0));
  CHECK_WRITE("if !exactly 42 (atleast 42 (x1 = 0, 1), 1) then 1",
      MakeIf(MakeLogicalCount(ex::NOT_EXACTLY, lhs, count2), n1, n0));

  args2[0] = l0;
  auto count = MakeCount(args2);
  CHECK_WRITE("if atleast 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::ATLEAST, lhs, count), l0), n1, n0));
  CHECK_WRITE("if atmost 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::ATMOST, lhs, count), l0), n1, n0));
  CHECK_WRITE("if exactly 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::EXACTLY, lhs, count), l0), n1, n0));
  CHECK_WRITE("if !atleast 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::NOT_ATLEAST, lhs, count), l0),
             n1, n0));
  CHECK_WRITE("if !atmost 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::NOT_ATMOST, lhs, count), l0),
             n1, n0));
  CHECK_WRITE("if !exactly 42 (0, 1) || 0 then 1",
      MakeIf(MakeBinaryLogical(
               ex::OR, MakeLogicalCount(ex::NOT_EXACTLY, lhs, count), l0),
             n1, n0));
}

TEST_F(ExprWriterTest, IteratedLogicalExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  LogicalExpr args[] = {MakeBinaryLogical(ex::AND, l0, l0), l0};
  CHECK_WRITE("if /* forall */ ((0 && 0) && 0) then 1",
      MakeIf(MakeIteratedLogical(ex::FORALL, args), n1, n0));
  args[0] = MakeBinaryLogical(ex::OR, l0, l0);
  CHECK_WRITE("if /* exists */ ((0 || 0) || 0) then 1",
      MakeIf(MakeIteratedLogical(ex::EXISTS, args), n1, n0));
  args[0] = l0;
  LogicalExpr args2[] = {MakeIteratedLogical(ex::FORALL, args), l0};
  CHECK_WRITE("if /* forall */ (/* forall */ (0 && 0) && 0) then 1",
      MakeIf(MakeIteratedLogical(ex::FORALL, args2), n1, n0));
  args2[0] = MakeIteratedLogical(ex::EXISTS, args);
  CHECK_WRITE("if /* exists */ (/* exists */ (0 || 0) || 0) then 1",
      MakeIf(MakeIteratedLogical(ex::EXISTS, args2), n1, n0));
}

TEST_F(ExprWriterTest, ImplicationExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  CHECK_WRITE("if 0 ==> 1 ==> 0 then 1",
      MakeIf(MakeImplication(MakeImplication(l0, l1, l0), l0, l0), n1, n0));
  CHECK_WRITE("if 0 ==> 1 ==> 0 else 1 then 1",
      MakeIf(MakeImplication(MakeImplication(l0, l1, l0), l0, l1), n1, n0));
  CHECK_WRITE("if 0 ==> 1 else 0 ==> 1 then 1",
      MakeIf(MakeImplication(l0, l1, MakeImplication(l0, l1, l0)), n1, n0));
  CHECK_WRITE("if 0 ==> (1 ==> 0) else 1 then 1",
      MakeIf(MakeImplication(l0, MakeImplication(l1, l0, l0), l1), n1, n0));
  CHECK_WRITE("if 0 ==> (1 ==> 0 else 1) then 1",
      MakeIf(MakeImplication(l0, MakeImplication(l1, l0, l1), l0), n1, n0));
  CHECK_WRITE("if 0 ==> 1 || 0 else 1 then 1",
              MakeIf(MakeImplication(
                       l0, MakeBinaryLogical(ex::OR, l1, l0), l1), n1, n0));
  CHECK_WRITE("if 0 ==> (1 <==> 0) else 1 then 1",
              MakeIf(MakeImplication(
                       l0, MakeBinaryLogical(ex::IFF, l1, l0), l1), n1, n0));
}

TEST_F(ExprWriterTest, PairwiseExprPrecedence) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  NumericExpr args[] = {
    MakeBinary(ex::ADD, n0, n1), MakeBinary(ex::ADD, n0, n1)
  };
  CHECK_WRITE("if alldiff(0 + 1, 0 + 1) then 1",
      MakeIf(MakePairwise(ex::ALLDIFF, args), n1, n0));
}

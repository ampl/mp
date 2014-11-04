/*
 Expression tests

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

#include "mp/expr.h"
#include "gtest-extra.h"

using mp::ExprFactory;
namespace expr = mp::expr;

TEST(ExprTest, Expr) {
  mp::Expr e;
  EXPECT_TRUE(e == 0);
}

TEST(ExprTest, NumericExpr) {
  mp::NumericExpr e;
  EXPECT_TRUE(e == 0);
}

TEST(ExprTest, LogicalExpr) {
  mp::LogicalExpr e;
  EXPECT_TRUE(e == 0);
}

TEST(ExprTest, NumericConstant) {
  mp::NumericConstant e;
  EXPECT_TRUE(e == 0);
  ExprFactory factory;
  e = factory.MakeNumericConstant(1.23);
  EXPECT_EQ(expr::CONSTANT, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(1.23, e.value());
}

TEST(ExprTest, Variable) {
  mp::Variable e;
  EXPECT_TRUE(e == 0);
  ExprFactory factory;
  e = factory.MakeVariable(42);
  EXPECT_EQ(expr::VARIABLE, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(42, e.index());
}

TEST(ExprTest, UnaryExpr) {
  mp::UnaryExpr e;
  EXPECT_TRUE(e == 0);
  ExprFactory factory;
  auto arg = factory.MakeNumericConstant(42);
  e = factory.MakeUnary(expr::ABS, arg);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::ABS, e.kind());
  EXPECT_EQ(arg, e.arg());
  EXPECT_DEBUG_DEATH(factory.MakeUnary(expr::ADD, arg),
                     "invalid expression kind");
  EXPECT_DEBUG_DEATH(factory.MakeUnary(expr::ABS, mp::NumericExpr()),
                     "invalid argument");
}

TEST(ExprTest, BinaryExpr) {
  mp::BinaryExpr e;
  EXPECT_TRUE(e == 0);
  ExprFactory factory;
  auto lhs = factory.MakeNumericConstant(42);
  auto rhs = factory.MakeVariable(0);
  e = factory.MakeBinary(expr::MUL, lhs, rhs);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::MUL, e.kind());
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_DEBUG_DEATH(factory.MakeBinary(expr::IF, lhs, rhs),
                     "invalid expression kind");
  EXPECT_DEBUG_DEATH(factory.MakeBinary(expr::MUL, mp::NumericExpr(), rhs),
                     "invalid argument");
  EXPECT_DEBUG_DEATH(factory.MakeBinary(expr::MUL, lhs, mp::NumericExpr()),
                     "invalid argument");
}

TEST(ExprTest, IfExpr) {
  mp::IfExpr e;
  EXPECT_TRUE(e == 0);
  ExprFactory factory;
  auto condition = factory.MakeLogicalConstant(true);
  auto true_expr = factory.MakeNumericConstant(42);
  auto false_expr = factory.MakeVariable(0);
  e = factory.MakeIf(condition, true_expr, false_expr);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::IF, e.kind());
  EXPECT_EQ(condition, e.condition());
  EXPECT_EQ(true_expr, e.true_expr());
  EXPECT_EQ(false_expr, e.false_expr());
  EXPECT_DEBUG_DEATH(factory.MakeIf(mp::LogicalExpr(), true_expr, false_expr),
                     "invalid argument");
  EXPECT_DEBUG_DEATH(factory.MakeIf(condition, mp::NumericExpr(), false_expr),
                     "invalid argument");
  factory.MakeIf(condition, true_expr, mp::NumericExpr());
}

TEST(ExprTest, PLTerm) {
  mp::PLTerm e;
  EXPECT_TRUE(e == 0);
  ExprFactory factory;
  ExprFactory::PLTermBuilder builder = factory.BeginPLTerm(2);
  builder.AddSlope(11);
  builder.AddBreakpoint(111);
  builder.AddSlope(22);
  builder.AddBreakpoint(222);
  builder.AddSlope(33);
  e = factory.EndPLTerm(builder, factory.MakeVariable(42));
  EXPECT_EQ(expr::PLTERM, e.kind());
  EXPECT_EQ(2, e.num_breakpoints());
  EXPECT_EQ(3, e.num_slopes());
  EXPECT_EQ(11, e.slope(0));
  EXPECT_EQ(22, e.slope(1));
  EXPECT_EQ(33, e.slope(2));
  EXPECT_DEBUG_DEATH(e.slope(-1), "index out of bounds");
  EXPECT_DEBUG_DEATH(e.slope(3), "index out of bounds");
  EXPECT_EQ(111, e.breakpoint(0));
  EXPECT_EQ(222, e.breakpoint(1));
  EXPECT_DEBUG_DEATH(e.breakpoint(-1), "index out of bounds");
  EXPECT_DEBUG_DEATH(e.breakpoint(2), "index out of bounds");
  EXPECT_EQ(42, e.var_index());
  EXPECT_DEBUG_DEATH(factory.BeginPLTerm(0), "invalid number of breakpoints");
}

#ifndef NDEBUG

TEST(ExprTest, TooManyBreakpoints) {
  ExprFactory factory;
  auto builder = factory.BeginPLTerm(1);
  builder.AddBreakpoint(0);
  EXPECT_DEBUG_DEATH(builder.AddBreakpoint(1), "too many breakpoints");
}

TEST(ExprTest, TooManySlopes) {
  ExprFactory factory;
  auto builder = factory.BeginPLTerm(1);
  builder.AddSlope(0);
  builder.AddSlope(1);
  EXPECT_DEBUG_DEATH(builder.AddSlope(2), "too many slopes");
}

TEST(ExprTest, InvalidPLTermArgument) {
  ExprFactory factory;
  auto builder = factory.BeginPLTerm(1);
  builder.AddSlope(0);
  builder.AddBreakpoint(0);
  builder.AddSlope(1);
  EXPECT_DEBUG_DEATH(factory.EndPLTerm(builder, mp::Variable()),
                     "invalid argument");
}

#endif

TEST(ExprTest, TooFewBreakpoints) {
  ExprFactory factory;
  auto builder = factory.BeginPLTerm(1);
  builder.AddSlope(0);
  builder.AddSlope(1);
  EXPECT_DEBUG_DEATH(factory.EndPLTerm(builder, factory.MakeVariable(0)),
                     "too few breakpoints");
}

TEST(ExprTest, TooFewSlopes) {
  ExprFactory factory;
  auto builder = factory.BeginPLTerm(1);
  builder.AddBreakpoint(0);
  builder.AddSlope(0);
  EXPECT_DEBUG_DEATH(factory.EndPLTerm(builder, factory.MakeVariable(0)),
                     "too few slopes");
}

TEST(ExprTest, Function) {
  mp::Function f;
  EXPECT_TRUE(f == 0);
}

TEST(ExprTest, CallExpr) {
  mp::CallExpr e;
  EXPECT_TRUE(e == 0);
  ExprFactory factory;
  mp::Function f = factory.AddFunction("foo");
  ExprFactory::CallExprBuilder builder = factory.BeginCall(f, 3);
  mp::Expr args[] = {
    factory.MakeNumericConstant(11),
    factory.MakeVariable(0),
    factory.MakeNumericConstant(22)
  };
  for (int i = 0; i < 3; ++i)
    builder.AddArg(args[i]);
  e = factory.EndCall(builder);
  EXPECT_EQ(expr::CALL, e.kind());
  EXPECT_EQ(3, e.num_args());
  mp::CallExpr::iterator it = e.begin();
  for (int i = 0; i < 3; ++i, ++it) {
    EXPECT_EQ(args[i], e.arg(i));
    EXPECT_EQ(args[i], *it);
  }
  EXPECT_EQ(e.end(), it);
  EXPECT_DEBUG_DEATH(e.arg(-1), "index out of bounds");
  EXPECT_DEBUG_DEATH(e.arg(3), "index out of bounds");
  EXPECT_DEBUG_DEATH(factory.BeginCall(f, -1), "invalid number of arguments");
  factory.BeginCall(f, 0);
  EXPECT_DEBUG_DEATH(factory.BeginCall(mp::Function(), 0), "invalid function");
}

TEST(ExprTest, CallExprIterator) {
  ExprFactory factory;
  mp::Function f = factory.AddFunction("foo");
  ExprFactory::CallExprBuilder builder = factory.BeginCall(f, 3);
  mp::Expr args[] = {
    factory.MakeNumericConstant(11),
    factory.MakeVariable(0),
    factory.MakeNumericConstant(22)
  };
  for (int i = 0; i < 3; ++i)
    builder.AddArg(args[i]);
  auto e = factory.EndCall(builder);
  mp::CallExpr::iterator i = e.begin();
  EXPECT_EQ(args[0], *i);
  EXPECT_EQ(expr::CONSTANT, i->kind());
  EXPECT_EQ(mp::CallExpr::iterator(e.begin()), i);
  auto j = i;
  EXPECT_TRUE(i == j);
  j = i++;
  EXPECT_TRUE(i != j);
  EXPECT_EQ(args[0], *j);
  EXPECT_EQ(args[1], *i);
  j = ++i;
  EXPECT_EQ(j, i);
  EXPECT_EQ(args[2], *i);
}

TEST(ExprTest, VarArgExpr) {
  mp::VarArgExpr e;
  EXPECT_TRUE(e == 0);
  ExprFactory factory;
  ExprFactory::VarArgExprBuilder builder = factory.BeginVarArg(expr::MAX, 3);
  mp::NumericExpr args[] = {
    factory.MakeNumericConstant(11),
    factory.MakeVariable(0),
    factory.MakeNumericConstant(22)
  };
  for (int i = 0; i < 3; ++i)
    builder.AddArg(args[i]);
  e = factory.EndVarArg(builder);
  EXPECT_EQ(expr::MAX, e.kind());
  EXPECT_EQ(3, e.num_args());
  mp::VarArgExpr::iterator it = e.begin();
  for (int i = 0; i < 3; ++i, ++it) {
    mp::NumericExpr arg = e.arg(i);
    EXPECT_EQ(args[i], arg);
    EXPECT_EQ(args[i], *it);
  }
  EXPECT_EQ(e.end(), it);
  EXPECT_DEBUG_DEATH(e.arg(-1), "index out of bounds");
  EXPECT_DEBUG_DEATH(e.arg(3), "index out of bounds");
  EXPECT_DEBUG_DEATH(factory.BeginVarArg(expr::MAX, 0),
                     "invalid number of arguments");
  factory.BeginVarArg(expr::MIN, 1);
  EXPECT_DEBUG_DEATH(factory.BeginVarArg(expr::SUM, 1),
                     "invalid expression kind");
}

TEST(ExprTest, SumExpr) {
  mp::SumExpr e;
  EXPECT_TRUE(e == 0);
  ExprFactory factory;
  ExprFactory::SumExprBuilder builder = factory.BeginSum(3);
  mp::NumericExpr args[] = {
    factory.MakeNumericConstant(11),
    factory.MakeVariable(0),
    factory.MakeNumericConstant(22)
  };
  for (int i = 0; i < 3; ++i)
    builder.AddArg(args[i]);
  e = factory.EndSum(builder);
  EXPECT_EQ(expr::SUM, e.kind());
  EXPECT_EQ(3, e.num_args());
  mp::SumExpr::iterator it = e.begin();
  for (int i = 0; i < 3; ++i, ++it) {
    mp::NumericExpr arg = e.arg(i);
    EXPECT_EQ(args[i], arg);
    EXPECT_EQ(args[i], *it);
  }
  EXPECT_EQ(e.end(), it);
  EXPECT_DEBUG_DEATH(e.arg(-1), "index out of bounds");
  EXPECT_DEBUG_DEATH(e.arg(3), "index out of bounds");
  EXPECT_DEBUG_DEATH(factory.BeginSum(-1), "invalid number of arguments");
  factory.BeginSum(0);
}

// TODO

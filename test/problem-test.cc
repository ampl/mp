/*
 Problem tests

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

#include "gtest/gtest.h"
#include "test-assert.h"

// Define MP_MAX_PROBLEM_ITEMS to a small value before including problem.h
// to test size overflow checks.
#define MP_MAX_PROBLEM_ITEMS 100

#include "mp/nl.h"
#include "mp/problem.h"

TEST(ProblemTest, ReadNL) {
  mp::Problem p;
  // Just make sure that ReadNLFile and ReadNLString are accessible.
  // They are tested elsewhere.
  EXPECT_THROW(ReadNLFile("nonexistent", p), fmt::SystemError);
  EXPECT_THROW(ReadNLString("", p), mp::Error);
}

// Include the source file to test the implementation.
#include "../src/problem.cc"

using mp::Problem;

#define EXPECT_LINEAR_EXPR(expr, indices, coefs) { \
  int num_terms = sizeof(indices) / sizeof(int); \
  EXPECT_EQ(num_terms, expr.num_terms()); \
  mp::LinearExpr::iterator it = expr.begin(); \
  for (int i = 0; i < num_terms; ++i, ++it) { \
    EXPECT_EQ(indices[i], it->var_index()); \
    EXPECT_EQ(coefs[i], it->coef()); \
  } \
  EXPECT_EQ(expr.end(),it); \
}

TEST(ProblemTest, LinearExpr) {
  mp::LinearExpr e;
  EXPECT_EQ(0, e.num_terms());
  EXPECT_GE(0, e.capacity());
  EXPECT_EQ(e.begin(), e.end());
  e.AddTerm(11, 2.2);
  EXPECT_EQ(1, e.num_terms());
  EXPECT_GE(1, e.capacity());
  auto i = e.begin();
  EXPECT_EQ(11, i->var_index());
  EXPECT_EQ(2.2, i->coef());
  e.Reserve(10);
  EXPECT_EQ(10, e.capacity());
}

TEST(ProblemTest, LinearExprIterator) {
  mp::LinearExpr e;
  e.AddTerm(11, 2.2);
  e.AddTerm(33, 4.4);
  EXPECT_EQ(2, e.num_terms());
  mp::LinearExpr::iterator i = e.begin();
  // Test dereference.
  EXPECT_EQ(11, (*i).var_index());
  // Test the arrow operator.
  EXPECT_EQ(11, i->var_index());
  EXPECT_EQ(2.2, i->coef());
  // Test postincrement.
  mp::LinearExpr::iterator j = i++;
  EXPECT_EQ(11, j->var_index());
  EXPECT_EQ(33, i->var_index());
  EXPECT_TRUE(i != j);
  // Test preincrement.
  i = ++j;
  EXPECT_EQ(33, j->var_index());
  EXPECT_EQ(33, i->var_index());
  EXPECT_TRUE(i == j);
  // Test end.
  EXPECT_NE(i, e.end());
  EXPECT_EQ(++i, e.end());
}

TEST(ProblemTest, AddVar) {
  Problem p;
  EXPECT_EQ(0, p.num_vars());
  p.AddVar(1.1, 2.2);
  EXPECT_EQ(1, p.num_vars());
  Problem::Variable var = p.var(0);
  const Problem::Variable cvar = var;
  EXPECT_EQ(0, cvar.index());
  EXPECT_EQ(1.1, cvar.lb());
  EXPECT_EQ(2.2, cvar.ub());
  EXPECT_EQ(mp::var::CONTINUOUS, cvar.type());
  p.AddVar(3.3, 4.4, mp::var::INTEGER);
  EXPECT_EQ(2, p.num_vars());
  var = p.var(1);
  EXPECT_EQ(1, var.index());
  EXPECT_EQ(3.3, var.lb());
  EXPECT_EQ(4.4, var.ub());
  EXPECT_EQ(mp::var::INTEGER, var.type());
}

TEST(ProblemTest, CompareVars) {
  Problem p;
  p.AddVar(0, 0);
  p.AddVar(0, 0);
  EXPECT_TRUE(p.var(0) == p.var(0));
  EXPECT_TRUE(p.var(0) != p.var(1));
  EXPECT_FALSE(p.var(0) != p.var(0));
  EXPECT_FALSE(p.var(0) == p.var(1));
}

TEST(ProblemTest, InvalidVarIndex) {
  Problem p;
  const int num_vars = 3;
  for (int i = 0; i < num_vars; ++i)
    p.AddVar(0, 0);
  EXPECT_ASSERT(p.var(-1), "invalid index");
  EXPECT_ASSERT(p.var(num_vars), "invalid index");
}

TEST(ProblemTest, MaxVars) {
  Problem p;
  for (int i = 0; i < MP_MAX_PROBLEM_ITEMS; ++i)
    p.AddVar(0, 0);
  EXPECT_EQ(MP_MAX_PROBLEM_ITEMS, p.num_vars());
  EXPECT_ASSERT(p.AddVar(0, 0), "too many variables");
}

TEST(ProblemTest, Vars) {
  Problem p;
  p.AddVar(11, 22);
  p.AddVar(33, 44, mp::var::INTEGER);
  Problem::VarRange vars = p.vars();
  Problem::VarRange::iterator i = vars.begin();
  // Test dereference.
  EXPECT_EQ(p.var(0), *i);
  // Test the arrow operator.
  EXPECT_EQ(11, i->lb());
  EXPECT_EQ(mp::var::CONTINUOUS, i->type());
  // Test postincrement.
  Problem::VarRange::iterator j = i++;
  EXPECT_EQ(p.var(0), *j);
  EXPECT_EQ(p.var(1), *i);
  EXPECT_TRUE(i != j);
  // Test preincrement.
  i = ++j;
  EXPECT_EQ(p.var(1), *j);
  EXPECT_EQ(p.var(1), *i);
  EXPECT_TRUE(i == j);
  // Test end.
  EXPECT_NE(i, vars.end());
  EXPECT_EQ(++i, vars.end());
  // Test invalid access.
  EXPECT_ASSERT(i->lb(), "invalid access");
  EXPECT_ASSERT(*i, "invalid access");
}

TEST(ProblemTest, MutVariable) {
  Problem p;
  p.AddVar(0, 1);
  Problem::MutVariable var = p.var(0);
  EXPECT_EQ(0, var.value());
  var.set_value(4.2);
  EXPECT_EQ(4.2, var.value());
  const Problem &cp = p;
  Problem::Variable cvar = cp.var(0);
  cvar = var;
}

TEST(ProblemTest, AddObj) {
  Problem p;
  EXPECT_EQ(0, p.num_objs());

  p.AddObj(mp::obj::MIN);
  EXPECT_EQ(1, p.num_objs());
  Problem::Objective obj = p.obj(0);
  EXPECT_EQ(mp::obj::MIN, obj.type());
  EXPECT_TRUE(!obj.nonlinear_expr());
  EXPECT_EQ(0, obj.linear_expr().num_terms());

  p.AddObj(mp::obj::MAX, mp::NumericExpr());
  EXPECT_EQ(2, p.num_objs());
  obj = p.obj(1);
  EXPECT_EQ(mp::obj::MAX, obj.type());
  EXPECT_TRUE(!obj.nonlinear_expr());
  EXPECT_EQ(0, obj.linear_expr().num_terms());

  auto nl_expr = p.MakeNumericConstant(42);
  Problem::LinearObjBuilder builder = p.AddObj(mp::obj::MIN, nl_expr);
  builder.AddTerm(0, 1.1);
  builder.AddTerm(3, 2.2);
  EXPECT_EQ(3, p.num_objs());
  obj = p.obj(2);
  EXPECT_EQ(mp::obj::MIN, obj.type());
  const int indices[] = {0, 3};
  const double coefs[] = {1.1, 2.2};
  EXPECT_LINEAR_EXPR(obj.linear_expr(), indices, coefs);
}

// Test adding linear after nonlinear objective and then accessing
// the nonlinear part of the linear objective.
TEST(ProblemTest, IncompleteNonlinearObj) {
  Problem p;
  p.AddObj(mp::obj::MIN, p.MakeNumericConstant(42));
  p.AddObj(mp::obj::MIN);
  EXPECT_EQ(mp::NumericExpr(), p.obj(1).nonlinear_expr());
}

TEST(ProblemTest, CompareObjs) {
  Problem p;
  p.AddObj(mp::obj::MIN);
  p.AddObj(mp::obj::MIN);
  EXPECT_TRUE(p.obj(0) == p.obj(0));
  EXPECT_TRUE(p.obj(0) != p.obj(1));
  EXPECT_FALSE(p.obj(0) != p.obj(0));
  EXPECT_FALSE(p.obj(0) == p.obj(1));
}

TEST(ProblemTest, InvalidObjIndex) {
  Problem p;
  const int num_objs = 3;
  for (int i = 0; i < num_objs; ++i)
    p.AddObj(mp::obj::MIN);
  EXPECT_ASSERT(p.obj(-1), "invalid index");
  EXPECT_ASSERT(p.obj(num_objs), "invalid index");
}

TEST(ProblemTest, MaxObjs) {
  Problem p;
  for (int i = 0; i < MP_MAX_PROBLEM_ITEMS; ++i)
    p.AddObj(mp::obj::MIN);
  EXPECT_EQ(MP_MAX_PROBLEM_ITEMS, p.num_objs());
  EXPECT_ASSERT(p.AddObj(mp::obj::MIN), "too many objectives");
}

TEST(ProblemTest, Objs) {
  Problem p;
  p.AddObj(mp::obj::MAX);
  p.AddObj(mp::obj::MIN);
  Problem::ObjRange objs = p.objs();
  Problem::ObjRange::iterator i = objs.begin();
  // Test dereference.
  EXPECT_EQ(p.obj(0), *i);
  // Test the arrow operator.
  EXPECT_EQ(mp::obj::MAX, i->type());
  // Test postincrement.
  Problem::ObjRange::iterator j = i++;
  EXPECT_EQ(p.obj(0), *j);
  EXPECT_EQ(p.obj(1), *i);
  EXPECT_TRUE(i != j);
  // Test preincrement.
  i = ++j;
  EXPECT_EQ(p.obj(1), *j);
  EXPECT_EQ(p.obj(1), *i);
  EXPECT_TRUE(i == j);
  // Test end.
  EXPECT_NE(i, objs.end());
  EXPECT_EQ(++i, objs.end());
  // Test invalid access.
  EXPECT_ASSERT(i->type(), "invalid access");
  EXPECT_ASSERT(*i, "invalid access");
}

TEST(ProblemTest, MutObjective) {
  Problem p;
  p.AddObj(mp::obj::MIN);
  Problem::MutObjective obj = p.obj(0);
  auto expr = p.MakeNumericConstant(42);
  obj.set_nonlinear_expr(expr);
  EXPECT_EQ(expr, obj.nonlinear_expr());
  obj.linear_expr().AddTerm(11, 2.2);
  const int indices[] = {11};
  const double coefs[] = {2.2};
  EXPECT_LINEAR_EXPR(obj.linear_expr(), indices, coefs);
  EXPECT_EQ(1, obj.linear_expr().num_terms());
  const Problem &cp = p;
  Problem::Objective cobj = cp.obj(0);
  cobj = obj;
}

TEST(ProblemTest, AddAlgebraicCon) {
  Problem p;
  EXPECT_EQ(0, p.num_algebraic_cons());

  p.AddCon(1.1, 2.2);
  EXPECT_EQ(1, p.num_algebraic_cons());
  Problem::AlgebraicCon con = p.algebraic_con(0);
  const Problem::AlgebraicCon ccon = con;
  EXPECT_EQ(1.1, ccon.lb());
  EXPECT_EQ(2.2, ccon.ub());
  EXPECT_TRUE(!ccon.nonlinear_expr());
  EXPECT_EQ(0, ccon.linear_expr().num_terms());

  p.AddCon(3.3, 4.4);
  EXPECT_EQ(2, p.num_algebraic_cons());
  con = p.algebraic_con(1);
  EXPECT_EQ(3.3, con.lb());
  EXPECT_EQ(4.4, con.ub());
  EXPECT_TRUE(!con.nonlinear_expr());
  EXPECT_EQ(0, con.linear_expr().num_terms());

  auto nl_expr = p.MakeNumericConstant(42);
  Problem::MutAlgebraicCon con1 = p.AddCon(5.5, 6.6);
  con1.set_nonlinear_expr(nl_expr);
  Problem::LinearConBuilder builder = con1.set_linear_expr(2);
  builder.AddTerm(0, 1.1);
  builder.AddTerm(3, 2.2);
  EXPECT_EQ(3, p.num_algebraic_cons());
  con = p.algebraic_con(2);
  EXPECT_EQ(5.5, con.lb());
  EXPECT_EQ(6.6, con.ub());
  EXPECT_EQ(nl_expr, con.nonlinear_expr());
  const int indices[] = {0, 3};
  const double coefs[] = {1.1, 2.2};
  EXPECT_LINEAR_EXPR(con.linear_expr(), indices, coefs);
}

// Test adding linear after nonlinear constraint and then accessing
// the nonlinear part of the linear constraint.
TEST(ProblemTest, IncompleteNonlinearCon) {
  Problem p;
  p.AddCon(0, 1).set_nonlinear_expr(p.MakeNumericConstant(42));
  p.AddCon(0, 1);
  EXPECT_EQ(mp::NumericExpr(), p.algebraic_con(1).nonlinear_expr());
}

TEST(ProblemTest, CompareAlgebraicCons) {
  Problem p;
  p.AddCon(1, 2);
  p.AddCon(1, 2);
  EXPECT_TRUE(p.algebraic_con(0) == p.algebraic_con(0));
  EXPECT_TRUE(p.algebraic_con(0) != p.algebraic_con(1));
  EXPECT_FALSE(p.algebraic_con(0) != p.algebraic_con(0));
  EXPECT_FALSE(p.algebraic_con(0) == p.algebraic_con(1));
}

TEST(ProblemTest, InvalidAlgebraicConIndex) {
  Problem p;
  const int num_cons = 3;
  for (int i = 0; i < num_cons; ++i)
    p.AddCon(1, 2);
  EXPECT_ASSERT(p.algebraic_con(-1), "invalid index");
  EXPECT_ASSERT(p.algebraic_con(num_cons), "invalid index");
}

TEST(ProblemTest, MaxAlgebraicCons) {
  Problem p;
  for (int i = 0; i < MP_MAX_PROBLEM_ITEMS; ++i)
    p.AddCon(1, 2);
  EXPECT_EQ(MP_MAX_PROBLEM_ITEMS, p.num_algebraic_cons());
  EXPECT_ASSERT(p.AddCon(1, 2), "too many algebraic constraints");
}

TEST(ProblemTest, AlgebraicCons) {
  Problem p;
  p.AddCon(1, 2);
  p.AddCon(3, 4);
  Problem::AlgebraicConRange cons = p.algebraic_cons();
  Problem::AlgebraicConRange::iterator i = cons.begin();
  // Test dereference.
  EXPECT_EQ(p.algebraic_con(0), *i);
  // Test the arrow operator.
  EXPECT_EQ(1, i->lb());
  // Test postincrement.
  Problem::AlgebraicConRange::iterator j = i++;
  EXPECT_EQ(p.algebraic_con(0), *j);
  EXPECT_EQ(p.algebraic_con(1), *i);
  EXPECT_TRUE(i != j);
  // Test preincrement.
  i = ++j;
  EXPECT_EQ(p.algebraic_con(1), *j);
  EXPECT_EQ(p.algebraic_con(1), *i);
  EXPECT_TRUE(i == j);
  // Test end.
  EXPECT_NE(i, cons.end());
  EXPECT_EQ(++i, cons.end());
  // Test invalid access.
  EXPECT_ASSERT(i->lb(), "invalid access");
  EXPECT_ASSERT(*i, "invalid access");
}

TEST(ProblemTest, MutAlgebraicCon) {
  Problem p;
  p.AddCon(0, 0);
  EXPECT_EQ(1, p.num_algebraic_cons());
  Problem::MutAlgebraicCon con = p.algebraic_con(0);
  EXPECT_EQ(0, con.dual());
  con.set_dual(4.2);
  EXPECT_EQ(4.2 , con.dual());
  auto expr = p.MakeNumericConstant(42);
  con.set_nonlinear_expr(expr);
  EXPECT_EQ(expr, con.nonlinear_expr());
  con.linear_expr().AddTerm(11, 2.2);
  const int indices[] = {11};
  const double coefs[] = {2.2};
  con.set_lb(3.3);
  con.set_ub(4.4);
  EXPECT_EQ(3.3, con.lb());
  EXPECT_EQ(4.4, con.ub());
  EXPECT_LINEAR_EXPR(con.linear_expr(), indices, coefs);
  EXPECT_EQ(1, con.linear_expr().num_terms());
  const Problem &cp = p;
  Problem::AlgebraicCon ccon = cp.algebraic_con(0);
  ccon = con;
}

TEST(ProblemTest, AddLogicalCon) {
  Problem p;
  EXPECT_EQ(0, p.num_logical_cons());

  mp::LogicalExpr expr = p.MakeLogicalConstant(true);
  p.AddCon(expr);
  EXPECT_EQ(1, p.num_logical_cons());
  Problem::LogicalCon con = p.logical_con(0);
  const Problem::LogicalCon ccon = con;
  EXPECT_EQ(expr, ccon.expr());

  mp::LogicalExpr expr2 = p.MakeNot(expr);
  p.AddCon(expr2);
  EXPECT_EQ(2, p.num_logical_cons());
  con = p.logical_con(1);
  EXPECT_EQ(expr2, con.expr());
}

TEST(ProblemTest, CompareLogicalCons) {
  Problem p;
  p.AddCon(p.MakeLogicalConstant(true));
  p.AddCon(p.MakeLogicalConstant(true));
  EXPECT_TRUE(p.logical_con(0) == p.logical_con(0));
  EXPECT_TRUE(p.logical_con(0) != p.logical_con(1));
  EXPECT_FALSE(p.logical_con(0) != p.logical_con(0));
  EXPECT_FALSE(p.logical_con(0) == p.logical_con(1));
}

TEST(ProblemTest, InvalidLogicalConIndex) {
  Problem p;
  const int num_cons = 3;
  for (int i = 0; i < num_cons; ++i)
    p.AddCon(p.MakeLogicalConstant(true));
  EXPECT_ASSERT(p.logical_con(-1), "invalid index");
  EXPECT_ASSERT(p.logical_con(num_cons), "invalid index");
}

TEST(ProblemTest, MaxLogicalCons) {
  Problem p;
  for (int i = 0; i < MP_MAX_PROBLEM_ITEMS; ++i)
    p.AddCon(p.MakeLogicalConstant(true));
  EXPECT_EQ(MP_MAX_PROBLEM_ITEMS, p.num_logical_cons());
  EXPECT_ASSERT(p.AddCon(p.MakeLogicalConstant(true)),
                "too many logical constraints");
}

TEST(ProblemTest, LogicalCons) {
  Problem p;
  mp::LogicalExpr expr = p.MakeLogicalConstant(false);
  p.AddCon(expr);
  p.AddCon(p.MakeLogicalConstant(true));
  Problem::LogicalConRange cons = p.logical_cons();
  Problem::LogicalConRange::iterator i = cons.begin();
  // Test dereference.
  EXPECT_EQ(p.logical_con(0), *i);
  // Test the arrow operator.
  EXPECT_EQ(expr, i->expr());
  // Test postincrement.
  Problem::LogicalConRange::iterator j = i++;
  EXPECT_EQ(p.logical_con(0), *j);
  EXPECT_EQ(p.logical_con(1), *i);
  EXPECT_TRUE(i != j);
  // Test preincrement.
  i = ++j;
  EXPECT_EQ(p.logical_con(1), *j);
  EXPECT_EQ(p.logical_con(1), *i);
  EXPECT_TRUE(i == j);
  // Test end.
  EXPECT_NE(i, cons.end());
  EXPECT_EQ(++i, cons.end());
  // Test invalid access.
  EXPECT_ASSERT(i->expr(), "invalid access");
  EXPECT_ASSERT(*i, "invalid access");
}

TEST(ProblemTest, HasNonlinearCons) {
  Problem p;
  p.AddVar(0, 1);
  EXPECT_FALSE(p.has_nonlinear_cons());
  p.AddCon(0, 0).set_linear_expr(0).AddTerm(0, 1);
  EXPECT_FALSE(p.has_nonlinear_cons());
  p.AddCon(0, 0).set_nonlinear_expr(p.MakeNumericConstant(42));
  EXPECT_TRUE(p.has_nonlinear_cons());
}

TEST(ProblemTest, AddCommonExpr) {
  Problem p;
  EXPECT_EQ(0, p.num_common_exprs());
  Problem::LinearExprBuilder builder =
      p.AddCommonExpr(p.MakeNumericConstant(42)).set_linear_expr(2);
  builder.AddTerm(0, 1.1);
  builder.AddTerm(3, 2.2);
  EXPECT_EQ(1, p.num_common_exprs());
  auto expr = p.common_expr(0);
  const int indices[] = {0, 3};
  const double coefs[] = {1.1, 2.2};
  EXPECT_LINEAR_EXPR(expr.linear_expr(), indices, coefs);
  p.AddCommonExpr(mp::NumericExpr());
  EXPECT_EQ(2, p.num_common_exprs());
  expr = p.common_expr(1);
  EXPECT_TRUE(!expr.nonlinear_expr());
  EXPECT_EQ(0, expr.linear_expr().num_terms());
}

// Test adding linear after nonlinear expression and then accessing
// the nonlinear part of the linear expression.
TEST(ProblemTest, IncompleteNonlinearCommonExpr) {
  Problem p;
  p.AddCommonExpr(p.MakeNumericConstant(42));
  p.AddCommonExpr(mp::NumericExpr());
  EXPECT_EQ(mp::NumericExpr(), p.common_expr(1).nonlinear_expr());
}

TEST(ProblemTest, CompareCommonExprs) {
  Problem p;
  p.AddCommonExpr(mp::NumericExpr());
  p.AddCommonExpr(mp::NumericExpr());
  EXPECT_TRUE(p.common_expr(0) == p.common_expr(0));
  EXPECT_TRUE(p.common_expr(0) != p.common_expr(1));
  EXPECT_FALSE(p.common_expr(0) != p.common_expr(0));
  EXPECT_FALSE(p.common_expr(0) == p.common_expr(1));
}

TEST(ProblemTest, InvalidCommonExprIndex) {
  Problem p;
  const int num_exprs = 3;
  for (int i = 0; i < num_exprs; ++i)
    p.AddCommonExpr(mp::NumericExpr());
  EXPECT_ASSERT(p.common_expr(-1), "invalid index");
  EXPECT_ASSERT(p.common_expr(num_exprs), "invalid index");
}

TEST(ProblemTest, MaxCommonExprs) {
  Problem p;
  for (int i = 0; i < MP_MAX_PROBLEM_ITEMS; ++i)
    p.AddCommonExpr(mp::NumericExpr());
  EXPECT_EQ(MP_MAX_PROBLEM_ITEMS, p.num_common_exprs());
  EXPECT_ASSERT(p.AddCommonExpr(mp::NumericExpr()), "too many expressions");
}

TEST(ProblemTest, Complementarity) {
  Problem p;
  EXPECT_FALSE(p.HasComplementarity());
  p.AddVar(0, 1);
  p.AddCon(0, 1);
  using mp::ComplInfo;
  p.SetComplementarity(0, 0, ComplInfo(0));
  EXPECT_TRUE(p.HasComplementarity());
  EXPECT_ASSERT(p.SetComplementarity(-1, 0, ComplInfo(0)), "invalid index");
  EXPECT_ASSERT(p.SetComplementarity(1, 0, ComplInfo(0)), "invalid index");
  EXPECT_ASSERT(p.SetComplementarity(0, -1, ComplInfo(0)), "invalid index");
  EXPECT_ASSERT(p.SetComplementarity(0, 1, ComplInfo(0)), "invalid index");
}

TEST(ProblemTest, RangeIteratorHasCategory) {
  Problem::VarRange::iterator::iterator_category();
}

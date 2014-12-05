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

#include "mp/problem.h"

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

TEST(ProblemTest, AddVar) {
  Problem p;
  EXPECT_EQ(0, p.num_vars());
  p.AddVar(1.1, 2.2);
  EXPECT_EQ(1, p.num_vars());
  Problem::Variable var = p.var(0);
  EXPECT_EQ(1.1, var.lb());
  EXPECT_EQ(2.2, var.ub());
  EXPECT_EQ(mp::var::CONTINUOUS, var.type());
  p.AddVar(3.3, 4.4, mp::var::INTEGER);
  EXPECT_EQ(2, p.num_vars());
  var = p.var(1);
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
  Problem::VarList vars = p.vars();
  Problem::VarList::iterator i = vars.begin();
  // Test dereference.
  EXPECT_EQ(p.var(0), *i);
  // Test the arrow operator.
  EXPECT_EQ(11, i->lb());
  EXPECT_EQ(mp::var::CONTINUOUS, i->type());
  // Test postincrement.
  Problem::VarList::iterator j = i++;
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
  Problem::ObjList objs = p.objs();
  Problem::ObjList::iterator i = objs.begin();
  // Test dereference.
  EXPECT_EQ(p.obj(0), *i);
  // Test the arrow operator.
  EXPECT_EQ(mp::obj::MAX, i->type());
  // Test postincrement.
  Problem::ObjList::iterator j = i++;
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

TEST(ProblemTest, AddAlgebraicCon) {
  Problem p;
  EXPECT_EQ(0, p.num_algebraic_cons());

  p.AddCon(1.1, 2.2);
  EXPECT_EQ(1, p.num_algebraic_cons());
  Problem::AlgebraicCon con = p.algebraic_con(0);
  EXPECT_EQ(1.1, con.lb());
  EXPECT_EQ(2.2, con.ub());
  EXPECT_TRUE(!con.nonlinear_expr());
  EXPECT_EQ(0, con.linear_expr().num_terms());

  p.AddCon(3.3, 4.4, mp::NumericExpr());
  EXPECT_EQ(2, p.num_algebraic_cons());
  con = p.algebraic_con(1);
  EXPECT_EQ(3.3, con.lb());
  EXPECT_EQ(4.4, con.ub());
  EXPECT_TRUE(!con.nonlinear_expr());
  EXPECT_EQ(0, con.linear_expr().num_terms());

  auto nl_expr = p.MakeNumericConstant(42);
  Problem::LinearObjBuilder builder = p.AddCon(5.5, 6.6, nl_expr);
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

TEST(ProblemTest, CompareCons) {
  Problem p;
  p.AddCon(1, 2);
  p.AddCon(1, 2);
  EXPECT_TRUE(p.algebraic_con(0) == p.algebraic_con(0));
  EXPECT_TRUE(p.algebraic_con(0) != p.algebraic_con(1));
  EXPECT_FALSE(p.algebraic_con(0) != p.algebraic_con(0));
  EXPECT_FALSE(p.algebraic_con(0) == p.algebraic_con(1));
}

TEST(ProblemTest, InvalidConIndex) {
  Problem p;
  const int num_cons = 3;
  for (int i = 0; i < num_cons; ++i)
    p.AddCon(1, 2);
  EXPECT_ASSERT(p.algebraic_con(-1), "invalid index");
  EXPECT_ASSERT(p.algebraic_con(num_cons), "invalid index");
}

TEST(ProblemTest, MaxCons) {
  Problem p;
  for (int i = 0; i < MP_MAX_PROBLEM_ITEMS; ++i)
    p.AddCon(1, 2);
  EXPECT_EQ(MP_MAX_PROBLEM_ITEMS, p.num_algebraic_cons());
  EXPECT_ASSERT(p.AddCon(1, 2), "too many algebraic constraints");
}

TEST(ProblemTest, Cons) {
  Problem p;
  p.AddCon(1, 2);
  p.AddCon(3, 4);
  Problem::AlgebraicConList cons = p.algebraic_cons();
  Problem::AlgebraicConList::iterator i = cons.begin();
  // Test dereference.
  EXPECT_EQ(p.algebraic_con(0), *i);
  // Test the arrow operator.
  EXPECT_EQ(1, i->lb());
  // Test postincrement.
  Problem::AlgebraicConList::iterator j = i++;
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

TEST(ProblemTest, SetInfo) {
  Problem p;
  auto info = mp::ProblemInfo();
  info.num_vars = 1;
  p.SetInfo(info);
  p.AddVar(0, 0);
  // TODO: test that there is no reallocation
}

// TODO: check the default definition of MP_MAX_PROBLEM_ITEMS
// TODO: more tests

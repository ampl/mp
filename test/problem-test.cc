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

#define MP_MAX_INDEX 100
#include "mp/problem.h"

TEST(ProblemTest, Ctor) {
  mp::Problem p;
  EXPECT_EQ(0, p.num_vars());
  EXPECT_EQ(0, p.num_objs());
  EXPECT_EQ(0, p.num_algebraic_cons());
  EXPECT_EQ(0, p.num_logical_cons());
}

TEST(ProblemTest, AddVar) {
  mp::Problem p;
  p.AddVar(11, 22);
  EXPECT_EQ(1, p.num_vars());
  EXPECT_EQ(11, p.var_lb(0));
  EXPECT_EQ(22, p.var_ub(0));
  EXPECT_EQ(mp::var::CONTINUOUS, p.var_type(0));
  p.AddVar(33, 44, mp::var::INTEGER);
  EXPECT_EQ(2, p.num_vars());
  EXPECT_EQ(33, p.var_lb(1));
  EXPECT_EQ(44, p.var_ub(1));
  EXPECT_EQ(mp::var::INTEGER, p.var_type(1));
}

TEST(ProblemTest, MaxVarIndex) {
  mp::Problem p;
  for (int i = 0; i < MP_MAX_INDEX; ++i)
    p.AddVar(0, 0);
  EXPECT_ASSERT(p.AddVar(0, 0), "too many variables");
}

TEST(ProblemTest, InvalidVarIndex) {
  mp::Problem p;
  const int num_vars = 3;
  for (int i = 0; i < num_vars; ++i)
    p.AddVar(0, 0);
  EXPECT_ASSERT(p.var_lb(-1), "invalid index");
  EXPECT_ASSERT(p.var_lb(num_vars), "invalid index");
  EXPECT_ASSERT(p.var_ub(-1), "invalid index");
  EXPECT_ASSERT(p.var_ub(num_vars), "invalid index");
  EXPECT_ASSERT(p.var_type(-1), "invalid index");
  EXPECT_ASSERT(p.var_type(num_vars), "invalid index");
}

// TODO: more tests

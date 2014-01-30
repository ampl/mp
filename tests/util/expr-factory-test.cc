/*
 Tests of the the AMPL expression factory.

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
#include "solvers/util/expr-factory.h"
#include "solvers/util/nl.h"

using ampl::ExprFactory;
using ampl::NLHeader;

namespace {

TEST(ExprFactoryTest, Ctor) {
  ExprFactory ef((NLHeader()));
}

TEST(ExprFactoryTest, CreateNumericConstant) {
  ExprFactory ef((NLHeader()));
  EXPECT_EQ(42.0, ef.CreateNumericConstant(42).value());
}

TEST(ExprFactoryTest, CreateVariable) {
  NLHeader header = {};
  header.num_vars = 10;
  ExprFactory ef(header);
  EXPECT_EQ(0, ef.CreateVariable(0).index());
  EXPECT_EQ(9, ef.CreateVariable(9).index());
  EXPECT_DEBUG_DEATH(ef.CreateVariable(-1);, "Assertion");  // NOLINT(*)
  EXPECT_DEBUG_DEATH(ef.CreateVariable(10);, "Assertion");  // NOLINT(*)
}

TEST(ExprFactoryTest, CompatibleWithFGRead) {
  NLHeader header = {};
  ExprFactory ef(header);
  ASL *asl = ASL_alloc(ASL_read_fg);
  // TODO: read asl and compare to the one produced by fg_read
  ASL_free(&asl);
}

// TODO: check if the asl produced by ExprFactory is binary compatible with
//       the one produced by reading an .nl file with fg_read.

// TODO: more tests
}

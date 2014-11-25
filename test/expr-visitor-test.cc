/*
 Expression visitor tests

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
#include "mp/expr-visitor.h"

using mp::internal::ExprTypes;

// IsSame<T, U>::VALUE is true iff T and U are the same type.
template <typename T, typename U>
struct IsSame {
  static const bool VALUE = false;
};

template <typename T>
struct IsSame<T, T> {
  static const bool VALUE = true;
};

// Checks if ExprTypes::T is a typedef of mp::T.
#define CHECK_EXPR_TYPE(T) \
  EXPECT_TRUE((IsSame<ExprTypes::T, mp::T>::VALUE))

TEST(ExprVisitorTest, ExprTypes) {
  CHECK_EXPR_TYPE(NumericExpr);
  CHECK_EXPR_TYPE(LogicalExpr);
  CHECK_EXPR_TYPE(NumericConstant);
  CHECK_EXPR_TYPE(Variable);
  CHECK_EXPR_TYPE(UnaryExpr);
  CHECK_EXPR_TYPE(BinaryExpr);
  CHECK_EXPR_TYPE(IfExpr);
  CHECK_EXPR_TYPE(PLTerm);
  CHECK_EXPR_TYPE(CallExpr);
  EXPECT_TRUE((IsSame<ExprTypes::VarArgExpr, mp::IteratedExpr>::VALUE));
  EXPECT_TRUE((IsSame<ExprTypes::SumExpr, mp::IteratedExpr>::VALUE));
  EXPECT_TRUE((IsSame<ExprTypes::NumberOfExpr, mp::IteratedExpr>::VALUE));
  CHECK_EXPR_TYPE(CountExpr);
  CHECK_EXPR_TYPE(LogicalConstant);
  CHECK_EXPR_TYPE(NotExpr);
  CHECK_EXPR_TYPE(BinaryLogicalExpr);
  CHECK_EXPR_TYPE(RelationalExpr);
  CHECK_EXPR_TYPE(LogicalCountExpr);
  CHECK_EXPR_TYPE(ImplicationExpr);
  CHECK_EXPR_TYPE(IteratedLogicalExpr);
  CHECK_EXPR_TYPE(PairwiseExpr);
}

class TestVisitor : public mp::ExprVisitor<TestVisitor, void, void> {};

// TODO

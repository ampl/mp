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

#include <stdexcept>
#include "gtest-extra.h"
#include "mock-allocator.h"
#include "test-assert.h"
#include "mp/expr.h"

using ::testing::_;
using ::testing::Return;

using mp::ExprFactory;
namespace expr = mp::expr;

class ExprTest : public ::testing::Test {
 protected:
  ExprFactory factory_;
};

TEST_F(ExprTest, Expr) {
  mp::Expr e;
  EXPECT_TRUE(e == 0);
}

TEST_F(ExprTest, NumericExpr) {
  mp::NumericExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::Expr(e);
}

TEST_F(ExprTest, LogicalExpr) {
  mp::LogicalExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::Expr(e);
}

TEST_F(ExprTest, NumericConstant) {
  mp::NumericConstant e;
  EXPECT_TRUE(e == 0);
  (void)mp::NumericExpr(e);
  e = factory_.MakeNumericConstant(1.23);
  EXPECT_EQ(expr::CONSTANT, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(1.23, e.value());
}

TEST_F(ExprTest, Variable) {
  mp::Reference e;
  EXPECT_TRUE(e == 0);
  (void)mp::NumericExpr(e);
  e = factory_.MakeVariable(42);
  EXPECT_EQ(expr::VARIABLE, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(42, e.index());
}

TEST_F(ExprTest, CommonExpr) {
  mp::Reference e;
  EXPECT_TRUE(e == 0);
  (void)mp::NumericExpr(e);
  e = factory_.MakeCommonExpr(42);
  EXPECT_EQ(expr::COMMON_EXPR, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(42, e.index());
}

TEST_F(ExprTest, UnaryExpr) {
  mp::UnaryExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::NumericExpr(e);
  auto arg = factory_.MakeNumericConstant(42);
  e = factory_.MakeUnary(expr::ABS, arg);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::ABS, e.kind());
  EXPECT_EQ(arg, e.arg());
  EXPECT_ASSERT(factory_.MakeUnary(expr::ADD, arg), "invalid expression kind");
  EXPECT_ASSERT(factory_.MakeUnary(expr::ABS, mp::NumericExpr()),
                "invalid argument");
}

TEST_F(ExprTest, BinaryExpr) {
  mp::BinaryExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::NumericExpr(e);
  auto lhs = factory_.MakeNumericConstant(42);
  auto rhs = factory_.MakeVariable(0);
  e = factory_.MakeBinary(expr::MUL, lhs, rhs);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::MUL, e.kind());
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_ASSERT(factory_.MakeBinary(expr::IF, lhs, rhs),
                "invalid expression kind");
  EXPECT_ASSERT(factory_.MakeBinary(expr::MUL, mp::NumericExpr(), rhs),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeBinary(expr::MUL, lhs, mp::NumericExpr()),
                "invalid argument");
}

TEST_F(ExprTest, IfExpr) {
  mp::IfExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::NumericExpr(e);
  auto condition = factory_.MakeLogicalConstant(true);
  auto true_expr = factory_.MakeNumericConstant(42);
  auto false_expr = factory_.MakeVariable(0);
  e = factory_.MakeIf(condition, true_expr, false_expr);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::IF, e.kind());
  EXPECT_EQ(condition, e.condition());
  EXPECT_EQ(true_expr, e.true_expr());
  EXPECT_EQ(false_expr, e.false_expr());
  EXPECT_ASSERT(factory_.MakeIf(mp::LogicalExpr(), true_expr, false_expr),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeIf(condition, mp::NumericExpr(), false_expr),
                "invalid argument");
  factory_.MakeIf(condition, true_expr, mp::NumericExpr());
}

TEST_F(ExprTest, PLTerm) {
  mp::PLTerm e;
  EXPECT_TRUE(e == 0);
  (void)mp::NumericExpr(e);
  ExprFactory::PLTermBuilder builder = factory_.BeginPLTerm(2);
  builder.AddSlope(11);
  builder.AddBreakpoint(111);
  builder.AddSlope(22);
  builder.AddBreakpoint(222);
  builder.AddSlope(33);
  auto var = factory_.MakeVariable(42);
  e = factory_.EndPLTerm(builder, var);
  EXPECT_EQ(expr::PLTERM, e.kind());
  EXPECT_EQ(2, e.num_breakpoints());
  EXPECT_EQ(3, e.num_slopes());
  EXPECT_EQ(11, e.slope(0));
  EXPECT_EQ(22, e.slope(1));
  EXPECT_EQ(33, e.slope(2));
  EXPECT_ASSERT(e.slope(-1), "index out of bounds");
  EXPECT_ASSERT(e.slope(3), "index out of bounds");
  EXPECT_EQ(111, e.breakpoint(0));
  EXPECT_EQ(222, e.breakpoint(1));
  EXPECT_ASSERT(e.breakpoint(-1), "index out of bounds");
  EXPECT_ASSERT(e.breakpoint(2), "index out of bounds");
  EXPECT_EQ(var, e.arg());
  EXPECT_ASSERT(factory_.BeginPLTerm(0), "invalid number of breakpoints");
  builder = factory_.BeginPLTerm(1);
  builder.AddSlope(0);
  builder.AddBreakpoint(0);
  builder.AddSlope(1);
  auto arg = factory_.MakeCommonExpr(0);
  e = factory_.EndPLTerm(builder, arg);
  EXPECT_EQ(arg, e.arg());
}

TEST_F(ExprTest, TooManyBreakpoints) {
  auto builder = factory_.BeginPLTerm(1);
  builder.AddBreakpoint(0);
  EXPECT_ASSERT(builder.AddBreakpoint(1), "too many breakpoints");
}

TEST_F(ExprTest, TooManySlopes) {
  auto builder = factory_.BeginPLTerm(1);
  builder.AddSlope(0);
  builder.AddSlope(1);
  EXPECT_ASSERT(builder.AddSlope(2), "too many slopes");
}

TEST_F(ExprTest, InvalidPLTermArgument) {
  auto builder = factory_.BeginPLTerm(1);
  builder.AddSlope(0);
  builder.AddBreakpoint(0);
  builder.AddSlope(1);
  EXPECT_ASSERT(factory_.EndPLTerm(builder, mp::Reference()),
                "invalid argument");
}

TEST_F(ExprTest, TooFewBreakpoints) {
  auto builder = factory_.BeginPLTerm(1);
  builder.AddSlope(0);
  builder.AddSlope(1);
  EXPECT_ASSERT(factory_.EndPLTerm(builder, factory_.MakeVariable(0)),
                "too few breakpoints");
}

TEST_F(ExprTest, TooFewSlopes) {
  auto builder = factory_.BeginPLTerm(1);
  builder.AddBreakpoint(0);
  builder.AddSlope(0);
  EXPECT_ASSERT(factory_.EndPLTerm(builder, factory_.MakeVariable(0)),
                "too few slopes");
}

TEST_F(ExprTest, Function) {
  mp::Function f;
  EXPECT_TRUE(f == 0);
}

TEST_F(ExprTest, AddFunction) {
  mp::Function f = factory_.AddFunction("foo", 42, mp::func::SYMBOLIC);
  EXPECT_STREQ("foo", f.name());
  EXPECT_EQ(42, f.num_args());
  EXPECT_EQ(mp::func::SYMBOLIC, f.type());
  f = factory_.AddFunction("bar", 0, mp::func::NUMERIC);
  EXPECT_STREQ("bar", f.name());
  EXPECT_EQ(0, f.num_args());
  EXPECT_EQ(mp::func::NUMERIC, f.type());
}

// Iterated expressions share the same builder so it is enough to test
// one, CallExprBuilder.

TEST_F(ExprTest, TooManyCallArgs) {
  mp::Function f = factory_.AddFunction("foo", 1);
  auto builder = factory_.BeginCall(f, 1);
  auto arg = factory_.MakeNumericConstant(0);
  builder.AddArg(arg);
  EXPECT_ASSERT(builder.AddArg(arg), "too many arguments");
}

TEST_F(ExprTest, InvalidCallArg) {
  mp::Function f = factory_.AddFunction("foo", 1);
  auto builder = factory_.BeginCall(f, 1);
  EXPECT_ASSERT(builder.AddArg(mp::NumericExpr()), "invalid argument");
}

TEST_F(ExprTest, TooFewCallArgs) {
  mp::Function f = factory_.AddFunction("foo", 1);
  auto builder = factory_.BeginCall(f, 1);
  EXPECT_ASSERT(factory_.EndCall(builder), "too few arguments");
}

// Expression iterators share the same implementation so it is enough to
// test CallExpr::iterator.
TEST_F(ExprTest, ExprIterator) {
  enum {NUM_ARGS = 3};
  mp::Function f = factory_.AddFunction("foo", NUM_ARGS);
  ExprFactory::CallExprBuilder builder = factory_.BeginCall(f, NUM_ARGS);
  mp::Expr args[NUM_ARGS] = {
    factory_.MakeNumericConstant(11),
    factory_.MakeVariable(0),
    factory_.MakeNumericConstant(22)
  };
  for (int i = 0; i < NUM_ARGS; ++i)
    builder.AddArg(args[i]);
  auto e = factory_.EndCall(builder);
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

TEST_F(ExprTest, LogicalConstant) {
  mp::LogicalConstant e;
  EXPECT_TRUE(e == 0);
  (void)mp::LogicalExpr(e);
  e = factory_.MakeLogicalConstant(false);
  EXPECT_EQ(expr::CONSTANT, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_FALSE(e.value());
  EXPECT_TRUE(factory_.MakeLogicalConstant(true).value());
}

TEST_F(ExprTest, NotExpr) {
  mp::NotExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::LogicalExpr(e);
  auto arg = factory_.MakeLogicalConstant(false);
  e = factory_.MakeNot(arg);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::NOT, e.kind());
  EXPECT_EQ(arg, e.arg());
  EXPECT_ASSERT(factory_.MakeNot(mp::LogicalExpr()), "invalid argument");
}

TEST_F(ExprTest, BinaryLogicalExpr) {
  mp::BinaryLogicalExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::LogicalExpr(e);
  auto lhs = factory_.MakeLogicalConstant(true);
  auto rhs = factory_.MakeLogicalConstant(false);
  e = factory_.MakeBinaryLogical(expr::AND, lhs, rhs);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::AND, e.kind());
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_ASSERT(factory_.MakeBinaryLogical(expr::LT, lhs, rhs),
                "invalid expression kind");
  EXPECT_ASSERT(factory_.MakeBinaryLogical(expr::AND, mp::LogicalExpr(), rhs),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeBinaryLogical(expr::AND, lhs, mp::LogicalExpr()),
                "invalid argument");
}

TEST_F(ExprTest, RelationalExpr) {
  mp::RelationalExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::LogicalExpr(e);
  auto lhs = factory_.MakeNumericConstant(42);
  auto rhs = factory_.MakeVariable(0);
  e = factory_.MakeRelational(expr::EQ, lhs, rhs);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::EQ, e.kind());
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_ASSERT(factory_.MakeRelational(expr::ATLEAST, lhs, rhs),
                "invalid expression kind");
  EXPECT_ASSERT(factory_.MakeRelational(expr::EQ, mp::NumericExpr(), rhs),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeRelational(expr::EQ, lhs, mp::NumericExpr()),
                "invalid argument");
}

TEST_F(ExprTest, LogicalCountExpr) {
  mp::LogicalCountExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::LogicalExpr(e);
  auto lhs = factory_.MakeNumericConstant(42);
  auto rhs = factory_.EndCount(factory_.BeginCount(0));
  e = factory_.MakeLogicalCount(expr::ATMOST, lhs, rhs);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::ATMOST, e.kind());
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_ASSERT(factory_.MakeLogicalCount(expr::IMPLICATION, lhs, rhs),
                "invalid expression kind");
  EXPECT_ASSERT(factory_.MakeLogicalCount(expr::ATMOST, mp::NumericExpr(), rhs),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeLogicalCount(expr::ATMOST, lhs, mp::CountExpr()),
                "invalid argument");
}

TEST_F(ExprTest, ImplicationExpr) {
  mp::ImplicationExpr e;
  (void)mp::LogicalExpr(e);
  EXPECT_TRUE(e == 0);
  auto condition = factory_.MakeLogicalConstant(true);
  auto true_expr = factory_.MakeLogicalConstant(false);
  auto false_expr = factory_.MakeLogicalConstant(true);
  e = factory_.MakeImplication(condition, true_expr, false_expr);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::IMPLICATION, e.kind());
  EXPECT_EQ(condition, e.condition());
  EXPECT_EQ(true_expr, e.true_expr());
  EXPECT_EQ(false_expr, e.false_expr());
  EXPECT_ASSERT(factory_.MakeImplication(
                  mp::LogicalExpr(), true_expr, false_expr),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeImplication(
                  condition, mp::LogicalExpr(), false_expr),
                "invalid argument");
  factory_.MakeImplication(condition, true_expr, mp::LogicalExpr());
}

template <typename ExprInfo>
class IteratedExprTest : public ::testing::Test {
 public:
  ExprInfo info_;
};

template <typename ExprType, expr::Kind K,
          typename BaseType = mp::NumericExpr>
struct ExprInfo {
  typedef BaseType Base;
  typedef ExprType Expr;
  typedef mp::NumericExpr Arg;
  static expr::Kind kind() { return K; }
  static int min_args() { return 0; }
};

struct CallInfo : ExprInfo<mp::CallExpr, expr::CALL> {
  typedef mp::Expr Arg;
  typedef ExprFactory::CallExprBuilder Builder;
  Builder BeginBuild(ExprFactory &f, int n) {
    return f.BeginCall(f.AddFunction("foo", n), n);
  }
  Expr EndBuild(ExprFactory &f, Builder b) { return f.EndCall(b); }
};

struct IteratedInfo : ExprInfo<mp::IteratedExpr, expr::SUM> {
  typedef ExprFactory::IteratedExprBuilder Builder;
  Builder BeginBuild(ExprFactory &f, int n) {
    return f.BeginIterated(expr::SUM, n);
  }
  Expr EndBuild(ExprFactory &f, Builder b) { return f.EndIterated(b); }
};

struct CountInfo : ExprInfo<mp::CountExpr, expr::COUNT> {
  typedef mp::LogicalExpr Arg;
  typedef ExprFactory::CountExprBuilder Builder;
  Builder BeginBuild(ExprFactory &f, int n) { return f.BeginCount(n); }
  Expr EndBuild(ExprFactory &f, Builder b) { return f.EndCount(b); }
};

struct IteratedLogicalInfo :
    ExprInfo<mp::IteratedLogicalExpr, expr::EXISTS, mp::LogicalExpr> {
  typedef mp::LogicalExpr Arg;
  typedef ExprFactory::IteratedLogicalExprBuilder Builder;
  Builder BeginBuild(ExprFactory &f, int n) {
    return f.BeginIteratedLogical(expr::EXISTS, n);
  }
  Expr EndBuild(ExprFactory &f, Builder b) { return f.EndIteratedLogical(b); }
};

struct PairwiseInfo :
    ExprInfo<mp::PairwiseExpr, expr::ALLDIFF, mp::LogicalExpr> {
  typedef ExprFactory::PairwiseExprBuilder Builder;
  Builder BeginBuild(ExprFactory &f, int n) {
    return f.BeginPairwise(expr::ALLDIFF, n);
  }
  Expr EndBuild(ExprFactory &f, Builder b) { return f.EndPairwise(b); }
};

typedef ::testing::Types<CallInfo, IteratedInfo, CountInfo,
                         IteratedLogicalInfo, PairwiseInfo> IteratedExprTypes;
TYPED_TEST_CASE(IteratedExprTest, IteratedExprTypes);

template <typename ExprType>
struct TestArgs {
  mp::NumericExpr args[3];
  explicit TestArgs(ExprFactory &f) {
    args[0] = f.MakeNumericConstant(11);
    args[1] = f.MakeVariable(0);
    args[2] = f.MakeNumericConstant(22);
  }
  mp::NumericExpr operator[](int i) const { return args[i]; }
};

template <>
struct TestArgs<mp::LogicalExpr> {
  mp::LogicalExpr args[3];
  explicit TestArgs(ExprFactory &f) {
    args[0] = f.MakeLogicalConstant(false);
    args[1] = f.MakeLogicalConstant(true);
    args[2] = f.MakeLogicalConstant(false);
  }
  mp::LogicalExpr operator[](int i) const { return args[i]; }
};

TYPED_TEST(IteratedExprTest, Test) {
  typename TypeParam::Expr e;
  EXPECT_TRUE(e == 0);
  (void)typename TypeParam::Base(e);
  ExprFactory factory;
  enum {NUM_ARGS = 3};
  TestArgs<typename TypeParam::Arg> args(factory);
  auto info = this->info_;
  typename TypeParam::Builder builder = info.BeginBuild(factory, NUM_ARGS);
  for (int i = info.min_args(); i < NUM_ARGS; ++i)
    builder.AddArg(args[i]);
  e = info.EndBuild(factory, builder);
  EXPECT_EQ(info.kind(), e.kind());
  EXPECT_EQ(3, e.num_args());
  typename TypeParam::Expr::iterator it = e.begin();
  for (int i = 0; i < NUM_ARGS; ++i, ++it) {
    typename TypeParam::Arg arg = e.arg(i);
    EXPECT_EQ(args[i], arg);
    EXPECT_EQ(args[i], *it);
  }
  EXPECT_EQ(e.end(), it);
  EXPECT_ASSERT(e.arg(-1), "index out of bounds");
  EXPECT_ASSERT(e.arg(NUM_ARGS), "index out of bounds");
  EXPECT_ASSERT(info.BeginBuild(factory, info.min_args() - 1),
                "invalid number of arguments");
  info.BeginBuild(factory, info.min_args());
}

TEST_F(ExprTest, StringLiteral) {
  mp::StringLiteral e;
  EXPECT_TRUE(e == 0);
  (void)mp::Expr(e);
  const char STR[] = "abc\0def";
  e = factory_.MakeStringLiteral(fmt::StringRef(STR, sizeof(STR)));
  EXPECT_EQ(expr::STRING, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(std::string(STR, sizeof(STR)), std::string(e.value(), sizeof(STR)));
}

TEST_F(ExprTest, InternalCast) {
  mp::Expr e = factory_.MakeNumericConstant(42);
  mp::internal::Cast<mp::NumericExpr>(e);
  mp::NumericConstant n = mp::internal::Cast<mp::NumericConstant>(e);
  EXPECT_EQ(42, n.value());
  EXPECT_ASSERT(mp::internal::Cast<mp::UnaryExpr>(e), "invalid cast");
}

TEST_F(ExprTest, Cast) {
  mp::Expr e = factory_.MakeNumericConstant(42);
  EXPECT_EQ(e, mp::Cast<mp::NumericExpr>(e));
  mp::NumericConstant n = mp::Cast<mp::NumericConstant>(e);
  EXPECT_EQ(42, n.value());
  EXPECT_EQ(mp::UnaryExpr(), mp::Cast<mp::UnaryExpr>(e));
}

TEST_F(ExprTest, InvalidCallExprFunction) {
  EXPECT_ASSERT(factory_.BeginCall(mp::Function(), 0), "invalid function");
}

TEST_F(ExprTest, InvalidIteratedExprKind) {
  EXPECT_ASSERT(factory_.BeginIterated(expr::COUNT, 1),
                "invalid expression kind");
}

TEST_F(ExprTest, InvalidNumberOfExprArg) {
  EXPECT_ASSERT(factory_.BeginNumberOf(1, mp::NumericExpr()),
      "invalid argument");
}

TEST_F(ExprTest, InvalidIteratedLogicalExprKind) {
  EXPECT_ASSERT(factory_.BeginIteratedLogical(expr::ALLDIFF, 1),
                "invalid expression kind");
}

TEST(ExprFactoryTest, ExprMemoryAllocation) {
  typedef AllocatorRef< MockAllocator<char> > Allocator;
  MockAllocator<char> alloc;
  mp::BasicExprFactory<Allocator> f((Allocator(&alloc)));
  char buffer[100];
  EXPECT_CALL(alloc, allocate(_)).WillOnce(Return(buffer));
  f.MakeNumericConstant(42);
  EXPECT_CALL(alloc, deallocate(buffer, _));
}

TEST(ExprFactoryTest, FuncMemoryAllocation) {
  typedef AllocatorRef< MockAllocator<char> > Allocator;
  MockAllocator<char> alloc;
  mp::BasicExprFactory<Allocator> f((Allocator(&alloc)));
  char buffer[100];
  EXPECT_CALL(alloc, allocate(_)).WillOnce(Return(buffer));
  f.AddFunction("f", 0);
  EXPECT_CALL(alloc, deallocate(buffer, _));
}

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

using mp::NumericExpr;
using mp::LogicalExpr;
using mp::ExprFactory;

namespace expr = mp::expr;

class ExprTest : public ::testing::Test {
 protected:
  ExprFactory factory_;
  LogicalExpr l0, l1;

  ExprTest() {
    l0 = factory_.MakeLogicalConstant(false);
    l1 = factory_.MakeLogicalConstant(true);
  }

  mp::NumericConstant MakeConst(double value) {
    return factory_.MakeNumericConstant(value);
  }

  mp::Reference MakeVariable(int index) { return factory_.MakeVariable(index); }

  mp::UnaryExpr MakeUnary(expr::Kind kind, NumericExpr arg) {
    return factory_.MakeUnary(kind, arg);
  }
  mp::BinaryExpr MakeBinary(expr::Kind kind, NumericExpr lhs, NumericExpr rhs) {
    return factory_.MakeBinary(kind, lhs, rhs);
  }

  template <int N>
  mp::CallExpr MakeCall(mp::Function f, NumericExpr (&args)[N]) {
    auto builder = factory_.BeginCall(f, N);
    for (int i = 0; i < N; ++i)
      builder.AddArg(args[i]);
    return factory_.EndCall(builder);
  }

  template <int N>
  mp::IteratedExpr MakeIterated(expr::Kind kind, NumericExpr (&args)[N]) {
    auto builder = factory_.BeginIterated(kind, N);
    for (int i = 0; i < N; ++i)
      builder.AddArg(args[i]);
    return factory_.EndIterated(builder);
  }

  template <int N>
  mp::SymbolicNumberOfExpr MakeSymbolicNumberOf(mp::Expr (&args)[N]) {
    auto builder = factory_.BeginSymbolicNumberOf(N, args[0]);
    for (int i = 1; i < N; ++i)
      builder.AddArg(args[i]);
    return factory_.EndSymbolicNumberOf(builder);
  }

  template <int N>
  mp::CountExpr MakeCount(LogicalExpr (&args)[N]) {
    auto builder = factory_.BeginCount(N);
    for (int i = 0; i < N; ++i)
      builder.AddArg(args[i]);
    return factory_.EndCount(builder);
  }

  mp::PLTerm MakePLTerm(int num_breakpoints, const double *breakpoints,
                        const double *slopes, mp::Reference arg) {
    auto builder = factory_.BeginPLTerm(num_breakpoints);
    for (int i = 0; i < num_breakpoints; ++i) {
      builder.AddSlope(slopes[i]);
      builder.AddBreakpoint(breakpoints[i]);
    }
    builder.AddSlope(slopes[num_breakpoints]);
    return factory_.EndPLTerm(builder, arg);
  }

  template <int N>
  mp::IteratedLogicalExpr MakeIteratedLogical(
      expr::Kind kind, LogicalExpr (&args)[N]) {
    auto builder = factory_.BeginIteratedLogical(kind, N);
    for (int i = 0; i < N; ++i)
      builder.AddArg(args[i]);
    return factory_.EndIteratedLogical(builder);
  }

  template <int N>
  mp::PairwiseExpr MakePairwise(expr::Kind kind, NumericExpr (&args)[N]) {
    auto builder = factory_.BeginPairwise(kind, N);
    for (int i = 0; i < N; ++i)
      builder.AddArg(args[i]);
    return factory_.EndPairwise(builder);
  }
};

TEST_F(ExprTest, Expr) {
  mp::Expr e;
  EXPECT_TRUE(e == 0);
}

TEST_F(ExprTest, NumericExpr) {
  NumericExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::Expr(e);
}

TEST_F(ExprTest, LogicalExpr) {
  LogicalExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::Expr(e);
}

TEST_F(ExprTest, NumericConstant) {
  mp::NumericConstant e;
  EXPECT_TRUE(e == 0);
  (void)NumericExpr(e);
  e = factory_.MakeNumericConstant(1.23);
  EXPECT_EQ(expr::NUMBER, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(1.23, e.value());
}

TEST_F(ExprTest, Variable) {
  mp::Reference e;
  EXPECT_TRUE(e == 0);
  (void)NumericExpr(e);
  e = factory_.MakeVariable(42);
  EXPECT_EQ(expr::VARIABLE, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(42, e.index());
}

TEST_F(ExprTest, CommonExpr) {
  mp::Reference e;
  EXPECT_TRUE(e == 0);
  (void)NumericExpr(e);
  e = factory_.MakeCommonExpr(42);
  EXPECT_EQ(expr::COMMON_EXPR, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(42, e.index());
}

TEST_F(ExprTest, UnaryExpr) {
  mp::UnaryExpr e;
  EXPECT_TRUE(e == 0);
  (void)NumericExpr(e);
  auto arg = factory_.MakeNumericConstant(42);
  e = MakeUnary(expr::ABS, arg);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::ABS, e.kind());
  EXPECT_EQ(arg, e.arg());
  EXPECT_ASSERT(MakeUnary(expr::ADD, arg), "invalid expression kind");
  EXPECT_ASSERT(MakeUnary(expr::ABS, NumericExpr()),
                "invalid argument");
}

TEST_F(ExprTest, BinaryExpr) {
  mp::BinaryExpr e;
  EXPECT_TRUE(e == 0);
  (void)NumericExpr(e);
  auto lhs = factory_.MakeNumericConstant(42);
  auto rhs = factory_.MakeVariable(0);
  e = MakeBinary(expr::MUL, lhs, rhs);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::MUL, e.kind());
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_ASSERT(MakeBinary(expr::IF, lhs, rhs), "invalid expression kind");
  EXPECT_ASSERT(MakeBinary(expr::MUL, NumericExpr(), rhs), "invalid argument");
  EXPECT_ASSERT(MakeBinary(expr::MUL, lhs, NumericExpr()), "invalid argument");
}

TEST_F(ExprTest, IfExpr) {
  mp::IfExpr e;
  EXPECT_TRUE(e == 0);
  (void)NumericExpr(e);
  auto condition = factory_.MakeLogicalConstant(true);
  auto true_expr = factory_.MakeNumericConstant(42);
  auto false_expr = factory_.MakeVariable(0);
  e = factory_.MakeIf(condition, true_expr, false_expr);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::IF, e.kind());
  EXPECT_EQ(condition, e.condition());
  EXPECT_EQ(true_expr, e.true_expr());
  EXPECT_EQ(false_expr, e.false_expr());
  EXPECT_ASSERT(factory_.MakeIf(LogicalExpr(), true_expr, false_expr),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeIf(condition, NumericExpr(), false_expr),
                "invalid argument");
  factory_.MakeIf(condition, true_expr, NumericExpr());
}

TEST_F(ExprTest, PLTerm) {
  mp::PLTerm e;
  EXPECT_TRUE(e == 0);
  (void)NumericExpr(e);
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
  mp::Function foo = factory_.AddFunction("foo", 42, mp::func::SYMBOLIC);
  f = foo;
  EXPECT_STREQ("foo", f.name());
  EXPECT_EQ(42, f.num_args());
  EXPECT_EQ(mp::func::SYMBOLIC, f.type());
  mp::Function bar = factory_.AddFunction("bar", 0, mp::func::NUMERIC);
  f = bar;
  EXPECT_STREQ("bar", f.name());
  EXPECT_EQ(0, f.num_args());
  EXPECT_EQ(mp::func::NUMERIC, f.type());
  EXPECT_EQ(f, bar);
  EXPECT_NE(f, foo);
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
  EXPECT_ASSERT(builder.AddArg(NumericExpr()), "invalid argument");
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
  EXPECT_EQ(expr::NUMBER, i->kind());
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
  (void)LogicalExpr(e);
  e = factory_.MakeLogicalConstant(false);
  EXPECT_EQ(expr::BOOL, e.kind());
  EXPECT_TRUE(e != 0);
  EXPECT_FALSE(e.value());
  EXPECT_TRUE(factory_.MakeLogicalConstant(true).value());
}

TEST_F(ExprTest, NotExpr) {
  mp::NotExpr e;
  EXPECT_TRUE(e == 0);
  (void)LogicalExpr(e);
  auto arg = factory_.MakeLogicalConstant(false);
  e = factory_.MakeNot(arg);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::NOT, e.kind());
  EXPECT_EQ(arg, e.arg());
  EXPECT_ASSERT(factory_.MakeNot(LogicalExpr()), "invalid argument");
}

TEST_F(ExprTest, BinaryLogicalExpr) {
  mp::BinaryLogicalExpr e;
  EXPECT_TRUE(e == 0);
  (void)LogicalExpr(e);
  auto lhs = factory_.MakeLogicalConstant(true);
  auto rhs = factory_.MakeLogicalConstant(false);
  e = factory_.MakeBinaryLogical(expr::AND, lhs, rhs);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::AND, e.kind());
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_ASSERT(factory_.MakeBinaryLogical(expr::LT, lhs, rhs),
                "invalid expression kind");
  EXPECT_ASSERT(factory_.MakeBinaryLogical(expr::AND, LogicalExpr(), rhs),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeBinaryLogical(expr::AND, lhs, LogicalExpr()),
                "invalid argument");
}

TEST_F(ExprTest, RelationalExpr) {
  mp::RelationalExpr e;
  EXPECT_TRUE(e == 0);
  (void)LogicalExpr(e);
  auto lhs = factory_.MakeNumericConstant(42);
  auto rhs = factory_.MakeVariable(0);
  e = factory_.MakeRelational(expr::EQ, lhs, rhs);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::EQ, e.kind());
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_ASSERT(factory_.MakeRelational(expr::ATLEAST, lhs, rhs),
                "invalid expression kind");
  EXPECT_ASSERT(factory_.MakeRelational(expr::EQ, NumericExpr(), rhs),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeRelational(expr::EQ, lhs, NumericExpr()),
                "invalid argument");
}

TEST_F(ExprTest, LogicalCountExpr) {
  mp::LogicalCountExpr e;
  EXPECT_TRUE(e == 0);
  (void)LogicalExpr(e);
  auto lhs = factory_.MakeNumericConstant(42);
  auto rhs = factory_.EndCount(factory_.BeginCount(0));
  e = factory_.MakeLogicalCount(expr::ATMOST, lhs, rhs);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::ATMOST, e.kind());
  EXPECT_EQ(lhs, e.lhs());
  EXPECT_EQ(rhs, e.rhs());
  EXPECT_ASSERT(factory_.MakeLogicalCount(expr::IMPLICATION, lhs, rhs),
                "invalid expression kind");
  EXPECT_ASSERT(factory_.MakeLogicalCount(expr::ATMOST, NumericExpr(), rhs),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeLogicalCount(expr::ATMOST, lhs, mp::CountExpr()),
                "invalid argument");
}

TEST_F(ExprTest, ImplicationExpr) {
  mp::ImplicationExpr e;
  (void)LogicalExpr(e);
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
                  LogicalExpr(), true_expr, false_expr),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeImplication(
                  condition, LogicalExpr(), false_expr),
                "invalid argument");
  factory_.MakeImplication(condition, true_expr, LogicalExpr());
}

template <typename ExprInfo>
class IteratedExprTest : public ::testing::Test {
 public:
  ExprInfo info_;
};

template <typename ExprType, expr::Kind K,
          typename BaseType = NumericExpr>
struct ExprInfo {
  typedef BaseType Base;
  typedef ExprType Expr;
  typedef NumericExpr Arg;
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
  typedef LogicalExpr Arg;
  typedef ExprFactory::CountExprBuilder Builder;
  Builder BeginBuild(ExprFactory &f, int n) { return f.BeginCount(n); }
  Expr EndBuild(ExprFactory &f, Builder b) { return f.EndCount(b); }
};

struct IteratedLogicalInfo :
    ExprInfo<mp::IteratedLogicalExpr, expr::EXISTS, LogicalExpr> {
  typedef LogicalExpr Arg;
  typedef ExprFactory::IteratedLogicalExprBuilder Builder;
  Builder BeginBuild(ExprFactory &f, int n) {
    return f.BeginIteratedLogical(expr::EXISTS, n);
  }
  Expr EndBuild(ExprFactory &f, Builder b) { return f.EndIteratedLogical(b); }
};

struct PairwiseInfo :
    ExprInfo<mp::PairwiseExpr, expr::ALLDIFF, LogicalExpr> {
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
  NumericExpr args[3];
  explicit TestArgs(ExprFactory &f) {
    args[0] = f.MakeNumericConstant(11);
    args[1] = f.MakeVariable(0);
    args[2] = f.MakeNumericConstant(22);
  }
  NumericExpr operator[](int i) const { return args[i]; }
};

template <>
struct TestArgs<LogicalExpr> {
  LogicalExpr args[3];
  explicit TestArgs(ExprFactory &f) {
    args[0] = f.MakeLogicalConstant(false);
    args[1] = f.MakeLogicalConstant(true);
    args[2] = f.MakeLogicalConstant(false);
  }
  LogicalExpr operator[](int i) const { return args[i]; }
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

TEST_F(ExprTest, SymbolicIfExpr) {
  mp::SymbolicIfExpr e;
  EXPECT_TRUE(e == 0);
  (void)mp::Expr(e);
  auto condition = factory_.MakeLogicalConstant(true);
  mp::Expr true_expr = factory_.MakeStringLiteral("a");
  mp::Expr false_expr = factory_.MakeVariable(0);
  e = factory_.MakeSymbolicIf(condition, true_expr, false_expr);
  EXPECT_TRUE(e != 0);
  EXPECT_EQ(expr::IFSYM, e.kind());
  EXPECT_EQ(condition, e.condition());
  EXPECT_EQ(true_expr, e.true_expr());
  EXPECT_EQ(false_expr, e.false_expr());
  EXPECT_ASSERT(factory_.MakeSymbolicIf(LogicalExpr(),
                                        true_expr, false_expr),
                "invalid argument");
  EXPECT_ASSERT(factory_.MakeSymbolicIf(condition, mp::Expr(), false_expr),
                "invalid argument");
  factory_.MakeSymbolicIf(condition, true_expr, mp::Expr());
}

TEST_F(ExprTest, UncheckedCast) {
  mp::Expr e = factory_.MakeNumericConstant(42);
  mp::internal::UncheckedCast<NumericExpr>(e);
  mp::NumericConstant n = mp::internal::UncheckedCast<mp::NumericConstant>(e);
  EXPECT_EQ(42, n.value());
  EXPECT_ASSERT(mp::internal::UncheckedCast<mp::UnaryExpr>(e), "invalid cast");
}

TEST_F(ExprTest, Cast) {
  mp::Expr e = factory_.MakeNumericConstant(42);
  EXPECT_EQ(e, mp::Cast<NumericExpr>(e));
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

TEST_F(ExprTest, NumberOfExpr) {
  auto arg0 = factory_.MakeVariable(0), arg1 = factory_.MakeVariable(1);
  auto builder = factory_.BeginNumberOf(2, arg0);
  builder.AddArg(arg1);
  auto expr = factory_.EndNumberOf(builder);
  EXPECT_EQ(2, expr.num_args());
  EXPECT_EQ(arg0, expr.arg(0));
  EXPECT_EQ(arg1, expr.arg(1));
}

TEST_F(ExprTest, InvalidNumberOfExprArg) {
  EXPECT_ASSERT(factory_.BeginNumberOf(1, NumericExpr()),
      "invalid argument");
}

TEST_F(ExprTest, SymbolicNumberOfExpr) {
  mp::Expr arg0 = factory_.MakeStringLiteral("abc");
  mp::Expr arg1 = factory_.MakeStringLiteral("def");
  auto builder = factory_.BeginSymbolicNumberOf(2, arg0);
  builder.AddArg(arg1);
  auto expr = factory_.EndSymbolicNumberOf(builder);
  EXPECT_EQ(2, expr.num_args());
  EXPECT_EQ(arg0, expr.arg(0));
  EXPECT_EQ(arg1, expr.arg(1));
}

TEST_F(ExprTest, InvalidSymbolicNumberOfExprArg) {
  EXPECT_ASSERT(factory_.BeginSymbolicNumberOf(1, mp::Expr()),
      "invalid argument");
}

TEST_F(ExprTest, InvalidIteratedLogicalExprKind) {
  EXPECT_ASSERT(factory_.BeginIteratedLogical(expr::ALLDIFF, 1),
                "invalid expression kind");
}

TEST_F(ExprTest, ConversionToExpr) {
  // Test that NumericConstant is not convertible to Expr&. If it was
  // convertible there would be an error because of an ambiguous call.
  // The conversion is forbidden because it compromises type safety
  // as illustrated in the following example:
  //   auto n = factory_.MakeNumericConstant(42);
  //   auto u = MakeUnary(expr::ABS, n);
  //   mp::Expr &e = n;
  //   e = u;
  struct Test {
    static void f(mp::Expr) {}
    static void f(mp::Expr &) {}
  };
  auto n = factory_.MakeNumericConstant(42);
  Test::f(n);
}

TEST_F(ExprTest, ArgumentDependentLookup) {
  NumericExpr e = factory_.MakeNumericConstant(42);
  // IsZero should be found by ADL:
  IsZero(e);
}

TEST_F(ExprTest, Format) {
  auto e = MakeUnary(expr::ABS, factory_.MakeNumericConstant(-42));
  EXPECT_EQ("abs(-42)", fmt::format("{}", e));
}

TEST_F(ExprTest, EqualNumericConstant) {
  EXPECT_TRUE(Equal(MakeConst(0.42), MakeConst(0.42)));
  EXPECT_FALSE(Equal(MakeConst(0.42), MakeConst(42)));
}

TEST_F(ExprTest, EqualVariable) {
  EXPECT_TRUE(Equal(MakeVariable(0), MakeVariable(0)));
  EXPECT_FALSE(Equal(MakeVariable(0), MakeVariable(1)));
  EXPECT_FALSE(Equal(MakeVariable(0), MakeConst(0)));
}

TEST_F(ExprTest, EqualCommonExpr) {
  EXPECT_TRUE(Equal(factory_.MakeCommonExpr(0), factory_.MakeCommonExpr(0)));
  EXPECT_FALSE(Equal(factory_.MakeCommonExpr(0), factory_.MakeCommonExpr(1)));
  EXPECT_FALSE(Equal(factory_.MakeCommonExpr(0), factory_.MakeVariable(0)));
}

TEST_F(ExprTest, EqualUnaryExpr) {
  NumericExpr e = MakeUnary(expr::MINUS, MakeVariable(0));
  EXPECT_TRUE(Equal(e, MakeUnary(expr::MINUS, MakeVariable(0))));
  EXPECT_FALSE(Equal(e, MakeVariable(0)));
  EXPECT_FALSE(Equal(e, MakeUnary(expr::FLOOR, MakeVariable(0))));
  EXPECT_FALSE(Equal(e, MakeUnary(expr::MINUS, MakeVariable(1))));
}

TEST_F(ExprTest, EqualBinaryExpr) {
  NumericExpr e = MakeBinary(expr::ADD, MakeVariable(0), MakeConst(42));
  EXPECT_TRUE(Equal(e, MakeBinary(expr::ADD, MakeVariable(0), MakeConst(42))));
  EXPECT_FALSE(Equal(e, MakeBinary(expr::SUB, MakeVariable(0), MakeConst(42))));
  EXPECT_FALSE(Equal(e, MakeBinary(expr::ADD, MakeConst(42), MakeVariable(0))));
  EXPECT_FALSE(Equal(e, MakeBinary(expr::ADD, MakeVariable(0), MakeConst(0))));
  EXPECT_FALSE(Equal(MakeConst(42), e));
}

TEST_F(ExprTest, EqualIfExpr) {
  NumericExpr e = factory_.MakeIf(
        factory_.MakeLogicalConstant(false), MakeVariable(1), MakeConst(42));
  EXPECT_TRUE(Equal(e, factory_.MakeIf(factory_.MakeLogicalConstant(false),
                                       MakeVariable(1), MakeConst(42))));
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  EXPECT_FALSE(Equal(e, MakeIterated(expr::SUM, args)));
  EXPECT_FALSE(Equal(e, factory_.MakeIf(factory_.MakeLogicalConstant(false),
                                        MakeVariable(1), MakeConst(0))));
  EXPECT_FALSE(Equal(e, MakeConst(42)));
}

TEST_F(ExprTest, EqualPLTerm) {
  double breaks[] = {5, 10};
  double slopes[] = {-1, 0, 1};
  mp::Reference x = MakeVariable(0), y = MakeVariable(1);
  NumericExpr e = MakePLTerm(2, breaks, slopes, x);
  EXPECT_TRUE(Equal(e, MakePLTerm(2, breaks, slopes, x)));
  EXPECT_FALSE(Equal(e, MakePLTerm(1, breaks, slopes, x)));
  EXPECT_FALSE(Equal(e, MakePLTerm(2, breaks, slopes, y)));
  double breaks2[] = {5, 11};
  EXPECT_FALSE(Equal(e, MakePLTerm(2, breaks2, slopes, x)));
  EXPECT_FALSE(Equal(e, MakeConst(42)));
}

TEST_F(ExprTest, EqualCallExpr) {
  NumericExpr args1[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  // args2 is used to make sure that Equal compares expressions structurally
  // instead of comparing pointers; don't replace with args1.
  NumericExpr args2[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  mp::Function f1 = factory_.AddFunction("f1", 0);
  NumericExpr e = MakeCall(f1, args1);
  EXPECT_TRUE(Equal(e, MakeCall(f1, args2)));
  NumericExpr args3[] = {MakeVariable(0), MakeVariable(1)};
  EXPECT_FALSE(Equal(e, MakeCall(f1, args3)));
  EXPECT_FALSE(Equal(MakeCall(f1, args3), MakeCall(f1, args1)));
  EXPECT_FALSE(Equal(e, MakeCall(factory_.AddFunction("f2", 0), args1)));
  NumericExpr args4[] = {MakeVariable(0), MakeVariable(1), MakeConst(0)};
  EXPECT_FALSE(Equal(e, MakeCall(f1, args4)));
  EXPECT_FALSE(Equal(e, MakeConst(42)));
}

TEST_F(ExprTest, EqualVarArgExpr) {
  NumericExpr args1[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  // args2 is used to make sure that Equal compares expressions structurally
  // instead of comparing pointers; don't replace with args1.
  NumericExpr args2[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  NumericExpr e = MakeIterated(expr::MIN, args1);
  EXPECT_TRUE(Equal(e, MakeIterated(expr::MIN, args2)));
  NumericExpr args3[] = {MakeVariable(0), MakeVariable(1)};
  EXPECT_FALSE(Equal(e, MakeIterated(expr::MIN, args3)));
  EXPECT_FALSE(Equal(MakeIterated(expr::MIN, args3),
                     MakeIterated(expr::MIN, args1)));
  EXPECT_FALSE(Equal(e, MakeIterated(expr::MAX, args1)));
  NumericExpr args4[] = {MakeVariable(0), MakeVariable(1), MakeConst(0)};
  EXPECT_FALSE(Equal(e, MakeIterated(expr::MIN, args4)));
  EXPECT_FALSE(Equal(e, MakeConst(42)));
}

TEST_F(ExprTest, EqualSumExpr) {
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  NumericExpr e = MakeIterated(expr::SUM, args);
  // args2 is used to make sure that Equal compares expressions structurally
  // instead of comparing pointers; don't replace with args.
  NumericExpr args2[] = {MakeVariable(0), MakeVariable(1), MakeConst(42)};
  EXPECT_TRUE(Equal(e, MakeIterated(expr::SUM, args2)));
  NumericExpr args3[] = {MakeVariable(0), MakeVariable(1)};
  EXPECT_FALSE(Equal(e, MakeIterated(expr::SUM, args3)));
  EXPECT_FALSE(Equal(MakeIterated(expr::SUM, args3), e));
  LogicalExpr args4[] = {l0, l1, l1};
  EXPECT_FALSE(Equal(e, MakeCount(args4)));
  NumericExpr args5[] = {MakeVariable(0), MakeVariable(1), MakeConst(0)};
  EXPECT_FALSE(Equal(e, MakeIterated(expr::SUM, args5)));
  EXPECT_FALSE(Equal(e, MakeConst(42)));
}

TEST_F(ExprTest, EqualNumberOfExpr) {
  NumericExpr args[] = {MakeVariable(0), MakeVariable(1)};
  EXPECT_TRUE(Equal(MakeIterated(expr::NUMBEROF, args),
                    MakeIterated(expr::NUMBEROF, args)));
}

// TODO
//TEST_F(ExprTest, EqualSymbolicNumberOfExpr) {
//  mp::Expr args[] = {factory_.MakeStringLiteral("test")};
//  EXPECT_TRUE(Equal(MakeSymbolicNumberOf(args), MakeSymbolicNumberOf(args)));
//}

TEST_F(ExprTest, EqualCountExpr) {
  LogicalExpr args[] = {l0, l1, l1};
  NumericExpr e = MakeCount(args);
  // args2 is used to make sure that Equal compares expressions structurally
  // instead of comparing pointers; don't replace with args.
  LogicalExpr args2[] = {
    factory_.MakeLogicalConstant(false),
    factory_.MakeLogicalConstant(true),
    factory_.MakeLogicalConstant(true)
  };
  EXPECT_TRUE(Equal(e, MakeCount(args2)));
  LogicalExpr args3[] = {l0, l1};
  EXPECT_FALSE(Equal(e, MakeCount(args3)));
  EXPECT_FALSE(Equal(MakeCount(args3), e));
  NumericExpr args4[] = {MakeConst(0), MakeConst(1), MakeConst(1)};
  EXPECT_FALSE(Equal(e, MakeIterated(expr::SUM, args4)));
  LogicalExpr args5[] = {l0, l1, l0};
  EXPECT_FALSE(Equal(e, MakeCount(args5)));
}

TEST_F(ExprTest, EqualLogicalConstant) {
  EXPECT_TRUE(Equal(l0, factory_.MakeLogicalConstant(false)));
  EXPECT_FALSE(Equal(l0, factory_.MakeLogicalConstant(true)));
}

TEST_F(ExprTest, EqualNotExpr) {
  EXPECT_TRUE(Equal(factory_.MakeNot(l0), factory_.MakeNot(l0)));
  EXPECT_FALSE(Equal(factory_.MakeNot(l0), factory_.MakeNot(l1)));
}

TEST_F(ExprTest, EqualBinaryLogicalExpr) {
  EXPECT_TRUE(Equal(factory_.MakeBinaryLogical(expr::OR, l0, l1),
                    factory_.MakeBinaryLogical(expr::OR, l0, l1)));
  EXPECT_FALSE(Equal(factory_.MakeBinaryLogical(expr::OR, l0, l1),
                     factory_.MakeBinaryLogical(expr::OR, l0, l0)));
}

TEST_F(ExprTest, EqualRelationalExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  EXPECT_TRUE(Equal(factory_.MakeRelational(expr::LT, n0, n1),
                    factory_.MakeRelational(expr::LT, n0, n1)));
  EXPECT_FALSE(Equal(factory_.MakeRelational(expr::LT, n0, n1),
                     factory_.MakeRelational(expr::LT, n0, n0)));
}

TEST_F(ExprTest, EqualLogicalCountExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  LogicalExpr args[] = {l0};
  auto count = MakeCount(args);
  EXPECT_TRUE(Equal(factory_.MakeLogicalCount(expr::ATMOST, n0, count),
                    factory_.MakeLogicalCount(expr::ATMOST, n0, count)));
  EXPECT_FALSE(Equal(factory_.MakeLogicalCount(expr::ATMOST, n0, count),
                     factory_.MakeLogicalCount(expr::ATMOST, n1, count)));
}

TEST_F(ExprTest, EqualImplicationExpr) {
  EXPECT_TRUE(Equal(factory_.MakeImplication(l0, l1, l0),
                    factory_.MakeImplication(l0, l1, l0)));
  EXPECT_FALSE(Equal(factory_.MakeImplication(l0, l1, l0),
                     factory_.MakeImplication(l0, l1, l1)));
}

TEST_F(ExprTest, EqualIteratedLogicalExpr) {
  LogicalExpr args[] = {l0}, args2[] = {l0}, args3[] = {l1};
  EXPECT_TRUE(Equal(MakeIteratedLogical(expr::EXISTS, args),
                    MakeIteratedLogical(expr::EXISTS, args2)));
  EXPECT_FALSE(Equal(MakeIteratedLogical(expr::EXISTS, args),
                     MakeIteratedLogical(expr::EXISTS, args3)));
}

TEST_F(ExprTest, EqualPairwiseExpr) {
  auto n0 = MakeConst(0), n1 = MakeConst(1);
  NumericExpr args[] = {n0}, args2[] = {n0}, args3[] = {n1};
  EXPECT_TRUE(Equal(MakePairwise(expr::ALLDIFF, args),
                    MakePairwise(expr::ALLDIFF, args2)));
  EXPECT_FALSE(Equal(MakePairwise(expr::ALLDIFF, args),
                     MakePairwise(expr::ALLDIFF, args3)));
}

// TODO: test Equal with StringExpr

TEST(ExprFactoryTest, ExprMemoryAllocation) {
  MockAllocator alloc;
  mp::BasicExprFactory< AllocatorRef<> > f((AllocatorRef<>(&alloc)));
  char buffer[100];
  EXPECT_CALL(alloc, allocate(_)).WillOnce(Return(buffer));
  f.MakeNumericConstant(42);
  EXPECT_CALL(alloc, deallocate(buffer, _));
}

TEST(ExprFactoryTest, FuncMemoryAllocation) {
  MockAllocator alloc;
  mp::BasicExprFactory< AllocatorRef<> > f((AllocatorRef<>(&alloc)));
  char buffer[100];
  EXPECT_CALL(alloc, allocate(_)).WillOnce(Return(buffer));
  f.AddFunction("f", 0);
  EXPECT_CALL(alloc, deallocate(buffer, _));
}

TEST(ExprFactoryTest, IntOverflow) {
  ExprFactory f;
  int int_max = std::numeric_limits<int>::max();
  EXPECT_THROW(f.BeginIterated(expr::SUM, int_max / sizeof(void*) + 2),
               mp::OverflowError);
  EXPECT_THROW(f.BeginPLTerm(int_max / (sizeof(double) * 2) + 1),
               mp::OverflowError);
  std::size_t max_size = std::numeric_limits<std::size_t>::max();
  EXPECT_THROW(f.AddFunction(fmt::StringRef("f", max_size), 0),
               mp::OverflowError);
  EXPECT_THROW(f.MakeStringLiteral(fmt::StringRef("s", int_max + 1)),
               mp::OverflowError);
}

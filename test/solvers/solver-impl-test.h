/*
 Solver implementation test

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

#ifndef TESTS_SOLVER_IMPL_TEST_H_
#define TESTS_SOLVER_IMPL_TEST_H_

#include "../gtest-extra.h"
#include "../solution-handler.h"

#ifdef HAVE_THREADS
# include <thread>
#endif

// Solver implementation test
class SolverImplTest : public ::testing::Test {
 protected:
  Solver solver_;

  typedef Solver::ProblemBuilder ProblemBuilder;
  ProblemBuilder pb_;

  typedef ProblemBuilder::NumericExpr NumericExpr;
  typedef ProblemBuilder::LogicalExpr LogicalExpr;
  typedef ProblemBuilder::Variable Variable;

  Variable x;
  Variable y;
  Variable z;

  FMT_DISALLOW_COPY_AND_ASSIGN(SolverImplTest);

  NumericExpr MakeConst(double value) {
    return pb_.MakeNumericConstant(value);
  }

  bool HasFeature(feature::Feature f) const {
    return (FEATURES & f) != 0;
  }

  class EvalResult {
   private:
    bool has_value_;
    double value_;
    double obj_value_;
    int solve_code_;

   public:
    explicit EvalResult(int solve_code = mp::sol::UNKNOWN)
    : has_value_(false), value_(), obj_value_(), solve_code_(solve_code) {}

    EvalResult(double value, double obj_value)
    : has_value_(true), value_(value), obj_value_(obj_value),
      solve_code_(mp::sol::UNKNOWN) {}

    bool has_value() const { return has_value_; }

    friend bool operator==(double lhs, const EvalResult &rhs) {
      if (!rhs.has_value_)
        throw std::runtime_error("no value");
      return lhs == rhs.value_;
    }

    double obj_value() const {
      if (!has_value_)
        throw std::runtime_error("no value");
      return obj_value_;
    }
    int solve_code() const { return solve_code_; }
    void set_solve_code(int code) { solve_code_ = code; }
  };

  EvalResult Solve(ProblemBuilder &pb);

  // Solves a problem containing a single constraint with the given
  // expression and returns the value of the variable with index 0.
  EvalResult Solve(LogicalExpr e,
      int var1, int var2, int var3 = 0, bool need_result = false);

  // Evaluates a numeric expression by constructing and solving a problem.
  EvalResult Eval(NumericExpr e,
      int var1 = 0, int var2 = 0, int var3 = 0) {
    return Solve(pb_.MakeRelational(mp::expr::EQ, pb_.MakeVariable(0), e),
                 var1, var2, var3, true);
  }

  // Evaluates a logical expression by constructing and solving a problem.
  EvalResult Eval(LogicalExpr e,
      int var1 = 0, int var2 = 0, int var3 = 0) {
    return Eval(pb_.MakeIf(e, MakeConst(1), MakeConst(0)), var1, var2, var3);
  }

  template <typename ArgHandler, typename Arg>
  void AddArgs(ArgHandler &handler, mp::ArrayRef<Arg> args) {
    for (std::size_t i = 0, n = args.size(); i < n; ++i)
      handler.AddArg(args[i]);
  }

  NumericExpr MakeSum(mp::ArrayRef<NumericExpr> args) {
    ProblemBuilder::NumericArgHandler handler = pb_.BeginSum(args.size());
    AddArgs(handler, args);
    return pb_.EndSum(handler);
  }

  ProblemBuilder::CountExpr MakeCount(mp::ArrayRef<LogicalExpr> args) {
    ProblemBuilder::LogicalArgHandler handler = pb_.BeginCount(args.size());
    AddArgs(handler, args);
    return pb_.EndCount(handler);
  }

  NumericExpr MakeNumberOf(NumericExpr value, mp::ArrayRef<NumericExpr> args) {
    ProblemBuilder::NumberOfArgHandler handler =
        pb_.BeginNumberOf(value, args.size());
    AddArgs(handler, args);
    return pb_.EndNumberOf(handler);
  }

  LogicalExpr MakeAllDiff(mp::ArrayRef<NumericExpr> args) {
    ProblemBuilder::AllDiffArgHandler handler = pb_.BeginAllDiff(args.size());
    AddArgs(handler, args);
    return pb_.EndAllDiff(handler);
  }

 public:
  SolverImplTest();
};

void Interrupt();

namespace var = mp::var;
namespace obj = mp::obj;

SolverImplTest::SolverImplTest() {
  auto info = mp::ProblemInfo();
  info.num_vars = 4;
  info.num_objs = 1;
  info.num_funcs = 2;
  pb_.SetInfo(info);
  x = pb_.MakeVariable(1);
  y = pb_.MakeVariable(2);
  z = pb_.MakeVariable(3);
}

SolverImplTest::EvalResult SolverImplTest::Solve(ProblemBuilder &pb) {
  struct TestSolutionHandler : mp::BasicSolutionHandler {
    EvalResult result;
    virtual ~TestSolutionHandler() {}
    void HandleSolution(int status, fmt::StringRef,
          const double *values, const double *, double obj_value) {
      result = values ? EvalResult(values[0], obj_value) : EvalResult();
      result.set_solve_code(status);
    }
  };
  TestSolutionHandler sh;
  solver_.Solve(pb, sh);
  return sh.result;
}

SolverImplTest::EvalResult SolverImplTest::Solve(
    LogicalExpr e, int var1, int var2, int var3, bool need_result) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_nl_integer_vars_in_cons = 4;
  info.num_logical_cons = 1;
  pb.SetInfo(info);
  double inf = std::numeric_limits<double>::infinity();
  pb.AddVar(need_result ? -inf : 0, need_result ? inf : 0, var::INTEGER);
  pb.AddVar(var1, var1, var::INTEGER);
  pb.AddVar(var2, var2, var::INTEGER);
  pb.AddVar(var3, var3, var::INTEGER);
  pb.AddCon(e);
  return Solve(pb);
}

TEST_F(SolverImplTest, Add) {
  NumericExpr e = pb_.MakeBinary(mp::expr::ADD, x, y);
  EXPECT_EQ(25, Eval(e, 10, 15));
  EXPECT_EQ(12, Eval(e, 19, -7));
}

TEST_F(SolverImplTest, Sub) {
  NumericExpr e = pb_.MakeBinary(mp::expr::SUB, x, y);
  EXPECT_EQ(-5, Eval(e, 10, 15));
  EXPECT_EQ(26, Eval(e, 19, -7));
}

TEST_F(SolverImplTest, Mul) {
  NumericExpr e = pb_.MakeBinary(mp::expr::MUL, x, y);
  EXPECT_EQ(150, Eval(e, 10, 15));
  EXPECT_EQ(-133, Eval(e, 19, -7));
}

TEST_F(SolverImplTest, Div) {
  NumericExpr e = pb_.MakeBinary(mp::expr::DIV, x, y);
  if (!HasFeature(feature::DIV)) {
    EXPECT_THROW_MSG(Eval(e, 150, 15), mp::Error, "unsupported: /");
    return;
  }
  EXPECT_EQ(10, Eval(e, 150, 15));
  EXPECT_EQ(-7, Eval(e, -133, 19));
}

TEST_F(SolverImplTest, Mod) {
  NumericExpr e = pb_.MakeBinary(mp::expr::MOD, x, y);
  EXPECT_EQ(0, Eval(e, 9, 3));
  EXPECT_EQ(2, Eval(e, 8, 3));
  EXPECT_EQ(-2, Eval(e, -8, 3));
  EXPECT_EQ(2, Eval(e, 8, -3));
  EXPECT_EQ(-2, Eval(e, -8, -3));
}

TEST_F(SolverImplTest, Pow) {
  NumericExpr e = pb_.MakeBinary(mp::expr::POW, x, y);
  if (!HasFeature(feature::POW)) {
    EXPECT_THROW_MSG(Eval(e, 2, 3), mp::Error, "unsupported: ^");
    return;
  }
  EXPECT_EQ(8, Eval(e, 2, 3));
  EXPECT_EQ(81, Eval(e, 3, 4));
}

TEST_F(SolverImplTest, Less) {
  NumericExpr e = pb_.MakeBinary(mp::expr::LESS, x, y);
  EXPECT_EQ(0, Eval(e, 10, 15));
  EXPECT_EQ(26, Eval(e, 19, -7));
}

TEST_F(SolverImplTest, Min) {
  typename ProblemBuilder::NumericArgHandler handler =
      pb_.BeginVarArg(mp::expr::MIN, 3);
  handler.AddArg(x);
  handler.AddArg(y);
  handler.AddArg(z);
  NumericExpr e = pb_.EndVarArg(handler);
  EXPECT_EQ(-7, Eval(e, 3, -7, 5));
  EXPECT_EQ(10, Eval(e, 10, 20, 30));
}

TEST_F(SolverImplTest, Max) {
  NumericExpr args[] = {x, y, z};
  NumericExpr e = pb_.MakeVarArg(mp::expr::MAX, args);
  EXPECT_EQ(5, Eval(e, 3, -7, 5));
  EXPECT_EQ(30, Eval(e, 30, 20, 10));
}

TEST_F(SolverImplTest, Floor) {
  NumericExpr e = pb_.MakeUnary(mp::expr::FLOOR, x);
  EXPECT_EQ(-42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
  if (!HasFeature(feature::FLOAT_CONST)) return;
  EXPECT_EQ(4, Eval(pb_.MakeUnary(mp::expr::FLOOR, MakeConst(4.9))));
  EXPECT_EQ(-5, Eval(pb_.MakeUnary(mp::expr::FLOOR, MakeConst(-4.1))));
}

TEST_F(SolverImplTest, Ceil) {
  NumericExpr e = pb_.MakeUnary(mp::expr::CEIL, x);
  EXPECT_EQ(-42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
  if (!HasFeature(feature::FLOAT_CONST)) return;
  EXPECT_EQ(5, Eval(pb_.MakeUnary(mp::expr::CEIL, MakeConst(4.1))));
  EXPECT_EQ(-4, Eval(pb_.MakeUnary(mp::expr::CEIL, MakeConst(-4.9))));
}

TEST_F(SolverImplTest, Abs) {
  NumericExpr e = pb_.MakeUnary(mp::expr::ABS, x);
  EXPECT_EQ(42, Eval(e, -42));
  EXPECT_EQ(42, Eval(e, 42));
}

TEST_F(SolverImplTest, Minus) {
  NumericExpr e = pb_.MakeUnary(mp::expr::MINUS, x);
  EXPECT_EQ(42, Eval(e, -42));
  EXPECT_EQ(-42, Eval(e, 42));
}

TEST_F(SolverImplTest, If) {
  NumericExpr e = pb_.MakeIf(pb_.MakeRelational(
                               mp::expr::EQ, x, MakeConst(1)), y, z);
  EXPECT_EQ(42, Eval(e, 1, 42, 10));
  EXPECT_EQ(10, Eval(e, 0, 42, 10));
  EXPECT_EQ(42, Eval(e, 1, 42, 42));
}

TEST_F(SolverImplTest, Tanh) {
  double arg = (std::log(1.5) - std::log(0.5)) / 2;
  NumericExpr e = pb_.MakeBinary(mp::expr::MUL,
      MakeConst(2), pb_.MakeUnary(mp::expr::TANH, MakeConst(arg)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW_MSG(Eval(e), mp::Error, "unsupported: tanh");
    return;
  }
  EXPECT_EQ(1, Eval(e));
}

TEST_F(SolverImplTest, Tan) {
  EXPECT_THROW_MSG(Eval(pb_.MakeUnary(mp::expr::TAN, x)),
                   mp::Error, "unsupported: tan");
}

TEST_F(SolverImplTest, Sqrt) {
  NumericExpr e = pb_.MakeUnary(mp::expr::SQRT, x);
  if (!HasFeature(feature::SQRT)) {
    EXPECT_THROW_MSG(Eval(e, 64), mp::Error, "unsupported: sqrt");
    return;
  }
  EXPECT_EQ(8, Eval(e, 64));
}

TEST_F(SolverImplTest, Sinh) {
  NumericExpr e = pb_.MakeUnary(
        mp::expr::SINH, MakeConst(std::log(2 + std::sqrt(5.0))));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW_MSG(Eval(e), mp::Error, "unsupported: sinh");
    return;
  }
  EXPECT_EQ(2, Eval(e));
}

TEST_F(SolverImplTest, Sin) {
  EXPECT_THROW_MSG(Eval(pb_.MakeUnary(mp::expr::SIN, x)),
                   mp::Error, "unsupported: sin");
}

TEST_F(SolverImplTest, Log10) {
  NumericExpr e = pb_.MakeUnary(mp::expr::LOG10, MakeConst(1000));
  if (!HasFeature(feature::LOG)) {
    EXPECT_THROW_MSG(Eval(e), mp::Error, "unsupported: log10");
    return;
  }
  EXPECT_EQ(3, Eval(e));
}

TEST_F(SolverImplTest, Log) {
  NumericExpr e = pb_.MakeUnary(mp::expr::LOG, MakeConst(std::exp(5.0)));
  if (!HasFeature(feature::LOG)) {
    EXPECT_THROW_MSG(Eval(e), mp::Error, "unsupported: log");
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_F(SolverImplTest, Exp) {
  NumericExpr e = pb_.MakeUnary(mp::expr::EXP, MakeConst(std::log(5.0)));
  if (!HasFeature(feature::EXP)) {
    EXPECT_THROW_MSG(Eval(e), mp::Error, "unsupported: exp");
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_F(SolverImplTest, Cosh) {
  double x = 5;
  NumericExpr e = pb_.MakeUnary(mp::expr::COSH,
      MakeConst(std::log(x + std::sqrt(x + 1) * std::sqrt(x - 1))));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW_MSG(Eval(e), mp::Error, "unsupported: cosh");
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_F(SolverImplTest, Cos) {
  EXPECT_THROW_MSG(Eval(pb_.MakeUnary(mp::expr::COS, x)),
                   mp::Error, "unsupported: cos");
}

TEST_F(SolverImplTest, Atanh) {
  NumericExpr x = pb_.MakeUnary(mp::expr::ATANH, MakeConst(std::tanh(5.0)));
  NumericExpr e = pb_.MakeUnary(mp::expr::FLOOR, pb_.MakeBinary(mp::expr::ADD,
      MakeConst(0.5), pb_.MakeBinary(mp::expr::MUL, MakeConst(1000000), x)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW_MSG(Eval(e), mp::Error, "unsupported: atanh");
    return;
  }
  EXPECT_EQ(5000000, Eval(e));
}

TEST_F(SolverImplTest, Atan2) {
  EXPECT_THROW_MSG(Eval(pb_.MakeBinary(mp::expr::ATAN2, x, y)),
                   mp::Error, "unsupported: atan2");
}

TEST_F(SolverImplTest, Atan) {
  EXPECT_THROW_MSG(Eval(pb_.MakeUnary(mp::expr::ATAN, x)),
                   mp::Error, "unsupported: atan");
}

TEST_F(SolverImplTest, Asinh) {
  NumericExpr e = pb_.MakeUnary(mp::expr::ASINH, MakeConst(std::sinh(5.0)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW_MSG(Eval(e), mp::Error, "unsupported: asinh");
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_F(SolverImplTest, Asin) {
  EXPECT_THROW_MSG(Eval(pb_.MakeUnary(mp::expr::ASIN, x)),
                   mp::Error, "unsupported: asin");
}

TEST_F(SolverImplTest, Acosh) {
  NumericExpr e = pb_.MakeUnary(mp::expr::ACOSH, MakeConst(std::cosh(5.0)));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW_MSG(Eval(e), mp::Error, "unsupported: acosh");
    return;
  }
  EXPECT_EQ(5, Eval(e));
}

TEST_F(SolverImplTest, Acos) {
  EXPECT_THROW_MSG(Eval(pb_.MakeUnary(mp::expr::ACOS, x)),
                   mp::Error, "unsupported: acos");
}

TEST_F(SolverImplTest, Sum) {
  NumericExpr args[] = {x, y, z};
  using mp::MakeArrayRef;
  EXPECT_EQ(0, Eval(MakeSum(MakeArrayRef(args, 0))));
  EXPECT_EQ(42, Eval(MakeSum(MakeArrayRef(args, 1)), 42));
  EXPECT_EQ(123, Eval(MakeSum(args), 100, 20, 3));
}

TEST_F(SolverImplTest, IntDiv) {
  NumericExpr e = pb_.MakeBinary(mp::expr::INT_DIV, x, y);
  EXPECT_EQ(3, Eval(e, 9, 3));
  EXPECT_EQ(2, Eval(e, 8, 3));
  EXPECT_EQ(-2, Eval(e, -8, 3));
  EXPECT_EQ(-2, Eval(e, 8, -3));
  EXPECT_EQ(2, Eval(e, -8, -3));
}

TEST_F(SolverImplTest, Precision) {
  EXPECT_THROW_MSG(Eval(pb_.MakeBinary(mp::expr::PRECISION, x, y)),
                   mp::Error, "unsupported: precision");
}

TEST_F(SolverImplTest, Round) {
  mp::expr::Kind round = mp::expr::ROUND;
  EXPECT_EQ(42, Eval(pb_.MakeBinary(round, x, MakeConst(0)), 42));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(4, Eval(pb_.MakeBinary(round, MakeConst(4.4), MakeConst(0))));
    EXPECT_EQ(5, Eval(pb_.MakeBinary(round, MakeConst(4.6), MakeConst(0))));
    EXPECT_EQ(-4, Eval(pb_.MakeBinary(round, MakeConst(-4.4), MakeConst(0))));
    EXPECT_EQ(-5, Eval(pb_.MakeBinary(round, MakeConst(-4.6), MakeConst(0))));
  }
  EXPECT_THROW_MSG(Eval(pb_.MakeBinary(round, x, MakeConst(1))),
                   mp::Error, "unsupported: round");
  EXPECT_THROW_MSG(Eval(pb_.MakeBinary(round, x, y)),
                   mp::Error, "unsupported: round");
}

TEST_F(SolverImplTest, Trunc) {
  mp::expr::Kind trunc = mp::expr::TRUNC;
  EXPECT_EQ(42, Eval(pb_.MakeBinary(trunc, x, MakeConst(0)), 42));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(4, Eval(pb_.MakeBinary(trunc, MakeConst(4.4), MakeConst(0))));
    EXPECT_EQ(4, Eval(pb_.MakeBinary(trunc, MakeConst(4.6), MakeConst(0))));
    EXPECT_EQ(-4, Eval(pb_.MakeBinary(trunc, MakeConst(-4.4), MakeConst(0))));
    EXPECT_EQ(-4, Eval(pb_.MakeBinary(trunc, MakeConst(-4.6), MakeConst(0))));
  }
  EXPECT_THROW_MSG(Eval(pb_.MakeBinary(trunc, x, MakeConst(1))),
                   mp::Error, "unsupported: trunc");
  EXPECT_THROW_MSG(Eval(pb_.MakeBinary(trunc, x, y)),
                   mp::Error, "unsupported: trunc");
}

TEST_F(SolverImplTest, Count) {
  LogicalExpr args[] = {
    pb_.MakeRelational(mp::expr::NE, x, MakeConst(0)),
    pb_.MakeRelational(mp::expr::NE, y, MakeConst(0)),
    pb_.MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  EXPECT_EQ(0, Eval(MakeCount(args)));
  EXPECT_EQ(1, Eval(MakeCount(args), 1));
  EXPECT_EQ(2, Eval(MakeCount(args), 0, 1, 1));
  EXPECT_EQ(3, Eval(MakeCount(args), 1, 1, 1));
}

TEST_F(SolverImplTest, NumberOf) {
  NumericExpr value = MakeConst(42);
  NumericExpr args[] = {x};
  EXPECT_EQ(0, Eval(MakeNumberOf(value, args)));
  EXPECT_EQ(1, Eval(MakeNumberOf(value, args), 42));
  NumericExpr args2[] = {x, y};
  EXPECT_EQ(0, Eval(MakeNumberOf(value, args2)));
  EXPECT_EQ(1, Eval(MakeNumberOf(value, args2), 0, 42));
  EXPECT_EQ(2, Eval(MakeNumberOf(value, args2), 42, 42));
}

TEST_F(SolverImplTest, PiecewiseLinear) {
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
  ProblemBuilder::PLTermHandler handler = pb_.BeginPLTerm(2);
  handler.AddSlope(-1);
  handler.AddBreakpoint(3);
  handler.AddSlope(0);
  handler.AddBreakpoint(6);
  handler.AddSlope(1);
  NumericExpr e = pb_.EndPLTerm(handler, pb_.MakeVariable(1));
  if (!HasFeature(feature::PLTERM)) {
    EXPECT_THROW_MSG(Eval(e, 42), mp::Error,
                     "unsupported: piecewise-linear term");
    return;
  }
  EXPECT_EQ(33, Eval(e, 42));
  EXPECT_EQ(-3, Eval(e, 4));
  EXPECT_EQ(1, Eval(e, -1));
}

TEST_F(SolverImplTest, UnsupportedFunctionCall) {
  EXPECT_THROW_MSG(pb_.SetFunction(0, "foo", 2, mp::func::NUMERIC),
                   mp::Error, "unsupported: function");
  EXPECT_THROW_MSG(pb_.BeginCall(0, 2),
                   mp::Error, "unsupported: function call");
}

TEST_F(SolverImplTest, PowConstExp) {
  EXPECT_EQ(16, Eval(pb_.MakeBinary(
                       mp::expr::POW_CONST_EXP, x, MakeConst(4)), 2));
}

TEST_F(SolverImplTest, Pow2) {
  EXPECT_EQ(49, Eval(pb_.MakeUnary(mp::expr::POW2, x), 7));
}

TEST_F(SolverImplTest, PowConstBase) {
  NumericExpr e = pb_.MakeBinary(mp::expr::POW_CONST_BASE, MakeConst(5), x);
  if (!HasFeature(feature::POW)) {
    EXPECT_THROW_MSG(Eval(e, 3), mp::Error, "unsupported: ^");
    return;
  }
  EXPECT_EQ(125, Eval(e, 3));
}

TEST_F(SolverImplTest, NumericConstant) {
  EXPECT_EQ(42, Eval(MakeConst(42)));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(42, Eval(pb_.MakeBinary(
                         mp::expr::MUL, MakeConst(0.42), MakeConst(100))));
    return;
  }
  EXPECT_THROW_MSG(Eval(MakeConst(0.42)), mp::Error,
    "value 0.42 can't be represented as int");
}

TEST_F(SolverImplTest, Var) {
  EXPECT_EQ(11, Eval(x, 11, 22));
  EXPECT_EQ(22, Eval(y, 11, 22));
  EXPECT_EQ(33, Eval(x, 33));
}

TEST_F(SolverImplTest, Or) {
  NumericExpr one = MakeConst(1);
  LogicalExpr e = pb_.MakeBinaryLogical(
      mp::expr::OR, pb_.MakeRelational(mp::expr::EQ, x, one),
        pb_.MakeRelational(mp::expr::EQ, y, one));
  EXPECT_EQ(0, Eval(e, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_F(SolverImplTest, And) {
  NumericExpr one = MakeConst(1);
  LogicalExpr e = pb_.MakeBinaryLogical(
      mp::expr::AND, pb_.MakeRelational(mp::expr::EQ, x, one),
        pb_.MakeRelational(mp::expr::EQ, y, one));
  EXPECT_EQ(0, Eval(e, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_F(SolverImplTest, LT) {
  LogicalExpr e = pb_.MakeRelational(mp::expr::LT, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_F(SolverImplTest, LE) {
  LogicalExpr e = pb_.MakeRelational(mp::expr::LE, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_F(SolverImplTest, EQ) {
  LogicalExpr e = pb_.MakeRelational(mp::expr::EQ, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(0, Eval(e, 5, 3));
}

TEST_F(SolverImplTest, GE) {
  LogicalExpr e = pb_.MakeRelational(mp::expr::GE, x, y);
  EXPECT_EQ(1, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_F(SolverImplTest, GT) {
  LogicalExpr e = pb_.MakeRelational(mp::expr::GT, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(0, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_F(SolverImplTest, NE) {
  LogicalExpr e = pb_.MakeRelational(mp::expr::NE, x, y);
  EXPECT_EQ(0, Eval(e, 3, 3));
  EXPECT_EQ(1, Eval(e, 3, 5));
  EXPECT_EQ(1, Eval(e, 5, 3));
}

TEST_F(SolverImplTest, Not) {
  LogicalExpr e = pb_.MakeNot(
        pb_.MakeRelational(mp::expr::EQ, x, MakeConst(1)));
  EXPECT_EQ(1, Eval(e, 0));
  EXPECT_EQ(0, Eval(e, 1));
}

TEST_F(SolverImplTest, AtLeast) {
  LogicalExpr args[] = {
    pb_.MakeRelational(mp::expr::NE, y, MakeConst(0)),
    pb_.MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = pb_.MakeLogicalCount(mp::expr::ATLEAST, x, MakeCount(args));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_F(SolverImplTest, AtMost) {
  LogicalExpr args[] = {
    pb_.MakeRelational(mp::expr::NE, y, MakeConst(0)),
    pb_.MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = pb_.MakeLogicalCount(mp::expr::ATMOST, x, MakeCount(args));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_F(SolverImplTest, Exactly) {
  LogicalExpr args[] = {
    pb_.MakeRelational(mp::expr::NE, y, MakeConst(0)),
    pb_.MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = pb_.MakeLogicalCount(mp::expr::EXACTLY, x, MakeCount(args));
  EXPECT_EQ(1, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(1, Eval(e, 2, 1, 1));
}

TEST_F(SolverImplTest, NotAtLeast) {
  LogicalExpr args[] = {
    pb_.MakeRelational(mp::expr::NE, y, MakeConst(0)),
    pb_.MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = pb_.MakeLogicalCount(
        mp::expr::NOT_ATLEAST, x, MakeCount(args));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_F(SolverImplTest, NotAtMost) {
  LogicalExpr args[] = {
    pb_.MakeRelational(mp::expr::NE, y, MakeConst(0)),
    pb_.MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = pb_.MakeLogicalCount(
        mp::expr::NOT_ATMOST, x, MakeCount(args));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(0, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_F(SolverImplTest, NotExactly) {
  LogicalExpr args[] = {
    pb_.MakeRelational(mp::expr::NE, y, MakeConst(0)),
    pb_.MakeRelational(mp::expr::NE, z, MakeConst(0))
  };
  LogicalExpr e = pb_.MakeLogicalCount(
        mp::expr::NOT_EXACTLY, x, MakeCount(args));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
  EXPECT_EQ(1, Eval(e, 2, 0, 1));
  EXPECT_EQ(0, Eval(e, 2, 1, 1));
}

TEST_F(SolverImplTest, ForAll) {
  LogicalExpr args[] = {
    pb_.MakeRelational(mp::expr::EQ, x, MakeConst(1)),
    pb_.MakeRelational(mp::expr::EQ, y, MakeConst(1)),
    pb_.MakeRelational(mp::expr::EQ, z, MakeConst(1))
  };
  LogicalExpr e = pb_.MakeIteratedLogical(mp::expr::FORALL, args);
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 0, 1));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(0, Eval(e, 0, 1, 1));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
}

TEST_F(SolverImplTest, Exists) {
  LogicalExpr args[] = {
    pb_.MakeRelational(mp::expr::EQ, x, MakeConst(1)),
    pb_.MakeRelational(mp::expr::EQ, y, MakeConst(1)),
    pb_.MakeRelational(mp::expr::EQ, z, MakeConst(1))
  };
  LogicalExpr e = pb_.MakeIteratedLogical(mp::expr::EXISTS, args);
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 0, 1));
  EXPECT_EQ(1, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 1));
  EXPECT_EQ(1, Eval(e, 1, 0, 0));
  EXPECT_EQ(1, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
}

TEST_F(SolverImplTest, Implication) {
  LogicalExpr e = pb_.MakeImplication(
      pb_.MakeRelational(mp::expr::EQ, x, MakeConst(1)),
      pb_.MakeRelational(mp::expr::EQ, y, MakeConst(1)),
      pb_.MakeRelational(mp::expr::EQ, z, MakeConst(1)));
  EXPECT_EQ(0, Eval(e, 0, 0, 0));
  EXPECT_EQ(1, Eval(e, 0, 0, 1));
  EXPECT_EQ(0, Eval(e, 0, 1, 0));
  EXPECT_EQ(1, Eval(e, 0, 1, 1));
  EXPECT_EQ(0, Eval(e, 1, 0, 0));
  EXPECT_EQ(0, Eval(e, 1, 0, 1));
  EXPECT_EQ(1, Eval(e, 1, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1, 1));
}

TEST_F(SolverImplTest, Iff) {
  LogicalExpr e = pb_.MakeBinaryLogical(mp::expr::IFF,
      pb_.MakeRelational(mp::expr::EQ, x, MakeConst(1)),
      pb_.MakeRelational(mp::expr::EQ, y, MakeConst(1)));
  EXPECT_EQ(1, Eval(e, 0, 0));
  EXPECT_EQ(0, Eval(e, 0, 1));
  EXPECT_EQ(0, Eval(e, 1, 0));
  EXPECT_EQ(1, Eval(e, 1, 1));
}

TEST_F(SolverImplTest, AllDiff) {
  NumericExpr args[] = {MakeConst(1), x, y};
  LogicalExpr e = MakeAllDiff(args);
  EXPECT_TRUE(Solve(e, 2, 3).has_value());
  EXPECT_FALSE(Solve(e, 2, 1).has_value());
  EXPECT_FALSE(Solve(e, 1, 1).has_value());
}

TEST_F(SolverImplTest, NestedAllDiff) {
  NumericExpr args[] = {MakeConst(1), x, y};
  EXPECT_THROW_MSG(Eval(pb_.MakeNot(MakeAllDiff(args)), 1, 2),
                   mp::Error, "unsupported: alldiff");
}

TEST_F(SolverImplTest, LogicalConstant) {
  EXPECT_EQ(0, Eval(pb_.MakeLogicalConstant(false)));
  EXPECT_EQ(1, Eval(pb_.MakeLogicalConstant(true)));
}

TEST_F(SolverImplTest, NonlinearObj) {
  ProblemBuilder pb;
  auto info = mp::ProblemInfo();
  info.num_objs = info.num_vars = info.num_common_exprs_in_objs = 1;
  pb.SetInfo(info);
  pb.AddVar(2, 2, var::INTEGER);
  NumericExpr x = pb.MakeVariable(0);
  pb.AddObj(obj::MIN, pb.MakeBinary(mp::expr::MUL, x, x), 0);
  EXPECT_EQ(4, Solve(pb).obj_value());
}

TEST_F(SolverImplTest, ObjConst) {
  ProblemBuilder pb;
  pb.AddVar(0, 0, var::INTEGER);
  pb.AddObj(obj::MIN, MakeConst(42), 0);
  EXPECT_EQ(42, Solve(pb).obj_value());
}

TEST_F(SolverImplTest, Minimize) {
  ProblemBuilder pb;
  pb.AddVar(42, 100, var::INTEGER);
  pb.AddObj(obj::MIN, pb.MakeVariable(0), 0);
  EXPECT_EQ(42, Solve(pb).obj_value());
}

TEST_F(SolverImplTest, Maximize) {
  ProblemBuilder pb;
  pb.AddVar(0, 42, var::INTEGER);
  pb.AddObj(obj::MAX, pb.MakeVariable(0), 0);
  EXPECT_EQ(42, Solve(pb).obj_value());
}

TEST_F(SolverImplTest, TimingOption) {
  struct TestOutputHandler : mp::OutputHandler {
    std::string output;

    virtual ~TestOutputHandler() {}
    void HandleOutput(fmt::StringRef output) {
      this->output += output;
    }
  };
  TestOutputHandler oh;
  solver_.set_output_handler(&oh);

  ProblemBuilder pb;
  pb.AddVar(42, 100, var::INTEGER);
  pb.AddObj(obj::MIN, pb.MakeVariable(0), 0);

  solver_.SetIntOption("timing", 0);
  mp::BasicSolutionHandler sol_handler;
  solver_.Solve(pb, sol_handler);
  EXPECT_TRUE(oh.output.find("Setup time = ") == std::string::npos);
  EXPECT_TRUE(oh.output.find("Solution time = ") == std::string::npos);
  EXPECT_TRUE(oh.output.find("Output time = ") == std::string::npos);

  solver_.SetIntOption("timing", 1);
  solver_.Solve(pb, sol_handler);
  EXPECT_TRUE(oh.output.find("Setup time = ") != std::string::npos);
  EXPECT_TRUE(oh.output.find("Solution time = ") != std::string::npos);
  EXPECT_TRUE(oh.output.find("Output time = ") != std::string::npos);
}

// ----------------------------------------------------------------------------
// Solve code tests

TEST_F(SolverImplTest, OptimalSolveCode) {
  pb_.AddObj(obj::MIN, MakeConst(42), 0);
  EXPECT_EQ(mp::sol::SOLVED, Solve(pb_).solve_code());
}

TEST_F(SolverImplTest, FeasibleSolveCode) {
  pb_.AddCon(MakeConst(42), 0, 100, 0);
  EXPECT_EQ(mp::sol::SOLVED, Solve(pb_).solve_code());
}

TEST_F(SolverImplTest, InfeasibleSolveCode) {
  pb_.AddCon(MakeConst(0), 1, 1, 0);
  EXPECT_EQ(mp::sol::INFEASIBLE, Solve(pb_).solve_code());
}

// ----------------------------------------------------------------------------
// Interrupt tests

#ifdef HAVE_THREADS
void Interrupt() {
  // Wait until started.
  while (mp::SignalHandler::stop())
    std::this_thread::yield();
  std::raise(SIGINT);
}

TEST_F(SolverImplTest, InterruptSolution) {
  std::thread t(Interrupt);
  Problem p;
  string message = Solve(p, "miplib/assign1").message;
  t.join();
  EXPECT_EQ(600, p.solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}
#endif

struct SolutionCounter : mp::BasicSolutionHandler {
  int num_solutions;
  SolutionCounter() : num_solutions(0) {}
  void HandleFeasibleSolution(
      fmt::StringRef, const double *, const double *, double) {
    ++num_solutions;
  }
};

TEST_F(SolverImplTest, CountSolutions) {
  ProblemBuilder pb;
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {
    pb.MakeVariable(0), pb.MakeVariable(1), pb.MakeVariable(2)
  };
  pb.AddCon(MakeAllDiff(args));
  solver_.SetIntOption("solutionlimit", 10);
  solver_.SetIntOption("countsolutions", 1);
  SolutionCounter sc;
  solver_.Solve(pb, sc);
  EXPECT_EQ(6, sc.num_solutions);
}

TEST_F(SolverImplTest, SatisfactionSolutionLimit) {
  ProblemBuilder pb;
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {
    pb.MakeVariable(0), pb.MakeVariable(1), pb.MakeVariable(2)
  };
  pb.AddCon(MakeAllDiff(args));
  solver_.SetIntOption("solutionlimit", 5);
  TestSolutionHandler sh;
  solver_.Solve(pb, sh);
  EXPECT_EQ(0, sh.status());
}

TEST_F(SolverImplTest, OptimizationSolutionLimit) {
  ProblemBuilder pb;
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {
    pb.MakeVariable(0), pb.MakeVariable(1), pb.MakeVariable(2)
  };
  pb.AddCon(MakeAllDiff(args));
  solver_.SetIntOption("solutionlimit", 2);
  TestSolutionHandler sh;
  solver_.Solve(pb, sh);
  EXPECT_GE(sh.status(), 400);
  EXPECT_LT(sh.status(), 499);
}

TEST_F(SolverImplTest, MultipleSolutions) {
  ProblemBuilder pb;
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddVar(1, 3, var::INTEGER);
  NumericExpr args[] = {
    pb.MakeVariable(0), pb.MakeVariable(1), pb.MakeVariable(2)
  };
  pb.AddCon(MakeAllDiff(args));
  solver_.SetIntOption("solutionlimit", 3);
  solver_.SetStrOption("solutionstub", "test");
  SolutionCounter sc;
  solver_.Solve(pb, sc);
  EXPECT_EQ(3, sc.num_solutions);
}

TEST_F(SolverImplTest, OptionValues) {
  for (auto i = solver_.option_begin(), e = solver_.option_end(); i != e; ++i) {
    for (mp::ValueArrayRef::iterator j = i->values().begin(),
        value_end = i->values().end(); j != value_end; ++j) {
      EXPECT_TRUE(j->value != 0);
    }
  }
}

TEST_F(SolverImplTest, CreateSolver) {
  EXPECT_STREQ(solver_.name(), mp::CreateSolver(0)->name());
}

#ifndef MP_TEST_DATA_DIR
# define MP_TEST_DATA_DIR "../data"
#endif

#endif  // TESTS_SOLVER_IMPL_TEST_H_

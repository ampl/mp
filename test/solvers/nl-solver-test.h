/*
 Nonliner solver test

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

#ifndef TESTS_NL_SOLVER_TEST_H_
#define TESTS_NL_SOLVER_TEST_H_

#include <cmath>
#include <string>

#if MP_THREAD && MP_THREAD_INTERRUPT
// Workaround the error "'yield' is not a member of 'std::this_thread'"
// on older gcc.
# if FMT_GCC_VERSION < 408
#  define _GLIBCXX_USE_SCHED_YIELD 1
# endif
# include <thread>
#endif

#include "../gtest-extra.h"
#include "../solution-handler.h"

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

template <typename Solver>
EvalResult Solve(Solver &solver, typename Solver::ProblemBuilder &pb) {
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
  solver.Solve(pb, sh);
  return sh.result;
}

// Solver implementation test
class NLSolverTest : public ::testing::Test {
 protected:
  Solver solver_;

  typedef Solver::ProblemBuilder ProblemBuilder;
  typedef ProblemBuilder::NumericExpr NumericExpr;
  typedef ProblemBuilder::LogicalExpr LogicalExpr;
  typedef ProblemBuilder::Variable Variable;

  FMT_DISALLOW_COPY_AND_ASSIGN(NLSolverTest);

  bool HasFeature(feature::Feature f) const {
    return (FEATURES & f) != 0;
  }

  void SetInfo(ProblemBuilder &pb) {
    auto info = mp::ProblemInfo();
    info.num_vars = info.num_nl_integer_vars_in_objs = 1;
    info.num_objs = info.num_nl_objs = 1;
    pb.SetInfo(info);
  }

  EvalResult Solve(ProblemBuilder &pb) { return ::Solve(solver_, pb); }

  class ExprFactory {
   protected:
    mutable Variable x, y, z;
    mutable NumericExpr one;

    void Init(ProblemBuilder &pb) const {
      x = pb.MakeVariable(1);
      y = pb.MakeVariable(2);
      z = pb.MakeVariable(3);
      one = pb.MakeNumericConstant(1);
    }

   public:
    virtual ~ExprFactory() {}
  };

  // An abstract factory for numeric expressions.
  class NumericExprFactory : public ExprFactory {
   protected:
    virtual NumericExpr Create(ProblemBuilder &pb) const = 0;

   public:
    NumericExpr operator()(ProblemBuilder &pb) const {
      Init(pb);
      return Create(pb);
    }
  };

  // A constant expression factory.
  class NumericConstantFactory : public NumericExprFactory {
   private:
    double value_;

   public:
    NumericConstantFactory(double value) : value_(value) {}
    virtual NumericExpr Create(ProblemBuilder &pb) const {
      return pb.MakeNumericConstant(value_);
    }
  };

  // A variable factory.
  class VariableFactory : public NumericExprFactory {
   private:
    int var_index_;

   public:
    VariableFactory(int var_index) : var_index_(var_index) {}

    virtual NumericExpr Create(ProblemBuilder &pb) const {
      return pb.MakeVariable(var_index_);
    }
  };

  // A unary expression factory.
  template <typename ArgFactory = VariableFactory>
  class UnaryExprFactory : public NumericExprFactory {
   private:
    mp::expr::Kind kind_;
    ArgFactory arg_factory_;

   public:
    explicit UnaryExprFactory(
        mp::expr::Kind kind, ArgFactory arg_factory = ArgFactory())
      : kind_(kind), arg_factory_(arg_factory) {}

    NumericExpr Create(ProblemBuilder &pb) const {
      return pb.MakeUnary(kind_, arg_factory_(pb));
    }
  };

  // A binary expression factory.
  template <typename LHSFactory = VariableFactory,
            typename RHSFactory = VariableFactory>
  class BinaryExprFactory : public NumericExprFactory {
   private:
    mp::expr::Kind kind_;
    LHSFactory lhs_factory_;
    RHSFactory rhs_factory_;

   public:
    BinaryExprFactory(mp::expr::Kind kind,
                      LHSFactory lhs_factory, RHSFactory rhs_factory)
      : kind_(kind), lhs_factory_(lhs_factory), rhs_factory_(rhs_factory) {}

    NumericExpr Create(ProblemBuilder &pb) const {
      return pb.MakeBinary(kind_, lhs_factory_(pb), rhs_factory_(pb));
    }
  };

  class VarArgFactory : public NumericExprFactory {
   private:
    mp::expr::Kind kind_;

   public:
    explicit VarArgFactory(mp::expr::Kind kind) : kind_(kind) {}

    NumericExpr Create(ProblemBuilder &pb) const {
      auto args = pb.BeginVarArg(kind_, 3);
      args.AddArg(x);
      args.AddArg(y);
      args.AddArg(z);
      return pb.EndVarArg(args);
    }
  };

  template <typename ArgFactory>
  UnaryExprFactory<ArgFactory>
      MakeUnaryExprFactory(mp::expr::Kind kind, ArgFactory arg) {
    return UnaryExprFactory<ArgFactory>(kind, arg);
  }

  UnaryExprFactory<NumericConstantFactory>
      MakeUnaryExprFactory(mp::expr::Kind kind, double arg) {
    return UnaryExprFactory<NumericConstantFactory>(
          kind, NumericConstantFactory(arg));
  }

  BinaryExprFactory<NumericConstantFactory, NumericConstantFactory>
      MakeBinaryExprFactory(mp::expr::Kind kind, double lhs, double rhs) {
    return BinaryExprFactory<NumericConstantFactory, NumericConstantFactory>(
          kind, NumericConstantFactory(lhs), NumericConstantFactory(rhs));
  }

  template <typename LHSFactory>
  BinaryExprFactory<LHSFactory, NumericConstantFactory>
      MakeBinaryExprFactory(mp::expr::Kind kind, LHSFactory lhs, double rhs) {
    return BinaryExprFactory<LHSFactory, NumericConstantFactory>(
          kind, lhs, NumericConstantFactory(rhs));
  }

  // An abstract factory for logical expressions.
  class LogicalExprFactory : public ExprFactory {
   protected:
    virtual LogicalExpr Create(ProblemBuilder &pb) const = 0;

   public:
    LogicalExpr operator()(ProblemBuilder &pb) const {
      Init(pb);
      return Create(pb);
    }
  };

  class BinaryLogicalExprFactory : public LogicalExprFactory {
   private:
    mp::expr::Kind kind_;

   public:
    explicit BinaryLogicalExprFactory(mp::expr::Kind kind) : kind_(kind) {}

    LogicalExpr Create(ProblemBuilder &pb) const {
      return pb.MakeBinaryLogical(kind_,
                                  pb.MakeRelational(mp::expr::EQ, x, one),
                                  pb.MakeRelational(mp::expr::EQ, y, one));
    }
  };

  class RelationalExprFactory : public LogicalExprFactory {
   private:
    mp::expr::Kind kind_;

   public:
    explicit RelationalExprFactory(mp::expr::Kind kind) : kind_(kind) {}

    LogicalExpr Create(ProblemBuilder &pb) const {
      return pb.MakeRelational(kind_, x, y);
    }
  };

  class LogicalCountExprFactory : public LogicalExprFactory {
   private:
    mp::expr::Kind kind_;

   public:
    explicit LogicalCountExprFactory(mp::expr::Kind kind) : kind_(kind) {}

    LogicalExpr Create(ProblemBuilder &pb) const {
      return pb.MakeLogicalCount(kind_, x, MakeCount(pb, 2));
    }
  };

  class IteratedLogicalExprFactory : public LogicalExprFactory {
   private:
    mp::expr::Kind kind_;

   public:
    explicit IteratedLogicalExprFactory(mp::expr::Kind kind) : kind_(kind) {}

    LogicalExpr Create(ProblemBuilder &pb) const {
      auto args = pb.BeginIteratedLogical(kind_, 3);
      for (int i = 1; i <= 3; ++i)
        args.AddArg(pb.MakeRelational(mp::expr::EQ, pb.MakeVariable(i), one));
      return pb.EndIteratedLogical(args);
    }
  };

  // Solves a problem containing a single constraint with the given
  // expression and returns the value of the variable with index 0.
  EvalResult Solve(const LogicalExprFactory &factory,
      int var1, int var2, int var3 = 0, bool need_result = false);

  // Constructs and evaluates a numeric expression.
  EvalResult Eval(const NumericExprFactory &factory,
                  int var1 = 0, int var2 = 0, int var3 = 0);

  // Constructs and evaluates a logical expression.
  EvalResult Eval(const LogicalExprFactory &factory,
                  int var1 = 0, int var2 = 0, int var3 = 0) {
    struct Factory : NumericExprFactory {
      const LogicalExprFactory &condition_;
      explicit Factory(const LogicalExprFactory &condition)
        : condition_(condition) {}
      NumericExpr Create(ProblemBuilder &pb) const {
        return pb.MakeIf(condition_(pb),
                         pb.MakeNumericConstant(1), pb.MakeNumericConstant(0));
      }
    } num_factory(factory);
    return Eval(num_factory, var1, var2, var3);
  }

  // Builds and evaluates an unary expression.
  // kind: expression kind
  // arg: value of the argument
  EvalResult EvalUnary(mp::expr::Kind kind, int arg) {
    return Eval(UnaryExprFactory<>(kind, VariableFactory(1)), arg);
  }

  // Builds and evaluates a binary expression.
  // kind: expression kind
  // lhs: value of the left-hand side (first argument)
  // rhs: value of the right-hand side (second argument)
  EvalResult EvalBinary(mp::expr::Kind kind, int lhs, int rhs) {
    return Eval(BinaryExprFactory<>(
                  kind, VariableFactory(1), VariableFactory(2)), lhs, rhs);
  }

  template <typename ArgHandler, typename Arg>
  void AddArgs(ArgHandler &handler, mp::ArrayRef<Arg> args) {
    for (std::size_t i = 0, n = args.size(); i < n; ++i)
      handler.AddArg(args[i]);
  }

  virtual void SetInfo(ProblemBuilder &pb, mp::ProblemInfo &info) {
    pb.SetInfo(info);
  }

 public:
  NLSolverTest() {}

  // Creates a test count expression with arguments var[i] != 0 for
  // i = start_var_index, ..., 3.
  static ProblemBuilder::CountExpr MakeCount(
      ProblemBuilder &pb, int start_var_index) {
    int num_args = 3 - start_var_index + 1;
    auto args = pb.BeginCount(num_args);
    auto zero = pb.MakeNumericConstant(0);
    for (int i = start_var_index; i <= 3; ++i)
      args.AddArg(pb.MakeRelational(mp::expr::NE, pb.MakeVariable(i), zero));
    return pb.EndCount(args);
  }

  static LogicalExpr MakeAllDiff(ProblemBuilder &pb, int start_var_index = 1) {
    int num_args = 3;
    auto args = pb.BeginPairwise(mp::expr::ALLDIFF, num_args);
    for (int i = 0; i < num_args; ++i)
      args.AddArg(pb.MakeVariable(i + start_var_index));
    return pb.EndPairwise(args);
  }
};

void Interrupt();

namespace var = mp::var;
namespace obj = mp::obj;

EvalResult NLSolverTest::Solve(
    const LogicalExprFactory &factory,
    int var1, int var2, int var3, bool need_result) {
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
  pb.AddCon(factory(pb));
  return Solve(pb);
}

// Constructs and evaluates a numeric expression expression.
EvalResult NLSolverTest::Eval(
    const NumericExprFactory &factory, int var1, int var2, int var3) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_nl_integer_vars_in_cons = 4;
  info.num_algebraic_cons = info.num_nl_cons = 1;
  info.num_con_nonzeros = 1;
  SetInfo(pb, info);
  auto inf = std::numeric_limits<double>::infinity();
  pb.AddVar(-inf, inf, mp::var::INTEGER);
  pb.AddVar(var1, var1, mp::var::INTEGER);
  pb.AddVar(var2, var2, mp::var::INTEGER);
  pb.AddVar(var3, var3, mp::var::INTEGER);
  auto cols = pb.GetColumnSizeHandler();
  cols.Add(1);
  for (int i = 1; i < info.num_vars - 1; ++i)
    cols.Add(0);
  auto con_builder = pb.AddCon(factory(pb), 0, 0, 0);
  con_builder.AddTerm(0, -1);
  return Solve(pb);
}

TEST_F(NLSolverTest, NonlinearObj) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  SetInfo(pb);
  pb.AddVar(2, 2, var::INTEGER);
  NumericExpr x = pb.MakeVariable(0);
  pb.AddObj(obj::MIN, pb.MakeBinary(mp::expr::MUL, x, x), 0);
  EXPECT_EQ(4, Solve(pb).obj_value());
}

TEST_F(NLSolverTest, ObjConst) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  SetInfo(pb);
  pb.AddVar(0, 0, var::INTEGER);
  auto obj = pb.AddObj(obj::MIN, pb.MakeNumericConstant(42), 0);
  obj.AddTerm(0, 1);
  EXPECT_EQ(42, Solve(pb).obj_value());
}

TEST_F(NLSolverTest, Minimize) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  SetInfo(pb);
  pb.AddVar(11, 22, var::INTEGER);
  pb.AddObj(obj::MIN, pb.MakeVariable(0), 0);
  EXPECT_EQ(11, Solve(pb).obj_value());
}

TEST_F(NLSolverTest, Maximize) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  SetInfo(pb);
  pb.AddVar(11, 22, var::INTEGER);
  pb.AddObj(obj::MAX, pb.MakeVariable(0), 0);
  EXPECT_EQ(22, Solve(pb).obj_value());
}

#ifndef MP_TSP_SIZE
# define MP_TSP_SIZE 10
#endif

// Creates a test travelling salesman problem.
template <typename ProblemBuilder>
void MakeTSP(ProblemBuilder &pb) {
  int n = MP_TSP_SIZE;
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_linear_binary_vars = n * n;
  info.num_objs = 1;
  info.num_algebraic_cons = 2 * n;
  pb.SetInfo(info);
  typedef typename ProblemBuilder::NumericExpr NumericExpr;
  // Add variables.
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      pb.AddVar(0, 1, mp::var::INTEGER);
  }
  // Add objective and constraints.
  auto obj = pb.AddObj(mp::obj::MIN, NumericExpr(), n * n);
  for (int i = 0; i < n; ++i) {
    auto in_con = pb.AddCon(NumericExpr(), 1, 1, n);   // exactly one incoming
    auto out_con = pb.AddCon(NumericExpr(), 1, 1, n);  // exactly one outgoing
    for (int j = 0; j < n; ++j) {
      // Distance is arbitrarily chosen to be i + j + 1.
      obj.AddTerm(i * n + j, i * j + 1);
      in_con.AddTerm(j * n + i, 1);
      out_con.AddTerm(i * n + j, 1);
    }
  }
}

// ----------------------------------------------------------------------------
// Expression tests

TEST_F(NLSolverTest, NumericConstant) {
  EXPECT_EQ(42, Eval(NumericConstantFactory(42)));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(42, Eval(MakeBinaryExprFactory(mp::expr::MUL, 0.42, 100.0)));
    return;
  }
  EXPECT_THROW_MSG(Eval(NumericConstantFactory(0.42)), mp::Error,
    "value 0.42 can't be represented as int");
}

TEST_F(NLSolverTest, Variable) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_nl_integer_vars_in_objs = 1;
  info.num_objs = info.num_nl_objs = 1;
  pb.SetInfo(info);
  pb.AddVar(42, 43, mp::var::INTEGER);
  pb.AddObj(mp::obj::MIN, pb.MakeVariable(0), 0);
  EXPECT_EQ(42, Solve(pb));
}

TEST_F(NLSolverTest, Floor) {
  auto kind = mp::expr::FLOOR;
  EXPECT_EQ(-42, EvalUnary(kind, -42));
  EXPECT_EQ(42, EvalUnary(kind, 42));
  if (!HasFeature(feature::FLOAT_CONST)) return;
  typedef UnaryExprFactory<NumericConstantFactory> Factory;
  EXPECT_EQ(4, Eval(Factory(kind, 4.9)));
  EXPECT_EQ(-5, Eval(Factory(kind, -4.1)));
}

TEST_F(NLSolverTest, Ceil) {
  auto kind = mp::expr::CEIL;
  EXPECT_EQ(-42, EvalUnary(kind, -42));
  EXPECT_EQ(42, EvalUnary(kind, 42));
  if (!HasFeature(feature::FLOAT_CONST)) return;
  typedef UnaryExprFactory<NumericConstantFactory> Factory;
  EXPECT_EQ(5, Eval(Factory(kind, 4.1)));
  EXPECT_EQ(-4, Eval(Factory(kind, -4.9)));
}

TEST_F(NLSolverTest, Abs) {
  auto kind = mp::expr::ABS;
  EXPECT_EQ(42, EvalUnary(kind, -42));
  EXPECT_EQ(42, EvalUnary(kind, 42));
}

TEST_F(NLSolverTest, Minus) {
  auto kind = mp::expr::MINUS;
  EXPECT_EQ(42, EvalUnary(kind, -42));
  EXPECT_EQ(-42, EvalUnary(kind, 42));
}

TEST_F(NLSolverTest, Tanh) {
  double arg = (std::log(1.5) - std::log(0.5)) / 2;
  auto factory = MakeBinaryExprFactory(
        mp::expr::MUL, MakeUnaryExprFactory(mp::expr::TANH, arg), 2);
  if (!HasFeature(feature::HYPERBOLIC))
    EXPECT_THROW_MSG(Eval(factory), mp::Error, "unsupported: tanh");
  else
    EXPECT_EQ(1, Eval(factory));
}

TEST_F(NLSolverTest, Tan) {
  auto factory = MakeUnaryExprFactory(mp::expr::TAN, std::atan(1.0));
  if (!HasFeature(feature::TRIGONOMETRIC))
    EXPECT_THROW(EvalUnary(mp::expr::TAN, 0), mp::UnsupportedError);
  else
    EXPECT_EQ(1, Eval(factory));
}

TEST_F(NLSolverTest, Sqrt) {
  auto kind = mp::expr::SQRT;
  if (!HasFeature(feature::SQRT))
    EXPECT_THROW_MSG(EvalUnary(kind, 64), mp::Error, "unsupported: sqrt");
  else
    EXPECT_EQ(8, EvalUnary(kind, 64));
}

TEST_F(NLSolverTest, Sinh) {
  auto factory = MakeUnaryExprFactory(
        mp::expr::SINH, std::log(2 + std::sqrt(5.0)));
  if (!HasFeature(feature::HYPERBOLIC))
    EXPECT_THROW_MSG(Eval(factory), mp::Error, "unsupported: sinh");
  else
    EXPECT_EQ(2, Eval(factory));
}

TEST_F(NLSolverTest, Sin) {
  auto factory = MakeUnaryExprFactory(mp::expr::SIN, std::asin(1.0));
  if (!HasFeature(feature::TRIGONOMETRIC))
    EXPECT_THROW(EvalUnary(mp::expr::SIN, 0), mp::UnsupportedError);
  else
    EXPECT_EQ(1, Eval(factory));
}

TEST_F(NLSolverTest, Log10) {
  auto factory = MakeUnaryExprFactory(mp::expr::LOG10, 1000.0);
  if (!HasFeature(feature::LOG))
    EXPECT_THROW_MSG(Eval(factory), mp::Error, "unsupported: log10");
  else
    EXPECT_EQ(3, Eval(factory));
}

TEST_F(NLSolverTest, Log) {
  auto factory = MakeUnaryExprFactory(mp::expr::LOG, std::exp(5.0));
  if (!HasFeature(feature::LOG))
    EXPECT_THROW_MSG(Eval(factory), mp::Error, "unsupported: log");
  else
    EXPECT_EQ(5, Eval(factory));
}

TEST_F(NLSolverTest, Exp) {
  auto factory = MakeUnaryExprFactory(mp::expr::EXP, std::log(5.0));
  if (!HasFeature(feature::EXP))
    EXPECT_THROW_MSG(Eval(factory), mp::Error, "unsupported: exp");
  else
    EXPECT_EQ(5, Eval(factory));
}

TEST_F(NLSolverTest, Cosh) {
  double x = 5;
  auto factory = MakeUnaryExprFactory(
        mp::expr::COSH, std::log(x + std::sqrt(x + 1) * std::sqrt(x - 1)));
  if (!HasFeature(feature::HYPERBOLIC))
    EXPECT_THROW_MSG(Eval(factory), mp::Error, "unsupported: cosh");
  else
    EXPECT_EQ(5, Eval(factory));
}

TEST_F(NLSolverTest, Cos) {
  auto factory = MakeUnaryExprFactory(mp::expr::COS, 0.0);
  if (!HasFeature(feature::TRIGONOMETRIC))
    EXPECT_THROW(EvalUnary(mp::expr::COS, 0), mp::UnsupportedError);
  else
    EXPECT_EQ(1, Eval(factory));
}

TEST_F(NLSolverTest, Atanh) {
  auto atanh = MakeUnaryExprFactory(mp::expr::ATANH, std::tanh(5.0));
  if (!HasFeature(feature::HYPERBOLIC)) {
    EXPECT_THROW_MSG(Eval(atanh), mp::Error, "unsupported: atanh");
    return;
  }
  auto factory = MakeUnaryExprFactory(
        mp::expr::FLOOR, MakeBinaryExprFactory(
          mp::expr::ADD,
          MakeBinaryExprFactory(mp::expr::MUL, atanh, 1000000.0), 0.5));
  EXPECT_EQ(5000000, Eval(factory));
}

TEST_F(NLSolverTest, Atan) {
  EXPECT_THROW(EvalUnary(mp::expr::ATAN, 0), mp::UnsupportedError);
}

TEST_F(NLSolverTest, Asinh) {
  auto factory = MakeUnaryExprFactory(mp::expr::ASINH, std::sinh(5.0));
  if (!HasFeature(feature::HYPERBOLIC))
    EXPECT_THROW_MSG(Eval(factory), mp::Error, "unsupported: asinh");
  else
    EXPECT_EQ(5, Eval(factory));
}

TEST_F(NLSolverTest, Asin) {
  EXPECT_THROW(EvalUnary(mp::expr::ASIN, 0), mp::UnsupportedError);
}

TEST_F(NLSolverTest, Acosh) {
  auto factory = MakeUnaryExprFactory(mp::expr::ACOSH, std::cosh(5.0));
  if (!HasFeature(feature::HYPERBOLIC))
    EXPECT_THROW_MSG(Eval(factory), mp::Error, "unsupported: acosh");
  else
    EXPECT_EQ(5, Eval(factory));
}

TEST_F(NLSolverTest, Acos) {
  EXPECT_THROW(EvalUnary(mp::expr::ACOS, 0), mp::UnsupportedError);
}

TEST_F(NLSolverTest, Pow2) {
  EXPECT_EQ(49, EvalUnary(mp::expr::POW2, 7));
}

TEST_F(NLSolverTest, Add) {
  auto kind = mp::expr::ADD;
  EXPECT_EQ(25, EvalBinary(kind, 10, 15));
  EXPECT_EQ(12, EvalBinary(kind, 19, -7));
}

TEST_F(NLSolverTest, Sub) {
  auto kind = mp::expr::SUB;
  EXPECT_EQ(-5, EvalBinary(kind, 10, 15));
  EXPECT_EQ(26, EvalBinary(kind, 19, -7));
}

TEST_F(NLSolverTest, Mul) {
  auto kind = mp::expr::MUL;
  EXPECT_EQ(150, EvalBinary(kind, 10, 15));
  EXPECT_EQ(-133, EvalBinary(kind, 19, -7));
}

TEST_F(NLSolverTest, Div) {
  auto kind = mp::expr::DIV;
  if (!HasFeature(feature::DIV)) {
    EXPECT_THROW_MSG(EvalBinary(kind, 150, 15), mp::Error, "unsupported: /");
    return;
  }
  EXPECT_EQ(10, EvalBinary(kind, 150, 15));
  EXPECT_EQ(-7, EvalBinary(kind, -133, 19));
}

TEST_F(NLSolverTest, IntDiv) {
  auto kind = mp::expr::INT_DIV;
  EXPECT_EQ( 3, EvalBinary(kind,  9,  3));
  EXPECT_EQ( 2, EvalBinary(kind,  8,  3));
  EXPECT_EQ(-2, EvalBinary(kind, -8,  3));
  EXPECT_EQ(-2, EvalBinary(kind,  8, -3));
  EXPECT_EQ( 2, EvalBinary(kind, -8, -3));
}

TEST_F(NLSolverTest, Mod) {
  auto kind = mp::expr::MOD;
  EXPECT_EQ( 0, EvalBinary(kind,  9,  3));
  EXPECT_EQ( 2, EvalBinary(kind,  8,  3));
  EXPECT_EQ(-2, EvalBinary(kind, -8,  3));
  EXPECT_EQ( 2, EvalBinary(kind,  8, -3));
  EXPECT_EQ(-2, EvalBinary(kind, -8, -3));
}

TEST_F(NLSolverTest, Pow) {
  auto kind = mp::expr::POW;
  if (!HasFeature(feature::POW)) {
    EXPECT_THROW_MSG(EvalBinary(kind, 2, 3), mp::Error, "unsupported: ^");
    return;
  }
  EXPECT_EQ(8, EvalBinary(kind, 2, 3));
  EXPECT_EQ(81, EvalBinary(kind, 3, 4));
}

TEST_F(NLSolverTest, PowConstBase) {
  struct Factory : NumericExprFactory {
    NumericExpr Create(ProblemBuilder &pb) const {
      return pb.MakeBinary(mp::expr::POW_CONST_BASE,
                           pb.MakeNumericConstant(5), x);
    }
  } factory;
  if (!HasFeature(feature::POW))
    EXPECT_THROW_MSG(Eval(factory, 3), mp::Error, "unsupported: ^");
  else
    EXPECT_EQ(125, Eval(factory, 3));
}

TEST_F(NLSolverTest, PowConstExp) {
  auto factory = MakeBinaryExprFactory(
        mp::expr::POW_CONST_EXP, VariableFactory(1), 4.0);
  EXPECT_EQ(16, Eval(factory, 2));
}

TEST_F(NLSolverTest, Less) {
  auto kind = mp::expr::LESS;
  EXPECT_EQ(0, EvalBinary(kind, 10, 15));
  EXPECT_EQ(26, EvalBinary(kind, 19, -7));
}

TEST_F(NLSolverTest, Atan2) {
  EXPECT_THROW(EvalBinary(mp::expr::ATAN2, 0, 0), mp::UnsupportedError);
}

TEST_F(NLSolverTest, Precision) {
  EXPECT_THROW(EvalBinary(mp::expr::PRECISION, 0, 0), mp::UnsupportedError);
}

TEST_F(NLSolverTest, Round) {
  auto kind = mp::expr::ROUND;
  EXPECT_EQ(42, Eval(MakeBinaryExprFactory(kind, VariableFactory(1), 0.0), 42));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(4, Eval(MakeBinaryExprFactory(kind, 4.4, 0.0)));
    EXPECT_EQ(5, Eval(MakeBinaryExprFactory(kind, 4.6, 0.0)));
    EXPECT_EQ(-4, Eval(MakeBinaryExprFactory(kind, -4.4, 0.0)));
    EXPECT_EQ(-5, Eval(MakeBinaryExprFactory(kind, -4.6, 0.0)));
  }
  const char *message = "unsupported: round with nonzero second parameter";
  EXPECT_THROW_MSG(
        Eval(MakeBinaryExprFactory(kind, VariableFactory(1), 1.0), 0),
        mp::Error, message);
  EXPECT_THROW_MSG(EvalBinary(kind, 0, 0), mp::Error, message);
}

TEST_F(NLSolverTest, Trunc) {
  auto kind = mp::expr::TRUNC;
  EXPECT_EQ(42, Eval(MakeBinaryExprFactory(kind, VariableFactory(1), 0.0), 42));
  if (HasFeature(feature::FLOAT_CONST)) {
    EXPECT_EQ(4, Eval(MakeBinaryExprFactory(kind, 4.4, 0.0)));
    EXPECT_EQ(4, Eval(MakeBinaryExprFactory(kind, 4.6, 0.0)));
    EXPECT_EQ(-4, Eval(MakeBinaryExprFactory(kind, -4.4, 0.0)));
    EXPECT_EQ(-4, Eval(MakeBinaryExprFactory(kind, -4.6, 0.0)));
  }
  const char *message = "unsupported: trunc with nonzero second parameter";
  EXPECT_THROW_MSG(
        Eval(MakeBinaryExprFactory(kind, VariableFactory(1), 1.0), 0),
        mp::Error, message);
  EXPECT_THROW_MSG(EvalBinary(kind, 0, 0), mp::Error, message);
}

TEST_F(NLSolverTest, If) {
  struct IfFactory : NumericExprFactory {
    NumericExpr Create(ProblemBuilder &pb) const {
      return pb.MakeIf(pb.MakeRelational(mp::expr::EQ, x, one), y, z);
    }
  } factory;
  EXPECT_EQ(42, Eval(factory, 1, 42, 10));
  EXPECT_EQ(10, Eval(factory, 0, 42, 10));
  EXPECT_EQ(42, Eval(factory, 1, 42, 42));
}

TEST_F(NLSolverTest, PLTerm) {
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
  struct IfFactory : NumericExprFactory {
    NumericExpr Create(ProblemBuilder &pb) const {
      auto plterm = pb.BeginPLTerm(2);
      plterm.AddSlope(-1);
      plterm.AddBreakpoint(3);
      plterm.AddSlope(0);
      plterm.AddBreakpoint(6);
      plterm.AddSlope(1);
      return pb.EndPLTerm(plterm, x);
    }
  } factory;
  if (!HasFeature(feature::PLTERM)) {
    EXPECT_THROW_MSG(Eval(factory, 42), mp::Error,
                     "unsupported: piecewise-linear term");
    return;
  }
  EXPECT_EQ(33, Eval(factory, 42));
  EXPECT_EQ(-3, Eval(factory, 4));
  EXPECT_EQ(1, Eval(factory, -1));
}

TEST_F(NLSolverTest, Call) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_linear_integer_vars = 1;
  info.num_objs = info.num_nl_objs = 1;
  info.num_funcs = 1;
  pb.SetInfo(info);
  pb.AddVar(0, 0, mp::var::INTEGER);
  try {
    pb.SetFunction(0, "foo", 2, mp::func::NUMERIC);
  } catch (const mp::UnsupportedError &) {
    // SetFunction may throw UnsupportedError.
    return;
  } catch (const mp::Error &e) {
    EXPECT_STREQ("function foo not available", e.what());
    return;
  }
  EXPECT_THROW({
      pb.SetFunction(0, "foo", 2, mp::func::NUMERIC);
      auto call = pb.BeginCall(0, 2);
      pb.AddObj(mp::obj::MIN, pb.EndCall(call), 0);
      Solve(pb);
  }, mp::UnsupportedError);
}

TEST_F(NLSolverTest, Min) {
  auto factory = VarArgFactory(mp::expr::MIN);
  EXPECT_EQ(-7, Eval(factory, 3, -7, 5));
  EXPECT_EQ(10, Eval(factory, 10, 20, 30));
}

TEST_F(NLSolverTest, Max) {
  auto factory = VarArgFactory(mp::expr::MAX);
  EXPECT_EQ(5, Eval(factory, 3, -7, 5));
  EXPECT_EQ(30, Eval(factory, 30, 20, 10));
}

TEST_F(NLSolverTest, Sum) {
  struct SumFactory : NumericExprFactory {
    int num_args_;
    explicit SumFactory(int num_args) : num_args_(num_args) {}
    NumericExpr Create(ProblemBuilder &pb) const {
      auto args = pb.BeginSum(num_args_);
      for (int i = 1; i <= num_args_; ++i)
        args.AddArg(pb.MakeVariable(i));
      return pb.EndSum(args);
    }
  };
  EXPECT_EQ(0, Eval(SumFactory(0), 100, 20, 3));
  EXPECT_EQ(100, Eval(SumFactory(1), 100, 20, 3));
  EXPECT_EQ(123, Eval(SumFactory(3), 100, 20, 3));
}

TEST_F(NLSolverTest, Count) {
  class CountFactory : public NumericExprFactory {
    NumericExpr Create(ProblemBuilder &pb) const { return MakeCount(pb, 1); }
  } factory;
  EXPECT_EQ(0, Eval(factory));
  EXPECT_EQ(1, Eval(factory, 1));
  EXPECT_EQ(2, Eval(factory, 0, 1, 1));
  EXPECT_EQ(3, Eval(factory, 1, 1, 1));
}

TEST_F(NLSolverTest, NumberOf) {
  struct NumberOfFactory : NumericExprFactory {
    int num_args_;
    explicit NumberOfFactory(int num_args) : num_args_(num_args) {}
    NumericExpr Create(ProblemBuilder &pb) const {
      auto value = pb.MakeNumericConstant(42);
      auto args = pb.BeginNumberOf(value, num_args_);
      for (int i = 1; i <= num_args_; ++i)
        args.AddArg(pb.MakeVariable(i));
      return pb.EndNumberOf(args);
    }
  } factory(1);
  EXPECT_EQ(0, Eval(factory));
  EXPECT_EQ(1, Eval(factory, 42));
  factory = NumberOfFactory(2);
  EXPECT_EQ(0, Eval(factory));
  EXPECT_EQ(1, Eval(factory, 0, 42));
  EXPECT_EQ(2, Eval(factory, 42, 42));
}

TEST_F(NLSolverTest, LogicalConstant) {
  struct LogicalConstantFactory : LogicalExprFactory {
    bool value_;
    explicit LogicalConstantFactory(bool value) : value_(value) {}
    LogicalExpr Create(ProblemBuilder &pb) const {
      return pb.MakeLogicalConstant(value_);
    }
  };
  EXPECT_EQ(0, Eval(LogicalConstantFactory(false)));
  EXPECT_EQ(1, Eval(LogicalConstantFactory(true)));
}

TEST_F(NLSolverTest, Not) {
  struct NotFactory : LogicalExprFactory {
    LogicalExpr Create(ProblemBuilder &pb) const {
      return pb.MakeNot(pb.MakeRelational(mp::expr::EQ, x, one));
    }
  } factory;
  EXPECT_EQ(1, Eval(factory, 0));
  EXPECT_EQ(0, Eval(factory, 1));
}

TEST_F(NLSolverTest, Or) {
  BinaryLogicalExprFactory factory(mp::expr::OR);
  EXPECT_EQ(0, Eval(factory, 0, 0));
  EXPECT_EQ(1, Eval(factory, 0, 1));
  EXPECT_EQ(1, Eval(factory, 1, 0));
  EXPECT_EQ(1, Eval(factory, 1, 1));
}

TEST_F(NLSolverTest, And) {
  BinaryLogicalExprFactory factory(mp::expr::AND);
  EXPECT_EQ(0, Eval(factory, 0, 0));
  EXPECT_EQ(0, Eval(factory, 0, 1));
  EXPECT_EQ(0, Eval(factory, 1, 0));
  EXPECT_EQ(1, Eval(factory, 1, 1));
}

TEST_F(NLSolverTest, Iff) {
  BinaryLogicalExprFactory factory(mp::expr::IFF);
  EXPECT_EQ(1, Eval(factory, 0, 0));
  EXPECT_EQ(0, Eval(factory, 0, 1));
  EXPECT_EQ(0, Eval(factory, 1, 0));
  EXPECT_EQ(1, Eval(factory, 1, 1));
}

TEST_F(NLSolverTest, LT) {
  RelationalExprFactory factory(mp::expr::LT);
  EXPECT_EQ(0, Eval(factory, 3, 3));
  EXPECT_EQ(1, Eval(factory, 3, 5));
  EXPECT_EQ(0, Eval(factory, 5, 3));
}

TEST_F(NLSolverTest, LE) {
  RelationalExprFactory factory(mp::expr::LE);
  EXPECT_EQ(1, Eval(factory, 3, 3));
  EXPECT_EQ(1, Eval(factory, 3, 5));
  EXPECT_EQ(0, Eval(factory, 5, 3));
}

TEST_F(NLSolverTest, EQ) {
  RelationalExprFactory factory(mp::expr::EQ);
  EXPECT_EQ(1, Eval(factory, 3, 3));
  EXPECT_EQ(0, Eval(factory, 3, 5));
  EXPECT_EQ(0, Eval(factory, 5, 3));
}

TEST_F(NLSolverTest, GE) {
  RelationalExprFactory factory(mp::expr::GE);
  EXPECT_EQ(1, Eval(factory, 3, 3));
  EXPECT_EQ(0, Eval(factory, 3, 5));
  EXPECT_EQ(1, Eval(factory, 5, 3));
}

TEST_F(NLSolverTest, GT) {
  RelationalExprFactory factory(mp::expr::GT);
  EXPECT_EQ(0, Eval(factory, 3, 3));
  EXPECT_EQ(0, Eval(factory, 3, 5));
  EXPECT_EQ(1, Eval(factory, 5, 3));
}

TEST_F(NLSolverTest, NE) {
  RelationalExprFactory factory(mp::expr::NE);
  EXPECT_EQ(0, Eval(factory, 3, 3));
  EXPECT_EQ(1, Eval(factory, 3, 5));
  EXPECT_EQ(1, Eval(factory, 5, 3));
}

TEST_F(NLSolverTest, AtLeast) {
  LogicalCountExprFactory factory(mp::expr::ATLEAST);
  EXPECT_EQ(1, Eval(factory, 0, 0, 0));
  EXPECT_EQ(1, Eval(factory, 0, 1, 0));
  EXPECT_EQ(0, Eval(factory, 1, 0, 0));
  EXPECT_EQ(1, Eval(factory, 1, 0, 1));
  EXPECT_EQ(1, Eval(factory, 1, 1, 1));
  EXPECT_EQ(0, Eval(factory, 2, 0, 1));
  EXPECT_EQ(1, Eval(factory, 2, 1, 1));
}

TEST_F(NLSolverTest, AtMost) {
  LogicalCountExprFactory factory(mp::expr::ATMOST);
  EXPECT_EQ(1, Eval(factory, 0, 0, 0));
  EXPECT_EQ(0, Eval(factory, 0, 1, 0));
  EXPECT_EQ(1, Eval(factory, 1, 0, 0));
  EXPECT_EQ(1, Eval(factory, 1, 0, 1));
  EXPECT_EQ(0, Eval(factory, 1, 1, 1));
  EXPECT_EQ(1, Eval(factory, 2, 0, 1));
  EXPECT_EQ(1, Eval(factory, 2, 1, 1));
}

TEST_F(NLSolverTest, Exactly) {
  LogicalCountExprFactory factory(mp::expr::EXACTLY);
  EXPECT_EQ(1, Eval(factory, 0, 0, 0));
  EXPECT_EQ(0, Eval(factory, 0, 1, 0));
  EXPECT_EQ(0, Eval(factory, 1, 0, 0));
  EXPECT_EQ(1, Eval(factory, 1, 0, 1));
  EXPECT_EQ(0, Eval(factory, 1, 1, 1));
  EXPECT_EQ(0, Eval(factory, 2, 0, 1));
  EXPECT_EQ(1, Eval(factory, 2, 1, 1));
}

TEST_F(NLSolverTest, NotAtLeast) {
  LogicalCountExprFactory factory(mp::expr::NOT_ATLEAST);
  EXPECT_EQ(0, Eval(factory, 0, 0, 0));
  EXPECT_EQ(0, Eval(factory, 0, 1, 0));
  EXPECT_EQ(1, Eval(factory, 1, 0, 0));
  EXPECT_EQ(0, Eval(factory, 1, 0, 1));
  EXPECT_EQ(0, Eval(factory, 1, 1, 1));
  EXPECT_EQ(1, Eval(factory, 2, 0, 1));
  EXPECT_EQ(0, Eval(factory, 2, 1, 1));
}

TEST_F(NLSolverTest, NotAtMost) {
  LogicalCountExprFactory factory(mp::expr::NOT_ATMOST);
  EXPECT_EQ(0, Eval(factory, 0, 0, 0));
  EXPECT_EQ(1, Eval(factory, 0, 1, 0));
  EXPECT_EQ(0, Eval(factory, 1, 0, 0));
  EXPECT_EQ(0, Eval(factory, 1, 0, 1));
  EXPECT_EQ(1, Eval(factory, 1, 1, 1));
  EXPECT_EQ(0, Eval(factory, 2, 0, 1));
  EXPECT_EQ(0, Eval(factory, 2, 1, 1));
}

TEST_F(NLSolverTest, NotExactly) {
  LogicalCountExprFactory factory(mp::expr::NOT_EXACTLY);
  EXPECT_EQ(0, Eval(factory, 0, 0, 0));
  EXPECT_EQ(1, Eval(factory, 0, 1, 0));
  EXPECT_EQ(1, Eval(factory, 1, 0, 0));
  EXPECT_EQ(0, Eval(factory, 1, 0, 1));
  EXPECT_EQ(1, Eval(factory, 1, 1, 1));
  EXPECT_EQ(1, Eval(factory, 2, 0, 1));
  EXPECT_EQ(0, Eval(factory, 2, 1, 1));
}

TEST_F(NLSolverTest, ForAll) {
  IteratedLogicalExprFactory factory(mp::expr::FORALL);
  EXPECT_EQ(0, Eval(factory, 0, 0, 0));
  EXPECT_EQ(0, Eval(factory, 0, 0, 1));
  EXPECT_EQ(0, Eval(factory, 0, 1, 0));
  EXPECT_EQ(0, Eval(factory, 0, 1, 1));
  EXPECT_EQ(0, Eval(factory, 1, 0, 0));
  EXPECT_EQ(0, Eval(factory, 1, 0, 1));
  EXPECT_EQ(0, Eval(factory, 1, 1, 0));
  EXPECT_EQ(1, Eval(factory, 1, 1, 1));
}

TEST_F(NLSolverTest, Exists) {
  IteratedLogicalExprFactory factory(mp::expr::EXISTS);
  EXPECT_EQ(0, Eval(factory, 0, 0, 0));
  EXPECT_EQ(1, Eval(factory, 0, 0, 1));
  EXPECT_EQ(1, Eval(factory, 0, 1, 0));
  EXPECT_EQ(1, Eval(factory, 0, 1, 1));
  EXPECT_EQ(1, Eval(factory, 1, 0, 0));
  EXPECT_EQ(1, Eval(factory, 1, 0, 1));
  EXPECT_EQ(1, Eval(factory, 1, 1, 0));
  EXPECT_EQ(1, Eval(factory, 1, 1, 1));
}

TEST_F(NLSolverTest, Implication) {
  struct ImplicationFactory : LogicalExprFactory {
    LogicalExpr Create(ProblemBuilder &pb) const {
      return pb.MakeImplication(
            pb.MakeRelational(mp::expr::EQ, x, one),
            pb.MakeRelational(mp::expr::EQ, y, one),
            pb.MakeRelational(mp::expr::EQ, z, one));
    }
  } factory;
  EXPECT_EQ(0, Eval(factory, 0, 0, 0));
  EXPECT_EQ(1, Eval(factory, 0, 0, 1));
  EXPECT_EQ(0, Eval(factory, 0, 1, 0));
  EXPECT_EQ(1, Eval(factory, 0, 1, 1));
  EXPECT_EQ(0, Eval(factory, 1, 0, 0));
  EXPECT_EQ(0, Eval(factory, 1, 0, 1));
  EXPECT_EQ(1, Eval(factory, 1, 1, 0));
  EXPECT_EQ(1, Eval(factory, 1, 1, 1));
}

TEST_F(NLSolverTest, AllDiff) {
  struct AllDiffFactory : LogicalExprFactory {
    LogicalExpr Create(ProblemBuilder &pb) const { return MakeAllDiff(pb); }
  } alldiff;
  EXPECT_TRUE(Solve(alldiff, 2, 1, 3).has_value());
  EXPECT_FALSE(Solve(alldiff, 1, 2, 1).has_value());
  EXPECT_FALSE(Solve(alldiff, 1, 1, 1).has_value());
}

TEST_F(NLSolverTest, NestedAllDiff) {
  struct Factory : LogicalExprFactory {
    LogicalExpr Create(ProblemBuilder &pb) const {
      return pb.MakeNot(MakeAllDiff(pb));
    }
  } factory;
  EXPECT_FALSE(Solve(factory, 2, 1, 3).has_value());
  EXPECT_TRUE(Solve(factory, 2, 1, 1).has_value());
}

// ----------------------------------------------------------------------------
// Solve code tests

TEST_F(NLSolverTest, OptimalSolveCode) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_linear_integer_vars = info.num_objs = 1;
  info.num_obj_nonzeros = 1;
  pb.SetInfo(info);
  pb.GetColumnSizeHandler();
  pb.AddVar(0, 0, mp::var::INTEGER);
  pb.AddObj(obj::MIN, NumericExpr(), 1).AddTerm(0, 1);
  EXPECT_EQ(mp::sol::SOLVED, Solve(pb).solve_code());
}

TEST_F(NLSolverTest, FeasibleSolveCode) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_linear_integer_vars = info.num_algebraic_cons = 1;
  pb.SetInfo(info);
  pb.AddVar(0, 0, mp::var::INTEGER);
  pb.AddCon(NumericExpr(), 0, 1, 1).AddTerm(0, 1);
  EXPECT_EQ(mp::sol::SOLVED, Solve(pb).solve_code());
}

TEST_F(NLSolverTest, InfeasibleSolveCode) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_linear_integer_vars = info.num_algebraic_cons = 1;
  pb.SetInfo(info);
  pb.AddVar(0, 0, mp::var::INTEGER);
  pb.AddCon(NumericExpr(), 1, 1, 1).AddTerm(0, 1);
  EXPECT_EQ(mp::sol::INFEASIBLE, Solve(pb).solve_code());
}

// ----------------------------------------------------------------------------
// Interrupt tests

#if MP_THREAD || !MP_THREAD_INTERRUPT

class TestInterrupterBase : public mp::Interrupter {
 private:
  Solver &solver_;

 public:
  explicit TestInterrupterBase(Solver &s) : solver_(s) {
    s.set_interrupter(this);
  }

  ~TestInterrupterBase() { solver_.set_interrupter(0); }

  bool Stop() const { return true; }

  void SetHandler(mp::InterruptHandler handler, void *data) { handler(data); }
};

#if !MP_THREAD_INTERRUPT
typedef TestInterrupterBase TestInterrupter;
#else
class TestInterrupter : public TestInterrupterBase {
 private:
  std::thread thread_;

  static void Run(mp::InterruptHandler handler, void *data) {
    while (!handler(data))
      std::this_thread::yield();
  }

 public:
  explicit TestInterrupter(Solver &s) : TestInterrupterBase(s) {}

  ~TestInterrupter() {
    if (thread_.joinable())
      thread_.join();
  }

  void SetHandler(mp::InterruptHandler handler, void *data) {
    thread_ = std::thread(Run, handler, data);
  }
};
#endif

TEST_F(NLSolverTest, Interrupt) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeTSP(pb);
  pb.EndBuild();
  TestInterrupter interrupter(solver_);
  TestSolutionHandler sh;
  solver_.Solve(pb, sh);
  EXPECT_EQ(600, sh.status());
  EXPECT_TRUE(sh.message().find("interrupted") != std::string::npos);
}

TEST_F(NLSolverTest, InitialValues) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeTSP(pb);
  for (int i = 0; i < MP_TSP_SIZE; ++i) {
    for (int j = 0; j < MP_TSP_SIZE; ++j)
      pb.SetInitialValue(i * MP_TSP_SIZE + j, i + j == MP_TSP_SIZE - 1 ? 1 : 0);
  }
  pb.EndBuild();
  TestInterrupter interrupter(solver_);
  TestSolutionHandler sh(MP_TSP_SIZE * MP_TSP_SIZE);
  solver_.Solve(pb, sh);
  const double *primal = sh.primal();
  if (!HasFeature(feature::INITIAL_VALUES)) {
    ASSERT_TRUE(primal == 0);
    return;
  }
  ASSERT_TRUE(primal != 0);
  for (int i = 0; i < MP_TSP_SIZE; ++i) {
    for (int j = 0; j < MP_TSP_SIZE; ++j)
      EXPECT_EQ(i + j == MP_TSP_SIZE - 1 ? 1 : 0, primal[i * MP_TSP_SIZE + j]);
  }
}

#endif

struct SolutionCounter : TestSolutionHandler {
  int num_solutions;
  SolutionCounter() : num_solutions(0) {}
  void HandleFeasibleSolution(
      fmt::StringRef, const double *, const double *, double) {
    ++num_solutions;
  }
};

// Makes a test problem with a single constraint alldiff(x0, x1, x2).
template <typename ProblemBuilder>
void MakeAllDiffProblem(ProblemBuilder &pb) {
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_nl_integer_vars_in_cons = 3;
  info.num_logical_cons = 1;
  pb.SetInfo(info);
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddVar(1, 3, var::INTEGER);
  pb.AddCon(NLSolverTest::MakeAllDiff(pb, 0));
}

TEST_F(NLSolverTest, CountSolutions) {
  if (!solver_.FindOption("solutionlimit"))
    return;
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeAllDiffProblem(pb);
  solver_.SetIntOption("solutionlimit", 10);
  solver_.SetIntOption("countsolutions", 1);
  SolutionCounter sc;
  solver_.Solve(pb, sc);
  EXPECT_EQ(6, sc.num_solutions);
}

TEST_F(NLSolverTest, SolutionLimit) {
  if (!solver_.FindOption("solutionlimit"))
    return;
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeAllDiffProblem(pb);
  solver_.SetIntOption("solutionlimit", 5);
  solver_.SetIntOption("countsolutions", 1);
  SolutionCounter sc;
  solver_.Solve(pb, sc);
  EXPECT_EQ(0, sc.status());
  EXPECT_EQ(5, sc.num_solutions);
}

TEST_F(NLSolverTest, MultipleSolutions) {
  if (!solver_.FindOption("solutionlimit"))
    return;
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  MakeAllDiffProblem(pb);
  solver_.SetIntOption("solutionlimit", 3);
  solver_.SetStrOption("solutionstub", "test");
  SolutionCounter sc;
  solver_.Solve(pb, sc);
  EXPECT_EQ(3, sc.num_solutions);
}

TEST_F(NLSolverTest, TimingOption) {
  struct TestOutputHandler : mp::OutputHandler {
    std::string output;

    virtual ~TestOutputHandler() {}
    void HandleOutput(fmt::StringRef output) {
      this->output += output;
    }
  };

  mp::BasicSolutionHandler sol_handler;

  {
    ProblemBuilder pb(solver_.GetProblemBuilder(""));
    MakeAllDiffProblem(pb);

    TestOutputHandler oh;
    solver_.set_output_handler(&oh);

    solver_.SetIntOption("timing", 0);
    solver_.Solve(pb, sol_handler);
    EXPECT_TRUE(oh.output.find("Setup time = ") == std::string::npos);
    EXPECT_TRUE(oh.output.find("Solution time = ") == std::string::npos);
    EXPECT_TRUE(oh.output.find("Output time = ") == std::string::npos);
  }

  {
    ProblemBuilder pb(solver_.GetProblemBuilder(""));
    MakeAllDiffProblem(pb);

    TestOutputHandler oh;
    solver_.set_output_handler(&oh);

    solver_.SetIntOption("timing", 1);
    solver_.Solve(pb, sol_handler);
    EXPECT_TRUE(oh.output.find("Setup time = ") != std::string::npos);
    EXPECT_TRUE(oh.output.find("Solution time = ") != std::string::npos);
    EXPECT_TRUE(oh.output.find("Output time = ") != std::string::npos);
  }
}

TEST_F(NLSolverTest, OptionValues) {
  for (auto i = solver_.option_begin(), e = solver_.option_end(); i != e; ++i) {
    for (mp::ValueArrayRef::iterator j = i->values().begin(),
        value_end = i->values().end(); j != value_end; ++j) {
      EXPECT_TRUE(j->value != 0);
    }
  }
}

// Solver should not use objno.
TEST_F(NLSolverTest, IgnoreObjNo) {
  solver_.SetIntOption("objno", 10);
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  SetInfo(pb);
  pb.AddVar(2, 2, var::INTEGER);
  NumericExpr x = pb.MakeVariable(0);
  pb.AddObj(obj::MIN, pb.MakeBinary(mp::expr::MUL, x, x), 0);
  EXPECT_EQ(4, Solve(pb).obj_value());
}

TEST_F(NLSolverTest, CreateSolver) {
  EXPECT_STREQ(solver_.name(), mp::CreateSolver(0)->name());
}

// Makes a problem for testing multiple objectives support and solves it.
template <typename Solver>
EvalResult SolveMultiObjTestProblem(Solver &solver, bool multiobj = false) {
  typename Solver::ProblemBuilder pb(solver.GetProblemBuilder(""));
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_nl_integer_vars_in_objs = 1;
  info.num_objs = info.num_nl_objs = multiobj ? 2 : 1;
  pb.SetInfo(info);
  pb.AddVar(0, 2, var::INTEGER);
  auto var = pb.MakeVariable(0);
  pb.AddObj(mp::obj::MIN, pb.MakeBinary(
              mp::expr::MOD, var, pb.MakeNumericConstant(2)), 0);
  if (multiobj)
    pb.AddObj(mp::obj::MAX, var, 0);
  return Solve(solver, pb);
}

TEST_F(NLSolverTest, MultiObjOption) {
  if (!HasFeature(feature::MULTIOBJ)) {
    EXPECT_FALSE(solver_.FindOption("multiobj"));
    return;
  }
  EXPECT_EQ(0, SolveMultiObjTestProblem(solver_));
  EXPECT_EQ(2, SolveMultiObjTestProblem(solver_, true));
}

struct TestOutputHandler : public mp::OutputHandler {
  std::string output;
  void HandleOutput(fmt::StringRef output) { this->output += output; }
};

// Test that providing an initial dual value doesn't cause an error.
TEST_F(NLSolverTest, InitialDualValue) {
  ProblemBuilder pb(solver_.GetProblemBuilder(""));
  auto info = mp::ProblemInfo();
  info.num_vars = info.num_algebraic_cons = 1;
  pb.SetInfo(info);
  pb.SetInitialDualValue(0, 0);
}

#endif  // TESTS_NL_SOLVER_TEST_H_

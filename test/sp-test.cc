/*
 Stochastic programming support tests

 Copyright (C) 2016 AMPL Optimization Inc

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

#include "sp.h"
#include "gtest-extra.h"
#include <gmock/gmock.h>

namespace expr = mp::expr;

// A test column-wise problem with specified number of variables.
class TestProblem : public mp::ColProblem {
 private:
  mp::ColProblemBuilder builder_;

 public:
  TestProblem(const mp::NLHeader &header) : builder_(*this) {
    builder_.OnHeader(header);
  }

  // Begins a random variable/vector.
  CallExprBuilder BeginRandom(int num_args) {
    auto random = AddFunction("random", 0);
    return BeginCall(random, num_args);
  }

  // Ends a random variable/vector.
  void EndRandom(CallExprBuilder builder) {
    auto zero = MakeNumericConstant(0);
    AddCon(MakeRelational(expr::NE, EndCall(builder), zero));
  }

  mp::ColProblemBuilder::ColumnSizeHandler OnColumnSizes() {
    return builder_.OnColumnSizes();
  }

  mp::ColProblemBuilder::LinearConHandler OnLinearConExpr(int con_index) {
    return builder_.OnLinearConExpr(con_index, 0);
  }

  Reference MakeTestRV() {
    auto rv = MakeVariable(0);
    auto builder = BeginRandom(2);
    builder.AddArg(rv);
    builder.AddArg(MakeNumericConstant(11));
    EndRandom(builder);
    return rv;
  }
};

mp::NLHeader MakeHeader(int num_vars, int num_integer_vars = 0) {
  mp::NLHeader header = mp::NLHeader();
  header.num_vars = num_vars;
  header.num_linear_integer_vars = num_integer_vars;
  return header;
}

// A test column-wise problem with specified number of variables and
// an empty constraint matrix.
class TestBasicProblem : public TestProblem {
 public:
  explicit TestBasicProblem(int num_vars, int num_integer_vars = 0)
    : TestProblem(MakeHeader(num_vars, num_integer_vars)) {
    auto col_sizes = OnColumnSizes();
    for (int i = 0; i < num_vars; ++i)
      col_sizes.Add(0);
  }
};

TEST(SPTest, EmptyProblem) {
  TestBasicProblem p(0);
  mp::SPAdapter sp(p);
  EXPECT_EQ(0, sp.num_vars());
  EXPECT_EQ(0, sp.num_cons());
  EXPECT_EQ(1, sp.num_stages());
  EXPECT_EQ(0, sp.num_rvs());
}

TEST(SPTest, LogicalConstraint) {
  TestBasicProblem p(1);
  p.AddCon(p.MakeRelational(
             expr::EQ, p.MakeVariable(0), p.MakeNumericConstant(0)));
  EXPECT_THROW_MSG(mp::SPAdapter sp(p);, mp::UnsupportedError,
                   "unsupported: logical constraint");
}

TEST(SPTest, NoFuncInLogicalConstraint) {
  TestBasicProblem p(1);
  p.AddCon(p.MakeRelational(
             expr::NE, p.MakeVariable(0), p.MakeNumericConstant(0)));
  EXPECT_THROW_MSG(mp::SPAdapter sp(p);, mp::UnsupportedError,
                   "unsupported: logical constraint");
}

TEST(SPTest, FuncInLogicalConstraint) {
  TestBasicProblem p(1);
  auto f = p.AddFunction("f", 0);
  auto builder = p.BeginCall(f, 0);
  p.AddCon(p.MakeRelational(
             expr::NE, p.EndCall(builder), p.MakeNumericConstant(0)));
  EXPECT_THROW_MSG(mp::SPAdapter sp(p);, mp::UnsupportedError,
                   "unsupported: logical constraint");
}

TEST(SPTest, NonzeroRHSInLogicalConstraint) {
  TestBasicProblem p(1);
  auto random = p.BeginRandom(0);
  p.AddCon(p.MakeRelational(
             expr::NE, p.EndCall(random), p.MakeNumericConstant(1)));
  EXPECT_THROW_MSG(mp::SPAdapter sp(p);, mp::UnsupportedError,
                   "unsupported: logical constraint");
}

TEST(SPTest, EmptyRandom) {
  TestBasicProblem p(1);
  auto random = p.BeginRandom(0);
  p.EndRandom(random);
  mp::SPAdapter sp(p);
  EXPECT_EQ(1, sp.num_rvs());
}

TEST(SPTest, RandomVar) {
  TestBasicProblem p(2);
  auto random = p.BeginRandom(1);
  random.AddArg(p.MakeVariable(0));
  p.EndRandom(random);
  mp::SPAdapter sp(p);
  EXPECT_EQ(1, sp.num_rvs());
  EXPECT_EQ(0, sp.rv(0).num_realizations());
  EXPECT_EQ(0, sp.rv(0).num_elements());
}

TEST(SPTest, EmptyRVWithProbability) {
  TestBasicProblem p(2);
  auto random = p.BeginRandom(1);
  random.AddArg(p.MakeNumericConstant(1));
  p.EndRandom(random);
  mp::SPAdapter sp(p);
  EXPECT_EQ(1, sp.num_rvs());
  EXPECT_EQ(1, sp.rv(0).num_realizations());
  EXPECT_EQ(0, sp.rv(0).num_elements());
}

TEST(SPTest, Probability) {
  TestBasicProblem p(2);
  auto random = p.BeginRandom(2);
  random.AddArg(p.MakeNumericConstant(0.3));
  random.AddArg(p.MakeNumericConstant(0.7));
  p.EndRandom(random);
  mp::SPAdapter sp(p);
  auto rv = sp.rv(0);
  EXPECT_EQ(2, rv.num_realizations());
  EXPECT_EQ(0.3, rv.probability(0));
  EXPECT_EQ(0.7, rv.probability(1));
}

TEST(SPTest, ProbabilityWithRandVar) {
  TestBasicProblem p(2);
  auto random = p.BeginRandom(5);
  random.AddArg(p.MakeNumericConstant(0.3));
  random.AddArg(p.MakeNumericConstant(0.7));
  random.AddArg(p.MakeVariable(0));
  random.AddArg(p.MakeNumericConstant(1));
  random.AddArg(p.MakeNumericConstant(2));
  p.EndRandom(random);
  mp::SPAdapter sp(p);
  auto rv = sp.rv(0);
  EXPECT_EQ(2, rv.num_realizations());
  EXPECT_EQ(0.3, rv.probability(0));
  EXPECT_EQ(0.7, rv.probability(1));
}

TEST(SPTest, InvalidProbability) {
  const double PROB[] = {-0.001, 1.001};
  for (auto prob: PROB) {
    TestBasicProblem p(2);
    auto random = p.BeginRandom(2);
    random.AddArg(p.MakeNumericConstant(prob));
    random.AddArg(p.MakeVariable(0));
    p.EndRandom(random);
    EXPECT_THROW_MSG(mp::SPAdapter sp(p);, mp::Error,
                     fmt::format("_slogcon[1]: invalid probability {}", prob));
  }
}

TEST(SPTest, InconsistentNumberOfRealizations) {
  TestBasicProblem p(2);
  auto random = p.BeginRandom(4);
  random.AddArg(p.MakeNumericConstant(0.3));
  random.AddArg(p.MakeNumericConstant(0.7));
  random.AddArg(p.MakeVariable(0));
  random.AddArg(p.MakeNumericConstant(1));
  p.EndRandom(random);
  EXPECT_THROW_MSG(mp::SPAdapter sp(p);, mp::Error,
                   "_slogcon[1]: inconsistent number of realizations");
}

TEST(SPTest, InvalidRandomArg) {
  TestBasicProblem p(2);
  auto random = p.BeginRandom(1);
  random.AddArg(p.MakeUnary(expr::MINUS, p.MakeNumericConstant(1)));
  p.EndRandom(random);
  EXPECT_THROW_MSG(mp::SPAdapter sp(p);, mp::Error,
                   fmt::format("_slogcon[1]: expected variable or constant"));
}

// TODO: test realizations

TEST(SPTest, SecondStageVar) {
  TestBasicProblem p(2);
  p.AddIntSuffix("stage", mp::suf::VAR, 1).SetValue(0, 2);
  mp::SPAdapter sp(p);
  EXPECT_EQ(2, sp.num_vars());
  EXPECT_EQ(0, sp.num_cons());
  EXPECT_EQ(2, sp.num_stages());
  EXPECT_EQ(1, sp.stage(0).num_vars());
  EXPECT_EQ(0, sp.stage(0).num_cons());
  EXPECT_EQ(1, sp.stage(1).num_vars());
  EXPECT_EQ(0, sp.stage(1).num_cons());
}

TEST(SPTest, OrderVarsByStage) {
  TestBasicProblem p(2, 1);
  EXPECT_EQ(mp::var::INTEGER, p.var(1).type());
  p.var(0).set_lb(42);
  p.AddIntSuffix("stage", mp::suf::VAR, 1).SetValue(0, 2);
  mp::SPAdapter sp(p);
  EXPECT_EQ(mp::var::INTEGER, sp.var(0).type());
  EXPECT_EQ(mp::var::CONTINUOUS, sp.var(1).type());
  EXPECT_EQ(42, sp.var(1).lb());
}

TEST(SPTest, SingleRealization) {
  TestBasicProblem p(2);
  auto random = p.BeginRandom(2);
  random.AddArg(p.MakeVariable(0));
  random.AddArg(p.MakeNumericConstant(42));
  p.EndRandom(random);
  mp::SPAdapter sp(p);
  EXPECT_EQ(1, sp.num_rvs());
  EXPECT_EQ(1, sp.rv(0).num_realizations());
  EXPECT_EQ(1, sp.rv(0).num_elements());
}

TEST(SPTest, MultipleRealizations) {
  TestBasicProblem p(2);
  auto random = p.BeginRandom(4);
  random.AddArg(p.MakeVariable(0));
  random.AddArg(p.MakeNumericConstant(11));
  random.AddArg(p.MakeNumericConstant(22));
  random.AddArg(p.MakeNumericConstant(33));
  p.EndRandom(random);
  mp::SPAdapter sp(p);
  EXPECT_EQ(1, sp.num_rvs());
  EXPECT_EQ(3, sp.rv(0).num_realizations());
  EXPECT_EQ(1, sp.rv(0).num_elements());
}

TEST(SPTest, Expectation) {
  TestBasicProblem p(1);
  p.AddObj(mp::obj::MIN, 0);
  auto call = p.BeginCall(p.AddFunction("expectation", 0), 1);
  auto expr = p.MakeNumericConstant(0);
  call.AddArg(expr);
  p.obj(0).set_nonlinear_expr(p.EndCall(call));
  mp::SPAdapter sp(p);
  EXPECT_EQ(expr, sp.obj(0).nonlinear_expr());
}

TEST(SPTest, Stage2VarOutsideOfExpectation) {
  TestBasicProblem p(1);
  p.AddIntSuffix("stage", mp::suf::VAR, 1).SetValue(0, 2);
  p.AddObj(mp::obj::MIN, 1).AddTerm(0, 1);
  EXPECT_THROW_MSG(mp::SPAdapter sp(p), mp::Error,
                   "second-stage variable outside of expectation in objective");
}

TEST(SPTest, TooFewArgsToExpectation) {
  TestBasicProblem p(1);
  p.AddObj(mp::obj::MIN, 0);
  auto call = p.BeginCall(p.AddFunction("expectation", 0), 0);
  p.obj(0).set_nonlinear_expr(p.EndCall(call));
  EXPECT_THROW_MSG(mp::SPAdapter sp(p), mp::Error,
                   "invalid arguments to expectation");
}

TEST(SPTest, TooManyArgsToExpectation) {
  TestBasicProblem p(1);
  p.AddObj(mp::obj::MIN, 0);
  auto call = p.BeginCall(p.AddFunction("expectation", 0), 2);
  call.AddArg(p.MakeNumericConstant(0));
  call.AddArg(p.MakeNumericConstant(0));
  p.obj(0).set_nonlinear_expr(p.EndCall(call));
  EXPECT_THROW_MSG(mp::SPAdapter sp(p), mp::Error,
                   "invalid arguments to expectation");
}

TEST(SPTest, StringArgToExpectation) {
  TestBasicProblem p(1);
  p.AddObj(mp::obj::MIN, 0);
  auto call = p.BeginCall(p.AddFunction("expectation", 0), 1);
  call.AddArg(p.MakeStringLiteral("foo"));
  p.obj(0).set_nonlinear_expr(p.EndCall(call));
  EXPECT_THROW_MSG(mp::SPAdapter sp(p), mp::Error,
                   "invalid arguments to expectation");
}

TEST(SPTest, RandomConMatrix) {
  auto header = MakeHeader(2);
  header.num_con_nonzeros = 1;
  TestProblem p(header);
  auto rv = p.MakeTestRV();
  auto con = p.AddCon(0, 0);
  con.set_nonlinear_expr(p.MakeBinary(expr::MUL, rv, p.MakeVariable(1)));
  auto cols = p.OnColumnSizes();
  cols.Add(0);
  cols.Add(1);
  p.OnLinearConExpr(0).AddTerm(1, 42);
  mp::SPAdapter sp(p);
  auto col = sp.column(0);
  auto it = col.begin();
  EXPECT_EQ(42, it->coef());
  EXPECT_EQ(0, it->con_index());
  EXPECT_EQ(col.end(), ++it);
}

TEST(SPTest, RandomRHS) {
  auto header = MakeHeader(2);
  header.num_con_nonzeros = 1;
  TestProblem p(header);
  p.MakeTestRV();
  p.AddCon(0, 0);
  auto cols = p.OnColumnSizes();
  cols.Add(1);
  cols.Add(0);
  p.OnLinearConExpr(0).AddTerm(0, 42);
  mp::SPAdapter sp(p);
  auto col = sp.column(0);
  EXPECT_EQ(col.begin(), col.end());
}

TEST(SPTest, MoreThan2StagesNotSupported) {
  TestBasicProblem p(1);
  p.AddIntSuffix("stage", mp::suf::VAR, 1).SetValue(0, 3);
  EXPECT_THROW_MSG(mp::SPAdapter sp(p), mp::Error,
                   "SP problems with more than 2 stages are not supported");
}

TEST(SPTest, NonlinearNotSupported) {
  TestBasicProblem p(1);
  p.AddCon(0, 0);
  p.algebraic_con(0).set_nonlinear_expr(
        p.MakeUnary(expr::POW2, p.MakeVariable(0)));
  mp::SPAdapter sp(p);
  mp::SPAdapter::Scenario scenario;
  EXPECT_THROW_MSG(sp.GetScenario(scenario, 0), mp::UnsupportedError,
                   "unsupported: ^2");
}

TEST(SPTest, GetScenario) {
  auto header = MakeHeader(2);
  header.num_con_nonzeros = 1;
  TestProblem p(header);
  p.MakeTestRV();
  p.AddCon(0, 0);
  auto cols = p.OnColumnSizes();
  cols.Add(1);
  cols.Add(0);
  p.OnLinearConExpr(0).AddTerm(0, -3);
  mp::SPAdapter sp(p);
  mp::SPAdapter::Scenario scenario;
  sp.GetScenario(scenario, 0);
  EXPECT_THAT(scenario.rhs_offsets(), testing::ElementsAre(33));
}

TEST(SPTest, RandoRHSInNonlinear) {
  auto header = MakeHeader(2);
  header.num_con_nonzeros = 1;
  TestProblem p(header);
  p.MakeTestRV();
  p.AddCon(0, 0).set_nonlinear_expr(p.MakeVariable(0));
  auto cols = p.OnColumnSizes();
  cols.Add(1);
  cols.Add(0);
  p.OnLinearConExpr(0).AddTerm(0, -3);
  mp::SPAdapter sp(p);
  mp::SPAdapter::Scenario scenario;
  sp.GetScenario(scenario, 0);
  EXPECT_THAT(scenario.rhs_offsets(), testing::ElementsAre(44));
}

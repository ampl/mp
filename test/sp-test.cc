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
#include "gtest/gtest.h"

// A test column-wise problem with specified number of variables.
class TestProblem : public mp::ColProblem {
 public:
  explicit TestProblem(int num_vars, int num_integer_vars = 0) {
    mp::ColProblemBuilder builder(*this);
    mp::NLHeader header = mp::NLHeader();
    header.num_vars = num_vars;
    header.num_linear_integer_vars = num_integer_vars;
    builder.OnHeader(header);
    auto col_sizes = builder.OnColumnSizes();
    for (int i = 0; i < num_vars; ++i)
      col_sizes.Add(0);
  }

  // Begins a random variable/vector.
  CallExprBuilder BeginRandom(int num_args) {
    auto random = AddFunction("random", 0);
    return BeginCall(random, num_args);
  }

  // Ends a random variable/vector.
  void EndRandom(CallExprBuilder builder) {
    auto zero = MakeNumericConstant(0);
    AddCon(MakeRelational(mp::expr::NE, EndCall(builder), zero));
  }
};

TEST(SPTest, EmptyProblem) {
  TestProblem p(0);
  mp::SPAdapter sp(p);
  EXPECT_EQ(0, sp.num_vars());
  EXPECT_EQ(0, sp.num_cons());
  EXPECT_EQ(1, sp.num_stages());
  EXPECT_EQ(0, sp.num_rvs());
}

TEST(SPTest, SecondStageVar) {
  TestProblem p(2);
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
  TestProblem p(2, 1);
  EXPECT_EQ(mp::var::INTEGER, p.var(1).type());
  p.var(0).set_lb(42);
  p.AddIntSuffix("stage", mp::suf::VAR, 1).SetValue(0, 2);
  mp::SPAdapter sp(p);
  EXPECT_EQ(mp::var::INTEGER, sp.var(0).type());
  EXPECT_EQ(mp::var::CONTINUOUS, sp.var(1).type());
  EXPECT_EQ(42, sp.var(1).lb());
}

TEST(SPTest, RandomVar) {
  TestProblem p(2);
  auto random = p.BeginRandom(1);
  random.AddArg(p.MakeVariable(0));
  p.EndRandom(random);
  mp::SPAdapter sp(p);
  EXPECT_EQ(1, sp.num_rvs());
  EXPECT_EQ(0, sp.rv(0).num_realizations());
  EXPECT_EQ(0, sp.rv(0).num_elements());
}

TEST(SPTest, SingleRealization) {
  TestProblem p(2);
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
  TestProblem p(2);
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

// TODO: more tests: "random-con-matrix", "random-con-matrix2", "random-rhs"

/*
TEST(SMPSWriterTest, NonlinearNotSupported) {
  WriteFile("test.nl", ReadFile(MP_TEST_DATA_DIR "/smps/nonlinear.nl"));
  EXPECT_THROW(Solve("test"), mp::Error);
}

TEST(SMPSWriterTest, MoreThan2StagesNotSupported) {
  WriteFile("test.nl", ReadFile(MP_TEST_DATA_DIR "/smps/three-stage.nl"));
  EXPECT_THROW(Solve("test"), mp::Error);
}

TEST(SMPSWriterTest, RangesNotSupported) {
  TempSMPSFiles files("smps/range-con");
  EXPECT_THROW(Solve("test"), mp::Error);
}

TEST(SMPSWriterTest, InconsistentProbabilities) {
  std::string filename = MP_TEST_DATA_DIR "/smps/inconsistent-probabilities";
  WriteFile("test.nl", ReadFile(filename + ".nl"));
  WriteFile("test.row", ReadFile(filename + ".row"));
  WriteFile("test.col", ReadFile(filename + ".col"));
  EXPECT_THROW(Solve("test"), mp::Error);
}
*/

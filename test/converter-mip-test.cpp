/*
 Tests of a MIP interface

 Copyright (C) 2020 AMPL Optimization Inc

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

 Author: Gleb Belov <Gleb.Belov@monash.edu>
 */

#include <cmath>
#include <set>

#include "gtest/gtest.h"

#include "converter-mip-test.h"

// TODO adapt for interface-to-end test of any concrete backend


namespace {

using namespace mip_converter_test;

TEST(RedefsMIPTest, PureMILP__01__seemsPassedToBackendCorrectly) {
  MIPInstance milp {
    { { minimize_,
        { { 1, 1, 2.5 },
          { 0, 1, 2 } } } },
    { {  0.0,  0.0, -1.2 },
      {  1.0,  3.0,  5.2 },
      {  I_,   F_,   F_ } },
    {
      { { { 2, 0.0, 4 },
          { 0, 1,   2 } }, -INFINITY, 56.4 }
    }
  };
  mp::Env e;
  MIPConverterTester tester(e);
  feedInstance(tester, milp);
  tester.ConvertModel();
  ASSERT_TRUE(tester.ObjsEqual(milp));
  ASSERT_TRUE(tester.VarBoundsEqual(milp));
  ASSERT_TRUE(tester.NConstrEqual(milp));
}



////////////// Absorbing disjunctions into indicators (cvt:pre:disj2ind) //////////////

using mp::expr::OR;
using mp::expr::AND;
using mp::expr::NE;
using mp::expr::EQ;
using mp::expr::LE;
using mp::expr::LT;
using mp::expr::GE;
using mp::var::INTEGER;
using mp::var::CONTINUOUS;

/// b==1 ==> x<=3, as AMPL sends it: b!=1 || x<=3
template <class Tester>
void AddImplLE(Tester& t, int b, int x, double rhs) {
  t.Add( t.Bin(OR, t.Cmp(NE, b, 1), t.Cmp(LE, x, rhs)) );
}

/// Check that the only indicator is \a ind,
/// on literal (b==bval), body x (cmp) rhs
template <class Ind>
void ExpectSingleIndicator(const std::vector<Ind>& inds,
                           int b, int bval, int x, double rhs) {
  ASSERT_EQ(1u, inds.size());
  EXPECT_EQ(b, inds[0].get_binary_var());
  EXPECT_EQ(bval, inds[0].get_binary_value());
  const auto& con = inds[0].get_constraint();
  ASSERT_EQ(1u, con.size());
  EXPECT_EQ(x, con.var(0));
  EXPECT_DOUBLE_EQ(rhs, con.rhs());
}

TEST(Disj2IndTest, ImplicationLE__becomesIndicatorOnB) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(-10, 10);
  AddImplLE(t, b, x, 3);
  t.ConvertModel();
  ExpectSingleIndicator(t.Backend().GetIndicatorConstraintLinLE(), b, 1, x, 3);
  EXPECT_EQ(0, t.Backend().NumLinCons());
  EXPECT_EQ(0, t.NumUsed<mp::CondLinConLE>());
  EXPECT_EQ(0, t.NumUsed<mp::NotConstraint>());
}

TEST(Disj2IndTest, ImplicationLE__optionOff__indicatorOnAuxVar) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "cvt:pre:disj2ind=0");
  int b = t.Var(0, 1, INTEGER), x = t.Var(-10, 10);
  AddImplLE(t, b, x, 3);
  t.ConvertModel();
  const auto& inds = t.Backend().GetIndicatorConstraintLinLE();
  ASSERT_EQ(1u, inds.size());
  EXPECT_NE(b, inds[0].get_binary_var());
  EXPECT_LE(1, t.Backend().NumLinCons());           // the OR row
  EXPECT_EQ(1, t.NumUsed<mp::CondLinConLE>());
}

TEST(Disj2IndTest, ImplicationLE__presolveOff__untouched) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "cvt:pre:all=0");
  int b = t.Var(0, 1, INTEGER), x = t.Var(-10, 10);
  AddImplLE(t, b, x, 3);
  t.ConvertModel();
  EXPECT_EQ(1, t.NumUsed<mp::CondLinConLE>());
}

TEST(Disj2IndTest, ImplicationLE__bigM__singleRowNoAuxBinary) {
  mp::Env e;
  Disj2IndTester<false, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(-10, 10);
  AddImplLE(t, b, x, 3);
  t.ConvertModel();
  EXPECT_EQ(1, t.Backend().NumLinCons());           // x + 7b <= 10
  EXPECT_EQ(0, t.NumUsed<mp::CondLinConLE>());
}

TEST(Disj2IndTest, ImplicationWithZero__literalB0) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(-10, 10);
  t.Add( t.Bin(OR, t.Cmp(NE, b, 0), t.Cmp(GE, x, 2)) );   // b==0 ==> x>=2
  t.ConvertModel();
  ExpectSingleIndicator(t.Backend().GetIndicatorConstraintLinGE(), b, 0, x, 2);
}

TEST(Disj2IndTest, StrictInequality__usesComparisonEps) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(-10, 10);
  int c = t.Var(0, 1, INTEGER), y = t.Var(-10, 10, INTEGER);
  t.Add( t.Bin(OR, t.Cmp(NE, b, 1), t.Cmp(LT, x, 3)) );
  t.Add( t.Bin(OR, t.Cmp(NE, c, 1), t.Cmp(LT, y, 3)) );
  t.ConvertModel();
  const auto& inds = t.Backend().GetIndicatorConstraintLinLE();
  ASSERT_EQ(2u, inds.size());
  for (const auto& ind: inds) {
    EXPECT_TRUE(ind.get_binary_var()==b || ind.get_binary_var()==c);
    const double rhs_expected = ind.get_binary_var()==b ? 3-1e-4 : 2;
    EXPECT_DOUBLE_EQ(rhs_expected, ind.get_constraint().rhs());
  }
}

TEST(Disj2IndTest, ConjunctionConsequent__twoIndicatorsOnB) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(-10, 10), y = t.Var(-10, 10);
  t.Add( t.Bin(OR, t.Cmp(NE, b, 1),
               t.Bin(AND, t.Cmp(LE, x, 3), t.Cmp(GE, y, 1))) );
  t.ConvertModel();
  ExpectSingleIndicator(t.Backend().GetIndicatorConstraintLinLE(), b, 1, x, 3);
  ExpectSingleIndicator(t.Backend().GetIndicatorConstraintLinGE(), b, 1, y, 1);
  EXPECT_EQ(0, t.NumUsed<mp::AndConstraint>());
  EXPECT_EQ(0, t.Backend().NumLinCons());
}

TEST(Disj2IndTest, DisequalityConsequent__literalInDisjunctionRow) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(-10, 10);
  t.Add( t.Bin(OR, t.Cmp(NE, b, 1), t.Cmp(NE, x, 3)) );   // b==1 ==> x!=3
  t.ConvertModel();
  EXPECT_EQ(1u, t.Backend().GetIndicatorConstraintLinLE().size());
  EXPECT_EQ(1u, t.Backend().GetIndicatorConstraintLinGE().size());
  const auto& rows = t.Backend().GetLinConGE();
  ASSERT_EQ(1u, rows.size());                    // -b + q1 + q2 >= 0
  EXPECT_EQ(3u, rows[0].size());
  EXPECT_DOUBLE_EQ(0.0, rows[0].rhs());
  EXPECT_EQ(0, t.NumUsed<mp::CondLinConEQ>());
}

TEST(Disj2IndTest, SharedConsequent__notAbsorbed) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int a = t.Var(0, 1, INTEGER), b = t.Var(0, 1, INTEGER);
  int x = t.Var(-10, 10);
  AddImplLE(t, a, x, 3);
  AddImplLE(t, b, x, 3);
  t.ConvertModel();
  const auto& inds = t.Backend().GetIndicatorConstraintLinLE();
  ASSERT_EQ(1u, inds.size());                    // on the shared result
  EXPECT_NE(a, inds[0].get_binary_var());
  EXPECT_NE(b, inds[0].get_binary_var());
}

TEST(Disj2IndTest, NativeOr__absorbedAnyway) {
  mp::Env e;
  Disj2IndTester<true, true> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(-10, 10);
  AddImplLE(t, b, x, 3);
  t.ConvertModel();
  ExpectSingleIndicator(t.Backend().GetIndicatorConstraintLinLE(), b, 1, x, 3);
  EXPECT_EQ(0u, t.Backend().GetOrConstraint().size());
}

TEST(Disj2IndTest, NativeConditionals__untouched) {
  mp::Env e;
  Disj2IndTester<true, true, true> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(-10, 10);
  AddImplLE(t, b, x, 3);
  t.ConvertModel();
  EXPECT_EQ(1u, t.Backend().GetCondLinConLE().size());
  EXPECT_EQ(1u, t.Backend().GetOrConstraint().size());
  EXPECT_EQ(0, t.Backend().NumIndicators());
}

/// b==1 ==> x==5 with integer x in 1..9;
/// if \a fEnc, also (x!=2 && x!=7) || y>=4, forcing x's unary encoding
template <class Tester>
void AddUEncCase(Tester& t, int b, int x, bool fEnc) {
  t.Add( t.Bin(OR, t.Cmp(NE, b, 1), t.Cmp(EQ, x, 5)) );
  if (fEnc) {
    int y = t.Var(1, 200);
    t.Add( t.Bin(OR, t.Bin(AND, t.Cmp(NE, x, 2), t.Cmp(NE, x, 7)),
                 t.Cmp(GE, y, 4)) );
  }
}

TEST(Disj2IndTest, UEncCandidate__noEncoding__absorbed) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(1, 9, INTEGER);
  AddUEncCase(t, b, x, false);
  t.ConvertModel();
  ExpectSingleIndicator(t.Backend().GetIndicatorConstraintLinEQ(), b, 1, x, 5);
  EXPECT_EQ(0, t.NumUsed<mp::CondLinConEQ>());
}

TEST(Disj2IndTest, UEncCandidate__encoding__flagRowKept) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(1, 9, INTEGER);
  AddUEncCase(t, b, x, true);
  t.ConvertModel();
  EXPECT_EQ(0u, t.Backend().GetIndicatorConstraintLinEQ().size());
  EXPECT_EQ(3, t.NumUsed<mp::CondLinConEQ>());   // bridged by the encoding
}

TEST(Disj2IndTest, UEncCandidate__bitOff__notAbsorbed) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "cvt:pre:disj2ind=3");
  int b = t.Var(0, 1, INTEGER), x = t.Var(1, 9, INTEGER);
  AddUEncCase(t, b, x, false);
  t.ConvertModel();
  const auto& inds = t.Backend().GetIndicatorConstraintLinEQ();
  ASSERT_EQ(1u, inds.size());
  EXPECT_NE(b, inds[0].get_binary_var());
}

TEST(Disj2IndTest, TwoConditionals__bit4__oneBinary) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "cvt:pre:disj2ind=15");
  int x = t.Var(-10, 10), y = t.Var(-10, 10);
  t.Add( t.Bin(OR, t.Cmp(LE, x, 0), t.Cmp(LE, y, 0)) );   // x<=0 || y<=0
  t.ConvertModel();
  const auto& inds = t.Backend().GetIndicatorConstraintLinLE();
  ASSERT_EQ(2u, inds.size());
  EXPECT_EQ(inds[0].get_binary_var(), inds[1].get_binary_var());
  EXPECT_NE(inds[0].get_binary_value(), inds[1].get_binary_value());
  EXPECT_EQ(0, t.Backend().NumLinCons());
}

TEST(Disj2IndTest, TwoConditionals__default__notAbsorbed) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int x = t.Var(-10, 10), y = t.Var(-10, 10);
  t.Add( t.Bin(OR, t.Cmp(LE, x, 0), t.Cmp(LE, y, 0)) );
  t.ConvertModel();
  const auto& inds = t.Backend().GetIndicatorConstraintLinLE();
  ASSERT_EQ(2u, inds.size());
  EXPECT_NE(inds[0].get_binary_var(), inds[1].get_binary_var());
}


/// x in {0} union [2, 10] encoded as x >= 2b, x <= 10b (b binary),
/// as AMPL does for semicontinuous variables
template <class Tester>
void AddSemicontinuous(Tester& t, int x, int b) {
  t.Lin2(1.0, x, -2.0, b, 0.0, INFINITY);
  t.Lin2(1.0, x, -10.0, b, -INFINITY, 0.0);
}

TEST(VarBndTest, ComparisonInGap__reusesBinary) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int x = t.Var(0, 10), b = t.Var(0, 1, INTEGER), a = t.Var(0, 1, INTEGER);
  AddSemicontinuous(t, x, b);
  t.Add( t.Bin(OR, t.Cmp(NE, a, 1), t.Cmp(GE, x, 1)) );   // a==1 ==> x>=1
  t.Add( t.Bin(OR, t.Cmp(LE, x, 0.5), t.Cmp(NE, a, 0)) );  // a==0 ==> x<=0.5
  t.ConvertModel();
  EXPECT_EQ(0, t.Backend().NumIndicators());     // clauses on a, b only
  EXPECT_EQ(0, t.template NumUsed<mp::CondLinConGE>());
  EXPECT_EQ(0, t.template NumUsed<mp::CondLinConLE>());
}

TEST(VarBndTest, ThresholdAboveLowerBound__notDecided) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int x = t.Var(0, 10), b = t.Var(0, 1, INTEGER), a = t.Var(0, 1, INTEGER);
  AddSemicontinuous(t, x, b);
  t.Add( t.Bin(OR, t.Cmp(NE, a, 1), t.Cmp(GE, x, 3)) );   // 3 > L=2
  t.ConvertModel();
  ExpectSingleIndicator(t.Backend().GetIndicatorConstraintLinGE(), a, 1, x, 3);
}

TEST(VarBndTest, OptionOff__notDecided) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "cvt:pre:vbnd=0 cvt:pre:disj2ind=0");
  int x = t.Var(0, 10), b = t.Var(0, 1, INTEGER), a = t.Var(0, 1, INTEGER);
  AddSemicontinuous(t, x, b);
  t.Add( t.Bin(OR, t.Cmp(NE, a, 1), t.Cmp(GE, x, 1)) );
  t.ConvertModel();
  EXPECT_EQ(1, t.template NumUsed<mp::CondLinConGE>());
}

/// (x<=1 && y>=2) || (x>=5 && y<=0) || z>=3
template <class Tester>
void AddThreeTermGdp(Tester& t, int x, int y, int z) {
  t.Add( t.Bin(OR,
               t.Bin(OR, t.And2(t.Cmp(LE, x, 1), t.Cmp(GE, y, 2)),
                         t.And2(t.Cmp(GE, x, 5), t.Cmp(LE, y, 0))),
               t.Cmp(GE, z, 3)) );
}

TEST(Disj2IndTest, ThreeTermGdp__bit16__indicatorsOnDisjuncts) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "cvt:pre:disj2ind=27");
  int x = t.Var(-10, 10), y = t.Var(-10, 10), z = t.Var(-10, 10);
  AddThreeTermGdp(t, x, y, z);
  t.ConvertModel();
  EXPECT_EQ(5, t.Backend().NumIndicators());
  EXPECT_EQ(0, t.template NumUsed<mp::CondLinConLE>());     // under the Ands
  EXPECT_EQ(1, t.template NumUsed<mp::CondLinConGE>());     // z>=3 only
  std::set<int> bins;
  for (const auto& ind: t.Backend().GetIndicatorConstraintLinLE())
    bins.insert(ind.get_binary_var());
  for (const auto& ind: t.Backend().GetIndicatorConstraintLinGE())
    bins.insert(ind.get_binary_var());
  EXPECT_EQ(3u, bins.size());                   // one per disjunct
}

TEST(Disj2IndTest, ThreeTermGdp__default__conjunctsReified) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int x = t.Var(-10, 10), y = t.Var(-10, 10), z = t.Var(-10, 10);
  AddThreeTermGdp(t, x, y, z);
  t.ConvertModel();
  EXPECT_EQ(2, t.template NumUsed<mp::CondLinConLE>());
}


/// b==1 ==> (x1==1 && x2==1) on binaries: Or(Not(b), And(x1, x2))
template <class Tester>
void AddImplOfConjunction(Tester& t, int b, int x1, int x2) {
  t.Add( t.Bin(OR, t.Cmp(NE, b, 1),
               t.And2(t.Cmp(EQ, x1, 1), t.Cmp(EQ, x2, 1))) );
}

TEST(Disj2IndTest, OrOverAndOfBinaries__bit32__clauses) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "cvt:pre:disj2ind=43");
  int b = t.Var(0, 1, INTEGER), x1 = t.Var(0, 1, INTEGER), x2 = t.Var(0, 1, INTEGER);
  AddImplOfConjunction(t, b, x1, x2);
  t.ConvertModel();
  EXPECT_EQ(0, t.template NumUsed<mp::AndConstraint>());
  EXPECT_EQ(2, t.Backend().NumLinCons());      // -b + xi >= 0
}

TEST(Disj2IndTest, OrOverAndOfBinaries__default__andKept) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x1 = t.Var(0, 1, INTEGER), x2 = t.Var(0, 1, INTEGER);
  AddImplOfConjunction(t, b, x1, x2);
  t.ConvertModel();
  EXPECT_EQ(1, t.template NumUsed<mp::AndConstraint>());
}


/// x == 0 || x >= 2  (x in [0, 10]): domain hole (0, 2);
/// a==1 ==> x >= 1  and  a==0 ==> x <= 0.5
template <class Tester>
void AddHoleCase(Tester& t, int x, int a) {
  t.Add( t.Bin(OR, t.Cmp(EQ, x, 0), t.Cmp(GE, x, 2)) );
  t.Add( t.Bin(OR, t.Cmp(NE, a, 1), t.Cmp(GE, x, 1)) );
  t.Add( t.Bin(OR, t.Cmp(LE, x, 0.5), t.Cmp(NE, a, 0)) );
}

TEST(VarBndTest, DomainHole__comparisonsShareOneBinary) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int x = t.Var(0, 10), a = t.Var(0, 1, INTEGER);
  AddHoleCase(t, x, a);
  t.ConvertModel();
  EXPECT_EQ(1, t.template NumUsed<mp::CondLinConGE>());    // x >= 2 only
  EXPECT_EQ(0, t.template NumUsed<mp::CondLinConLE>());
}

TEST(VarBndTest, DomainHole__optionBitOff__separateComparisons) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "cvt:pre:vbnd=1 cvt:pre:disj2ind=0");
  int x = t.Var(0, 10), a = t.Var(0, 1, INTEGER);
  AddHoleCase(t, x, a);
  t.ConvertModel();
  EXPECT_EQ(2, t.template NumUsed<mp::CondLinConGE>());    // x >= 2, x >= 1
}

}


/// cvt:pre:auxint: implied-integer auxiliaries submitted as continuous
TEST(AuxIntTest, IfThenOfIntegers__continuousSubmitted_integerInside) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(0, 10, INTEGER),
      y = t.Var(0, 10, INTEGER);
  auto ite = t.If(t.Cmp(EQ, b, 1), x, y);
  t.Minimize(ite);
  t.Add( t.CmpE(GE, ite, 4.5) );        // integer result: >= 5
  t.ConvertModel();
  // The if-then result, and the root comparison's result fixed to 1
  EXPECT_EQ(2, t.NumAuxVars(INTEGER, CONTINUOUS));
  EXPECT_EQ(0, t.NumAuxVars(INTEGER, INTEGER));
  int r = y+1;                           // the if-then result
  EXPECT_EQ(INTEGER, t.InternalType(r));
  EXPECT_EQ(CONTINUOUS, t.Backend().type(r));
  int nrows = 0;
  for (const auto& lc: t.Backend().GetLinConGE())
    if (1 == lc.size() && r == lc.pvars()[0]) {
      EXPECT_EQ(5.0, lc.lb() / lc.pcoefs()[0]);  // rounded as for an integer
      ++nrows;
    }
  EXPECT_EQ(1, nrows);
}

TEST(AuxIntTest, OptionOn__integerSubmitted) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "cvt:pre:auxint=1");
  int b = t.Var(0, 1, INTEGER), x = t.Var(0, 10, INTEGER),
      y = t.Var(0, 10, INTEGER);
  t.Minimize( t.If(t.Cmp(EQ, b, 1), x, y) );
  t.ConvertModel();
  EXPECT_EQ(0, t.NumAuxVars(INTEGER, CONTINUOUS));
  EXPECT_EQ(1, t.NumAuxVars(INTEGER, INTEGER));
}

TEST(AuxIntTest, IfThenOfContinuous__continuous) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(0, 10), y = t.Var(0, 10, INTEGER);
  t.Minimize( t.If(t.Cmp(EQ, b, 1), x, y) );
  t.ConvertModel();
  EXPECT_EQ(0, t.NumAuxVars(INTEGER, CONTINUOUS));
  EXPECT_EQ(0, t.NumAuxVars(INTEGER, INTEGER));
}

/// redefvar_01: y <= 4 || (if b==0 then 1 else 0),
/// the if-then result becomes an indicator binary (cvt:pre:disj2ind)
TEST(AuxIntTest, IfThenAsIndicatorBinary__binary) {
  mp::Env e;
  Disj2IndTester<true, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), y = t.Var(0, 10);
  t.Add( t.Bin(OR, t.Cmp(LE, y, 4), t.CmpE(NE, t.IfC(t.Cmp(EQ, b, 0), 1, 0), 0)) );
  t.ConvertModel();
  ASSERT_LE(1, t.Backend().NumIndicators());
  for (const auto& ic: t.Backend().GetIndicatorConstraintLinLE()) {
    auto bv = ic.get_binary_var();
    EXPECT_EQ(INTEGER, t.Backend().type(bv));
    EXPECT_EQ(0.0, t.Backend().lb(bv));
    EXPECT_EQ(1.0, t.Backend().ub(bv));
  }
}

TEST(AuxIntTest, ComparisonResults__integer) {
  mp::Env e;
  Disj2IndTester<false, false> t(e, "cvt:pre:disj2ind=0");
  int x = t.Var(-10, 10);
  t.Add( t.Bin(OR, t.Cmp(GE, x, 1), t.Cmp(LE, x, -1)) );   // big-M
  t.ConvertModel();
  EXPECT_EQ(2, t.NumAuxVars(INTEGER, INTEGER));       // the 2 results
  EXPECT_EQ(1, t.NumAuxVars(INTEGER, CONTINUOUS));    // OR result, fixed
}

TEST(AuxIntTest, CountResult__continuous) {
  mp::Env e;
  Disj2IndTester<false, false> t(e, "");
  int x = t.Var(-10, 10);
  t.Minimize( t.Count2(t.Cmp(GE, x, 1), t.Cmp(LE, x, -1)) );
  t.ConvertModel();
  EXPECT_EQ(2, t.NumAuxVars(INTEGER, INTEGER));       // comparison results
  EXPECT_EQ(1, t.NumAuxVars(INTEGER, CONTINUOUS));    // count
}

TEST(AuxIntTest, FixedIntVar__continuousSubmittedOnly) {
  mp::Env e;
  Disj2IndTester<false, false> t(e, "");
  int x = t.Var(3, 3, INTEGER), y = t.Var(0, 10, INTEGER);
  t.Lin2(1.0, x, 1.0, y, -INFINITY, 8.0);
  t.ConvertModel();
  EXPECT_EQ(CONTINUOUS, t.Backend().type(x));
  EXPECT_EQ(INTEGER, t.InternalType(x));
  EXPECT_EQ(INTEGER, t.Backend().type(y));
}


/// Big-M for  b==0 ==> y >= 0  where y's best-known lower bound is 0
/// (unsubmitted): the row must stay, else y is free in the solver
TEST(BigMTest, ImpliedByBestBoundsOnly__rowKept) {
  mp::Env e;
  Disj2IndTester<false, false> t(e, "");
  int b = t.Var(0, 1, INTEGER), x = t.Var(50, 270, INTEGER);
  auto& m = t.GetModel();
  auto ite = m.MakeIf(t.Cmp(EQ, b, 1), m.MakeVariable(x),
                      m.MakeNumericConstant(0.0));
  t.Minimize(ite);
  t.ConvertModel();
  int r = -1;                                   // the if-then result:
  for (int v=t.Backend().NumVars(); --v > x; )  // integer inside, not fixed
    if (INTEGER == t.InternalType(v) && t.Backend().lb(v) < t.Backend().ub(v))
      r = v;
  ASSERT_GE(r, 0);
  // A row bounding y from below by a constant, i.e.,
  // its other variables (the else constant) are fixed
  auto IsLowerRowOnR = [&t, r](const auto& lc, double sign) {
    bool has_r = false;
    for (size_t i=0; i<lc.size(); ++i) {
      if (r == lc.pvars()[i]) {
        if (sign * lc.pcoefs()[i] <= 0.0)
          return false;
        has_r = true;
      } else if (t.Backend().lb(lc.pvars()[i]) != t.Backend().ub(lc.pvars()[i]))
        return false;
    }
    return has_r;
  };
  int nrows = 0;
  for (const auto& lc: t.Backend().GetLinConGE())
    nrows += IsLowerRowOnR(lc, 1.0);
  for (const auto& lc: t.Backend().GetLinConLE())
    nrows += IsLowerRowOnR(lc, -1.0);
  for (const auto& lc: t.Backend().GetLinConRange())
    nrows += IsLowerRowOnR(lc, 1.0) && lc.lb() > -1e20;
  EXPECT_EQ(1, nrows);
}

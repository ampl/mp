/*
 JaCoP solver tests.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include <algorithm>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <csignal>
#include <cstdlib>

#include "gtest/gtest.h"

#include "solvers/jacop/jacop.h"
#include "solvers/util/expr.h"

extern "C" {
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
}

#include "tests/args.h"
#include "tests/expr_builder.h"
#include "tests/solution_handler.h"
#include "tests/util.h"
#include "tests/config.h"

#ifdef HAVE_THREADS
# include <thread>
#endif

#if defined(_MSC_VER)
# define isnan _isnan
#else
# define isnan std::isnan
#endif

using std::ifstream;
using std::size_t;
using std::string;
using std::vector;

using ampl::Class;
using ampl::Variable;
using ampl::ExprBuilder;
using ampl::LogicalExpr;
using ampl::NumericExpr;
using ampl::UnsupportedExprError;
using ampl::NLToJaCoPConverter;

#define DATA_DIR "../data/"

namespace {

class JaCoPConverterTest : public ::testing::Test, public ExprBuilder {
 protected:
  static ampl::JVM jvm_;
  int min_int_;
  int max_int_;
  Variable x;
  Variable y;
  Variable z;

  void InitVars(NLToJaCoPConverter &c, int var1, int var2, int var3);
  double Solve(NLToJaCoPConverter &c);

  // Converts a numeric expression from AMPL to Gecode form and evaluates it.
  // Returns the value of the Gecode expression.
  double ConvertAndEval(NumericExpr e,
      int var1 = 0, int var2 = 0, int var3 = 0);
  double ConvertAndEval(LogicalExpr e,
      int var1 = 0, int var2 = 0, int var3 = 0);

  JaCoPConverterTest() : x(AddVar(1)), y(AddVar(2)), z(AddVar(3)) {
    jclass domain_class = jvm_.FindClass("JaCoP/core/IntDomain");
    min_int_ = jvm_.GetStaticIntField(
        domain_class, jvm_.GetStaticFieldID(domain_class, "MinInt", "I"));
    max_int_ = jvm_.GetStaticIntField(
        domain_class, jvm_.GetStaticFieldID(domain_class, "MaxInt", "I"));
  }

  double min() const { return min_int_; }
  double max() const { return max_int_; }
};

ampl::JVM JaCoPConverterTest::jvm_;

void JaCoPConverterTest::InitVars(
    NLToJaCoPConverter &c, int var1, int var2, int var3) {
  Class<ampl::IntVar> &cls = c.var_class();
  jobjectArray vars = c.var_array(); // TODO: var array is null
  jvm_.SetObjectArrayElement(vars, 0, cls.NewObject(jvm_, var1, var1));
  jvm_.SetObjectArrayElement(vars, 1, cls.NewObject(jvm_, var2, var2));
  jvm_.SetObjectArrayElement(vars, 2, cls.NewObject(jvm_, var3, var3));
}

double JaCoPConverterTest::Solve(NLToJaCoPConverter &c) {
  Class<ampl::DepthFirstSearch> dfs_class;
  jobject search = dfs_class.NewObject(jvm_);
  jmethodID setPrintInfo =
      jvm_.GetMethod(dfs_class.get(), "setPrintInfo", "(Z)V");
  jvm_.CallVoidMethod(search, setPrintInfo, false);
  jmethodID labeling = jvm_.GetMethod(dfs_class.get(), "labeling",
      "(LJaCoP/core/Store;LJaCoP/search/SelectChoicePoint;)Z");
  jobject indomain = jvm_.NewObject("JaCoP/search/IndomainMin", "()V");
  jobject select = jvm_.NewObject("JaCoP/search/InputOrderSelect",
      "(LJaCoP/core/Store;[LJaCoP/core/Var;LJaCoP/search/Indomain;)V",
      c.store(), c.var_array(), indomain);
  jboolean found = jvm_.CallBooleanMethod(search, labeling, c.store(), select);
  if (!found) return std::numeric_limits<double>::quiet_NaN();
  jmethodID value = jvm_.GetMethod(c.var_class().get(), "value", "()I");
  return jvm_.CallIntMethod(c.vars()[0], value);
}

// Converts a numeric expression from AMPL to Gecode form and evaluates it.
// Returns the value of the Gecode expression.
double JaCoPConverterTest::ConvertAndEval(
    NumericExpr e, int var1, int var2, int var3) {
  ampl::Problem p;
  p.AddVar(min(), max(), ampl::INTEGER);
  p.AddVar(var1, var1, ampl::INTEGER);
  p.AddVar(var2, var2, ampl::INTEGER);
  p.AddVar(var3, var3, ampl::INTEGER);
  p.AddCon(AddRelational(EQ, AddVar(0), e));
  NLToJaCoPConverter converter(jvm_);
  converter.Convert(p);
  return Solve(converter);
}

double JaCoPConverterTest::ConvertAndEval(
    LogicalExpr e, int var1, int var2, int var3) {
  // TODO: handle AllDiff specially
  return ConvertAndEval(AddIf(e, AddNum(1), AddNum(0)), var1, var2, var3);
}

TEST_F(JaCoPConverterTest, ConvertPlus) {
  NumericExpr e = AddBinary(OPPLUS, x, y);
  EXPECT_EQ(25, ConvertAndEval(e, 10, 15));
  EXPECT_EQ(12, ConvertAndEval(e, 19, -7));
  EXPECT_NE(0, isnan(ConvertAndEval(e, max(), 1)));
}

TEST_F(JaCoPConverterTest, ConvertMinus) {
  NumericExpr e = AddBinary(OPMINUS, x, y);
  EXPECT_EQ(-5, ConvertAndEval(e, 10, 15));
  EXPECT_EQ(26, ConvertAndEval(e, 19, -7));
  EXPECT_NE(0, isnan(ConvertAndEval(e, min(), 1)));
}

TEST_F(JaCoPConverterTest, ConvertMult) {
  NumericExpr e = AddBinary(OPMULT, x, y);
  EXPECT_EQ(150, ConvertAndEval(e, 10, 15));
  EXPECT_EQ(-133, ConvertAndEval(e, 19, -7));
  EXPECT_NE(0, isnan(ConvertAndEval(e, max(), 2)));
}

TEST_F(JaCoPConverterTest, ConvertDiv) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OPDIV, x, y)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertRem) {
  NumericExpr e = AddBinary(OPREM, x, y);
  EXPECT_EQ(0, ConvertAndEval(e, 9, 3));
  EXPECT_EQ(2, ConvertAndEval(e, 8, 3));
  EXPECT_EQ(-2, ConvertAndEval(e, -8, 3));
  EXPECT_EQ(2, ConvertAndEval(e, 8, -3));
  EXPECT_EQ(-2, ConvertAndEval(e, -8, -3));
}

TEST_F(JaCoPConverterTest, ConvertPow) {
  NumericExpr e = AddBinary(OPPOW, x, y);
  EXPECT_EQ(8, ConvertAndEval(e, 2, 3));
  EXPECT_EQ(81, ConvertAndEval(e, 3, 4));
}

TEST_F(JaCoPConverterTest, ConvertNumericLess) {
  NumericExpr e = AddBinary(OPLESS, x, y);
  EXPECT_EQ(0, ConvertAndEval(e, 10, 15));
  EXPECT_EQ(26, ConvertAndEval(e, 19, -7));
}

TEST_F(JaCoPConverterTest, ConvertMin) {
  NumericExpr e = AddVarArg(MINLIST, x, y, z);
  EXPECT_EQ(-7, ConvertAndEval(e, 3, -7, 5));
  EXPECT_EQ(10, ConvertAndEval(e, 10, 20, 30));
  EXPECT_THROW(ConvertAndEval(AddVarArg(MINLIST)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertMax) {
  NumericExpr e = AddVarArg(MAXLIST, x, y, z);
  EXPECT_EQ(5, ConvertAndEval(e, 3, -7, 5));
  EXPECT_EQ(30, ConvertAndEval(e, 30, 20, 10));
  EXPECT_THROW(ConvertAndEval(AddVarArg(MAXLIST)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertFloor) {
  NumericExpr e = AddUnary(FLOOR, x);
  EXPECT_EQ(-42, ConvertAndEval(e, -42));
  EXPECT_EQ(42, ConvertAndEval(e, 42));
  EXPECT_EQ(6, ConvertAndEval(AddUnary(FLOOR, AddUnary(OP_sqrt, x)), 42));
}

TEST_F(JaCoPConverterTest, ConvertCeil) {
  NumericExpr e = AddUnary(CEIL, x);
  EXPECT_EQ(-42, ConvertAndEval(e, -42));
  EXPECT_EQ(42, ConvertAndEval(e, 42));
}

TEST_F(JaCoPConverterTest, ConvertAbs) {
  NumericExpr e = AddUnary(ABS, x);
  EXPECT_EQ(42, ConvertAndEval(e, -42));
  EXPECT_EQ(42, ConvertAndEval(e, 42));
}

TEST_F(JaCoPConverterTest, ConvertUnaryMinus) {
  NumericExpr e = AddUnary(OPUMINUS, x);
  EXPECT_EQ(42, ConvertAndEval(e, -42));
  EXPECT_EQ(-42, ConvertAndEval(e, 42));
}

TEST_F(JaCoPConverterTest, ConvertIf) {
  NumericExpr e = AddIf(AddRelational(EQ, x, AddNum(1)), y, z);
  EXPECT_EQ(42, ConvertAndEval(e, 1, 42, 10));
  EXPECT_EQ(10, ConvertAndEval(e, 0, 42, 10));
}

TEST_F(JaCoPConverterTest, ConvertTanh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_tanh, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertTan) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_tan, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertSqrt) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_sqrt, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertSinh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_sinh, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertSin) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_sin, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertLog10) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_log10, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertLog) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_log, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertExp) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_exp, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertCosh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_cosh, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertCos) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_cos, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertAtanh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_atanh, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertAtan2) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OP_atan2, x, y)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertAtan) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_atan, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertAsinh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_asinh, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertAsin) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_asin, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertAcosh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_acosh, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertAcos) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_acos, x)), UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertSum) {
  EXPECT_EQ(0, ConvertAndEval(AddSum()));
  EXPECT_EQ(42, ConvertAndEval(AddSum(x), 42));
  EXPECT_EQ(123, ConvertAndEval(AddSum(x, y, z), 100, 20, 3));
}

TEST_F(JaCoPConverterTest, ConvertIntDiv) {
  NumericExpr e = AddBinary(OPintDIV, x, y);
  EXPECT_EQ(3, ConvertAndEval(e, 9, 3));
  EXPECT_EQ(2, ConvertAndEval(e, 8, 3));
  EXPECT_EQ(-2, ConvertAndEval(e, -8, 3));
  EXPECT_EQ(-2, ConvertAndEval(e, 8, -3));
  EXPECT_EQ(2, ConvertAndEval(e, -8, -3));
}

TEST_F(JaCoPConverterTest, ConvertPrecision) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OPprecision, x, y)),
      UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertRound) {
  EXPECT_EQ(42, ConvertAndEval(AddBinary(OPround, x, AddNum(0)), 42));
  EXPECT_THROW(ConvertAndEval(AddBinary(OPround, x, AddNum(1))),
      UnsupportedExprError);
  EXPECT_THROW(ConvertAndEval(AddBinary(OPround, x, y)),
        UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertTrunc) {
  EXPECT_EQ(42, ConvertAndEval(AddBinary(OPtrunc, x, AddNum(0)), 42));
  EXPECT_THROW(ConvertAndEval(AddBinary(OPtrunc, x, AddNum(1))),
      UnsupportedExprError);
  EXPECT_THROW(ConvertAndEval(AddBinary(OPtrunc, x, y)),
        UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertCount) {
  LogicalExpr a(AddRelational(NE, x, AddNum(0)));
  LogicalExpr b(AddRelational(NE, y, AddNum(0)));
  LogicalExpr c(AddRelational(NE, z, AddNum(0)));
  EXPECT_EQ(0, ConvertAndEval(AddCount(a, b, c)));
  EXPECT_EQ(1, ConvertAndEval(AddCount(a, b, c), 1));
  EXPECT_EQ(2, ConvertAndEval(AddCount(a, b, c), 0, 1, 1));
  EXPECT_EQ(3, ConvertAndEval(AddCount(a, b, c), 1, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertNumberOf) {
  ampl::NumericConstant val = AddNum(42);
  EXPECT_EQ(0, ConvertAndEval(AddNumberOf(val, x)));
  EXPECT_EQ(1, ConvertAndEval(AddNumberOf(val, x), 42));
  EXPECT_EQ(0, ConvertAndEval(AddNumberOf(val, x, y)));
  EXPECT_EQ(1, ConvertAndEval(AddNumberOf(val, x, y), 0, 42));
  EXPECT_EQ(2, ConvertAndEval(AddNumberOf(val, x, y), 42, 42));
  EXPECT_EQ(3, ConvertAndEval(AddBinary(OPPLUS,
      AddNumberOf(val, x, y), AddNumberOf(AddNum(11), y, z)),
      42, 42, 11));
}

TEST_F(JaCoPConverterTest, ConvertPLTerm) {
  double args[] = {-1, 5, 0, 10, 1};
  EXPECT_THROW(ConvertAndEval(AddPLTerm(5, args, 0)),
      UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertPowConstExp) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OP1POW, x, AddNum(42))),
      UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertPow2) {
  EXPECT_EQ(49, ConvertAndEval(AddUnary(OP2POW, x), 7));
}

TEST_F(JaCoPConverterTest, ConvertPowConstBase) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OPCPOW, AddNum(42), x)),
      UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertNum) {
  EXPECT_EQ(42, ConvertAndEval(AddNum(42)));
  std::string message;
  try {
    ConvertAndEval(AddNum(0.42));
  } catch (const ampl::Error &e) {
    message = e.what();
  }
  EXPECT_EQ("value 0.42 can't be represented as int", message);
  EXPECT_EQ(min(), ConvertAndEval(AddNum(min())));
  // TODO
  //EXPECT_THROW(ConvertAndEval(AddNum(min() - 1)), Gecode::Int::OutOfLimits);
  EXPECT_EQ(max(), ConvertAndEval(AddNum(max())));
  // TODO
  //EXPECT_THROW(ConvertAndEval(AddNum(max() + 1)), Gecode::Int::OutOfLimits);
}

TEST_F(JaCoPConverterTest, ConvertVar) {
  EXPECT_EQ(11, ConvertAndEval(x, 11, 22));
  EXPECT_EQ(22, ConvertAndEval(y, 11, 22));
  EXPECT_EQ(33, ConvertAndEval(x, 33));
}

TEST_F(JaCoPConverterTest, ConvertOr) {
  NumericExpr one = AddNum(1);
  LogicalExpr e = AddBinaryLogical(
      OPOR, AddRelational(EQ, x, one), AddRelational(EQ, y, one));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertAnd) {
  NumericExpr one = AddNum(1);
  LogicalExpr e = AddBinaryLogical(
      OPAND, AddRelational(EQ, x, one), AddRelational(EQ, y, one));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertLess) {
  LogicalExpr e = AddRelational(LT, x, y);
  EXPECT_EQ(0, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(1, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(0, ConvertAndEval(e, 5, 3));
}

TEST_F(JaCoPConverterTest, ConvertLessEqual) {
  LogicalExpr e = AddRelational(LE, x, y);
  EXPECT_EQ(1, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(1, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(0, ConvertAndEval(e, 5, 3));
}

TEST_F(JaCoPConverterTest, ConvertEqual) {
  LogicalExpr e = AddRelational(EQ, x, y);
  EXPECT_EQ(1, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(0, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(0, ConvertAndEval(e, 5, 3));
}

TEST_F(JaCoPConverterTest, ConvertGreaterEqual) {
  LogicalExpr e = AddRelational(GE, x, y);
  EXPECT_EQ(1, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(0, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(1, ConvertAndEval(e, 5, 3));
}

TEST_F(JaCoPConverterTest, ConvertGreater) {
  LogicalExpr e = AddRelational(GT, x, y);
  EXPECT_EQ(0, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(0, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(1, ConvertAndEval(e, 5, 3));
}

TEST_F(JaCoPConverterTest, ConvertNotEqual) {
  LogicalExpr e = AddRelational(NE, x, y);
  EXPECT_EQ(0, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(1, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(1, ConvertAndEval(e, 5, 3));
}

TEST_F(JaCoPConverterTest, ConvertNot) {
  LogicalExpr e = AddNot(AddRelational(EQ, x, AddNum(1)));
  EXPECT_EQ(1, ConvertAndEval(e, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 1));
}

TEST_F(JaCoPConverterTest, ConvertAtLeast) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPATLEAST, AddVar(1), AddCount(a, b));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 1, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 2, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 2, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertAtMost) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPATMOST, AddVar(1), AddCount(a, b));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 1, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 2, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 2, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertExactly) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPEXACTLY, AddVar(1), AddCount(a, b));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 1, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 1, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 2, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 2, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertNotAtLeast) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPNOTATLEAST, AddVar(1), AddCount(a, b));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 1, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 2, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 2, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertNotAtMost) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPNOTATMOST, AddVar(1), AddCount(a, b));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 1, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 2, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 2, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertNotExactly) {
  LogicalExpr a(AddRelational(NE, y, AddNum(0)));
  LogicalExpr b(AddRelational(NE, z, AddNum(0)));
  LogicalExpr e = AddLogicalCount(OPNOTEXACTLY, AddVar(1), AddCount(a, b));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 2, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 2, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertForAll) {
  LogicalExpr e = AddIteratedLogical(ANDLIST,
      AddRelational(EQ, x, AddNum(1)),
      AddRelational(EQ, y, AddNum(1)),
      AddRelational(EQ, z, AddNum(1)));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 1, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 1, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertExists) {
  LogicalExpr e = AddIteratedLogical(ORLIST,
      AddRelational(EQ, x, AddNum(1)),
      AddRelational(EQ, y, AddNum(1)),
      AddRelational(EQ, z, AddNum(1)));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 1, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertImplication) {
  LogicalExpr e = AddImplication(
      AddRelational(EQ, x, AddNum(1)),
      AddRelational(EQ, y, AddNum(1)),
      AddRelational(EQ, z, AddNum(1)));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 1, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1, 1));
  e = AddImplication(
      AddRelational(EQ, x, AddNum(1)),
      AddRelational(EQ, y, AddNum(1)),
      AddBool(false));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertIff) {
  LogicalExpr e = AddBinaryLogical(OP_IFF,
      AddRelational(EQ, x, AddNum(1)),
      AddRelational(EQ, y, AddNum(1)));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertAllDiff) {
  LogicalExpr e = AddAllDiff(AddNum(1), x, y);
  EXPECT_EQ(1, ConvertAndEval(e, 2, 3));
  EXPECT_EQ(0, ConvertAndEval(e, 2, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 1));
}

TEST_F(JaCoPConverterTest, ConvertNestedAllDiff) {
  EXPECT_THROW(ConvertAndEval(AddNot(AddAllDiff(AddNum(1), x, y)), 1, 2),
      UnsupportedExprError);
}

TEST_F(JaCoPConverterTest, ConvertLogicalConstant) {
  EXPECT_EQ(0, ConvertAndEval(AddBool(false)));
  EXPECT_EQ(1, ConvertAndEval(AddBool(1)));
}

// ----------------------------------------------------------------------------
// Solver tests

// TODO
/*class JaCoPSolverTest : public ::testing::Test, public ExprBuilder {
 protected:
  ampl::GecodeSolver solver_;

  class ParseResult {
   private:
    bool result_;
    std::string error_;

    void True() const {}
    typedef void (ParseResult::*SafeBool)() const;

   public:
    ParseResult(bool result, std::string error)
    : result_(result), error_(error) {}

    operator SafeBool() const { return result_ ? &ParseResult::True : 0; }

    std::string error() const {
      EXPECT_FALSE(result_);
      return error_;
    }
  };

  struct TestErrorHandler : ampl::ErrorHandler {
    std::string error;
    void HandleError(fmt::StringRef message) {
      error += message.c_str();
    }
  };

  ParseResult ParseOptions(const char *opt1, const char *opt2 = nullptr) {
    TestErrorHandler eh;
    solver_.set_error_handler(&eh);
    bool result = solver_.ParseOptions(
        Args(opt1, opt2), solver_, ampl::BasicSolver::NO_OPTION_ECHO);
    if (result)
      EXPECT_EQ("", eh.error);
    return ParseResult(result, eh.error);
  }

  int RunSolver(const char *stub = nullptr, const char *opt1 = nullptr,
      const char *opt2 = nullptr, const char *opt3 = nullptr) {
    return solver_.Run(Args("gecode", "-s", stub, opt1, opt2, opt3));
  }

  SolveResult Solve(const char *stub, const char *opt1 = nullptr,
      const char *opt2 = nullptr, const char *opt3 = nullptr);
};

SolveResult JaCoPSolverTest::Solve(const char *stub,
    const char *opt1, const char *opt2, const char *opt3) {
  TestSolutionHandler sh;
  solver_.set_solution_handler(&sh);
  RunSolver(stub, opt1, opt2, opt3);
  const string &message = sh.message();
  int solve_code = sh.solve_code();
  EXPECT_GE(solve_code, 0);
  bool solved = true;
  if (solve_code < 100)
    EXPECT_TRUE(message.find("optimal solution") != string::npos);
  else if (solve_code < 200)
    EXPECT_TRUE(message.find("feasible solution") != string::npos);
  else
    solved = false;
  return SolveResult(solved, sh.obj_value(), message);
}

TEST_F(JaCoPSolverTest, ObjConst) {
  EXPECT_EQ(42, Solve(DATA_DIR "objconstint").obj);
}

TEST_F(JaCoPSolverTest, ContinuousVarsNotSupported) {
  EXPECT_THROW(RunSolver(DATA_DIR "objconst"), std::runtime_error);
}

TEST_F(JaCoPSolverTest, SolveAssign0) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign0").obj);
}

TEST_F(JaCoPSolverTest, SolveAssign1) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign1").obj);
}

// Disabled because it takes too long to solve.
TEST_F(JaCoPSolverTest, DISABLED_SolveBalassign0) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign0").obj);
}

// Disabled because it takes too long to solve.
TEST_F(JaCoPSolverTest, DISABLED_SolveBalassign1) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign1").obj);
}

TEST_F(JaCoPSolverTest, SolveFlowshp0) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp0").obj);
}

TEST_F(JaCoPSolverTest, SolveFlowshp1) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp1").obj);
}

TEST_F(JaCoPSolverTest, SolveFlowshp2) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp2").obj);
}

TEST_F(JaCoPSolverTest, SolveGrpassign0) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign0").obj);
}

// Disabled because variables in subscripts are not yet allowed.
TEST_F(JaCoPSolverTest, DISABLED_SolveGrpassign1) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1").obj);
}

// Disabled because object-valued variables are not yet allowed.
TEST_F(JaCoPSolverTest, DISABLED_SolveGrpassign1a) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1a").obj);
}

TEST_F(JaCoPSolverTest, SolveMagic) {
  EXPECT_TRUE(Solve(DATA_DIR "magic").solved);
}

TEST_F(JaCoPSolverTest, SolveMapcoloring) {
  EXPECT_TRUE(Solve(DATA_DIR "mapcoloring").solved);
}

TEST_F(JaCoPSolverTest, SolveNQueens) {
  EXPECT_TRUE(Solve(DATA_DIR "nqueens").solved);
}

TEST_F(JaCoPSolverTest, SolveNQueens0) {
  EXPECT_TRUE(Solve(DATA_DIR "nqueens0").solved);
}

// Disabled because it takes somewhat long (compared to other tests).
TEST_F(JaCoPSolverTest, DISABLED_SolveOpenShop) {
  EXPECT_EQ(1955, Solve(DATA_DIR "openshop").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(JaCoPSolverTest, DISABLED_SolveParty1) {
  EXPECT_EQ(61, Solve(DATA_DIR "party1").obj);
}

// Disabled because Gecode doesn't support 'alldiff' as a subexpression.
TEST_F(JaCoPSolverTest, DISABLED_SolveParty2) {
  EXPECT_EQ(3, Solve(DATA_DIR "party2").obj);
}

// Disabled because it takes somewhat long (compared to other tests).
TEST_F(JaCoPSolverTest, DISABLED_SolvePhoto9) {
  EXPECT_EQ(10, Solve(DATA_DIR "photo9").obj);
}

// Disabled because it takes somewhat long (compared to other tests).
TEST_F(JaCoPSolverTest, DISABLED_SolvePhoto11) {
  EXPECT_EQ(12, Solve(DATA_DIR "photo11").obj);
}

TEST_F(JaCoPSolverTest, SolveSched0) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched0").obj);
}

TEST_F(JaCoPSolverTest, SolveSched1) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched1").obj);
}

TEST_F(JaCoPSolverTest, SolveSched2) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched2").obj);
}

TEST_F(JaCoPSolverTest, SolveSendMoreMoney) {
  EXPECT_TRUE(Solve(DATA_DIR "send-more-money").solved);
}

TEST_F(JaCoPSolverTest, SolveSendMostMoney) {
  EXPECT_NEAR(10876, Solve(DATA_DIR "send-most-money").obj, 1e-5);
}

TEST_F(JaCoPSolverTest, SolveSeq0) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0").obj, 1e-5);
}

TEST_F(JaCoPSolverTest, SolveSeq0a) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0a").obj, 1e-5);
}

TEST_F(JaCoPSolverTest, SolveSudokuHard) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuHard").solved);
}

TEST_F(JaCoPSolverTest, SolveSudokuVeryEasy) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuVeryEasy").solved);
}

// ----------------------------------------------------------------------------
// Solve code tests

TEST_F(JaCoPSolverTest, OptimalSolveCode) {
  Solve(DATA_DIR "objconstint");
  EXPECT_EQ(0, solver_.problem().solve_code());
}

TEST_F(JaCoPSolverTest, FeasibleSolveCode) {
  Solve(DATA_DIR "feasible");
  EXPECT_EQ(100, solver_.problem().solve_code());
}

TEST_F(JaCoPSolverTest, InfeasibleSolveCode) {
  Solve(DATA_DIR "infeasible");
  EXPECT_EQ(200, solver_.problem().solve_code());
}

// ----------------------------------------------------------------------------
// Interrupt tests

#ifdef HAVE_THREADS
void Interrupt() {
  // Wait until started.
  while (ampl::SignalHandler::stop())
    std::this_thread::yield();
  std::raise(SIGINT);
}

TEST_F(JaCoPSolverTest, InterruptSolution) {
  std::thread t(Interrupt);
  std::string message = Solve(DATA_DIR "miplib/assign1").message;
  t.join();
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find("interrupted") != string::npos);
}
#endif

// ----------------------------------------------------------------------------
// Option tests

TEST_F(JaCoPSolverTest, VersionOption) {
  EXPECT_FALSE((solver_.flags() & ASL_OI_show_version) != 0);
  EXPECT_TRUE(ParseOptions("version"));
  EXPECT_TRUE((solver_.flags() & ASL_OI_show_version) != 0);
}

TEST_F(JaCoPSolverTest, WantsolOption) {
  EXPECT_EQ(0, solver_.wantsol());
  EXPECT_TRUE(ParseOptions("wantsol=1"));
  EXPECT_EQ(1, solver_.wantsol());
  EXPECT_TRUE(ParseOptions("wantsol=5"));
  EXPECT_EQ(5, solver_.wantsol());
}

TEST_F(JaCoPSolverTest, ADOption) {
  EXPECT_EQ(Gecode::Search::Options().a_d, solver_.options().a_d);
  EXPECT_TRUE(ParseOptions("a_d=42"));
  EXPECT_EQ(42u, solver_.options().a_d);
  EXPECT_EQ("Invalid value -1 for option a_d", ParseOptions("a_d=-1").error());
}

TEST_F(JaCoPSolverTest, CDOption) {
  EXPECT_EQ(Gecode::Search::Options().c_d, solver_.options().c_d);
  EXPECT_TRUE(ParseOptions("c_d=42"));
  EXPECT_EQ(42u, solver_.options().c_d);
  EXPECT_EQ("Invalid value -1 for option c_d", ParseOptions("c_d=-1").error());
}

TEST_F(JaCoPSolverTest, FailLimitOption) {
  std::string message =
      Solve(DATA_DIR "miplib/assign1", "faillimit=10").message;
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find(" 11 fails") != string::npos);
  EXPECT_EQ("Invalid value -1 for option faillimit",
      ParseOptions("faillimit=-1").error());
}

TEST_F(JaCoPSolverTest, MemoryLimitOption) {
  Solve(DATA_DIR "miplib/assign1", "memorylimit=1000000");
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_EQ("Invalid value -1 for option memorylimit",
      ParseOptions("memorylimit=-1").error());
}

TEST_F(JaCoPSolverTest, NodeLimitOption) {
  std::string message =
      Solve(DATA_DIR "miplib/assign1", "nodelimit=10").message;
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_TRUE(message.find("11 nodes") != string::npos);
  EXPECT_EQ("Invalid value -1 for option nodelimit",
      ParseOptions("nodelimit=-1").error());
}

TEST_F(JaCoPSolverTest, TimeLimitOption) {
  Solve(DATA_DIR "miplib/assign1", "timelimit=0.1");
  EXPECT_EQ(600, solver_.problem().solve_code());
  EXPECT_EQ("Invalid value -1 for option timelimit",
      ParseOptions("timelimit=-1").error());
}

TEST_F(JaCoPSolverTest, ThreadsOption) {
  EXPECT_EQ(Gecode::Search::Options().threads, solver_.options().threads);
  EXPECT_TRUE(ParseOptions("threads=0.5"));
  EXPECT_EQ(0.5, solver_.options().threads);
  EXPECT_TRUE(ParseOptions("threads=-10"));
  EXPECT_EQ(-10.0, solver_.options().threads);
}

template <typename T>
struct OptionValue {
  const char *name;
  T value;
};

const OptionValue<Gecode::IntValBranch> VAL_BRANCHINGS[] = {
    {"min",        Gecode::INT_VAL_MIN},
    {"med",        Gecode::INT_VAL_MED},
    {"max",        Gecode::INT_VAL_MAX},
    {"rnd",        Gecode::INT_VAL_RND},
    {"split_min",  Gecode::INT_VAL_SPLIT_MIN},
    {"split_max",  Gecode::INT_VAL_SPLIT_MAX},
    {"range_min",  Gecode::INT_VAL_RANGE_MIN},
    {"range_max",  Gecode::INT_VAL_RANGE_MAX},
    {"values_min", Gecode::INT_VALUES_MIN},
    {"values_max", Gecode::INT_VALUES_MAX},
    {}
};

TEST_F(JaCoPSolverTest, ValBranchingOption) {
  EXPECT_EQ(Gecode::INT_VAL_MIN, solver_.val_branching());
  unsigned count = 0;
  for (const OptionValue<Gecode::IntValBranch>
      *p = VAL_BRANCHINGS; p->name; ++p, ++count) {
    EXPECT_TRUE(ParseOptions(
        c_str(fmt::Format("val_branching={}") << p->name)));
    EXPECT_EQ(p->value, solver_.val_branching());
  }
  EXPECT_EQ(10u, count);
}

const OptionValue<Gecode::IntVarBranch> VAR_BRANCHINGS[] = {
    {"none",            Gecode::INT_VAR_NONE},
    {"rnd",             Gecode::INT_VAR_RND},
    {"degree_min",      Gecode::INT_VAR_DEGREE_MIN},
    {"degree_max",      Gecode::INT_VAR_DEGREE_MAX},
    {"afc_min",         Gecode::INT_VAR_AFC_MIN},
    {"afc_max",         Gecode::INT_VAR_AFC_MAX},
    {"min_min",         Gecode::INT_VAR_MIN_MIN},
    {"min_max",         Gecode::INT_VAR_MIN_MAX},
    {"max_min",         Gecode::INT_VAR_MAX_MIN},
    {"max_max",         Gecode::INT_VAR_MAX_MAX},
    {"size_min",        Gecode::INT_VAR_SIZE_MIN},
    {"size_max",        Gecode::INT_VAR_SIZE_MAX},
    {"size_degree_min", Gecode::INT_VAR_SIZE_DEGREE_MIN},
    {"size_degree_max", Gecode::INT_VAR_SIZE_DEGREE_MAX},
    {"size_afc_min",    Gecode::INT_VAR_SIZE_AFC_MIN},
    {"size_afc_max",    Gecode::INT_VAR_SIZE_AFC_MAX},
    {"regret_min_min",  Gecode::INT_VAR_REGRET_MIN_MIN},
    {"regret_min_max",  Gecode::INT_VAR_REGRET_MIN_MAX},
    {"regret_max_min",  Gecode::INT_VAR_REGRET_MAX_MIN},
    {"regret_max_max",  Gecode::INT_VAR_REGRET_MAX_MAX},
    {}
};

TEST_F(JaCoPSolverTest, VarBranchingOption) {
  EXPECT_EQ(Gecode::INT_VAR_SIZE_MIN, solver_.var_branching());
  unsigned count = 0;
  for (const OptionValue<Gecode::IntVarBranch>
      *p = VAR_BRANCHINGS; p->name; ++p, ++count) {
    EXPECT_TRUE(ParseOptions(
        c_str(fmt::Format("var_branching={}") << p->name)));
    EXPECT_EQ(p->value, solver_.var_branching());
  }
  EXPECT_EQ(20u, count);
}

TEST_F(JaCoPSolverTest, OutLevOption) {
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve(DATA_DIR "objconstint");
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("", ReadFile("out"));

  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve(DATA_DIR "objconstint", "outlev=1");
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  EXPECT_EQ("outlev=1\n"
      " Max Depth      Nodes      Fails      Best Obj\n"
      "                                            42\n", ReadFile("out"));

  EXPECT_TRUE(ParseOptions("outlev=0"));
  EXPECT_TRUE(ParseOptions("outlev=1"));
  EXPECT_EQ("Invalid value -1 for option outlev",
      ParseOptions("outlev=-1").error());
  EXPECT_EQ("Invalid value 2 for option outlev",
      ParseOptions("outlev=2").error());
}

TEST_F(JaCoPSolverTest, OutFreqOption) {
  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve(DATA_DIR "party1", "outlev=1", "outfreq=1", "timelimit=2.5");
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  std::string out = ReadFile("out");
  EXPECT_EQ(6, std::count(out.begin(), out.end(), '\n'));

  EXPECT_EXIT({
    FILE *f = freopen("out", "w", stdout);
    Solve(DATA_DIR "party1", "outlev=1", "outfreq=2", "timelimit=2.5");
    fclose(f);
    exit(0);
  }, ::testing::ExitedWithCode(0), "");
  out = ReadFile("out");
  EXPECT_EQ(5, std::count(out.begin(), out.end(), '\n'));

  EXPECT_EQ("Invalid value -1 for option outfreq",
      ParseOptions("outfreq=-1").error());
  EXPECT_EQ("Invalid value 0 for option outfreq",
      ParseOptions("outfreq=0").error());
}*/
}

/*
 Gecode driver tests.

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

#include <gecode/search.hh>

#include <algorithm>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>

#include "gtest/gtest.h"

#include "solvers/gecode/gecode.h"
#include "solvers/util/expr.h"

extern "C" {
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
}

#include "tests/expr_builder.h"
#include "tests/config.h"

using std::ifstream;
using std::size_t;
using std::string;
using std::vector;

using ampl::Variable;
using ampl::ExprBuilder;
using ampl::GecodeProblem;
using ampl::LogicalExpr;
using ampl::NumericExpr;
using ampl::UnsupportedExprError;

#define DATA_DIR "../data/"

namespace {

bool AreBothSpaces(char lhs, char rhs) { return lhs == ' ' && rhs == ' '; }

// Replace all occurrences of string old_s in s with new_s.
void Replace(string &s, const string &old_s, const string &new_s) {
  size_t pos = 0;
  while ((pos = s.find(old_s, pos)) != string::npos) {
    s.replace(pos, old_s.length(), new_s);
    pos += new_s.length();
  }
}

struct SolveResult {
  bool solved;
  double obj;
  SolveResult(bool solved, double obj) : solved(solved), obj(obj) {}
};

// Helper class that copies arguments to comply with the IlogCPDriver::run
// function signature and avoid unwanted modification.
class Args {
 private:
  std::size_t argc_;
  vector<char> store_;
  vector<char*> argv_;

 public:
  Args() : argc_(0) {}

  char **get() {
    argv_.resize(argc_ + 1);
    for (std::size_t i = 0, j = 0; i < argc_; j += strlen(&store_[j]) + 1, ++i)
      argv_[i] = &store_[j];
    return &argv_[0];
  }

  Args& operator+(const char *arg) {
    if (arg) {
      ++argc_;
      store_.insert(store_.end(), arg, arg + strlen(arg) + 1);
    }
    return *this;
  }
};

static void InitVars(GecodeProblem &p, int var1, int var2, int var3) {
  Gecode::IntVarArray &vars = p.vars();
  vars[0] = Gecode::IntVar(p,
              Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  vars[1] = Gecode::IntVar(p, var1, var1);
  vars[2] = Gecode::IntVar(p, var2, var2);
  vars[3] = Gecode::IntVar(p, var3, var3);
}

static double Solve(GecodeProblem &p) {
  Gecode::DFS<GecodeProblem> engine(&p);
  std::auto_ptr<GecodeProblem> solution(engine.next());
  return solution.get() ?
      solution->vars()[0].val() : std::numeric_limits<double>::quiet_NaN();
}

// Converts a numeric expression from AMPL to Gecode form and evaluates it.
// Returns the value of the Gecode expression.
static double ConvertAndEval(NumericExpr e,
    int var1 = 0, int var2 = 0, int var3 = 0) {
  ampl::NLToGecodeConverter converter(4);
  GecodeProblem &p = converter.problem();
  InitVars(p, var1, var2, var3);
  Gecode::rel(p, converter.problem().vars()[0] == converter.Visit(e));
  return Solve(p);
}

static double ConvertAndEval(LogicalExpr e,
    int var1 = 0, int var2 = 0, int var3 = 0) {
  ampl::NLToGecodeConverter converter(4);
  GecodeProblem &p = converter.problem();
  InitVars(p, var1, var2, var3);
  Gecode::BoolVar result(p, 0, 1);
  Gecode::BoolExpr expr = converter.Visit(e);
  if (!ampl::Cast<ampl::AllDiffExpr>(e)) {
    Gecode::rel(p, result == expr);
    Gecode::channel(p, result, p.vars()[0]);
  } else {
    Gecode::DFS<GecodeProblem> engine(&p);
    std::auto_ptr<GecodeProblem> solution(engine.next());
    return solution.get() ? 1 : 0;
  }
  return Solve(p);
}

class GecodeTest : public ::testing::Test, public ExprBuilder {
 protected:
  Variable x;
  Variable y;
  Variable z;

  GecodeTest() : x(AddVar(1)), y(AddVar(2)), z(AddVar(3)) {}

  int RunDriver(const char *stub = nullptr, const char *opt = nullptr) {
    //return p.run((Args() + "gecode" + "-s" + stub + opt).get());
    // TODO
    return 0;
  }

  SolveResult Solve(const char *stub, const char *opt = nullptr);
};

SolveResult GecodeTest::Solve(const char *stub, const char *opt) {
  RunDriver(stub, opt);
  ifstream ifs((string(stub) + ".sol").c_str());
  string line;
  getline(ifs, line);
  bool solved = line.find("optimal solution") != string::npos;
  if (!solved) solved = line.find("feasible solution") != string::npos;
  getline(ifs, line);
  const char obj[] = "objective ";
  size_t pos = line.find(obj);
  return SolveResult(solved, pos != string::npos ?
      atof(line.c_str() + pos + sizeof(obj) - 1) :
      std::numeric_limits<double>::quiet_NaN());
}

TEST_F(GecodeTest, ConvertPlus) {
  NumericExpr e = AddBinary(OPPLUS, x, y);
  EXPECT_EQ(25, ConvertAndEval(e, 10, 15));
  EXPECT_EQ(12, ConvertAndEval(e, 19, -7));
  EXPECT_NE(0, isnan(ConvertAndEval(e, Gecode::Int::Limits::max, 1)));
}

TEST_F(GecodeTest, ConvertMinus) {
  NumericExpr e = AddBinary(OPMINUS, x, y);
  EXPECT_EQ(-5, ConvertAndEval(e, 10, 15));
  EXPECT_EQ(26, ConvertAndEval(e, 19, -7));
  EXPECT_NE(0, isnan(ConvertAndEval(e, Gecode::Int::Limits::min, 1)));
}

TEST_F(GecodeTest, ConvertMult) {
  NumericExpr e = AddBinary(OPMULT, x, y);
  EXPECT_EQ(150, ConvertAndEval(e, 10, 15));
  EXPECT_EQ(-133, ConvertAndEval(e, 19, -7));
  EXPECT_NE(0, isnan(ConvertAndEval(e, Gecode::Int::Limits::max, 2)));
}

TEST_F(GecodeTest, ConvertDiv) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OPDIV, x, y)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertRem) {
  NumericExpr e = AddBinary(OPREM, x, y);
  EXPECT_EQ(0, ConvertAndEval(e, 9, 3));
  EXPECT_EQ(2, ConvertAndEval(e, 8, 3));
  EXPECT_EQ(-2, ConvertAndEval(e, -8, 3));
  EXPECT_EQ(2, ConvertAndEval(e, 8, -3));
  EXPECT_EQ(-2, ConvertAndEval(e, -8, -3));
}

TEST_F(GecodeTest, ConvertPow) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OPPOW, x, y)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertNumericLess) {
  NumericExpr e = AddBinary(OPLESS, x, y);
  EXPECT_EQ(0, ConvertAndEval(e, 10, 15));
  EXPECT_EQ(26, ConvertAndEval(e, 19, -7));
}

TEST_F(GecodeTest, ConvertMin) {
  NumericExpr e = AddVarArg(MINLIST, x, y, z);
  EXPECT_EQ(-7, ConvertAndEval(e, 3, -7, 5));
  EXPECT_EQ(10, ConvertAndEval(e, 10, 20, 30));
  EXPECT_THROW(ConvertAndEval(AddVarArg(MINLIST)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertMax) {
  NumericExpr e = AddVarArg(MAXLIST, x, y, z);
  EXPECT_EQ(5, ConvertAndEval(e, 3, -7, 5));
  EXPECT_EQ(30, ConvertAndEval(e, 30, 20, 10));
  EXPECT_THROW(ConvertAndEval(AddVarArg(MAXLIST)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertFloor) {
  NumericExpr e = AddUnary(FLOOR, x);
  EXPECT_EQ(-42, ConvertAndEval(e, -42));
  EXPECT_EQ(42, ConvertAndEval(e, 42));
  EXPECT_EQ(6, ConvertAndEval(AddUnary(FLOOR, AddUnary(OP_sqrt, x)), 42));
}

TEST_F(GecodeTest, ConvertCeil) {
  NumericExpr e = AddUnary(CEIL, x);
  EXPECT_EQ(-42, ConvertAndEval(e, -42));
  EXPECT_EQ(42, ConvertAndEval(e, 42));
}

TEST_F(GecodeTest, ConvertAbs) {
  NumericExpr e = AddUnary(ABS, x);
  EXPECT_EQ(42, ConvertAndEval(e, -42));
  EXPECT_EQ(42, ConvertAndEval(e, 42));
}

TEST_F(GecodeTest, ConvertUnaryMinus) {
  NumericExpr e = AddUnary(OPUMINUS, x);
  EXPECT_EQ(42, ConvertAndEval(e, -42));
  EXPECT_EQ(-42, ConvertAndEval(e, 42));
}

TEST_F(GecodeTest, ConvertIf) {
  NumericExpr e = AddIf(AddRelational(EQ, x, AddNum(1)), y, z);
  EXPECT_EQ(42, ConvertAndEval(e, 1, 42, 10));
  EXPECT_EQ(10, ConvertAndEval(e, 0, 42, 10));
}

TEST_F(GecodeTest, ConvertTanh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_tanh, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertTan) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_tan, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertSqrt) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_sqrt, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertSinh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_sinh, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertSin) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_sin, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertLog10) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_log10, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertLog) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_log, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertExp) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_exp, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertCosh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_cosh, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertCos) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_cos, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertAtanh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_atanh, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertAtan2) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OP_atan2, x, y)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertAtan) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_atan, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertAsinh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_asinh, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertAsin) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_asin, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertAcosh) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_acosh, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertAcos) {
  EXPECT_THROW(ConvertAndEval(AddUnary(OP_acos, x)), UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertSum) {
  EXPECT_EQ(0, ConvertAndEval(AddSum()));
  EXPECT_EQ(42, ConvertAndEval(AddSum(x), 42));
  EXPECT_EQ(123, ConvertAndEval(AddSum(x, y, z), 100, 20, 3));
}

TEST_F(GecodeTest, ConvertIntDiv) {
  NumericExpr e = AddBinary(OPintDIV, x, y);
  EXPECT_EQ(3, ConvertAndEval(e, 9, 3));
  EXPECT_EQ(2, ConvertAndEval(e, 8, 3));
  EXPECT_EQ(-2, ConvertAndEval(e, -8, 3));
  EXPECT_EQ(-2, ConvertAndEval(e, 8, -3));
  EXPECT_EQ(2, ConvertAndEval(e, -8, -3));
}

TEST_F(GecodeTest, ConvertPrecision) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OPprecision, x, y)),
      UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertRound) {
  EXPECT_EQ(42, ConvertAndEval(AddBinary(OPround, x, AddNum(0)), 42));
  EXPECT_THROW(ConvertAndEval(AddBinary(OPround, x, AddNum(1))),
      UnsupportedExprError);
  EXPECT_THROW(ConvertAndEval(AddBinary(OPround, x, y)),
        UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertTrunc) {
  EXPECT_EQ(42, ConvertAndEval(AddBinary(OPtrunc, x, AddNum(0)), 42));
  EXPECT_THROW(ConvertAndEval(AddBinary(OPtrunc, x, AddNum(1))),
      UnsupportedExprError);
  EXPECT_THROW(ConvertAndEval(AddBinary(OPtrunc, x, y)),
        UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertCount) {
  LogicalExpr a(AddRelational(NE, x, AddNum(0)));
  LogicalExpr b(AddRelational(NE, y, AddNum(0)));
  LogicalExpr c(AddRelational(NE, z, AddNum(0)));
  EXPECT_EQ(0, ConvertAndEval(AddCount(a, b, c)));
  EXPECT_EQ(1, ConvertAndEval(AddCount(a, b, c), 1));
  EXPECT_EQ(2, ConvertAndEval(AddCount(a, b, c), 0, 1, 1));
  EXPECT_EQ(3, ConvertAndEval(AddCount(a, b, c), 1, 1, 1));
}

TEST_F(GecodeTest, ConvertNumberOf) {
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

TEST_F(GecodeTest, ConvertPLTerm) {
  double args[] = {-1, 5, 0, 10, 1};
  EXPECT_THROW(ConvertAndEval(AddPLTerm(5, args, 0)),
      UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertPowConstExp) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OP1POW, x, AddNum(42))),
      UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertPow2) {
  EXPECT_EQ(49, ConvertAndEval(AddUnary(OP2POW, x), 7));
}

TEST_F(GecodeTest, ConvertPowConstBase) {
  EXPECT_THROW(ConvertAndEval(AddBinary(OPCPOW, AddNum(42), x)),
      UnsupportedExprError);
}

TEST_F(GecodeTest, ConvertNum) {
  EXPECT_EQ(42, ConvertAndEval(AddNum(42)));
  EXPECT_THROW(ConvertAndEval(AddNum(0.42)), UnsupportedExprError);
  int min = Gecode::Int::Limits::min;
  EXPECT_EQ(min, ConvertAndEval(AddNum(min)));
  EXPECT_THROW(ConvertAndEval(AddNum(min - 1)), Gecode::Int::OutOfLimits);
  int max = Gecode::Int::Limits::max;
  EXPECT_EQ(max, ConvertAndEval(AddNum(max)));
  EXPECT_THROW(ConvertAndEval(AddNum(max + 1)), Gecode::Int::OutOfLimits);
}

TEST_F(GecodeTest, ConvertVar) {
  EXPECT_EQ(11, ConvertAndEval(x, 11, 22));
  EXPECT_EQ(22, ConvertAndEval(y, 11, 22));
  EXPECT_EQ(33, ConvertAndEval(x, 33));
}

TEST_F(GecodeTest, ConvertOr) {
  NumericExpr one = AddNum(1);
  LogicalExpr e = AddBinaryLogical(
      OPOR, AddRelational(EQ, x, one), AddRelational(EQ, y, one));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 1));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1));
}

TEST_F(GecodeTest, ConvertAnd) {
  NumericExpr one = AddNum(1);
  LogicalExpr e = AddBinaryLogical(
      OPAND, AddRelational(EQ, x, one), AddRelational(EQ, y, one));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1));
}

TEST_F(GecodeTest, ConvertLess) {
  LogicalExpr e = AddRelational(LT, x, y);
  EXPECT_EQ(0, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(1, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(0, ConvertAndEval(e, 5, 3));
}

TEST_F(GecodeTest, ConvertLessEqual) {
  LogicalExpr e = AddRelational(LE, x, y);
  EXPECT_EQ(1, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(1, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(0, ConvertAndEval(e, 5, 3));
}

TEST_F(GecodeTest, ConvertEqual) {
  LogicalExpr e = AddRelational(EQ, x, y);
  EXPECT_EQ(1, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(0, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(0, ConvertAndEval(e, 5, 3));
}

TEST_F(GecodeTest, ConvertGreaterEqual) {
  LogicalExpr e = AddRelational(GE, x, y);
  EXPECT_EQ(1, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(0, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(1, ConvertAndEval(e, 5, 3));
}

TEST_F(GecodeTest, ConvertGreater) {
  LogicalExpr e = AddRelational(GT, x, y);
  EXPECT_EQ(0, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(0, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(1, ConvertAndEval(e, 5, 3));
}

TEST_F(GecodeTest, ConvertNotEqual) {
  LogicalExpr e = AddRelational(NE, x, y);
  EXPECT_EQ(0, ConvertAndEval(e, 3, 3));
  EXPECT_EQ(1, ConvertAndEval(e, 3, 5));
  EXPECT_EQ(1, ConvertAndEval(e, 5, 3));
}

TEST_F(GecodeTest, ConvertNot) {
  LogicalExpr e = AddNot(AddRelational(EQ, x, AddNum(1)));
  EXPECT_EQ(1, ConvertAndEval(e, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 1));
}

TEST_F(GecodeTest, ConvertAtLeast) {
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

TEST_F(GecodeTest, ConvertAtMost) {
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

TEST_F(GecodeTest, ConvertExactly) {
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

TEST_F(GecodeTest, ConvertNotAtLeast) {
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

TEST_F(GecodeTest, ConvertNotAtMost) {
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

TEST_F(GecodeTest, ConvertNotExactly) {
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

TEST_F(GecodeTest, ConvertForAll) {
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

TEST_F(GecodeTest, ConvertExists) {
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

TEST_F(GecodeTest, ConvertImplication) {
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

TEST_F(GecodeTest, ConvertIff) {
  LogicalExpr e = AddBinaryLogical(OP_IFF,
      AddRelational(EQ, x, AddNum(1)),
      AddRelational(EQ, y, AddNum(1)));
  EXPECT_EQ(1, ConvertAndEval(e, 0, 0));
  EXPECT_EQ(0, ConvertAndEval(e, 0, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 0));
  EXPECT_EQ(1, ConvertAndEval(e, 1, 1));
}

TEST_F(GecodeTest, ConvertAllDiff) {
  LogicalExpr e = AddAllDiff(AddNum(1), x, y);
  EXPECT_EQ(1, ConvertAndEval(e, 2, 3));
  EXPECT_EQ(0, ConvertAndEval(e, 2, 1));
  EXPECT_EQ(0, ConvertAndEval(e, 1, 1));
}

TEST_F(GecodeTest, ConvertLogicalConstant) {
  EXPECT_EQ(0, ConvertAndEval(AddBool(false)));
  EXPECT_EQ(1, ConvertAndEval(AddBool(1)));
}

// ----------------------------------------------------------------------------
// Driver tests
// TODO
/*
TEST_F(GecodeTest, Usage) {
  FILE *saved_stderr = Stderr;
  Stderr = fopen("out", "w");
  p.run((Args() + "ilogcp").get());
  fclose(Stderr);
  Stderr = saved_stderr;

  ifstream ifs("out");
  enum { BUFFER_SIZE = 4096 };
  char buffer[BUFFER_SIZE];
  string text;
  while (ifs) {
    ifs.read(buffer, BUFFER_SIZE);
    text += string(buffer, static_cast<string::size_type>(ifs.gcount()));
  }
  EXPECT_TRUE(text.find("usage: ") != string::npos);
}

TEST_F(GecodeTest, ObjConst) {
  EXPECT_EQ(0, RunDriver(DATA_DIR "objconst"));
  IloModel::Iterator iter(mod_);
  ASSERT_NE(0, iter.ok());
  IloObjective obj = (*iter).asObjective();
  EXPECT_EQ(42, obj.getConstant());
}

TEST_F(GecodeTest, CPOptimizerDoesntSupportContinuousVars) {
  EXPECT_EQ(1, RunDriver(DATA_DIR "objconst", "optimizer=cp"));
}

TEST_F(GecodeTest, SolveNumberOfCplex) {
  p.use_numberof(false);
  RunDriver(DATA_DIR "numberof", "optimizer=cplex");
}

TEST_F(GecodeTest, SolveAssign0) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign0").obj);
}

TEST_F(GecodeTest, SolveAssign1) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign1").obj);
}

TEST_F(GecodeTest, SolveBalassign0) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign0").obj);
}

TEST_F(GecodeTest, SolveBalassign1) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign1").obj);
}

TEST_F(GecodeTest, SolveFlowshp0) {
  EXPECT_NEAR(22, Solve(DATA_DIR "flowshp0").obj, 1e-5);
}

TEST_F(GecodeTest, SolveFlowshp1) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp1").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(GecodeTest, DISABLED_SolveFlowshp2) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp2").obj);
}

TEST_F(GecodeTest, SolveGrpassign0) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign0").obj);
}

// Disabled because variables in subscripts are not yet allowed.
TEST_F(GecodeTest, DISABLED_SolveGrpassign1) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1").obj);
}

// Disabled because object-valued variables are not yet allowed.
TEST_F(GecodeTest, DISABLED_SolveGrpassign1a) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1a").obj);
}

TEST_F(GecodeTest, SolveMagic) {
  EXPECT_TRUE(Solve(DATA_DIR "magic").solved);
}

TEST_F(GecodeTest, SolveMapcoloring) {
  EXPECT_TRUE(Solve(DATA_DIR "mapcoloring").solved);
}

TEST_F(GecodeTest, SolveNQueens) {
  EXPECT_TRUE(Solve(DATA_DIR "nqueens").solved);
}

TEST_F(GecodeTest, SolveNQueens0) {
  EXPECT_EQ(0, Solve(DATA_DIR "nqueens0").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(GecodeTest, DISABLED_SolveParty1) {
  EXPECT_EQ(61, Solve(DATA_DIR "party1").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(GecodeTest, DISABLED_SolveParty2) {
  EXPECT_EQ(3, Solve(DATA_DIR "party2").obj);
}

TEST_F(GecodeTest, SolveSched0) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched0").obj);
}

TEST_F(GecodeTest, SolveSched1) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched1").obj);
}

TEST_F(GecodeTest, SolveSched2) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched2").obj);
}

TEST_F(GecodeTest, SolveSendMoreMoney) {
  EXPECT_TRUE(Solve(DATA_DIR "send-more-money").solved);
}

TEST_F(GecodeTest, SolveSendMostMoney) {
  EXPECT_NEAR(10876, Solve(DATA_DIR "send-most-money",
      "relativeoptimalitytolerance=1e-5").obj, 1e-5);
}

TEST_F(GecodeTest, SolveSeq0) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0").obj, 1e-5);
}

TEST_F(GecodeTest, SolveSeq0a) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0a").obj, 1e-5);
}

TEST_F(GecodeTest, SolveSudokuHard) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuHard").solved);
}

TEST_F(GecodeTest, SolveSudokuVeryEasy) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuVeryEasy").solved);
}

// ----------------------------------------------------------------------------
// Solve code tests

TEST_F(GecodeTest, OptimalSolveCode) {
  Solve(DATA_DIR "objconst");
  EXPECT_EQ(0, p.solve_code());
}

TEST_F(GecodeTest, FeasibleSolveCode) {
  Solve(DATA_DIR "feasible");
  EXPECT_EQ(100, p.solve_code());
}

TEST_F(GecodeTest, InfeasibleSolveCode) {
  Solve(DATA_DIR "infeasible");
  EXPECT_EQ(200, p.solve_code());
}

TEST_F(GecodeTest, InfeasibleOrUnboundedSolveCode) {
  Solve(DATA_DIR "unbounded");
  EXPECT_EQ(201, p.solve_code());
}*/
}

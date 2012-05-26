#include <ilconcert/ilomodel.h>

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <cstdlib>

#include "gtest/gtest.h"

#include "solvers/concert/concert.h"
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
#include "tests/config.h"

using std::size_t;
using std::string;

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

// Returns a string representation of the argument.
template <typename T>
string str(T t) {
  std::ostringstream ss;
  ss << t;
  string s = ss.str();

  // Replace adjacent duplicate spaces and possible trailing space.
  string::iterator end = std::unique(s.begin(), s.end(), AreBothSpaces);
  if (*(end - 1) == ' ') --end;
  s.erase(end, s.end());

  // Normalize representation of infinity.
  Replace(s, "1.#INF", "inf");
  return s;
}

// Evaluates a Concert expression.
double eval(IloExpr e) {
  return e.getImpl()->eval(IloAlgorithm());
}

// A functor for deleting ASL expressions recursively.
struct ExprDeleter {
  void operator()(expr *e) const;
};

void ExprDeleter::operator()(expr *e) const {
  if (!e) return;
  switch (optype[reinterpret_cast<size_t>(e->op)]) {
  case OPTYPE_UNARY:
    (*this)(e->L.e);
    delete e;
    break;
  case OPTYPE_BINARY:
    (*this)(e->L.e);
    (*this)(e->R.e);
    delete e;
    break;
  case OPTYPE_VARARG: {
    expr_va *eva = reinterpret_cast<expr_va*>(e);
    for (de *d = reinterpret_cast<expr_va*>(e)->L.d; d->e; ++d)
      (*this)(d->e);
    delete eva;
    break;
  }
  case OPTYPE_PLTERM:
    std::free(e->L.p);
    (*this)(e->R.e);
    delete e;
    break;
  case OPTYPE_IF: {
    expr_if *eif = reinterpret_cast<expr_if*>(e);
    (*this)(eif->e);
    (*this)(eif->T);
    (*this)(eif->F);
    delete eif;
    break;
  }
  case OPTYPE_SUM: case OPTYPE_COUNT:
    for (expr **i = e->L.ep, **end = e->R.ep; i < end; ++i)
      (*this)(*i);
    delete e;
    break;
  case OPTYPE_NUMBER:
    delete reinterpret_cast<expr_n*>(e);
    break;
  default:
    delete e;
    break;
  }
}

#if HAVE_UNIQUE_PTR
typedef std::unique_ptr<expr, ExprDeleter> ExprPtr;
#else
// Fall back to std::auto_ptr - leaks subexpressions.
typedef std::auto_ptr<expr> ExprPtr;
ExprPtr move(ExprPtr e) { return e; }
#endif

class ConcertTest : public ::testing::Test {
 protected:
  void SetUp() {
    env = IloEnv();
    mod = IloModel(env);
    Var = IloNumVarArray(env, 3);
    Var[0] = IloNumVar(env, 0, 1, "x");
    Var[1] = IloNumVar(env, 0, 1, "y");
    Var[2] = IloNumVar(env, 0, 1, "theta");
    usenumberof = 0;
  }

  void TearDown() {
    Var = IloNumVarArray();
  }

  // Creates an ASL expression representing a number.
  static ExprPtr NewNum(double n) {
    expr_n e = {reinterpret_cast<efunc_n*>(OPNUM), n};
    return ExprPtr(reinterpret_cast<expr*>(new expr_n(e)));
  }

  // Creates an ASL expression representing a variable.
  static ExprPtr NewVar(int var_index) {
    expr e = {reinterpret_cast<efunc*>(OPVARVAL), var_index, 0, {0}, {0}, 0};
    return ExprPtr(new expr(e));
  }

  // Creates an unary ASL expression.
  static ExprPtr NewUnary(int opcode, ExprPtr arg) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {arg.release()}, {0}, 0};
    return ExprPtr(new expr(e));
  }

  // Creates a binary ASL expression.
  static ExprPtr NewBinary(int opcode, ExprPtr lhs, ExprPtr rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.release()}, {rhs.release()}, 0};
    return ExprPtr(new expr(e));
  }

  static de MakeDE(ExprPtr e) {
    de result = {e.release(), 0, {0}};
    return result;
  }

  // Creates a variable-argument ASL expression with up to 3 arguments.
  static ExprPtr NewVarArg(int opcode, ExprPtr e1, ExprPtr e2,
      ExprPtr e3 = ExprPtr());

  static ExprPtr NewPLTerm(int size, const double *args, int var_index);

  // Creates an ASL expression representing if-then-else.
  static ExprPtr NewIf(int opcode, ExprPtr condition,
      ExprPtr true_expr, ExprPtr false_expr);

  // Creates an ASL expression representing a sum with up to 3 arguments.
  static ExprPtr NewSum(int opcode, ExprPtr arg1, ExprPtr arg2,
      ExprPtr arg3 = ExprPtr());

  static double EvalRem(double lhs, double rhs) {
    return eval(build_expr(NewBinary(OPREM, NewNum(lhs), NewNum(rhs)).get()));
  }
};

ExprPtr ConcertTest::NewVarArg(int opcode, ExprPtr e1, ExprPtr e2, ExprPtr e3) {
  expr_va e = {reinterpret_cast<efunc*>(opcode), 0,
               {0}, {0}, 0, 0, 0};
  expr_va *copy = new expr_va(e);
  ExprPtr result(reinterpret_cast<expr*>(copy));
  de *args = new de[4];
  args[0] = MakeDE(move(e1));
  args[1] = MakeDE(move(e2));
  args[2] = MakeDE(move(e3));
  args[3] = MakeDE(ExprPtr());
  copy->L.d = args;
  return result;
}

ExprPtr ConcertTest::NewPLTerm(int size, const double *args, int var_index) {
  expr e = {reinterpret_cast<efunc*>(OPPLTERM), 0, 0, {0}, {0}, 0};
  ExprPtr pl(new expr(e));
  pl->L.p = static_cast<plterm*>(
      std::calloc(1, sizeof(plterm) + sizeof(real) * (size - 1)));
  pl->L.p->n = 3;
  real *bs = pl->L.p->bs;
  for (int i = 0; i < size; i++)
    bs[i] = args[i];
  pl->R.e = NewVar(var_index).release();
  return pl;
}

ExprPtr ConcertTest::NewIf(int opcode, ExprPtr condition,
    ExprPtr true_expr, ExprPtr false_expr) {
  expr_if e = {reinterpret_cast<efunc*>(opcode), 0, condition.release(),
               true_expr.release(), false_expr.release(),
               0, 0, 0, 0, {0}, {0}, 0, 0};
  return ExprPtr(reinterpret_cast<expr*>(new expr_if(e)));
}

ExprPtr ConcertTest::NewSum(int opcode,
    ExprPtr arg1, ExprPtr arg2, ExprPtr arg3) {
  expr e = {reinterpret_cast<efunc*>(opcode), 0, 0, {0}, {0}, 0};
  ExprPtr sum(new expr(e));
  expr** args = sum->L.ep = new expr*[3];
  sum->R.ep = args + (arg3 ? 3 : 2);
  args[0] = arg1.release();
  args[1] = arg2.release();
  args[2] = arg3.release();
  return sum;
}

TEST_F(ConcertTest, ConvertNum) {
  EXPECT_EQ("0.42", str(build_expr(NewNum(0.42).get())));
}

TEST_F(ConcertTest, ConvertVar) {
  EXPECT_EQ("theta", str(build_expr(NewVar(2).get())));
}

TEST_F(ConcertTest, ConvertPlus) {
  EXPECT_EQ("x + 42", str(build_expr(
    NewBinary(OPPLUS, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("x + y", str(build_expr(
    NewBinary(OPPLUS, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertMinus) {
  EXPECT_EQ("x + -42", str(build_expr(
    NewBinary(OPMINUS, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("x + -1 * y", str(build_expr(
    NewBinary(OPMINUS, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertMult) {
  EXPECT_EQ("42 * x", str(build_expr(
    NewBinary(OPMULT, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("x * y", str(build_expr(
    NewBinary(OPMULT, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertDiv) {
  EXPECT_EQ("x / 42", str(build_expr(
    NewBinary(OPDIV, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("x / y", str(build_expr(
    NewBinary(OPDIV, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertRem) {
  EXPECT_EQ("x + trunc(x / y ) * y * -1", str(build_expr(
    NewBinary(OPREM, NewVar(0), NewVar(1)).get())));
  EXPECT_EQ(0, EvalRem(9, 3));
  EXPECT_EQ(2, EvalRem(8, 3));
  EXPECT_EQ(-2, EvalRem(-8, 3));
  EXPECT_EQ(2, EvalRem(8, -3));
  EXPECT_EQ(-2, EvalRem(-8, -3));
  EXPECT_EQ(1.5, EvalRem(7.5, 3));
}

TEST_F(ConcertTest, ConvertPow) {
  EXPECT_EQ("x ^ 42", str(build_expr(
    NewBinary(OPPOW, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("x ^ y", str(build_expr(
    NewBinary(OPPOW, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertLess) {
  EXPECT_EQ("max(x + -42 , 0)", str(build_expr(
    NewBinary(OPLESS, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("max(x + -1 * y , 0)", str(build_expr(
    NewBinary(OPLESS, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertMin) {
  EXPECT_EQ("min( [x , y , 42 ])", str(build_expr(
    NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertMax) {
  EXPECT_EQ("max([x , y , 42 ])", str(build_expr(
    NewVarArg(MAXLIST, NewVar(0), NewVar(1), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertFloor) {
  EXPECT_EQ("floor(x )", str(build_expr(NewUnary(FLOOR, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertCeil) {
  EXPECT_EQ("ceil(x )", str(build_expr(NewUnary(CEIL, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAbs) {
  EXPECT_EQ("abs(x )", str(build_expr(NewUnary(ABS, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertUMinus) {
  EXPECT_EQ("-1 * x", str(build_expr(NewUnary(OPUMINUS, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertUnsupportedOpThrows) {
  int ops[] = {
      OPOR, OPAND, LT, LE, EQ, GE, GT, NE, OPNOT,
      OPprecision, OPFUNCALL, OPIFSYM, OPHOL, OPNUMBEROFs,
      OPATLEAST, OPATMOST, OPEXACTLY, OPALLDIFF,
      OPNOTATLEAST, OPNOTATMOST, OPNOTEXACTLY,
      ANDLIST, ORLIST, OPIMPELSE, OP_IFF, N_OPS, -1, 500
  };
  size_t i = 0;
  for (size_t num_ops = sizeof(ops) / sizeof(*ops); i < num_ops; ++i) {
    EXPECT_THROW(build_expr(
      NewBinary(ops[i], NewVar(0), NewVar(1)).get()), UnsupportedExprError);
  }
  // Paranoid: make sure that the loop body has been executed enough times.
  EXPECT_EQ(28u, i);
}

TEST_F(ConcertTest, ConvertIf) {
  EXPECT_EQ("IloNumVar(7)[-inf..inf]", str(build_expr(NewIf(OPIFnl,
      NewBinary(EQ, NewVar(0), NewNum(0)), NewVar(1), NewNum(42)).get())));

  IloModel::Iterator iter(mod);
  ASSERT_TRUE(iter.ok());
  IloIfThenI *ifTrue = dynamic_cast<IloIfThenI*>((*iter).getImpl());
  ASSERT_TRUE(ifTrue != nullptr);
  EXPECT_EQ("x == 0", str(ifTrue->getLeft()));
  EXPECT_EQ("IloNumVar(7)[-inf..inf] == y", str(ifTrue->getRight()));

  ++iter;
  ASSERT_TRUE(iter.ok());
  IloIfThenI *ifFalse = dynamic_cast<IloIfThenI*>((*iter).getImpl());
  ASSERT_TRUE(ifFalse != nullptr);
  IloNotI *ifNot = dynamic_cast<IloNotI*>(ifFalse->getLeft().getImpl());
  EXPECT_EQ("x == 0", str(ifNot->getConstraint()));
  EXPECT_EQ("IloNumVar(7)[-inf..inf] == 42", str(ifFalse->getRight()));

  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(ConcertTest, ConvertTanh) {
  EXPECT_EQ("exp(2 * x ) + -1 / exp(2 * x ) + 1",
            str(build_expr(NewUnary(OP_tanh, NewVar(0)).get())));
  // Concert incorrectly omits brackets around the dividend and divisor
  // above, so test also by evaluating the expression at several points.
  IloExpr e(build_expr(NewUnary(OP_tanh, NewNum(1)).get()));
  EXPECT_NEAR(0.761594, eval(e), 1e-5);
  e = build_expr(NewUnary(OP_tanh, NewNum(0)).get());
  EXPECT_EQ(0, eval(e));
  e = build_expr(NewUnary(OP_tanh, NewNum(-2)).get());
  EXPECT_NEAR(-0.964027, eval(e), 1e-5);
}

TEST_F(ConcertTest, ConvertTan) {
  EXPECT_EQ("tan(x )", str(build_expr(NewUnary(OP_tan, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertSqrt) {
  EXPECT_EQ("x ^ 0.5",
            str(build_expr(NewUnary(OP_sqrt, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertSinh) {
  EXPECT_EQ("exp(x ) * 0.5 + exp(-1 * x ) * -0.5",
            str(build_expr(NewUnary(OP_sinh, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertSin) {
  EXPECT_EQ("sin(x )", str(build_expr(NewUnary(OP_sin, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertLog10) {
  EXPECT_EQ("log(x )/ 2.30259",
            str(build_expr(NewUnary(OP_log10, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertLog) {
  EXPECT_EQ("log(x )", str(build_expr(NewUnary(OP_log, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertExp) {
  EXPECT_EQ("exp(x )", str(build_expr(NewUnary(OP_exp, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertCosh) {
  EXPECT_EQ("exp(x ) * 0.5 + exp(-1 * x ) * 0.5",
            str(build_expr(NewUnary(OP_cosh, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertCos) {
  EXPECT_EQ("cos(x )", str(build_expr(NewUnary(OP_cos, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAtanh) {
  EXPECT_EQ("log(x + 1 ) * 0.5 + log(-1 * x + 1 ) * -0.5",
            str(build_expr(NewUnary(OP_atanh, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAtan2) {
  EXPECT_EQ("IloNumVar(8)[-inf..inf]",
            str(build_expr(NewBinary(OP_atan2, NewVar(1), NewVar(0)).get())));

  IloModel::Iterator iter(mod);
  ASSERT_TRUE(iter.ok());
  IloIfThenI *ifXNonnegative = dynamic_cast<IloIfThenI*>((*iter).getImpl());
  ASSERT_TRUE(ifXNonnegative != nullptr);
  EXPECT_EQ("0 <= x", str(ifXNonnegative->getLeft()));
  EXPECT_EQ("IloNumVar(8)[-inf..inf] == arc-tan(y / x )",  // (1)
            str(ifXNonnegative->getRight()));

  ++iter;
  ASSERT_TRUE(iter.ok());
  IloIfThenI *ifDiffSigns = dynamic_cast<IloIfThenI*>((*iter).getImpl());
  ASSERT_TRUE(ifDiffSigns != nullptr);
  EXPECT_EQ("(x <= 0 ) && (0 <= y )", str(ifDiffSigns->getLeft()));
  EXPECT_EQ("IloNumVar(8)[-inf..inf] == arc-tan(y / x ) + 3.14159",  // (2)
            str(ifDiffSigns->getRight()));

  ++iter;
  ASSERT_TRUE(iter.ok());
  IloIfThenI *ifSameSigns = dynamic_cast<IloIfThenI*>((*iter).getImpl());
  ASSERT_TRUE(ifSameSigns != nullptr);
  EXPECT_EQ("(x <= 0 ) && (y <= 0 )", str(ifSameSigns->getLeft()));
  EXPECT_EQ("IloNumVar(8)[-inf..inf] == arc-tan(y / x ) + -3.14159",
            str(ifSameSigns->getRight()));

  ++iter;
  EXPECT_FALSE(iter.ok());

  // Check that (1) and (2) both yield NaN when x == 0 and y == 0.
  double zero = 0, d = IloArcTan(zero / zero);
  EXPECT_TRUE(d != d);
  double d1 = d + 3.14;
  EXPECT_TRUE(d1 != d1);
}

TEST_F(ConcertTest, ConvertAtan) {
  EXPECT_EQ("arc-tan(x )",
            str(build_expr(NewUnary(OP_atan, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAsinh) {
  EXPECT_EQ("log(x + square(x ) + 1 ^ 0.5)",
            str(build_expr(NewUnary(OP_asinh, NewVar(0)).get())));
  // Concert incorrectly omits brackets around square(x) + 1
  // above, so test also by evaluating the expression at several points.
  IloExpr e(build_expr(NewUnary(OP_asinh, NewNum(1)).get()));
  EXPECT_NEAR(0.881373, eval(e), 1e-5);
  e = build_expr(NewUnary(OP_asinh, NewNum(0)).get());
  EXPECT_EQ(0, eval(e));
  e = build_expr(NewUnary(OP_asinh, NewNum(-2)).get());
  EXPECT_NEAR(-1.443635, eval(e), 1e-5);
}

TEST_F(ConcertTest, ConvertAsin) {
  EXPECT_EQ("arc-sin(x )",
            str(build_expr(NewUnary(OP_asin, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAcosh) {
  EXPECT_EQ("log(x + x + 1 ^ 0.5 * x + -1 ^ 0.5)",
            str(build_expr(NewUnary(OP_acosh, NewVar(0)).get())));
  // Concert incorrectly omits brackets around x + 1 and x + -1
  // above, so test also by evaluating the expression at several points.
  IloExpr e(build_expr(NewUnary(OP_acosh, NewNum(1)).get()));
  EXPECT_NEAR(0, eval(e), 1e-5);
  e = build_expr(NewUnary(OP_acosh, NewNum(10)).get());
  EXPECT_NEAR(2.993222, eval(e), 1e-5);
  e = build_expr(NewUnary(OP_acosh, NewNum(0)).get());
  double n = eval(e);
  EXPECT_TRUE(n != n);
}

TEST_F(ConcertTest, ConvertAcos) {
  EXPECT_EQ("arc-cos(x )",
            str(build_expr(NewUnary(OP_acos, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertSum) {
  EXPECT_EQ("x + y + 42", str(build_expr(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertIntDiv) {
  EXPECT_EQ("trunc(x / y )", str(build_expr(
    NewBinary(OPintDIV, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertRound) {
  EXPECT_EQ("round(x )", str(build_expr(
    NewBinary(OPround, NewVar(0), NewNum(0)).get())));

  EXPECT_EQ(1235, eval(build_expr(
    NewBinary(OPround, NewNum(1234.56), NewNum(0)).get())));
  EXPECT_EQ(3, eval(build_expr(
    NewBinary(OPround, NewNum(2.5), NewNum(0)).get())));
  EXPECT_EQ(-2, eval(build_expr(
    NewBinary(OPround, NewNum(-2.5), NewNum(0)).get())));

  EXPECT_THROW(build_expr(
    NewBinary(OPround, NewVar(0), NewVar(1)).get()), UnsupportedExprError);
}

TEST_F(ConcertTest, ConvertTrunc) {
  EXPECT_EQ("trunc(x )", str(build_expr(
    NewBinary(OPtrunc, NewVar(0), NewNum(0)).get())));
  EXPECT_EQ(1234, eval(build_expr(
    NewBinary(OPtrunc, NewNum(1234.56), NewNum(0)).get())));
  EXPECT_THROW(build_expr(
    NewBinary(OPtrunc, NewVar(0), NewVar(1)).get()), UnsupportedExprError);
}

TEST_F(ConcertTest, Convert1Pow) {
  EXPECT_EQ("x ^ 42", str(build_expr(
    NewBinary(OP1POW, NewVar(0), NewNum(42)).get())));
}

TEST_F(ConcertTest, Convert2Pow) {
  EXPECT_EQ("square(x )", str(build_expr(
    NewUnary(OP2POW, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertCPow) {
  EXPECT_EQ("42 ^ x", str(build_expr(
    NewBinary(OPCPOW, NewNum(42), NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertPLTerm) {
  double args[] = {-1, 5, 0, 10, 1};
  EXPECT_EQ("piecewiselinear(x[0..1] , [5, 10], [-1, 0, 1], 0, 0)",
            str(build_expr(NewPLTerm(5, args, 0).get())));
}

TEST_F(ConcertTest, ConvertCount) {
  ExprPtr a(NewBinary(EQ, NewVar(0), NewNum(0)));
  ExprPtr b(NewBinary(LE, NewVar(1), NewNum(42)));
  ExprPtr c(NewBinary(GE, NewVar(2), NewNum(0)));
  EXPECT_EQ("x == 0 + y <= 42 + 0 <= theta", str(build_expr(
      NewSum(OPCOUNT, move(a), move(b), move(c)).get())));
}

TEST_F(ConcertTest, ConvertNumberOf) {
  usenumberof = 1;
  EXPECT_EQ("x == theta + y == theta", str(build_expr(
      NewSum(OPNUMBEROF, NewVar(2), NewVar(0), NewVar(1)).get())));
  usenumberof = 0;
  EXPECT_EQ("x == 42 + y == 42", str(build_expr(
      NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)).get())));
  // TODO: test if usenumberof is 1
}

TEST_F(ConcertTest, ConvertVarSubVar) {
  EXPECT_EQ("num-exprs[IloIntVar(4)[0..2] ]", str(build_expr(
      NewUnary(OPVARSUBVAR, NewVar(0)).get())));

  IloModel::Iterator iter(mod);
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(4)[0..2] == x", str(*iter));

  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(ConcertTest, SameNum) {
  EXPECT_TRUE(same_expr(NewNum(0.42).get(), NewNum(0.42).get()));
  EXPECT_FALSE(same_expr(NewNum(0.42).get(), NewNum(42).get()));
}

TEST_F(ConcertTest, SameVar) {
  EXPECT_TRUE(same_expr(NewVar(0).get(), NewVar(0).get()));
  EXPECT_FALSE(same_expr(NewVar(0).get(), NewVar(1).get()));
  EXPECT_FALSE(same_expr(NewVar(0).get(), NewNum(0).get()));
}

TEST_F(ConcertTest, SameUnary) {
  EXPECT_TRUE(same_expr(NewUnary(OPUMINUS, NewVar(0)).get(),
                        NewUnary(OPUMINUS, NewVar(0)).get()));
  EXPECT_FALSE(same_expr(NewUnary(OPUMINUS, NewVar(0)).get(),
                         NewVar(0).get()));
  EXPECT_FALSE(same_expr(NewUnary(OPUMINUS, NewVar(0)).get(),
                         NewUnary(FLOOR, NewVar(0)).get()));
  EXPECT_FALSE(same_expr(NewUnary(OPUMINUS, NewVar(0)).get(),
                         NewUnary(OPUMINUS, NewVar(1)).get()));
}

TEST_F(ConcertTest, SameBinary) {
  EXPECT_TRUE(same_expr(NewBinary(OPPLUS, NewVar(0), NewNum(42)).get(),
                        NewBinary(OPPLUS, NewVar(0), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(NewBinary(OPPLUS, NewVar(0), NewNum(42)).get(),
                         NewBinary(OPMINUS, NewVar(0), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(NewBinary(OPPLUS, NewVar(0), NewNum(42)).get(),
                         NewBinary(OPPLUS, NewNum(42), NewVar(0)).get()));
  EXPECT_FALSE(same_expr(NewBinary(OPPLUS, NewVar(0), NewNum(42)).get(),
                         NewBinary(OPPLUS, NewVar(0), NewNum(0)).get()));
  EXPECT_FALSE(same_expr(NewNum(42).get(),
                         NewBinary(OPPLUS, NewVar(0), NewNum(42)).get()));
}

TEST_F(ConcertTest, SameVarArg) {
  EXPECT_TRUE(same_expr(
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewVarArg(MINLIST, NewVar(0), NewVar(1)).get()));
  EXPECT_FALSE(same_expr(
      NewVarArg(MINLIST, NewVar(0), NewVar(1)).get(),
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewVarArg(MAXLIST, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(0)).get()));
  EXPECT_FALSE(same_expr(
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewNum(42).get()));
}

TEST_F(ConcertTest, SamePLTerm) {
  double args[] = {-1, 5, 0, 10, 1};
  EXPECT_TRUE(same_expr(
      NewPLTerm(5, args, 0).get(),
      NewPLTerm(5, args, 0).get()));
  EXPECT_FALSE(same_expr(
      NewPLTerm(5, args, 0).get(),
      NewPLTerm(3, args, 0).get()));
  EXPECT_FALSE(same_expr(
      NewPLTerm(5, args, 0).get(),
      NewPLTerm(5, args, 1).get()));
  double args2[] = {-1, 5, 0, 11, 1};
  EXPECT_FALSE(same_expr(
      NewPLTerm(5, args, 0).get(),
      NewPLTerm(5, args2, 0).get()));
  EXPECT_FALSE(same_expr(
      NewPLTerm(5, args, 0).get(),
      NewNum(42).get()));
}

TEST_F(ConcertTest, SameIf) {
  EXPECT_TRUE(same_expr(
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewIf(OPIFSYM, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(0)).get()));
  EXPECT_FALSE(same_expr(
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewNum(42).get()));
}

TEST_F(ConcertTest, SameSum) {
  EXPECT_TRUE(same_expr(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1)).get()));
  EXPECT_FALSE(same_expr(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1)).get(),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(0)).get()));
  EXPECT_FALSE(same_expr(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewNum(42).get()));
}

TEST_F(ConcertTest, SameCount) {
  EXPECT_TRUE(same_expr(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewSum(OPCOUNT, NewVar(0), NewVar(1)).get()));
  EXPECT_FALSE(same_expr(
      NewSum(OPCOUNT, NewVar(0), NewVar(1)).get(),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)).get()));
  EXPECT_FALSE(same_expr(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(0)).get()));
  EXPECT_FALSE(same_expr(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)).get(),
      NewNum(42).get()));
}

TEST_F(ConcertTest, ConvertSameExprThrowsOnUnsupportedOp) {
  EXPECT_THROW(same_expr(
      NewUnary(OPFUNCALL, ExprPtr()).get(),
      NewUnary(OPFUNCALL, ExprPtr()).get()),
      UnsupportedExprError);
  EXPECT_THROW(same_expr(
      NewUnary(OPHOL, ExprPtr()).get(),
      NewUnary(OPHOL, ExprPtr()).get()),
      UnsupportedExprError);
}

TEST_F(ConcertTest, ConvertSameExprThrowsOnUnknownOp) {
  EXPECT_THROW(same_expr(
      NewUnary(7, ExprPtr()).get(),
      NewUnary(7, ExprPtr()).get()),
      Error);
}
}

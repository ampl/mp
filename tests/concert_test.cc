// Concert driver tests.

#include <ilconcert/ilodiffi.h>
#include <ilconcert/ilopathi.h>
#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

#include "gtest/gtest.h"

#include "solvers/concert/concert.h"
#include "solvers/concert/util.h"
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
#include "tests/config.h"

using std::ifstream;
using std::size_t;
using std::string;
using std::vector;

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

struct SolveResult {
  bool solved;
  double obj;
  SolveResult(bool solved, double obj) : solved(solved), obj(obj) {}
};

class ConcertTest : public ::testing::Test {
 protected:
  Driver d;
  IloEnv env;
  IloModel mod;

  void SetUp() {
    env = d.env();
    mod = d.mod();
    IloNumVarArray vars = IloNumVarArray(env, 3);
    vars[0] = IloIntVar(env, 0, 1, "x");
    vars[1] = IloNumVar(env, 0, 1, "y");
    vars[2] = IloNumVar(env, 0, 1, "theta");
    d.set_vars(vars);
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

  double EvalRem(double lhs, double rhs) {
    return eval(d.build_expr(NewBinary(OPREM, NewNum(lhs), NewNum(rhs)).get()));
  }

  int RunDriver(const char *stub, const char *opt);

  SolveResult Solve(const char *stub);
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

int ConcertTest::RunDriver(const char *stub = nullptr,
                           const char *opt = nullptr) {
  // Copy arguments to comply with the Driver::run function signature and
  // avoid unwanted modification.
  const char *args[] = {"concert", "-s", stub, opt};
  vector<char> store;
  size_t num_args = sizeof(args) / sizeof(*args);
  if (!stub) num_args -= 2;
  else if (!opt) --num_args;
  for (size_t i = 0; i < num_args; ++i) {
    const char *arg = args[i];
    store.insert(store.end(), arg, arg + strlen(arg) + 1);
  }
  vector<char*> argv(num_args + 1);
  for (size_t i = 0, j = 0; i < num_args; j += strlen(args[i]) + 1, ++i)
    argv[i] = &store[j];
  try {
    return d.run(num_args, &argv[0]);
  } catch (const IloException& e) {
    throw std::runtime_error(e.getMessage());
  }
  return 0;
}

SolveResult ConcertTest::Solve(const char *stub) {
  RunDriver(stub);
  ifstream ifs((string(stub) + ".sol").c_str());
  string line;
  getline(ifs, line);
  getline(ifs, line);
  bool solved =line.find("solution found") != string::npos;
  getline(ifs, line);
  const char obj[] = "objective ";
  size_t pos = line.find(obj);
  double zero = 0;
  return SolveResult(solved, pos != string::npos ?
      atof(line.c_str() + pos + sizeof(obj) - 1) : zero / zero);
}


TEST_F(ConcertTest, ConvertNum) {
  EXPECT_EQ("0.42", str(d.build_expr(NewNum(0.42).get())));
}

TEST_F(ConcertTest, ConvertVar) {
  EXPECT_EQ("theta", str(d.build_expr(NewVar(2).get())));
}

TEST_F(ConcertTest, ConvertPlus) {
  EXPECT_EQ("x + 42", str(d.build_expr(
    NewBinary(OPPLUS, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("x + y", str(d.build_expr(
    NewBinary(OPPLUS, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertMinus) {
  EXPECT_EQ("x + -42", str(d.build_expr(
    NewBinary(OPMINUS, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("x + -1 * y", str(d.build_expr(
    NewBinary(OPMINUS, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertMult) {
  EXPECT_EQ("42 * x", str(d.build_expr(
    NewBinary(OPMULT, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("x * y", str(d.build_expr(
    NewBinary(OPMULT, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertDiv) {
  EXPECT_EQ("x / 42", str(d.build_expr(
    NewBinary(OPDIV, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("x / y", str(d.build_expr(
    NewBinary(OPDIV, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertRem) {
  EXPECT_EQ("x + trunc(x / y ) * y * -1", str(d.build_expr(
    NewBinary(OPREM, NewVar(0), NewVar(1)).get())));
  EXPECT_EQ(0, EvalRem(9, 3));
  EXPECT_EQ(2, EvalRem(8, 3));
  EXPECT_EQ(-2, EvalRem(-8, 3));
  EXPECT_EQ(2, EvalRem(8, -3));
  EXPECT_EQ(-2, EvalRem(-8, -3));
  EXPECT_EQ(1.5, EvalRem(7.5, 3));
}

TEST_F(ConcertTest, ConvertPow) {
  EXPECT_EQ("x ^ 42", str(d.build_expr(
    NewBinary(OPPOW, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("x ^ y", str(d.build_expr(
    NewBinary(OPPOW, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertLess) {
  EXPECT_EQ("max(x + -42 , 0)", str(d.build_expr(
    NewBinary(OPLESS, NewVar(0), NewNum(42)).get())));
  EXPECT_EQ("max(x + -1 * y , 0)", str(d.build_expr(
    NewBinary(OPLESS, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertMin) {
  EXPECT_EQ("min( [x , y , 42 ])", str(d.build_expr(
    NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertMax) {
  EXPECT_EQ("max([x , y , 42 ])", str(d.build_expr(
    NewVarArg(MAXLIST, NewVar(0), NewVar(1), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertFloor) {
  EXPECT_EQ("floor(x )", str(d.build_expr(NewUnary(FLOOR, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertCeil) {
  EXPECT_EQ("ceil(x )", str(d.build_expr(NewUnary(CEIL, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAbs) {
  EXPECT_EQ("abs(x )", str(d.build_expr(NewUnary(ABS, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertUMinus) {
  EXPECT_EQ("-1 * x", str(d.build_expr(NewUnary(OPUMINUS, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertUnsupportedExprThrows) {
  int ops[] = {
      OPOR, OPAND, LT, LE, EQ, GE, GT, NE, OPNOT,
      OPprecision, OPFUNCALL, OPIFSYM, OPHOL,
      OPATLEAST, OPATMOST, OPEXACTLY,
      OPNOTATLEAST, OPNOTATMOST, OPNOTEXACTLY,
      OPIMPELSE, OP_IFF, N_OPS, -1, 500
  };
  size_t i = 0;
  for (size_t num_ops = sizeof(ops) / sizeof(*ops); i < num_ops; ++i) {
    EXPECT_THROW(d.build_expr(
      NewBinary(ops[i], NewVar(0), NewVar(1)).get()), UnsupportedExprError);
  }
  // Paranoid: make sure that the loop body has been executed enough times.
  EXPECT_EQ(24u, i);

  EXPECT_THROW(d.build_expr(
    NewSum(OPNUMBEROFs, NewVar(0), NewVar(1)).get()), UnsupportedExprError);
  EXPECT_THROW(d.build_expr(
    NewSum(OPALLDIFF, NewVar(0), NewVar(1)).get()), UnsupportedExprError);
  EXPECT_THROW(d.build_expr(
    NewSum(ANDLIST, NewVar(0), NewVar(1)).get()), UnsupportedExprError);
  EXPECT_THROW(d.build_expr(
    NewSum(ORLIST, NewVar(0), NewVar(1)).get()), UnsupportedExprError);
}

TEST_F(ConcertTest, ConvertIf) {
  EXPECT_EQ("IloNumVar(7)[-inf..inf]", str(d.build_expr(NewIf(OPIFnl,
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
            str(d.build_expr(NewUnary(OP_tanh, NewVar(0)).get())));
  // Concert incorrectly omits brackets around the dividend and divisor
  // above, so test also by evaluating the expression at several points.
  IloExpr e(d.build_expr(NewUnary(OP_tanh, NewNum(1)).get()));
  EXPECT_NEAR(0.761594, eval(e), 1e-5);
  e = d.build_expr(NewUnary(OP_tanh, NewNum(0)).get());
  EXPECT_EQ(0, eval(e));
  e = d.build_expr(NewUnary(OP_tanh, NewNum(-2)).get());
  EXPECT_NEAR(-0.964027, eval(e), 1e-5);
}

TEST_F(ConcertTest, ConvertTan) {
  EXPECT_EQ("tan(x )", str(d.build_expr(NewUnary(OP_tan, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertSqrt) {
  EXPECT_EQ("x ^ 0.5",
            str(d.build_expr(NewUnary(OP_sqrt, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertSinh) {
  EXPECT_EQ("exp(x ) * 0.5 + exp(-1 * x ) * -0.5",
            str(d.build_expr(NewUnary(OP_sinh, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertSin) {
  EXPECT_EQ("sin(x )", str(d.build_expr(NewUnary(OP_sin, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertLog10) {
  EXPECT_EQ("log(x )/ 2.30259",
            str(d.build_expr(NewUnary(OP_log10, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertLog) {
  EXPECT_EQ("log(x )", str(d.build_expr(NewUnary(OP_log, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertExp) {
  EXPECT_EQ("exp(x )", str(d.build_expr(NewUnary(OP_exp, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertCosh) {
  EXPECT_EQ("exp(x ) * 0.5 + exp(-1 * x ) * 0.5",
            str(d.build_expr(NewUnary(OP_cosh, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertCos) {
  EXPECT_EQ("cos(x )", str(d.build_expr(NewUnary(OP_cos, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAtanh) {
  EXPECT_EQ("log(x + 1 ) * 0.5 + log(-1 * x + 1 ) * -0.5",
            str(d.build_expr(NewUnary(OP_atanh, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAtan2) {
  EXPECT_EQ("IloNumVar(8)[-inf..inf]",
            str(d.build_expr(NewBinary(OP_atan2, NewVar(1), NewVar(0)).get())));

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
            str(d.build_expr(NewUnary(OP_atan, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAsinh) {
  EXPECT_EQ("log(x + square(x ) + 1 ^ 0.5)",
            str(d.build_expr(NewUnary(OP_asinh, NewVar(0)).get())));
  // Concert incorrectly omits brackets around square(x) + 1
  // above, so test also by evaluating the expression at several points.
  IloExpr e(d.build_expr(NewUnary(OP_asinh, NewNum(1)).get()));
  EXPECT_NEAR(0.881373, eval(e), 1e-5);
  e = d.build_expr(NewUnary(OP_asinh, NewNum(0)).get());
  EXPECT_EQ(0, eval(e));
  e = d.build_expr(NewUnary(OP_asinh, NewNum(-2)).get());
  EXPECT_NEAR(-1.443635, eval(e), 1e-5);
}

TEST_F(ConcertTest, ConvertAsin) {
  EXPECT_EQ("arc-sin(x )",
            str(d.build_expr(NewUnary(OP_asin, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAcosh) {
  EXPECT_EQ("log(x + x + 1 ^ 0.5 * x + -1 ^ 0.5)",
            str(d.build_expr(NewUnary(OP_acosh, NewVar(0)).get())));
  // Concert incorrectly omits brackets around x + 1 and x + -1
  // above, so test also by evaluating the expression at several points.
  IloExpr e(d.build_expr(NewUnary(OP_acosh, NewNum(1)).get()));
  EXPECT_NEAR(0, eval(e), 1e-5);
  e = d.build_expr(NewUnary(OP_acosh, NewNum(10)).get());
  EXPECT_NEAR(2.993222, eval(e), 1e-5);
  e = d.build_expr(NewUnary(OP_acosh, NewNum(0)).get());
  double n = eval(e);
  EXPECT_TRUE(n != n);
}

TEST_F(ConcertTest, ConvertAcos) {
  EXPECT_EQ("arc-cos(x )",
            str(d.build_expr(NewUnary(OP_acos, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertSum) {
  EXPECT_EQ("x + y + 42", str(d.build_expr(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertIntDiv) {
  EXPECT_EQ("trunc(x / y )", str(d.build_expr(
    NewBinary(OPintDIV, NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, ConvertRound) {
  EXPECT_EQ("round(x )", str(d.build_expr(
    NewBinary(OPround, NewVar(0), NewNum(0)).get())));

  EXPECT_EQ(1235, eval(d.build_expr(
    NewBinary(OPround, NewNum(1234.56), NewNum(0)).get())));
  EXPECT_EQ(3, eval(d.build_expr(
    NewBinary(OPround, NewNum(2.5), NewNum(0)).get())));
  EXPECT_EQ(-2, eval(d.build_expr(
    NewBinary(OPround, NewNum(-2.5), NewNum(0)).get())));

  EXPECT_THROW(d.build_expr(
    NewBinary(OPround, NewVar(0), NewVar(1)).get()), UnsupportedExprError);
}

TEST_F(ConcertTest, ConvertTrunc) {
  EXPECT_EQ("trunc(x )", str(d.build_expr(
    NewBinary(OPtrunc, NewVar(0), NewNum(0)).get())));
  EXPECT_EQ(1234, eval(d.build_expr(
    NewBinary(OPtrunc, NewNum(1234.56), NewNum(0)).get())));
  EXPECT_THROW(d.build_expr(
    NewBinary(OPtrunc, NewVar(0), NewVar(1)).get()), UnsupportedExprError);
}

TEST_F(ConcertTest, Convert1Pow) {
  EXPECT_EQ("x ^ 42", str(d.build_expr(
    NewBinary(OP1POW, NewVar(0), NewNum(42)).get())));
}

TEST_F(ConcertTest, Convert2Pow) {
  EXPECT_EQ("square(x )", str(d.build_expr(
    NewUnary(OP2POW, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertCPow) {
  EXPECT_EQ("42 ^ x", str(d.build_expr(
    NewBinary(OPCPOW, NewNum(42), NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertPLTerm) {
  double args[] = {-1, 5, 0, 10, 1};
  EXPECT_EQ("piecewiselinear(x[0..1] , [5, 10], [-1, 0, 1], 0, 0)",
            str(d.build_expr(NewPLTerm(5, args, 0).get())));
}

TEST_F(ConcertTest, ConvertCount) {
  ExprPtr a(NewBinary(EQ, NewVar(0), NewNum(0)));
  ExprPtr b(NewBinary(LE, NewVar(1), NewNum(42)));
  ExprPtr c(NewBinary(GE, NewVar(2), NewNum(0)));
  EXPECT_EQ("x == 0 + y <= 42 + 0 <= theta", str(d.build_expr(
      NewSum(OPCOUNT, move(a), move(b), move(c)).get())));
}

TEST_F(ConcertTest, ConvertNumberOf) {
  d.use_numberof();
  EXPECT_EQ("x == theta + y == theta", str(d.build_expr(
      NewSum(OPNUMBEROF, NewVar(2), NewVar(0), NewVar(1)).get())));
  d.use_numberof(false);
  EXPECT_EQ("x == 42 + y == 42", str(d.build_expr(
      NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)).get())));
}

TEST_F(ConcertTest, IloArrayCopyingIsCheap) {
  IloIntArray array(env);
  array.add(42);
  EXPECT_TRUE(array.getImpl() != nullptr);
  EXPECT_EQ(array.getImpl(), IloIntArray(array).getImpl());
}

TEST_F(ConcertTest, ConvertSingleNumberOfToIloDistribute) {
  d.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  EXPECT_EQ("IloIntVar(10)" + bounds, str(d.build_expr(
      NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)).get())));
  d.finish_building_numberof();
  IloModel::Iterator iter(mod);
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(4)" + bounds + " == x", str(*iter));
  ++iter;
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(7)" + bounds + " == y", str(*iter));
  ++iter;
  ASSERT_TRUE(iter.ok());
  IloDistributeI *dist = dynamic_cast<IloDistributeI*>((*iter).getImpl());
  ASSERT_TRUE(dist != nullptr);
  EXPECT_EQ("[IloIntVar(10)" + bounds + " ]", str(dist->getCardVarArray()));
  EXPECT_EQ("[IloIntVar(4)" + bounds + " , IloIntVar(7)" + bounds + " ]",
      str(dist->getVarArray()));
  EXPECT_EQ("[42]", str(dist->getValueArray()));
  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(ConcertTest, ConvertTwoNumberOfsWithSameValuesToIloDistribute) {
  d.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  ExprPtr expr(NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)));
  EXPECT_EQ("IloIntVar(10)" + bounds, str(d.build_expr(expr.get())));
  EXPECT_EQ("IloIntVar(10)" + bounds, str(d.build_expr(
      NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)).get())));
  d.finish_building_numberof();
  IloModel::Iterator iter(mod);
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(4)" + bounds + " == x", str(*iter));
  ++iter;
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(7)" + bounds + " == y", str(*iter));
  ++iter;
  ASSERT_TRUE(iter.ok());
  IloDistributeI *dist = dynamic_cast<IloDistributeI*>((*iter).getImpl());
  ASSERT_TRUE(dist != nullptr);
  EXPECT_EQ("[IloIntVar(10)" + bounds + " ]", str(dist->getCardVarArray()));
  EXPECT_EQ("[IloIntVar(4)" + bounds + " , IloIntVar(7)" + bounds + " ]",
      str(dist->getVarArray()));
  EXPECT_EQ("[42]", str(dist->getValueArray()));
  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(ConcertTest, ConvertTwoNumberOfsWithDiffValuesToIloDistribute) {
  d.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  ExprPtr expr(NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)));
  EXPECT_EQ("IloIntVar(10)" + bounds, str(d.build_expr(expr.get())));
  EXPECT_EQ("IloIntVar(12)" + bounds, str(d.build_expr(
      NewSum(OPNUMBEROF, NewNum(43), NewVar(0), NewVar(1)).get())));
  d.finish_building_numberof();
  IloModel::Iterator iter(mod);
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(4)" + bounds + " == x", str(*iter));
  ++iter;
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(7)" + bounds + " == y", str(*iter));
  ++iter;
  ASSERT_TRUE(iter.ok());
  IloDistributeI *dist = dynamic_cast<IloDistributeI*>((*iter).getImpl());
  ASSERT_TRUE(dist != nullptr);
  EXPECT_EQ("[IloIntVar(10)" + bounds + " , IloIntVar(12)" + bounds + " ]",
      str(dist->getCardVarArray()));
  EXPECT_EQ("[IloIntVar(4)" + bounds + " , IloIntVar(7)" + bounds + " ]",
      str(dist->getVarArray()));
  EXPECT_EQ("[42, 43]", str(dist->getValueArray()));
  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(ConcertTest, ConvertTwoNumberOfsWithDiffExprs) {
  d.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  ExprPtr expr(NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)));
  EXPECT_EQ("IloIntVar(10)" + bounds, str(d.build_expr(expr.get())));
  EXPECT_EQ("IloIntVar(15)" + bounds, str(d.build_expr(
      NewSum(OPNUMBEROF, NewNum(42), NewVar(2)).get())));
  d.finish_building_numberof();
  IloModel::Iterator iter(mod);
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(4)" + bounds + " == x", str(*iter));
  ++iter;
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(7)" + bounds + " == y", str(*iter));
  ++iter;
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(12)" + bounds + " == theta", str(*iter));
  ++iter;
  ASSERT_TRUE(iter.ok());
  IloDistributeI *dist = dynamic_cast<IloDistributeI*>((*iter).getImpl());
  ASSERT_TRUE(dist != nullptr);
  EXPECT_EQ("[IloIntVar(10)" + bounds + " ]", str(dist->getCardVarArray()));
  EXPECT_EQ("[IloIntVar(4)" + bounds + " , IloIntVar(7)" + bounds + " ]",
      str(dist->getVarArray()));
  EXPECT_EQ("[42]", str(dist->getValueArray()));
  ++iter;
  ASSERT_TRUE(iter.ok());
  dist = dynamic_cast<IloDistributeI*>((*iter).getImpl());
  ASSERT_TRUE(dist != nullptr);
  EXPECT_EQ("[IloIntVar(15)" + bounds + " ]", str(dist->getCardVarArray()));
  EXPECT_EQ("[IloIntVar(12)" + bounds + " ]",
      str(dist->getVarArray()));
  EXPECT_EQ("[42]", str(dist->getValueArray()));
  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(ConcertTest, ConvertVarSubVar) {
  EXPECT_EQ("num-exprs[IloIntVar(4)[0..2] ]", str(d.build_expr(
      NewUnary(OPVARSUBVAR, NewVar(0)).get())));

  IloModel::Iterator iter(mod);
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(4)[0..2] == x", str(*iter));

  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(ConcertTest, ConvertIncompleteConstraintExprThrows) {
  int ops[] = {
      OPPLUS, OPMINUS, OPMULT, OPDIV, OPREM,
      OPPOW, OPLESS, FLOOR, CEIL, ABS, OPUMINUS,
      OP_tanh, OP_tan, OP_sqrt, OP_sinh, OP_sin, OP_log10,
      OP_log, OP_exp, OP_cosh, OP_cos, OP_atanh, OP_atan2,
      OP_atan, OP_asinh, OP_asin, OP_acosh, OP_acos,
      OPintDIV, OPprecision, OPround, OPtrunc,
      OP1POW, OP2POW, OPCPOW, OPFUNCALL, OPHOL, OPVARVAL,
      N_OPS, -1, 500
  };
  size_t i = 0;
  for (size_t num_ops = sizeof(ops) / sizeof(*ops); i < num_ops; ++i) {
    EXPECT_THROW(d.build_constr(
      NewBinary(ops[i], NewVar(0), NewVar(1)).get()),
      IncompleteConstraintExprError);
  }
  // Paranoid: make sure that the loop body has been executed enough times.
  EXPECT_EQ(41u, i);

  EXPECT_THROW(d.build_constr(
    NewVarArg(MINLIST, NewVar(0), NewVar(1)).get()),
    IncompleteConstraintExprError);
  EXPECT_THROW(d.build_constr(
    NewVarArg(MAXLIST, NewVar(0), NewVar(1)).get()),
    IncompleteConstraintExprError);
  EXPECT_THROW(d.build_constr(
    NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(0)).get()),
    IncompleteConstraintExprError);
  EXPECT_THROW(d.build_constr(
    NewIf(OPIFSYM, NewVar(0), NewVar(1), NewNum(0)).get()),
    IncompleteConstraintExprError);
  EXPECT_THROW(d.build_constr(
    NewSum(OPSUMLIST, NewVar(0), NewVar(1)).get()),
    IncompleteConstraintExprError);
  EXPECT_THROW(d.build_constr(
    NewSum(OPCOUNT, NewVar(0), NewVar(1)).get()),
    IncompleteConstraintExprError);
  EXPECT_THROW(d.build_constr(
    NewSum(OPNUMBEROF, NewVar(0), NewVar(1)).get()),
    IncompleteConstraintExprError);
  EXPECT_THROW(d.build_constr(
    NewSum(OPNUMBEROFs, NewVar(0), NewVar(1)).get()),
    IncompleteConstraintExprError);
  EXPECT_THROW(d.build_constr(
    NewPLTerm(0, 0, 0).get()), IncompleteConstraintExprError);
}

TEST_F(ConcertTest, ConvertFalse) {
  EXPECT_EQ("IloNumVar(4)[1..1] == 0", str(d.build_constr(NewNum(0).get())));
}

TEST_F(ConcertTest, ConvertTrue) {
  EXPECT_EQ("IloNumVar(4)[1..1] == 1", str(d.build_constr(NewNum(1).get())));
}

TEST_F(ConcertTest, ThrowsOnNumInConstraint) {
  EXPECT_THROW(d.build_constr(NewNum(42).get()), Error);
}

TEST_F(ConcertTest, ConvertLT) {
  IloConstraint c(d.build_constr(NewBinary(LT, NewVar(0), NewNum(42)).get()));
  IloNotI *n = dynamic_cast<IloNotI*>(c.getImpl());
  ASSERT_TRUE(n != nullptr);
  EXPECT_EQ("42 <= x", str(n->getConstraint()));
}

TEST_F(ConcertTest, ConvertLE) {
  EXPECT_EQ("x <= 42", str(d.build_constr(
      NewBinary(LE, NewVar(0), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertEQ) {
  EXPECT_EQ("x == 42", str(d.build_constr(
      NewBinary(EQ, NewVar(0), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertGE) {
  EXPECT_EQ("42 <= x", str(d.build_constr(
      NewBinary(GE, NewVar(0), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertGT) {
  IloConstraint c(d.build_constr(NewBinary(GT, NewVar(0), NewNum(42)).get()));
  IloNotI *n = dynamic_cast<IloNotI*>(c.getImpl());
  ASSERT_TRUE(n != nullptr);
  EXPECT_EQ("x <= 42", str(n->getConstraint()));
}

TEST_F(ConcertTest, ConvertNE) {
  EXPECT_EQ("x != 42", str(d.build_constr(
      NewBinary(NE, NewVar(0), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertAtMost) {
  EXPECT_EQ("42 <= x", str(d.build_constr(
      NewBinary(OPATMOST, NewVar(0), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertNotAtMost) {
  IloConstraint c(d.build_constr(
      NewBinary(OPNOTATMOST, NewVar(0), NewNum(42)).get()));
  IloNotI *n = dynamic_cast<IloNotI*>(c.getImpl());
  ASSERT_TRUE(n != nullptr);
  EXPECT_EQ("42 <= x", str(n->getConstraint()));
}

TEST_F(ConcertTest, ConvertAtLeast) {
  EXPECT_EQ("x <= 42", str(d.build_constr(
      NewBinary(OPATLEAST, NewVar(0), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertNotAtLeast) {
  IloConstraint c(d.build_constr(
      NewBinary(OPNOTATLEAST, NewVar(0), NewNum(42)).get()));
  IloNotI *n = dynamic_cast<IloNotI*>(c.getImpl());
  ASSERT_TRUE(n != nullptr);
  EXPECT_EQ("x <= 42", str(n->getConstraint()));
}

TEST_F(ConcertTest, ConvertExactly) {
  EXPECT_EQ("x == 42", str(d.build_constr(
      NewBinary(OPEXACTLY, NewVar(0), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertNotExactly) {
  EXPECT_EQ("x != 42", str(d.build_constr(
      NewBinary(OPNOTEXACTLY, NewVar(0), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertOr) {
  IloConstraint c(d.build_constr(NewBinary(OPOR,
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2))).get()));
  IloIfThenI *ifThen = dynamic_cast<IloIfThenI*>(c.getImpl());
  ASSERT_TRUE(ifThen != nullptr);
  IloNotI *n = dynamic_cast<IloNotI*>(ifThen->getLeft().getImpl());
  ASSERT_TRUE(n != nullptr);
  EXPECT_EQ("x == 1", str(n->getConstraint()));
  EXPECT_EQ("x == 2", str(ifThen->getRight()));
}

TEST_F(ConcertTest, CheckOrTruthTable) {
  IloNumVarArray vars = d.vars();
  vars[0].setBounds(0, 0);
  vars[1].setBounds(0, 0);
  mod.add(d.build_constr(NewBinary(OPOR,
        NewBinary(EQ, NewVar(0), NewNum(1)),
        NewBinary(EQ, NewVar(1), NewNum(1))).get()));
  IloCP cp(mod);
  EXPECT_FALSE(cp.solve());
  vars[0].setBounds(0, 0);
  vars[1].setBounds(1, 1);
  EXPECT_TRUE(cp.solve());
  vars[0].setBounds(1, 1);
  vars[1].setBounds(0, 0);
  EXPECT_TRUE(cp.solve());
  vars[0].setBounds(1, 1);
  vars[1].setBounds(1, 1);
  EXPECT_TRUE(cp.solve());
}

TEST_F(ConcertTest, ConvertExists) {
  EXPECT_EQ("(x == 1 ) || (x == 2 ) || (x == 3 )",
      str(d.build_constr(NewSum(ORLIST,
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2)),
      NewBinary(EQ, NewVar(0), NewNum(3))).get())));
}

TEST_F(ConcertTest, ConvertAnd) {
  EXPECT_EQ("(x == 1 ) && (x == 2 )", str(d.build_constr(NewBinary(OPAND,
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2))).get())));
}

TEST_F(ConcertTest, ConvertForAll) {
  EXPECT_EQ("(x == 1 ) && (x == 2 ) && (x == 3 )",
      str(d.build_constr(NewSum(ANDLIST,
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2)),
      NewBinary(EQ, NewVar(0), NewNum(3))).get())));
}

TEST_F(ConcertTest, ConvertNot) {
  IloConstraint c(d.build_constr(
      NewUnary(OPNOT, NewBinary(LE, NewVar(0), NewNum(42))).get()));
  IloNotI *n = dynamic_cast<IloNotI*>(c.getImpl());
  ASSERT_TRUE(n != nullptr);
  EXPECT_EQ("x <= 42", str(n->getConstraint()));
}

TEST_F(ConcertTest, ConvertIff) {
  EXPECT_EQ("x == 1 == x == 2", str(d.build_constr(NewBinary(OP_IFF,
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2))).get())));
}

TEST_F(ConcertTest, ConvertImpElse) {
  IloConstraint con(d.build_constr(NewIf(OPIMPELSE,
      NewBinary(EQ, NewVar(0), NewNum(0)),
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2))).get()));
  IloAndI* conjunction = dynamic_cast<IloAndI*>(con.getImpl());
  ASSERT_TRUE(conjunction != nullptr);

  IloAndI::Iterator iter(conjunction);
  ASSERT_TRUE(iter.ok());
  IloIfThenI *ifTrue = dynamic_cast<IloIfThenI*>(*iter);
  ASSERT_TRUE(ifTrue != nullptr);
  EXPECT_EQ("x == 0", str(ifTrue->getLeft()));
  EXPECT_EQ("x == 1", str(ifTrue->getRight()));

  ++iter;
  ASSERT_TRUE(iter.ok());
  IloIfThenI *ifFalse = dynamic_cast<IloIfThenI*>(*iter);
  ASSERT_TRUE(ifFalse != nullptr);
  IloNotI *ifNot = dynamic_cast<IloNotI*>(ifFalse->getLeft().getImpl());
  EXPECT_EQ("x == 0", str(ifNot->getConstraint()));
  EXPECT_EQ("x == 2", str(ifFalse->getRight()));

  ++iter;
  EXPECT_FALSE(iter.ok());
}

TEST_F(ConcertTest, ConvertAllDiff) {
  IloConstraint con(d.build_constr(
      NewSum(OPALLDIFF, NewVar(0), NewNum(42)).get()));
  IloAllDiffI* diff = dynamic_cast<IloAllDiffI*>(con.getImpl());
  ASSERT_TRUE(diff != nullptr);
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  EXPECT_EQ("[x[0..1] , IloIntVar(4)" + bounds + " ]",
      str(diff->getExprArray()));

  IloModel::Iterator iter(mod);
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(4)" + bounds +" == 42", str(*iter));
  ++iter;
  EXPECT_FALSE(iter.ok());
}

// ----------------------------------------------------------------------------
// Util tests

TEST_F(ConcertTest, GetOpName) {
  EXPECT_STREQ("+", get_opname(OPPLUS));
  EXPECT_STREQ("-", get_opname(OPMINUS));
  EXPECT_STREQ("*", get_opname(OPMULT));
  EXPECT_STREQ("/", get_opname(OPDIV));
  EXPECT_STREQ("mod", get_opname(OPREM));
  EXPECT_STREQ("^", get_opname(OPPOW));
  EXPECT_STREQ("less", get_opname(OPLESS));
  EXPECT_STREQ("min", get_opname(MINLIST));
  EXPECT_STREQ("max", get_opname(MAXLIST));
  EXPECT_STREQ("floor", get_opname(FLOOR));
  EXPECT_STREQ("ceil", get_opname(CEIL));
  EXPECT_STREQ("abs", get_opname(ABS));
  EXPECT_STREQ("unary -", get_opname(OPUMINUS));
  EXPECT_STREQ("||", get_opname(OPOR));
  EXPECT_STREQ("&&", get_opname(OPAND));
  EXPECT_STREQ("<", get_opname(LT));
  EXPECT_STREQ("<=", get_opname(LE));
  EXPECT_STREQ("=", get_opname(EQ));
  EXPECT_STREQ(">=", get_opname(GE));
  EXPECT_STREQ(">", get_opname(GT));
  EXPECT_STREQ("!=", get_opname(NE));
  EXPECT_STREQ("!", get_opname(OPNOT));
  EXPECT_STREQ("if-then-else", get_opname(OPIFnl));
  EXPECT_STREQ("tanh", get_opname(OP_tanh));
  EXPECT_STREQ("tan", get_opname(OP_tan));
  EXPECT_STREQ("sqrt", get_opname(OP_sqrt));
  EXPECT_STREQ("sinh", get_opname(OP_sinh));
  EXPECT_STREQ("sin", get_opname(OP_sin));
  EXPECT_STREQ("log10", get_opname(OP_log10));
  EXPECT_STREQ("log", get_opname(OP_log));
  EXPECT_STREQ("exp", get_opname(OP_exp));
  EXPECT_STREQ("cosh", get_opname(OP_cosh));
  EXPECT_STREQ("cos", get_opname(OP_cos));
  EXPECT_STREQ("atanh", get_opname(OP_atanh));
  EXPECT_STREQ("atan2", get_opname(OP_atan2));
  EXPECT_STREQ("atan", get_opname(OP_atan));
  EXPECT_STREQ("asinh", get_opname(OP_asinh));
  EXPECT_STREQ("asin", get_opname(OP_asin));
  EXPECT_STREQ("acosh", get_opname(OP_acosh));
  EXPECT_STREQ("acos", get_opname(OP_acos));
  EXPECT_STREQ("sum", get_opname(OPSUMLIST));
  EXPECT_STREQ("div", get_opname(OPintDIV));
  EXPECT_STREQ("precision", get_opname(OPprecision));
  EXPECT_STREQ("round", get_opname(OPround));
  EXPECT_STREQ("trunc", get_opname(OPtrunc));
  EXPECT_STREQ("count", get_opname(OPCOUNT));
  EXPECT_STREQ("numberof", get_opname(OPNUMBEROF));
  EXPECT_STREQ("string numberof", get_opname(OPNUMBEROFs));
  EXPECT_STREQ("atleast", get_opname(OPATLEAST));
  EXPECT_STREQ("atmost", get_opname(OPATMOST));
  EXPECT_STREQ("pl term", get_opname(OPPLTERM));
  EXPECT_STREQ("string if-then-else", get_opname(OPIFSYM));
  EXPECT_STREQ("exactly", get_opname(OPEXACTLY));
  EXPECT_STREQ("not atleast", get_opname(OPNOTATLEAST));
  EXPECT_STREQ("not atmost", get_opname(OPNOTATMOST));
  EXPECT_STREQ("not exactly", get_opname(OPNOTEXACTLY));
  EXPECT_STREQ("forall", get_opname(ANDLIST));
  EXPECT_STREQ("exists", get_opname(ORLIST));
  EXPECT_STREQ("implies else", get_opname(OPIMPELSE));
  EXPECT_STREQ("iff", get_opname(OP_IFF));
  EXPECT_STREQ("alldiff", get_opname(OPALLDIFF));
  EXPECT_STREQ("1pow", get_opname(OP1POW));
  EXPECT_STREQ("^2", get_opname(OP2POW));
  EXPECT_STREQ("cpow", get_opname(OPCPOW));
  EXPECT_STREQ("function call", get_opname(OPFUNCALL));
  EXPECT_STREQ("number", get_opname(OPNUM));
  EXPECT_STREQ("string", get_opname(OPHOL));
  EXPECT_STREQ("variable", get_opname(OPVARVAL));
  EXPECT_STREQ("unknown", get_opname(N_OPS));
  EXPECT_STREQ("unknown", get_opname(-1));
  EXPECT_STREQ("unknown", get_opname(500));
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

TEST_F(ConcertTest, SameExprThrowsOnUnsupportedOp) {
  EXPECT_THROW(same_expr(
      NewUnary(OPFUNCALL, ExprPtr()).get(),
      NewUnary(OPFUNCALL, ExprPtr()).get()),
      UnsupportedExprError);
  EXPECT_THROW(same_expr(
      NewUnary(OPHOL, ExprPtr()).get(),
      NewUnary(OPHOL, ExprPtr()).get()),
      UnsupportedExprError);
}

// ----------------------------------------------------------------------------
// Driver tests

TEST_F(ConcertTest, Usage) {
  system("concert/concert 2> out");
  ifstream ifs("out");
  size_t BUFFER_SIZE = 4096;
  char buffer[BUFFER_SIZE];
  string text;
  while (ifs) {
    ifs.read(buffer, BUFFER_SIZE);
    text += string(buffer, ifs.gcount());
  }
  EXPECT_TRUE(text.find("usage: ") != string::npos);
}

TEST_F(ConcertTest, ObjConst) {
  RunDriver("data/objconst");
  IloModel::Iterator iter(mod);
  ASSERT_TRUE(iter.ok());
  IloObjective obj = (*iter).asObjective();
  EXPECT_EQ(42, obj.getConstant());
}

TEST_F(ConcertTest, SolveNumberOfCplex) {
  d.use_numberof(false);
  RunDriver("data/numberof", "ilogcplex");
}

TEST_F(ConcertTest, SolveAssign0) {
  EXPECT_EQ(61, Solve("data/assign0").obj);
}

// Disabled because variables in subscripts are not yet allowed.
TEST_F(ConcertTest, DISABLED_SolveAssign1) {
  EXPECT_EQ(61, Solve("data/assign1").obj);
}

// Disabled because of a syntax error in the model.
TEST_F(ConcertTest, DISABLED_SolveAssign1a) {
  EXPECT_EQ(61, Solve("data/assign1a").obj);
}

TEST_F(ConcertTest, SolveBalassign0) {
  EXPECT_EQ(14, Solve("data/balassign0").obj);
}

TEST_F(ConcertTest, SolveBalassign1) {
  EXPECT_EQ(14, Solve("data/balassign1").obj);
}

TEST_F(ConcertTest, SolveMagic) {
  EXPECT_TRUE(Solve("data/magic").solved);
}

TEST_F(ConcertTest, SolveNQueens) {
  EXPECT_TRUE(Solve("data/nqueens").solved);
}

TEST_F(ConcertTest, SolveNQueens0) {
  EXPECT_EQ(0, Solve("data/nqueens0").obj);
}

// Disabled because of an .nl input problem.
TEST_F(ConcertTest, DISABLED_SolveParty1) {
  EXPECT_EQ(61, Solve("data/party1").obj);
}

TEST_F(ConcertTest, SolveParty2) {
  EXPECT_EQ(3, Solve("data/party2").obj);
}

// ----------------------------------------------------------------------------
// Option tests

TEST_F(ConcertTest, DebugExpr0) {
  RunDriver("data/magic", "debugexpr=0");
  EXPECT_EQ(0, d.get_option(Driver::DEBUGEXPR));
}

TEST_F(ConcertTest, DebugExpr1) {
  RunDriver("data/magic", "debugexpr=1");
  EXPECT_EQ(1, d.get_option(Driver::DEBUGEXPR));
}

TEST_F(ConcertTest, DebugExpr42) {
  RunDriver("data/magic", "debugexpr=42");
  EXPECT_EQ(42, d.get_option(Driver::DEBUGEXPR));
}

TEST_F(ConcertTest, IlogSolver) {
  RunDriver("data/magic", "ilogsolver");
  EXPECT_EQ(0, d.get_option(Driver::ILOGOPTTYPE));
  EXPECT_TRUE(dynamic_cast<IloCplexI*>(d.alg().getImpl()) == nullptr);
}

TEST_F(ConcertTest, IlogCplex) {
  RunDriver("data/objconst", "ilogcplex");
  EXPECT_EQ(1, d.get_option(Driver::ILOGOPTTYPE));
  EXPECT_TRUE(dynamic_cast<IloCplexI*>(d.alg().getImpl()) != nullptr);
}
}

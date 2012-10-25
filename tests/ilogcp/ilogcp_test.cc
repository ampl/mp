// Ilogcp driver tests.

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

#include "solvers/ilogcp/ilogcp.h"
#include "solvers/util/expr.h"

extern "C" {
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
}

#include "tests/config.h"

using std::ifstream;
using std::size_t;
using std::string;
using std::vector;

using namespace ampl;

#define DATA_DIR "../data/"

namespace ampl {
class ExprBuilder {
 public:
  static expr *GetExprPtr(Expr e) {
    return e.expr_;
  }
};
}

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

struct SolveResult {
  bool solved;
  double obj;
  SolveResult(bool solved, double obj) : solved(solved), obj(obj) {}
};

// Helper class that copies arguments to comply with the IlogCPDriver::run
// function signature and avoid unwanted modification.
class Args {
 private:
  int argc_;
  vector<char> store_;
  vector<char*> argv_;

 public:
  Args() : argc_(0) {}

  char **get() {
    argv_.resize(argc_ + 1);
    for (int i = 0, j = 0; i < argc_; j += strlen(&store_[j]) + 1, ++i)
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

struct EnumValue {
  const char *name;
  IloCP::ParameterValues value;
};

expr *ptr(Expr e) {
  return ampl::ExprBuilder::GetExprPtr(e);
}

class IlogCPTest : public ::testing::Test {
 protected:
  IlogCPDriver d;
  IloEnv env_;
  IloModel mod_;
  std::vector<expr*> exprs_;
  
  void SetUp() {
    env_ = d.env();
    mod_ = d.mod();
    IloNumVarArray vars = IloNumVarArray(env_, 3);
    vars[0] = IloIntVar(env_, 0, 1, "x");
    vars[1] = IloNumVar(env_, 0, 1, "y");
    vars[2] = IloNumVar(env_, 0, 1, "theta");
    d.set_vars(vars);
  }
  
  void TearDown();

  NumericExpr NewExpr(expr *e) {
    exprs_.push_back(e);
    return NumericExpr(e);
  }
  
  // Creates an ASL expression representing a number.
  NumericExpr NewNum(double n) {
    expr_n e = {reinterpret_cast<efunc_n*>(OPNUM), n};
    return NewExpr(reinterpret_cast<expr*>(new expr_n(e)));
  }

  // Creates an ASL expression representing a variable.
  NumericExpr NewVar(int var_index) {
    expr e = {reinterpret_cast<efunc*>(OPVARVAL), var_index, 0, {0}, {0}, 0};
    return NewExpr(new expr(e));
  }

  // Creates an unary ASL expression.
  NumericExpr NewUnary(int opcode, NumericExpr arg) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0, {ptr(arg)}, {0}, 0};
    return NewExpr(new expr(e));
  }

  // Creates a binary ASL expression.
  NumericExpr NewBinary(int opcode, NumericExpr lhs, NumericExpr rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {ptr(lhs)}, {ptr(rhs)}, 0};
    return NewExpr(new expr(e));
  }

  static de MakeDE(NumericExpr e) {
    de result = {ptr(e), 0, {0}};
    return result;
  }

  // Creates a variable-argument ASL expression with up to 3 arguments.
  NumericExpr NewVarArg(int opcode, NumericExpr e1,
      NumericExpr e2, NumericExpr e3 = NumericExpr());

  NumericExpr NewPLTerm(int size, const double *args, int var_index);

  // Creates an ASL expression representing if-then-else.
  NumericExpr NewIf(int opcode, NumericExpr condition,
      NumericExpr true_expr, NumericExpr false_expr);

  // Creates an ASL expression representing a sum with up to 3 arguments.
  NumericExpr NewSum(int opcode, NumericExpr arg1,
      NumericExpr arg2, NumericExpr arg3 = NumericExpr());

  IloConstraint ConvertLogical(NumericExpr e) {
    return d.Visit(LogicalExpr(ptr(e)));
  }

  double EvalRem(double lhs, double rhs) {
    return eval(d.Visit(NewBinary(OPREM, NewNum(lhs), NewNum(rhs))));
  }

  int RunDriver(const char *stub = nullptr, const char *opt = nullptr) {
    try {
      return d.run((Args() + "ilogcp" + "-s" + stub + opt).get());
    } catch (const IloException &e) {
      throw std::runtime_error(e.getMessage());
    }
    return 0;
  }

  bool ParseOptions(const char *opt1, const char *opt2 = nullptr) {
    try {
      return d.parse_options((Args() + opt1 + opt2).get());
    } catch (const IloException &e) {
      throw std::runtime_error(e.getMessage());
    }
    return false;
  }

  SolveResult Solve(const char *stub);

  template <typename T>
  static string Option(const char *name, T value) {
    std::ostringstream os;
    os << name << "=" << value;
    return os.str();
  }

  void CheckIntCPOption(const char *option, IloCP::IntParam param,
      int start, int end, int offset = 0, bool accepts_auto = false,
      const EnumValue *values = 0);

  template <typename ParamT>
  void CheckIntCPLEXOption(const char *option,
      ParamT param, int start, int end) {
    for (int i = std::max(start, -9); i <= std::min(end, 9); ++i) {
      EXPECT_TRUE(ParseOptions("optimizer=cplex", Option(option, i).c_str()));
      CPLEXOptimizer *opt = dynamic_cast<CPLEXOptimizer*>(d.optimizer());
      ASSERT_TRUE(opt != nullptr);
      EXPECT_EQ(i, opt->cplex().getParam(param))
        << "Failed option: " << option;
    }
    if (end != INT_MAX) {
      EXPECT_FALSE(ParseOptions("optimizer=cplex",
          Option(option, end + 1).c_str()));
    }
    if (start != INT_MIN) {
      int small = start - 1;
      EXPECT_FALSE(ParseOptions("optimizer=cplex",
          Option(option, small).c_str()));
      EXPECT_FALSE(ParseOptions("optimizer=cp", Option(option, start).c_str()));
    }
  }

  void CheckDblCPOption(const char *option, IloCP::NumParam param,
      double good, double bad);
};
  
void IlogCPTest::TearDown() {
  for (vector<expr*>::const_iterator
       i = exprs_.begin(), end = exprs_.end(); i != end; ++i) {
    expr *e = *i;
    if (NumericExpr(e).opcode() >= N_OPS)
      continue;
    switch (NumericExpr(e).type()) {
      case OPTYPE_VARARG:
        delete reinterpret_cast<expr_va*>(e);
        break;
      case OPTYPE_PLTERM:
        std::free(e->L.p);
        delete e;
        break;
      case OPTYPE_IF:
        delete reinterpret_cast<expr_if*>(e);
        break;
      case OPTYPE_NUMBER:
        delete reinterpret_cast<expr_n*>(e);
        break;
      default:
        delete e;
        break;
    }
  }
  exprs_.clear();
}

NumericExpr IlogCPTest::NewVarArg(int opcode,
    NumericExpr e1, NumericExpr e2, NumericExpr e3) {
  expr_va e = {reinterpret_cast<efunc*>(opcode), 0, {0}, {0}, 0, 0, 0};
  expr_va *copy = new expr_va(e);
  expr *result(reinterpret_cast<expr*>(copy));
  de *args = new de[4];
  args[0] = MakeDE(e1);
  args[1] = MakeDE(e2);
  args[2] = MakeDE(e3);
  args[3] = MakeDE(NumericExpr());
  copy->L.d = args;
  return NewExpr(result);
}

NumericExpr IlogCPTest::NewPLTerm(int size, const double *args, int var_index) {
  expr e = {reinterpret_cast<efunc*>(OPPLTERM), 0, 0, {0}, {0}, 0};
  NumericExpr pl(NewExpr(new expr(e)));
  ptr(pl)->L.p = static_cast<plterm*>(
      std::calloc(1, sizeof(plterm) + sizeof(real) * (size - 1)));
  ptr(pl)->L.p->n = (size + 1) / 2;
  real *bs = ptr(pl)->L.p->bs;
  for (int i = 0; i < size; i++)
    bs[i] = args[i];
  ptr(pl)->R.e = ptr(NewVar(var_index));
  return pl;
}

NumericExpr IlogCPTest::NewIf(int opcode, NumericExpr condition,
    NumericExpr true_expr, NumericExpr false_expr) {
  expr_if e = {reinterpret_cast<efunc*>(opcode), 0, ptr(condition),
               ptr(true_expr), ptr(false_expr),
               0, 0, 0, 0, {0}, {0}, 0, 0};
  return NewExpr(reinterpret_cast<expr*>(new expr_if(e)));
}

NumericExpr IlogCPTest::NewSum(int opcode,
    NumericExpr arg1, NumericExpr arg2, NumericExpr arg3) {
  expr e = {reinterpret_cast<efunc*>(opcode), 0, 0, {0}, {0}, 0};
  NumericExpr sum(NewExpr(new expr(e)));
  expr** args = ptr(sum)->L.ep = new expr*[3];
  ptr(sum)->R.ep = args + (ptr(arg3) ? 3 : 2);
  args[0] = ptr(arg1);
  args[1] = ptr(arg2);
  args[2] = ptr(arg3);
  return sum;
}

SolveResult IlogCPTest::Solve(const char *stub) {
  RunDriver(stub);
  ifstream ifs((string(stub) + ".sol").c_str());
  string line;
  getline(ifs, line);
  bool solved = line.find("optimal solution") != string::npos;
  if (!solved) solved = line.find("feasible solution") != string::npos;
  getline(ifs, line);
  const char obj[] = "objective ";
  size_t pos = line.find(obj);
  double zero = 0;
  return SolveResult(solved, pos != string::npos ?
      atof(line.c_str() + pos + sizeof(obj) - 1) : zero / zero);
}

void IlogCPTest::CheckIntCPOption(const char *option,
    IloCP::IntParam param, int start, int end, int offset, bool accepts_auto,
    const EnumValue *values) {
  for (int i = start; i <= std::min(end, 9); ++i) {
    EXPECT_TRUE(ParseOptions("optimizer=cp", Option(option, i).c_str()));
    CPOptimizer *opt = dynamic_cast<CPOptimizer*>(d.optimizer());
    ASSERT_TRUE(opt != nullptr);
    EXPECT_EQ(offset + i, opt->solver().getParameter(param))
      << "Failed option: " << option;
  }
  if (end != INT_MAX)
    EXPECT_FALSE(ParseOptions("optimizer=cp", Option(option, end + 1).c_str()));
  if (accepts_auto) {
    EXPECT_TRUE(ParseOptions("optimizer=cp", Option(option, -1).c_str()));
    CPOptimizer *opt = dynamic_cast<CPOptimizer*>(d.optimizer());
    ASSERT_TRUE(opt != nullptr);
    EXPECT_EQ(IloCP::Auto, opt->solver().getParameter(param));

    EXPECT_TRUE(ParseOptions("optimizer=cp", Option(option, "auto").c_str()));
    opt = dynamic_cast<CPOptimizer*>(d.optimizer());
    ASSERT_TRUE(opt != nullptr);
    EXPECT_EQ(IloCP::Auto, opt->solver().getParameter(param));
  }
  int small = start - 1;
  if (accepts_auto && small == -1)
    --small;
  EXPECT_FALSE(ParseOptions("optimizer=cp", Option(option, small).c_str()));
  EXPECT_FALSE(ParseOptions("optimizer=cplex", Option(option, start).c_str()));
  if (values) {
    int count = 0;
    for (const EnumValue *v = values; v->name; ++v, ++count) {
      EXPECT_TRUE(ParseOptions("optimizer=cp", Option(option, v->name).c_str()));
      CPOptimizer *opt = dynamic_cast<CPOptimizer*>(d.optimizer());
      ASSERT_TRUE(opt != nullptr);
      EXPECT_EQ(v->value, opt->solver().getParameter(param))
        << "Failed option: " << option;
    }
    EXPECT_EQ(end - start + 1, count);
  }
}

void IlogCPTest::CheckDblCPOption(const char *option,
    IloCP::NumParam param, double good, double bad) {
  EXPECT_TRUE(ParseOptions("optimizer=cp", Option(option, good).c_str()));
  CPOptimizer *opt = dynamic_cast<CPOptimizer*>(d.optimizer());
  ASSERT_TRUE(opt != nullptr);
  EXPECT_EQ(good, opt->solver().getParameter(param))
    << "Failed option: " << option;

  EXPECT_FALSE(ParseOptions("optimizer=cp", Option(option, bad).c_str()));
  EXPECT_FALSE(ParseOptions("optimizer=cplex", Option(option, good).c_str()));
}

TEST_F(IlogCPTest, ConvertNum) {
  EXPECT_EQ("0.42", str(d.Visit(NewNum(0.42))));
}

TEST_F(IlogCPTest, ConvertVar) {
  EXPECT_EQ("theta", str(d.Visit(NewVar(2))));
}

TEST_F(IlogCPTest, ConvertPlus) {
  EXPECT_EQ("x + 42", str(d.Visit(NewBinary(OPPLUS, NewVar(0), NewNum(42)))));
  EXPECT_EQ("x + y", str(d.Visit(NewBinary(OPPLUS, NewVar(0), NewVar(1)))));
}

TEST_F(IlogCPTest, ConvertMinus) {
  EXPECT_EQ("x + -42", str(d.Visit(
    NewBinary(OPMINUS, NewVar(0), NewNum(42)))));
  EXPECT_EQ("x + -1 * y", str(d.Visit(
    NewBinary(OPMINUS, NewVar(0), NewVar(1)))));
}

TEST_F(IlogCPTest, ConvertMult) {
  EXPECT_EQ("42 * x", str(d.Visit(
    NewBinary(OPMULT, NewVar(0), NewNum(42)))));
  EXPECT_EQ("x * y", str(d.Visit(
    NewBinary(OPMULT, NewVar(0), NewVar(1)))));
}

TEST_F(IlogCPTest, ConvertDiv) {
  EXPECT_EQ("x / 42", str(d.Visit(
    NewBinary(OPDIV, NewVar(0), NewNum(42)))));
  EXPECT_EQ("x / y", str(d.Visit(
    NewBinary(OPDIV, NewVar(0), NewVar(1)))));
}

TEST_F(IlogCPTest, ConvertRem) {
  EXPECT_EQ("x + trunc(x / y ) * y * -1",
    str(d.Visit(NewBinary(OPREM, NewVar(0), NewVar(1)))));
  EXPECT_EQ(0, EvalRem(9, 3));
  EXPECT_EQ(2, EvalRem(8, 3));
  EXPECT_EQ(-2, EvalRem(-8, 3));
  EXPECT_EQ(2, EvalRem(8, -3));
  EXPECT_EQ(-2, EvalRem(-8, -3));
  EXPECT_EQ(1.5, EvalRem(7.5, 3));
}

TEST_F(IlogCPTest, ConvertPow) {
  EXPECT_EQ("x ^ 42", str(d.Visit(
    NewBinary(OPPOW, NewVar(0), NewNum(42)))));
  EXPECT_EQ("x ^ y", str(d.Visit(
    NewBinary(OPPOW, NewVar(0), NewVar(1)))));
}

TEST_F(IlogCPTest, ConvertLess) {
  EXPECT_EQ("max(x + -42 , 0)", str(d.Visit(
    NewBinary(OPLESS, NewVar(0), NewNum(42)))));
  EXPECT_EQ("max(x + -1 * y , 0)", str(d.Visit(
    NewBinary(OPLESS, NewVar(0), NewVar(1)))));
}

TEST_F(IlogCPTest, ConvertMin) {
  EXPECT_EQ("min( [x , y , 42 ])", str(d.Visit(
    NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertMax) {
  EXPECT_EQ("max([x , y , 42 ])", str(d.Visit(
    NewVarArg(MAXLIST, NewVar(0), NewVar(1), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertFloor) {
  EXPECT_EQ("floor(x )", str(d.Visit(NewUnary(FLOOR, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertCeil) {
  EXPECT_EQ("ceil(x )", str(d.Visit(NewUnary(CEIL, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertAbs) {
  EXPECT_EQ("abs(x )", str(d.Visit(NewUnary(ABS, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertUMinus) {
  EXPECT_EQ("-1 * x", str(d.Visit(NewUnary(OPUMINUS, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertInvalidExprThrows) {
  int ops[] = {
      OPOR, OPAND, LT, LE, EQ, GE, GT, NE, OPNOT,
      OPFUNCALL, OPIFSYM, OPHOL,
      OPATLEAST, OPATMOST, OPEXACTLY,
      OPNOTATLEAST, OPNOTATMOST, OPNOTEXACTLY,
      OPIMPELSE, OP_IFF, N_OPS, -1, 500
  };
  size_t i = 0;
  for (size_t num_ops = sizeof(ops) / sizeof(*ops); i < num_ops; ++i) {
    EXPECT_THROW(d.Visit(
      NewBinary(ops[i], NewVar(0), NewVar(1))), InvalidNumericExprError);
  }
  // Paranoid: make sure that the loop body has been executed enough times.
  EXPECT_EQ(23u, i);

  EXPECT_THROW(d.Visit(
    NewSum(OPNUMBEROFs, NewVar(0), NewVar(1))), InvalidNumericExprError);
  EXPECT_THROW(d.Visit(
    NewSum(OPALLDIFF, NewVar(0), NewVar(1))), InvalidNumericExprError);
  EXPECT_THROW(d.Visit(
    NewSum(ANDLIST, NewVar(0), NewVar(1))), InvalidNumericExprError);
  EXPECT_THROW(d.Visit(
    NewSum(ORLIST, NewVar(0), NewVar(1))), InvalidNumericExprError);
  
  EXPECT_THROW(d.Visit(
    NewBinary(OPprecision, NewVar(0), NewVar(1))), UnsupportedExprError);
}

TEST_F(IlogCPTest, ConvertIf) {
  EXPECT_EQ("IloNumVar(7)[-inf..inf]", str(d.Visit(NewIf(OPIFnl,
      NewBinary(EQ, NewVar(0), NewNum(0)), NewVar(1), NewNum(42)))));

  IloModel::Iterator iter(mod_);
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

TEST_F(IlogCPTest, ConvertTanh) {
  EXPECT_EQ("exp(2 * x ) + -1 / exp(2 * x ) + 1",
            str(d.Visit(NewUnary(OP_tanh, NewVar(0)))));
  // Concert incorrectly omits brackets around the dividend and divisor
  // above, so test also by evaluating the expression at several points.
  IloExpr e(d.Visit(NewUnary(OP_tanh, NewNum(1))));
  EXPECT_NEAR(0.761594, eval(e), 1e-5);
  e = d.Visit(NewUnary(OP_tanh, NewNum(0)));
  EXPECT_EQ(0, eval(e));
  e = d.Visit(NewUnary(OP_tanh, NewNum(-2)));
  EXPECT_NEAR(-0.964027, eval(e), 1e-5);
}

TEST_F(IlogCPTest, ConvertTan) {
  EXPECT_EQ("tan(x )", str(d.Visit(NewUnary(OP_tan, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertSqrt) {
  EXPECT_EQ("x ^ 0.5",
            str(d.Visit(NewUnary(OP_sqrt, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertSinh) {
  EXPECT_EQ("exp(x ) * 0.5 + exp(-1 * x ) * -0.5",
            str(d.Visit(NewUnary(OP_sinh, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertSin) {
  EXPECT_EQ("sin(x )", str(d.Visit(NewUnary(OP_sin, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertLog10) {
  EXPECT_EQ("log(x )/ 2.30259",
            str(d.Visit(NewUnary(OP_log10, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertLog) {
  EXPECT_EQ("log(x )", str(d.Visit(NewUnary(OP_log, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertExp) {
  EXPECT_EQ("exp(x )", str(d.Visit(NewUnary(OP_exp, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertCosh) {
  EXPECT_EQ("exp(x ) * 0.5 + exp(-1 * x ) * 0.5",
            str(d.Visit(NewUnary(OP_cosh, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertCos) {
  EXPECT_EQ("cos(x )", str(d.Visit(NewUnary(OP_cos, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertAtanh) {
  EXPECT_EQ("log(x + 1 ) * 0.5 + log(-1 * x + 1 ) * -0.5",
            str(d.Visit(NewUnary(OP_atanh, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertAtan2) {
  EXPECT_EQ("IloNumVar(8)[-inf..inf]",
            str(d.Visit(NewBinary(OP_atan2, NewVar(1), NewVar(0)))));

  IloModel::Iterator iter(mod_);
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

TEST_F(IlogCPTest, ConvertAtan) {
  EXPECT_EQ("arc-tan(x )",
            str(d.Visit(NewUnary(OP_atan, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertAsinh) {
  EXPECT_EQ("log(x + square(x ) + 1 ^ 0.5)",
            str(d.Visit(NewUnary(OP_asinh, NewVar(0)))));
  // Concert incorrectly omits brackets around square(x) + 1
  // above, so test also by evaluating the expression at several points.
  IloExpr e(d.Visit(NewUnary(OP_asinh, NewNum(1))));
  EXPECT_NEAR(0.881373, eval(e), 1e-5);
  e = d.Visit(NewUnary(OP_asinh, NewNum(0)));
  EXPECT_EQ(0, eval(e));
  e = d.Visit(NewUnary(OP_asinh, NewNum(-2)));
  EXPECT_NEAR(-1.443635, eval(e), 1e-5);
}

TEST_F(IlogCPTest, ConvertAsin) {
  EXPECT_EQ("arc-sin(x )",
            str(d.Visit(NewUnary(OP_asin, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertAcosh) {
  EXPECT_EQ("log(x + x + 1 ^ 0.5 * x + -1 ^ 0.5)",
            str(d.Visit(NewUnary(OP_acosh, NewVar(0)))));
  // Concert incorrectly omits brackets around x + 1 and x + -1
  // above, so test also by evaluating the expression at several points.
  IloExpr e(d.Visit(NewUnary(OP_acosh, NewNum(1))));
  EXPECT_NEAR(0, eval(e), 1e-5);
  e = d.Visit(NewUnary(OP_acosh, NewNum(10)));
  EXPECT_NEAR(2.993222, eval(e), 1e-5);
  e = d.Visit(NewUnary(OP_acosh, NewNum(0)));
  double n = eval(e);
  EXPECT_TRUE(n != n);
}

TEST_F(IlogCPTest, ConvertAcos) {
  EXPECT_EQ("arc-cos(x )", str(d.Visit(NewUnary(OP_acos, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertSum) {
  EXPECT_EQ("x + y + 42", str(d.Visit(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertIntDiv) {
  EXPECT_EQ("trunc(x / y )",
    str(d.Visit(NewBinary(OPintDIV, NewVar(0), NewVar(1)))));
}

TEST_F(IlogCPTest, ConvertRound) {
  EXPECT_EQ("round(x )",
    str(d.Visit(NewBinary(OPround, NewVar(0), NewNum(0)))));

  EXPECT_EQ(1235, eval(d.Visit(
    NewBinary(OPround, NewNum(1234.56), NewNum(0)))));
  EXPECT_EQ(3, eval(d.Visit(
    NewBinary(OPround, NewNum(2.5), NewNum(0)))));
  EXPECT_EQ(-2, eval(d.Visit(
    NewBinary(OPround, NewNum(-2.5), NewNum(0)))));

  EXPECT_THROW(d.Visit(
    NewBinary(OPround, NewVar(0), NewVar(1))), UnsupportedExprError);
}

TEST_F(IlogCPTest, ConvertTrunc) {
  EXPECT_EQ("trunc(x )", str(d.Visit(
    NewBinary(OPtrunc, NewVar(0), NewNum(0)))));
  EXPECT_EQ(1234, eval(d.Visit(
    NewBinary(OPtrunc, NewNum(1234.56), NewNum(0)))));
  EXPECT_THROW(d.Visit(
    NewBinary(OPtrunc, NewVar(0), NewVar(1))), UnsupportedExprError);
}

TEST_F(IlogCPTest, Convert1Pow) {
  EXPECT_EQ("x ^ 42", str(d.Visit(NewBinary(OP1POW, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, Convert2Pow) {
  EXPECT_EQ("square(x )", str(d.Visit(NewUnary(OP2POW, NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertCPow) {
  EXPECT_EQ("42 ^ x", str(d.Visit(
    NewBinary(OPCPOW, NewNum(42), NewVar(0)))));
}

TEST_F(IlogCPTest, ConvertPLTerm) {
  double args[] = {-1, 5, 0, 10, 1};
  EXPECT_EQ("piecewiselinear(x[0..1] , [5, 10], [-1, 0, 1], 0, 0)",
            str(d.Visit(NewPLTerm(5, args, 0))));
}

TEST_F(IlogCPTest, ConvertCount) {
  NumericExpr a(NewBinary(EQ, NewVar(0), NewNum(0)));
  NumericExpr b(NewBinary(LE, NewVar(1), NewNum(42)));
  NumericExpr c(NewBinary(GE, NewVar(2), NewNum(0)));
  EXPECT_EQ("x == 0 + y <= 42 + 0 <= theta",
      str(d.Visit(NewSum(OPCOUNT, a, b, c))));
}

TEST_F(IlogCPTest, ConvertNumberOf) {
  d.use_numberof();
  EXPECT_EQ("x == theta + y == theta",
      str(d.Visit(NewSum(OPNUMBEROF, NewVar(2), NewVar(0), NewVar(1)))));
  d.use_numberof(false);
  EXPECT_EQ("x == 42 + y == 42",
      str(d.Visit(NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)))));
}

TEST_F(IlogCPTest, IloArrayCopyingIsCheap) {
  IloIntArray array(env_);
  array.add(42);
  EXPECT_TRUE(array.getImpl() != nullptr);
  EXPECT_EQ(array.getImpl(), IloIntArray(array).getImpl());
}

TEST_F(IlogCPTest, ConvertSingleNumberOfToIloDistribute) {
  d.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  EXPECT_EQ("IloIntVar(10)" + bounds, str(d.Visit(
      NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)))));
  d.FinishBuildingNumberOf();
  IloModel::Iterator iter(mod_);
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

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithSameValuesToIloDistribute) {
  d.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  NumericExpr expr(NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)));
  EXPECT_EQ("IloIntVar(10)" + bounds, str(d.Visit(NumericExpr(ptr(expr)))));
  EXPECT_EQ("IloIntVar(10)" + bounds, str(d.Visit(
      NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)))));
  d.FinishBuildingNumberOf();
  IloModel::Iterator iter(mod_);
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

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithDiffValuesToIloDistribute) {
  d.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  EXPECT_EQ("IloIntVar(10)" + bounds,
      str(d.Visit(NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)))));
  EXPECT_EQ("IloIntVar(12)" + bounds,
      str(d.Visit(NewSum(OPNUMBEROF, NewNum(43), NewVar(0), NewVar(1)))));
  d.FinishBuildingNumberOf();
  IloModel::Iterator iter(mod_);
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

TEST_F(IlogCPTest, ConvertTwoNumberOfsWithDiffExprs) {
  d.use_numberof();
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  EXPECT_EQ("IloIntVar(10)" + bounds,
      str(d.Visit(NewSum(OPNUMBEROF, NewNum(42), NewVar(0), NewVar(1)))));
  EXPECT_EQ("IloIntVar(15)" + bounds,
      str(d.Visit(NewSum(OPNUMBEROF, NewNum(42), NewVar(2)))));
  d.FinishBuildingNumberOf();
  IloModel::Iterator iter(mod_);
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

TEST_F(IlogCPTest, ConvertInvalidLogicalExprThrows) {
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
    EXPECT_THROW(ConvertLogical(
      NewBinary(ops[i], NewVar(0), NewVar(1))),
      InvalidLogicalExprError);
  }
  // Paranoid: make sure that the loop body has been executed enough times.
  EXPECT_EQ(41u, i);

  EXPECT_THROW(ConvertLogical(
    NewVarArg(MINLIST, NewVar(0), NewVar(1))),
    InvalidLogicalExprError);
  EXPECT_THROW(ConvertLogical(
    NewVarArg(MAXLIST, NewVar(0), NewVar(1))),
    InvalidLogicalExprError);
  EXPECT_THROW(ConvertLogical(
    NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(0))),
    InvalidLogicalExprError);
  EXPECT_THROW(ConvertLogical(
    NewIf(OPIFSYM, NewVar(0), NewVar(1), NewNum(0))),
    InvalidLogicalExprError);
  EXPECT_THROW(ConvertLogical(
    NewSum(OPSUMLIST, NewVar(0), NewVar(1))),
    InvalidLogicalExprError);
  EXPECT_THROW(ConvertLogical(
    NewSum(OPCOUNT, NewVar(0), NewVar(1))),
    InvalidLogicalExprError);
  EXPECT_THROW(ConvertLogical(
    NewSum(OPNUMBEROF, NewVar(0), NewVar(1))),
    InvalidLogicalExprError);
  EXPECT_THROW(ConvertLogical(
    NewSum(OPNUMBEROFs, NewVar(0), NewVar(1))),
    InvalidLogicalExprError);
  EXPECT_THROW(ConvertLogical(
    NewPLTerm(0, 0, 0)), InvalidLogicalExprError);
}

TEST_F(IlogCPTest, ConvertFalse) {
  EXPECT_EQ("IloNumVar(4)[1..1] == 0", str(ConvertLogical(NewNum(0))));
}

TEST_F(IlogCPTest, ConvertTrue) {
  EXPECT_EQ("IloNumVar(4)[1..1] == 1", str(ConvertLogical(NewNum(1))));
  EXPECT_EQ("IloNumVar(7)[1..1] == 1", str(ConvertLogical(NewNum(42))));
}

TEST_F(IlogCPTest, ConvertLT) {
  EXPECT_EQ("x <= 41", str(ConvertLogical(
      NewBinary(LT, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertLE) {
  EXPECT_EQ("x <= 42", str(ConvertLogical(
      NewBinary(LE, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertEQ) {
  EXPECT_EQ("x == 42", str(ConvertLogical(
      NewBinary(EQ, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertGE) {
  EXPECT_EQ("42 <= x", str(ConvertLogical(
      NewBinary(GE, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertGT) {
  EXPECT_EQ("43 <= x", str(ConvertLogical(
      NewBinary(GT, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertNE) {
  EXPECT_EQ("x != 42", str(ConvertLogical(
      NewBinary(NE, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertAtMost) {
  EXPECT_EQ("42 <= x", str(ConvertLogical(
      NewBinary(OPATMOST, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertNotAtMost) {
  IloConstraint c(ConvertLogical(
      NewBinary(OPNOTATMOST, NewVar(0), NewNum(42))));
  IloNotI *n = dynamic_cast<IloNotI*>(c.getImpl());
  ASSERT_TRUE(n != nullptr);
  EXPECT_EQ("42 <= x", str(n->getConstraint()));
}

TEST_F(IlogCPTest, ConvertAtLeast) {
  EXPECT_EQ("x <= 42", str(ConvertLogical(
      NewBinary(OPATLEAST, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertNotAtLeast) {
  IloConstraint c(ConvertLogical(
      NewBinary(OPNOTATLEAST, NewVar(0), NewNum(42))));
  IloNotI *n = dynamic_cast<IloNotI*>(c.getImpl());
  ASSERT_TRUE(n != nullptr);
  EXPECT_EQ("x <= 42", str(n->getConstraint()));
}

TEST_F(IlogCPTest, ConvertExactly) {
  EXPECT_EQ("x == 42", str(ConvertLogical(
      NewBinary(OPEXACTLY, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertNotExactly) {
  EXPECT_EQ("x != 42", str(ConvertLogical(
      NewBinary(OPNOTEXACTLY, NewVar(0), NewNum(42)))));
}

TEST_F(IlogCPTest, ConvertOr) {
  IloConstraint c(ConvertLogical(NewBinary(OPOR,
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2)))));
  IloIfThenI *ifThen = dynamic_cast<IloIfThenI*>(c.getImpl());
  ASSERT_TRUE(ifThen != nullptr);
  IloNotI *n = dynamic_cast<IloNotI*>(ifThen->getLeft().getImpl());
  ASSERT_TRUE(n != nullptr);
  EXPECT_EQ("x == 1", str(n->getConstraint()));
  EXPECT_EQ("x == 2", str(ifThen->getRight()));
}

TEST_F(IlogCPTest, CheckOrTruthTable) {
  IloNumVarArray vars = d.vars();
  vars[0].setBounds(0, 0);
  vars[1].setBounds(0, 0);
  mod_.add(ConvertLogical(NewBinary(OPOR,
        NewBinary(EQ, NewVar(0), NewNum(1)),
        NewBinary(EQ, NewVar(1), NewNum(1)))));
  IloCP cp(mod_);
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

TEST_F(IlogCPTest, ConvertExists) {
  EXPECT_EQ("(x == 1 ) || (x == 2 ) || (x == 3 )",
      str(ConvertLogical(NewSum(ORLIST,
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2)),
      NewBinary(EQ, NewVar(0), NewNum(3))))));
}

TEST_F(IlogCPTest, ConvertAnd) {
  EXPECT_EQ("(x == 1 ) && (x == 2 )", str(ConvertLogical(NewBinary(OPAND,
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2))))));
}

TEST_F(IlogCPTest, ConvertForAll) {
  EXPECT_EQ("(x == 1 ) && (x == 2 ) && (x == 3 )",
      str(ConvertLogical(NewSum(ANDLIST,
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2)),
      NewBinary(EQ, NewVar(0), NewNum(3))))));
}

TEST_F(IlogCPTest, ConvertNot) {
  IloConstraint c(ConvertLogical(
      NewUnary(OPNOT, NewBinary(LE, NewVar(0), NewNum(42)))));
  IloNotI *n = dynamic_cast<IloNotI*>(c.getImpl());
  ASSERT_TRUE(n != nullptr);
  EXPECT_EQ("x <= 42", str(n->getConstraint()));
}

TEST_F(IlogCPTest, ConvertIff) {
  EXPECT_EQ("x == 1 == x == 2", str(ConvertLogical(NewBinary(OP_IFF,
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2))))));
}

TEST_F(IlogCPTest, ConvertImpElse) {
  IloConstraint con(ConvertLogical(NewIf(OPIMPELSE,
      NewBinary(EQ, NewVar(0), NewNum(0)),
      NewBinary(EQ, NewVar(0), NewNum(1)),
      NewBinary(EQ, NewVar(0), NewNum(2)))));
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

TEST_F(IlogCPTest, ConvertAllDiff) {
  IloConstraint con(ConvertLogical(
      NewSum(OPALLDIFF, NewVar(0), NewNum(42))));
  IloAllDiffI* diff = dynamic_cast<IloAllDiffI*>(con.getImpl());
  ASSERT_TRUE(diff != nullptr);
  std::ostringstream os;
  os << "[" << IloIntMin << ".." << IloIntMax << "]";
  string bounds = os.str();
  EXPECT_EQ("[x[0..1] , IloIntVar(4)" + bounds + " ]",
      str(diff->getExprArray()));

  IloModel::Iterator iter(mod_);
  ASSERT_TRUE(iter.ok());
  EXPECT_EQ("IloIntVar(4)" + bounds +" == 42", str(*iter));
  ++iter;
  EXPECT_FALSE(iter.ok());
}

// ----------------------------------------------------------------------------
// Util tests

TEST_F(IlogCPTest, EqualNum) {
  EXPECT_TRUE(Equal(NewNum(0.42), NewNum(0.42)));
  EXPECT_FALSE(Equal(NewNum(0.42), NewNum(42)));
}

TEST_F(IlogCPTest, EqualVar) {
  EXPECT_TRUE(Equal(NewVar(0), NewVar(0)));
  EXPECT_FALSE(Equal(NewVar(0), NewVar(1)));
  EXPECT_FALSE(Equal(NewVar(0), NewNum(0)));
}

TEST_F(IlogCPTest, EqualUnary) {
  EXPECT_TRUE(Equal(NewUnary(OPUMINUS, NewVar(0)),
                    NewUnary(OPUMINUS, NewVar(0))));
  EXPECT_FALSE(Equal(NewUnary(OPUMINUS, NewVar(0)),
                     NewVar(0)));
  EXPECT_FALSE(Equal(NewUnary(OPUMINUS, NewVar(0)),
                     NewUnary(FLOOR, NewVar(0))));
  EXPECT_FALSE(Equal(NewUnary(OPUMINUS, NewVar(0)),
                     NewUnary(OPUMINUS, NewVar(1))));
}

TEST_F(IlogCPTest, EqualBinary) {
  EXPECT_TRUE(Equal(NewBinary(OPPLUS, NewVar(0), NewNum(42)),
                        NewBinary(OPPLUS, NewVar(0), NewNum(42))));
  EXPECT_FALSE(Equal(NewBinary(OPPLUS, NewVar(0), NewNum(42)),
                         NewBinary(OPMINUS, NewVar(0), NewNum(42))));
  EXPECT_FALSE(Equal(NewBinary(OPPLUS, NewVar(0), NewNum(42)),
                         NewBinary(OPPLUS, NewNum(42), NewVar(0))));
  EXPECT_FALSE(Equal(NewBinary(OPPLUS, NewVar(0), NewNum(42)),
                         NewBinary(OPPLUS, NewVar(0), NewNum(0))));
  EXPECT_FALSE(Equal(NewNum(42),
                         NewBinary(OPPLUS, NewVar(0), NewNum(42))));
}

TEST_F(IlogCPTest, EqualVarArg) {
  EXPECT_TRUE(Equal(
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)),
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)),
      NewVarArg(MINLIST, NewVar(0), NewVar(1))));
  EXPECT_FALSE(Equal(
      NewVarArg(MINLIST, NewVar(0), NewVar(1)),
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)),
      NewVarArg(MAXLIST, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)),
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(0))));
  EXPECT_FALSE(Equal(
      NewVarArg(MINLIST, NewVar(0), NewVar(1), NewNum(42)),
      NewNum(42)));
}

TEST_F(IlogCPTest, EqualPLTerm) {
  double args[] = {-1, 5, 0, 10, 1};
  EXPECT_TRUE(Equal(
      NewPLTerm(5, args, 0),
      NewPLTerm(5, args, 0)));
  EXPECT_FALSE(Equal(
      NewPLTerm(5, args, 0),
      NewPLTerm(3, args, 0)));
  EXPECT_FALSE(Equal(
      NewPLTerm(5, args, 0),
      NewPLTerm(5, args, 1)));
  double args2[] = {-1, 5, 0, 11, 1};
  EXPECT_FALSE(Equal(
      NewPLTerm(5, args, 0),
      NewPLTerm(5, args2, 0)));
  EXPECT_FALSE(Equal(
      NewPLTerm(5, args, 0),
      NewNum(42)));
}

TEST_F(IlogCPTest, EqualIf) {
  EXPECT_TRUE(Equal(
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(42)),
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(42)),
      NewIf(OPIFSYM, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(42)),
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(0))));
  EXPECT_FALSE(Equal(
      NewIf(OPIFnl, NewVar(0), NewVar(1), NewNum(42)),
      NewNum(42)));
}

TEST_F(IlogCPTest, EqualSum) {
  EXPECT_TRUE(Equal(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1))));
  EXPECT_FALSE(Equal(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(0))));
  EXPECT_FALSE(Equal(
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42)),
      NewNum(42)));
}

TEST_F(IlogCPTest, EqualCount) {
  EXPECT_TRUE(Equal(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)),
      NewSum(OPCOUNT, NewVar(0), NewVar(1))));
  EXPECT_FALSE(Equal(
      NewSum(OPCOUNT, NewVar(0), NewVar(1)),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)),
      NewSum(OPSUMLIST, NewVar(0), NewVar(1), NewNum(42))));
  EXPECT_FALSE(Equal(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)),
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(0))));
  EXPECT_FALSE(Equal(
      NewSum(OPCOUNT, NewVar(0), NewVar(1), NewNum(42)),
      NewNum(42)));
}

TEST_F(IlogCPTest, EqualExprThrowsOnUnsupportedOp) {
  EXPECT_THROW(Equal(
      NewUnary(OPFUNCALL, NumericExpr()),
      NewUnary(OPFUNCALL, NumericExpr())),
      UnsupportedExprError);
  EXPECT_THROW(Equal(
      NewUnary(OPHOL, NumericExpr()),
      NewUnary(OPHOL, NumericExpr())),
      UnsupportedExprError);
}

// ----------------------------------------------------------------------------
// Driver tests

TEST_F(IlogCPTest, Usage) {
  FILE *saved_stderr = Stderr;
  Stderr = fopen("out", "w");
  d.run((Args() + "ilogcp").get());
  fclose(Stderr);
  Stderr = saved_stderr;

  ifstream ifs("out");
  enum { BUFFER_SIZE = 4096 };
  char buffer[BUFFER_SIZE];
  string text;
  while (ifs) {
    ifs.read(buffer, BUFFER_SIZE);
    text += string(buffer, ifs.gcount());
  }
  EXPECT_TRUE(text.find("usage: ") != string::npos);
}

TEST_F(IlogCPTest, ObjConst) {
  EXPECT_EQ(0, RunDriver(DATA_DIR "objconst"));
  IloModel::Iterator iter(mod_);
  ASSERT_TRUE(iter.ok());
  IloObjective obj = (*iter).asObjective();
  EXPECT_EQ(42, obj.getConstant());
}

TEST_F(IlogCPTest, CPOptimizerDoesntSupportContinuousVars) {
  EXPECT_EQ(1, RunDriver(DATA_DIR "objconst", "optimizer=cp"));
}

TEST_F(IlogCPTest, SolveNumberOfCplex) {
  d.use_numberof(false);
  RunDriver(DATA_DIR "numberof", "optimizer=cplex");
}

TEST_F(IlogCPTest, SolveAssign0) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign0").obj);
}

TEST_F(IlogCPTest, SolveAssign1) {
  EXPECT_EQ(6, Solve(DATA_DIR "assign1").obj);
}

TEST_F(IlogCPTest, SolveBalassign0) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign0").obj);
}

TEST_F(IlogCPTest, SolveBalassign1) {
  EXPECT_EQ(14, Solve(DATA_DIR "balassign1").obj);
}

TEST_F(IlogCPTest, SolveFlowshp0) {
  EXPECT_NEAR(22, Solve(DATA_DIR "flowshp0").obj, 1e-5);
}

TEST_F(IlogCPTest, SolveFlowshp1) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp1").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(IlogCPTest, DISABLED_SolveFlowshp2) {
  EXPECT_EQ(22, Solve(DATA_DIR "flowshp2").obj);
}

TEST_F(IlogCPTest, SolveGrpassign0) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign0").obj);
}

// Disabled because variables in subscripts are not yet allowed.
TEST_F(IlogCPTest, DISABLED_SolveGrpassign1) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1").obj);
}

// Disabled because object-valued variables are not yet allowed.
TEST_F(IlogCPTest, DISABLED_SolveGrpassign1a) {
  EXPECT_EQ(61, Solve(DATA_DIR "grpassign1a").obj);
}

TEST_F(IlogCPTest, SolveMagic) {
  EXPECT_TRUE(Solve(DATA_DIR "magic").solved);
}

TEST_F(IlogCPTest, SolveMapcoloring) {
  EXPECT_TRUE(Solve(DATA_DIR "mapcoloring").solved);
}

TEST_F(IlogCPTest, SolveMoney) {
  EXPECT_TRUE(Solve(DATA_DIR "money").solved);
}

TEST_F(IlogCPTest, SolveNQueens) {
  EXPECT_TRUE(Solve(DATA_DIR "nqueens").solved);
}

TEST_F(IlogCPTest, SolveNQueens0) {
  EXPECT_EQ(0, Solve(DATA_DIR "nqueens0").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(IlogCPTest, DISABLED_SolveParty1) {
  EXPECT_EQ(61, Solve(DATA_DIR "party1").obj);
}

// Disabled because it's too difficult to solve.
TEST_F(IlogCPTest, DISABLED_SolveParty2) {
  EXPECT_EQ(3, Solve(DATA_DIR "party2").obj);
}

TEST_F(IlogCPTest, SolveSched0) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched0").obj);
}

TEST_F(IlogCPTest, SolveSched1) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched1").obj);
}

TEST_F(IlogCPTest, SolveSched2) {
  EXPECT_EQ(5, Solve(DATA_DIR "sched2").obj);
}

TEST_F(IlogCPTest, SolveSeq0) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0").obj, 1e-5);
}

TEST_F(IlogCPTest, SolveSeq0a) {
  EXPECT_NEAR(332, Solve(DATA_DIR "seq0a").obj, 1e-5);
}

TEST_F(IlogCPTest, SolveSudokuHard) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuHard").solved);
}

TEST_F(IlogCPTest, SolveSudokuVeryEasy) {
  EXPECT_TRUE(Solve(DATA_DIR "sudokuVeryEasy").solved);
}

// ----------------------------------------------------------------------------
// Option tests

TEST_F(IlogCPTest, VersionOption) {
  EXPECT_EQ(0, d.show_version());
  EXPECT_TRUE(ParseOptions("version"));
  EXPECT_EQ(1, d.show_version());
}

TEST_F(IlogCPTest, WantsolOption) {
  EXPECT_EQ(0, d.wantsol());
  EXPECT_TRUE(ParseOptions("wantsol=1"));
  EXPECT_EQ(1, d.wantsol());
  EXPECT_TRUE(ParseOptions("wantsol=5"));
  EXPECT_EQ(5, d.wantsol());
}

TEST_F(IlogCPTest, DebugExprOption) {
  EXPECT_TRUE(ParseOptions("debugexpr=0"));
  EXPECT_EQ(0, d.get_option(IlogCPDriver::DEBUGEXPR));
  EXPECT_TRUE(ParseOptions("debugexpr=1"));
  EXPECT_EQ(1, d.get_option(IlogCPDriver::DEBUGEXPR));
  EXPECT_TRUE(ParseOptions("debugexpr=42"));
  EXPECT_EQ(42, d.get_option(IlogCPDriver::DEBUGEXPR));
  EXPECT_FALSE(ParseOptions("debugexpr=oops"));
}

TEST_F(IlogCPTest, OptimizerOption) {
  EXPECT_EQ(IlogCPDriver::AUTO, d.get_option(IlogCPDriver::OPTIMIZER));

  EXPECT_TRUE(ParseOptions("optimizer=cplex"));
  EXPECT_EQ(IlogCPDriver::CPLEX, d.get_option(IlogCPDriver::OPTIMIZER));
  EXPECT_TRUE(dynamic_cast<IloCplexI*>(d.alg().getImpl()) != nullptr);

  EXPECT_TRUE(ParseOptions("optimizer=cp"));
  EXPECT_EQ(IlogCPDriver::CP, d.get_option(IlogCPDriver::OPTIMIZER));
  EXPECT_TRUE(dynamic_cast<IloCplexI*>(d.alg().getImpl()) == nullptr);
}

TEST_F(IlogCPTest, TimingOption) {
  EXPECT_TRUE(ParseOptions("timing=0"));
  EXPECT_EQ(0, d.get_option(IlogCPDriver::TIMING));
  EXPECT_TRUE(ParseOptions("timing=1"));
  EXPECT_EQ(1, d.get_option(IlogCPDriver::TIMING));
  EXPECT_FALSE(ParseOptions("timing=42"));
  EXPECT_FALSE(ParseOptions("timing=oops"));
}

TEST_F(IlogCPTest, UseNumberOfOption) {
  EXPECT_TRUE(ParseOptions("usenumberof=0"));
  EXPECT_EQ(0, d.get_option(IlogCPDriver::USENUMBEROF));
  EXPECT_TRUE(ParseOptions("usenumberof=1"));
  EXPECT_EQ(1, d.get_option(IlogCPDriver::USENUMBEROF));
  EXPECT_FALSE(ParseOptions("usenumberof=42"));
  EXPECT_FALSE(ParseOptions("timing=oops"));
}

TEST_F(IlogCPTest, CPFlagOptions) {
  const EnumValue flags[] = {
      {"off", IloCP::Off},
      {"on",  IloCP::On},
      {nullptr}
  };
  CheckIntCPOption("constraintaggregation", IloCP::ConstraintAggregation,
      0, 1, IloCP::Off, false, flags);
  CheckIntCPOption("dynamicprobing", IloCP::DynamicProbing,
      0, 1, IloCP::Off, true, flags);
  CheckIntCPOption("temporalrelaxation", IloCP::TemporalRelaxation,
      0, 1, IloCP::Off, false, flags);
}

TEST_F(IlogCPTest, CPInferenceLevelOptions) {
  const EnumValue inf_levels[] = {
      {"default",  IloCP::Default},
      {"low",      IloCP::Low},
      {"basic",    IloCP::Basic},
      {"medium",   IloCP::Medium},
      {"extended", IloCP::Extended},
      {nullptr}
  };
  CheckIntCPOption("alldiffinferencelevel", IloCP::AllDiffInferenceLevel,
      0, 4, IloCP::Default, false, inf_levels);
  CheckIntCPOption("defaultinferencelevel", IloCP::DefaultInferenceLevel,
      1, 4, IloCP::Default, false, inf_levels + 1);
  CheckIntCPOption("distributeinferencelevel", IloCP::DistributeInferenceLevel,
      0, 4, IloCP::Default, false, inf_levels);
}

TEST_F(IlogCPTest, CPDefaultVerbosityQuiet) {
  EXPECT_TRUE(ParseOptions("optimizer=cp"));
  CPOptimizer *opt = dynamic_cast<CPOptimizer*>(d.optimizer());
  ASSERT_TRUE(opt != nullptr);
  EXPECT_EQ(IloCP::Quiet, opt->solver().getParameter(IloCP::LogVerbosity));
}

TEST_F(IlogCPTest, CPVerbosityOptions) {
  const EnumValue verbosities[] = {
      {"quiet",   IloCP::Quiet},
      {"terse",   IloCP::Terse},
      {"normal",  IloCP::Normal},
      {"verbose", IloCP::Verbose},
      {nullptr}
  };
  CheckIntCPOption("logverbosity", IloCP::LogVerbosity,
      0, 3, IloCP::Quiet, false, verbosities);
  CheckIntCPOption("outlev", IloCP::LogVerbosity,
      0, 3, IloCP::Quiet, false, verbosities);
  CheckIntCPOption("propagationlog", IloCP::PropagationLog,
      0, 3, IloCP::Quiet, false, verbosities);
}

TEST_F(IlogCPTest, CPSearchTypeOption) {
  const EnumValue types[] = {
      {"depthfirst", IloCP::DepthFirst},
      {"restart",    IloCP::Restart},
      {"multipoint", IloCP::MultiPoint},
      {nullptr}
  };
  CheckIntCPOption("searchtype", IloCP::SearchType,
      0, 2, IloCP::DepthFirst, CPX_VERSION > 1220, types);
}

TEST_F(IlogCPTest, CPTimeModeOption) {
  const EnumValue modes[] = {
      {"cputime",     IloCP::CPUTime},
      {"elapsedtime", IloCP::ElapsedTime},
      {nullptr}
  };
  CheckIntCPOption("timemode", IloCP::TimeMode,
      0, 1, IloCP::CPUTime, false, modes);
}

TEST_F(IlogCPTest, CPOptions) {
  CheckIntCPOption("branchlimit", IloCP::BranchLimit, 0, INT_MAX);
  CheckIntCPOption("choicepointlimit", IloCP::ChoicePointLimit, 0, INT_MAX);
  CheckDblCPOption("dynamicprobingstrength",
      IloCP::DynamicProbingStrength, 42, -1);
  CheckIntCPOption("faillimit", IloCP::FailLimit, 0, INT_MAX);
  CheckIntCPOption("logperiod", IloCP::LogPeriod, 1, INT_MAX);
  CheckIntCPOption("multipointnumberofsearchpoints",
      IloCP::MultiPointNumberOfSearchPoints, 2, INT_MAX);
  CheckDblCPOption("optimalitytolerance", IloCP::OptimalityTolerance, 42, -1);
  CheckIntCPOption("randomseed", IloCP::RandomSeed, 0, INT_MAX);
  CheckDblCPOption("relativeoptimalitytolerance",
      IloCP::RelativeOptimalityTolerance, 42, -1);
  CheckDblCPOption("restartgrowthfactor", IloCP::RestartGrowthFactor, 42, -1);
  CheckIntCPOption("restartfaillimit", IloCP::RestartFailLimit, 1, INT_MAX);
  CheckIntCPOption("solutionlimit", IloCP::SolutionLimit, 0, INT_MAX);
  CheckDblCPOption("timelimit", IloCP::TimeLimit, 42, -1);
  if (CPX_VERSION > 1220)
    CheckIntCPOption("workers", IloCP::Workers, 0, INT_MAX, 0, true);
  else CheckIntCPOption("workers", IloCP::Workers, 1, 4, 0, false);
}

TEST_F(IlogCPTest, CPLEXDefaultMIPDisplayZero) {
  EXPECT_TRUE(ParseOptions("optimizer=cplex"));
  CPLEXOptimizer *opt = dynamic_cast<CPLEXOptimizer*>(d.optimizer());
  ASSERT_TRUE(opt != nullptr);
  EXPECT_EQ(0, opt->cplex().getParam(IloCplex::MIPDisplay));
}

TEST_F(IlogCPTest, CPLEXOptions) {
  CheckIntCPLEXOption("mipdisplay", IloCplex::MIPDisplay, 0, 5);
  CheckIntCPLEXOption("mipinterval",
      IloCplex::MIPInterval, INT_MIN, INT_MAX);
}

// ----------------------------------------------------------------------------
// Solve code tests

TEST_F(IlogCPTest, OptimalSolveCode) {
  Solve(DATA_DIR "objconst");
  EXPECT_EQ(0, d.get_asl()->p.solve_code_);
}

TEST_F(IlogCPTest, FeasibleSolveCode) {
  Solve(DATA_DIR "feasible");
  EXPECT_EQ(100, d.get_asl()->p.solve_code_);
}

TEST_F(IlogCPTest, InfeasibleSolveCode) {
  Solve(DATA_DIR "infeasible");
  EXPECT_EQ(200, d.get_asl()->p.solve_code_);
}

TEST_F(IlogCPTest, InfeasibleOrUnboundedSolveCode) {
  Solve(DATA_DIR "unbounded");
  EXPECT_EQ(201, d.get_asl()->p.solve_code_);
}
}

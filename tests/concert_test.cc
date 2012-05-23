#include <ilconcert/ilomodel.h>

#include <algorithm>
#include <memory>
#include <sstream>

#include "gtest/gtest.h"

#include "solvers/concert/concert.h"
#include "solvers/asl.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"
#include "tests/config.h"

namespace {

bool AreBothSpaces(char lhs, char rhs) { return lhs == ' ' && rhs == ' '; }

// Returns a string representation of the argument.
template <typename T>
std::string str(T t) {
  std::ostringstream ss;
  ss << t;
  auto s = ss.str();

  // Replace adjacent duplicate spaces and possible trailing space.
  auto end = std::unique(s.begin(), s.end(), AreBothSpaces);
  if (*(end - 1) == ' ') --end;
  s.erase(end, s.end());

  return s;
}

// A functor for deleting ASL expressions recursively.
struct ExprDeleter {
  void operator()(expr *e) const;
};

void ExprDeleter::operator()(expr *e) const {
  if (!e) return;
  size_t op = reinterpret_cast<size_t>(e->op);
  switch (op) {
  case OPNUM:
    delete reinterpret_cast<expr_n*>(e);
    return;
  case MINLIST: case MAXLIST: {
    expr_va *eva = reinterpret_cast<expr_va*>(e);
    for (de *d = reinterpret_cast<expr_va*>(e)->L.d; d->e; ++d)
      (*this)(d->e);
    delete eva;
    return;
  }
  case OPIFnl: {
    expr_if *eif = reinterpret_cast<expr_if*>(e);
    (*this)(eif->e);
    (*this)(eif->T);
    (*this)(eif->F);
    delete eif;
    return;
  }
  }
  // Delete subexpressions recursively.
  (*this)(e->L.e);
  (*this)(e->R.e);
  delete e;
}

class ConcertTest : public ::testing::Test {
 protected:
  void SetUp() {
    env = IloEnv();
    mod = IloModel(env);
    Var = IloNumVarArray(env, 3);
    Var[0] = IloNumVar(env, 0, 1, "x");
    Var[1] = IloNumVar(env, 0, 1, "y");
    Var[2] = IloNumVar(env, 0, 1, "theta");
  }

  void TearDown() {
    Var = IloNumVarArray();
  }

#if HAVE_UNIQUE_PTR
  typedef std::unique_ptr<expr, ExprDeleter> ExprPtr;
#else
  // Fall back to std::auto_ptr - leaks subexpressions.
  typedef std::auto_ptr<expr> ExprPtr;
#endif

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

  // Creates a variable-argument ASL expression with 3 arguments.
  static ExprPtr NewExpr(int opcode, ExprPtr e1, ExprPtr e2, ExprPtr e3) {
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

  // Creates an if ASL expression.
  static ExprPtr NewIf(ExprPtr condition,
      ExprPtr true_expr, ExprPtr false_expr) {
    expr_if e = {reinterpret_cast<efunc*>(OPIFnl), 0, condition.release(),
                 true_expr.release(), false_expr.release(),
                 0, 0, 0, 0, {0}, {0}, 0, 0};
    return ExprPtr(reinterpret_cast<expr*>(new expr_if(e)));
  }
};

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
    NewExpr(MINLIST, NewVar(0), NewVar(1), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertMax) {
  EXPECT_EQ("max([x , y , 42 ])", str(build_expr(
    NewExpr(MAXLIST, NewVar(0), NewVar(1), NewNum(42)).get())));
}

TEST_F(ConcertTest, ConvertFloor) {
  EXPECT_EQ("floor(x )", str(build_expr(
    NewUnary(FLOOR, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertCeil) {
  EXPECT_EQ("ceil(x )", str(build_expr(
    NewUnary(CEIL, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertAbs) {
  EXPECT_EQ("abs(x )", str(build_expr(
    NewUnary(ABS, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertUMinus) {
  EXPECT_EQ("-1 * x", str(build_expr(
    NewUnary(OPUMINUS, NewVar(0)).get())));
}

TEST_F(ConcertTest, ConvertLogicalOrComparisonThrows) {
  int ops[] = {OPOR, OPAND, LT, LE, EQ, GE, GT, NE, OPNOT};
  size_t i = 0;
  for (size_t num_ops = sizeof(ops) / sizeof(*ops); i < num_ops; ++i) {
    EXPECT_THROW(build_expr(
      NewBinary(ops[i], NewVar(0), NewVar(1)).get()), Error);
  }
  // Paranoid: make sure that the loop body has been executed enough times.
  EXPECT_EQ(9u, i);
}

TEST_F(ConcertTest, ConvertIf) {
  build_expr(NewIf(NewBinary(EQ,
    NewVar(0), NewNum(0)), NewVar(1), NewNum(42)).get());

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
}


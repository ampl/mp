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
  if (reinterpret_cast<size_t>(e->op) == OPNUM) {
    delete reinterpret_cast<expr_n*>(e);
    return;
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

  // Creates a binary ASL expression.
  static ExprPtr NewBinary(int opcode, ExprPtr lhs, ExprPtr rhs) {
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.release()}, {rhs.release()}, 0};
    return ExprPtr(new expr(e));
  }

  // Creates a variable-argument ASL expression.
  static ExprPtr NewExpr(int opcode, ExprPtr lhs, ExprPtr rhs) {
    // TODO
    expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
              {lhs.release()}, {rhs.release()}, 0};
    return ExprPtr(new expr(e));
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
  build_expr(NewBinary(OPREM, NewVar(0), NewVar(1)).get());
  IloModel::Iterator iter(mod);

  ASSERT_TRUE(iter.ok());
  IloIfThenI *ifNonnegative = dynamic_cast<IloIfThenI*>((*iter).getImpl());
  ASSERT_TRUE(ifNonnegative != nullptr);
  EXPECT_EQ("0 <= x / y", str(ifNonnegative->getLeft()));
  EXPECT_EQ("IloNumVar(10)[-inf..inf] == floor(x / y )",
            str(ifNonnegative->getRight()));

  ++iter;
  ASSERT_TRUE(iter.ok());
  IloIfThenI *ifNegative = dynamic_cast<IloIfThenI*>((*iter).getImpl());
  ASSERT_TRUE(ifNegative != nullptr);
  IloNotI *ifNot = dynamic_cast<IloNotI*>(ifNegative->getLeft().getImpl());
  EXPECT_EQ("0 <= x / y", str(ifNot->getConstraint()));
  EXPECT_EQ("IloNumVar(10)[-inf..inf] == ceil(x / y )",
            str(ifNegative->getRight()));

  ++iter;
  EXPECT_FALSE(iter.ok());
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
}


#include "concert.h"

#include "gtest/gtest.h"
#include <ilconcert/ilomodel.h>
#include <algorithm>
#include <memory>
#include <sstream>

#include "asl.h"
#include "nlp.h"
#include "opcode.hd"

using namespace std;

namespace {

bool both_are_spaces(char lhs, char rhs) { return lhs == ' ' && rhs == ' '; }

// Returns a string representation of the specified expression.
string str(IloExpr e) {
  ostringstream ss;
  ss << e;
  string s = ss.str();

  // Replace adjacent duplicate spaces and possible trailing space.
  auto end = unique(s.begin(), s.end(), both_are_spaces);
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

typedef unique_ptr<expr, ExprDeleter> ExprPtr;

// Creates an ASL expression representing a number.
ExprPtr new_num(double n) {
  return ExprPtr(reinterpret_cast<expr*>(
    new expr_n{reinterpret_cast<efunc_n*>(OPNUM), n}));
}

ExprPtr new_var() {
  expr e = {reinterpret_cast<efunc*>(OPVARVAL), 0, 0, {0}, {0}, 0};
  return ExprPtr(new expr(e));
}

// Creates a binary ASL expression.
ExprPtr new_binary(int opcode, ExprPtr lhs, ExprPtr rhs) {
  expr e = {reinterpret_cast<efunc*>(opcode), 0, 0,
            {lhs.release()}, {rhs.release()}, 0};
  return ExprPtr(new expr(e));
}

TEST(CONCERTTest, ConvertNum) {
  EXPECT_EQ("0.42", str(build_expr(new_num(0.42).get())));
}

TEST(CONCERTTest, ConvertVar) {
  Var = IloNumVarArray(env, 1);
  Var[0] = IloNumVar(env, 0, 1, "theta");
  EXPECT_EQ("theta", str(build_expr(new_var().get())));
}

TEST(CONCERTTest, ConvertPlus) {
  Var = IloNumVarArray(env, 1);
  Var[0] = IloNumVar(env, 0, 1, "x");
  EXPECT_EQ("x + 42", str(build_expr(
    new_binary(OPPLUS, new_num(42), new_var()).get())));
}

}


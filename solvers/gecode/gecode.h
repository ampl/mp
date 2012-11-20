#ifndef AMPL_SOLVERS_GECODE_H
#define AMPL_SOLVERS_GECODE_H

#include <memory>
#include <gecode/minimodel.hh>

#include "solvers/util/driver.h"

struct Option_Info;
struct ASL_fg;

namespace ampl {

class GecodeProblem: public Gecode::Space,
  public ExprVisitor<GecodeProblem, Gecode::LinExpr, Gecode::BoolExpr> {
 public:
  typedef Gecode::LinExpr LinExpr;
  typedef Gecode::BoolExpr BoolExpr;

 private:
  Gecode::IntVarArray vars_;
  Gecode::IntVar obj_;
  Gecode::IntRelType obj_irt_; // IRT_NQ - no objective,
                               // IRT_LE - minimization, IRT_GR - maximization
  static const Gecode::BoolExpr DUMMY_EXPR;

  // Creates an integer variable to represent an expression.
  Gecode::IntVar CreateVar(LinExpr e) {
    Gecode::IntVar var(*this,
        Gecode::Int::Limits::min, Gecode::Int::Limits::max);
    rel(*this, var == e);
    return var;
  }

  // Creates a boolean variable to represent an expression.
  Gecode::BoolVar CreateVar(BoolExpr e) {
    Gecode::BoolVar var(*this, 0, 1);
    rel(*this, var == e);
    return var;
  }

  BoolExpr Convert(Gecode::BoolOpType op, IteratedLogicalExpr e);

  static void RequireNonzeroConstRHS(
      BinaryExpr e, const std::string &func_name);

 public:
  GecodeProblem(int num_vars) :
    vars_(*this, num_vars), obj_irt_(Gecode::IRT_NQ) {}

  GecodeProblem(bool share, GecodeProblem &s);

  Gecode::Space *copy(bool share);

  Gecode::IntVarArray &vars() { return vars_; }
  Gecode::IntVar &obj() { return obj_; }

  void SetObj(Driver::ObjType obj_type, const LinExpr &expr);

  virtual void constrain(const Gecode::Space &best);

  // The methods below perform conversion of AMPL expressions into
  // equivalent Gecode expressions. Gecode doesn't support the following
  // expressions/functions:
  // * division other than integer one
  // * trigonometric functions
  //   http://www.gecode.org/pipermail/users/2011-March/003177.html
  // * log, log10, exp, pow
  // * sqrt other than floor(sqrt())

  template <typename Grad>
  LinExpr ConvertExpr(Grad *grad, NumericExpr nonlinear);

  LinExpr VisitPlus(BinaryExpr e) {
    return Visit(e.lhs()) + Visit(e.rhs());
  }

  LinExpr VisitMinus(BinaryExpr e) {
    return Visit(e.lhs()) - Visit(e.rhs());
  }

  LinExpr VisitMult(BinaryExpr e) {
    return Visit(e.lhs()) * Visit(e.rhs());
  }

  LinExpr VisitRem(BinaryExpr e) {
    return Visit(e.lhs()) % Visit(e.rhs());
  }

  LinExpr VisitNumericLess(BinaryExpr e) {
    return max(Visit(e.lhs()) - Visit(e.rhs()), 0);
  }

  LinExpr VisitMin(VarArgExpr e);

  LinExpr VisitMax(VarArgExpr e);

  LinExpr VisitFloor(UnaryExpr e);

  LinExpr VisitCeil(UnaryExpr e) {
    // ceil does nothing because Gecode supports only integer expressions.
    return Visit(e.arg());
  }

  LinExpr VisitAbs(UnaryExpr e) {
    return abs(Visit(e.arg()));
  }

  LinExpr VisitUnaryMinus(UnaryExpr e) {
    return -Visit(e.arg());
  }

  LinExpr VisitIf(IfExpr e);

  LinExpr VisitSum(SumExpr e);

  LinExpr VisitIntDiv(BinaryExpr e) {
    return Visit(e.lhs()) / Visit(e.rhs());
  }

  LinExpr VisitRound(BinaryExpr e) {
    // round does nothing because Gecode supports only integer expressions.
    RequireNonzeroConstRHS(e, "round");
    return Visit(e.lhs());
  }

  LinExpr VisitTrunc(BinaryExpr e) {
    // trunc does nothing because Gecode supports only integer expressions.
    RequireNonzeroConstRHS(e, "trunc");
    return Visit(e.lhs());
  }

  LinExpr VisitCount(CountExpr e);

  LinExpr VisitNumberOf(NumberOfExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitPLTerm(PiecewiseLinearTerm t) {
    return VisitUnhandledNumericExpr(t); // TODO
  }

  LinExpr VisitConstExpPow(BinaryExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitPow2(UnaryExpr e) {
    return sqr(Visit(e.arg()));
  }

  LinExpr VisitConstBasePow(BinaryExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitNumericConstant(NumericConstant c) {
    double value = c.value();
    if (static_cast<int>(value) != value)
      throw UnsupportedExprError("non-integer constant");
    return c.value();
  }

  LinExpr VisitVariable(Variable v) {
    return vars_[v.index()];
  }

  BoolExpr VisitOr(BinaryLogicalExpr e) {
    return Visit(e.lhs()) || Visit(e.rhs());
  }

  BoolExpr VisitAnd(BinaryLogicalExpr e) {
    return Visit(e.lhs()) && Visit(e.rhs());
  }

  BoolExpr VisitLess(RelationalExpr e) {
    return Visit(e.lhs()) < Visit(e.rhs());
  }

  BoolExpr VisitLessEqual(RelationalExpr e) {
    return Visit(e.lhs()) <= Visit(e.rhs());
  }

  BoolExpr VisitEqual(RelationalExpr e) {
    return Visit(e.lhs()) == Visit(e.rhs());
  }

  BoolExpr VisitGreaterEqual(RelationalExpr e) {
    return Visit(e.lhs()) >= Visit(e.rhs());
  }

  BoolExpr VisitGreater(RelationalExpr e) {
    return Visit(e.lhs()) > Visit(e.rhs());
  }

  BoolExpr VisitNotEqual(RelationalExpr e) {
    return Visit(e.lhs()) != Visit(e.rhs());
  }

  BoolExpr VisitNot(NotExpr e) {
    return !Visit(e.arg());
  }

  BoolExpr VisitAtLeast(RelationalExpr e) {
    return Visit(e.lhs()) <= Visit(e.rhs());
  }

  BoolExpr VisitAtMost(RelationalExpr e) {
    return Visit(e.lhs()) >= Visit(e.rhs());
  }

  BoolExpr VisitExactly(RelationalExpr e) {
    return Visit(e.lhs()) == Visit(e.rhs());
  }

  BoolExpr VisitNotAtLeast(RelationalExpr e) {
    return Visit(e.lhs()) > Visit(e.rhs());
  }

  BoolExpr VisitNotAtMost(RelationalExpr e) {
    return Visit(e.lhs()) < Visit(e.rhs());
  }

  BoolExpr VisitNotExactly(RelationalExpr e) {
    return Visit(e.lhs()) != Visit(e.rhs());
  }

  BoolExpr VisitForAll(IteratedLogicalExpr e) {
    return Convert(Gecode::BOT_AND, e);
  }

  BoolExpr VisitExists(IteratedLogicalExpr e) {
    return Convert(Gecode::BOT_OR, e);
  }

  BoolExpr VisitImplication(ImplicationExpr e) {
    BoolExpr condition = Visit(e.condition());
    BoolExpr expr = condition >> Visit(e.true_expr());
    LogicalConstant c = Cast<LogicalConstant>(e.false_expr());
    if (c && !c.value())
      return expr;
    return expr || (condition && Visit(e.false_expr()));
  }

  BoolExpr VisitIff(BinaryLogicalExpr e) {
    return Visit(e.lhs()) == Visit(e.rhs());
  }

  BoolExpr VisitAllDiff(AllDiffExpr e);

  BoolExpr VisitLogicalConstant(LogicalConstant c) {
    return VisitUnhandledLogicalExpr(c); // TODO
  }
};

// The Gecode driver for AMPL.
class GecodeDriver : public Driver {
 private:
  std::auto_ptr<Option_Info> oinfo_;

 public:
  GecodeDriver();

  // Run the driver.
  int run(char **argv);
};

}

#endif // AMPL_SOLVERS_GECODE_H

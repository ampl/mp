#ifndef AMPL_SOLVERS_GECODE_H
#define AMPL_SOLVERS_GECODE_H

#include <memory>
#include <gecode/minimodel.hh>

#include "solvers/util/driver.h"

struct Option_Info;
struct ASL_fg;

namespace ampl {

class GecodeProblem: public Gecode::Space {
 private:
  Gecode::IntVarArray vars_;
  Gecode::IntVar obj_;
  Gecode::IntRelType obj_irt_; // IRT_NQ - no objective,
                               // IRT_LE - minimization, IRT_GR - maximization
 public:
  GecodeProblem(int num_vars) :
    vars_(*this, num_vars), obj_irt_(Gecode::IRT_NQ) {}

  GecodeProblem(bool share, GecodeProblem &s);

  Gecode::Space *copy(bool share);

  Gecode::IntVarArray &vars() { return vars_; }
  Gecode::IntVar &obj() { return obj_; }

  void SetObj(Problem::ObjType obj_type, const Gecode::LinExpr &expr);

  virtual void constrain(const Gecode::Space &best);
};

class NLToGecodeConverter :
  private ExprVisitor<NLToGecodeConverter, Gecode::LinExpr, Gecode::BoolExpr> {
 private:
   GecodeProblem problem_;

   friend class ExprVisitor;

   typedef Gecode::LinExpr LinExpr;
   typedef Gecode::BoolExpr BoolExpr;

   // Creates an integer variable to represent an expression.
   Gecode::IntVar CreateVar(LinExpr e) {
     Gecode::IntVar var(problem_,
         Gecode::Int::Limits::min, Gecode::Int::Limits::max);
     rel(problem_, var == e);
     return var;
   }

   // Creates a boolean variable to represent an expression.
   Gecode::BoolVar CreateVar(BoolExpr e) {
     Gecode::BoolVar var(problem_, 0, 1);
     rel(problem_, var == e);
     return var;
   }

   BoolExpr Convert(Gecode::BoolOpType op, IteratedLogicalExpr e);

   static void RequireNonzeroConstRHS(
       BinaryExpr e, const std::string &func_name);

   template <typename Grad>
   LinExpr ConvertExpr(Grad *grad, NumericExpr nonlinear);

   // The methods below perform conversion of AMPL NL expressions into
   // equivalent Gecode expressions. Gecode doesn't support the following
   // expressions/functions:
   // * division other than integer one
   // * trigonometric functions
   //   http://www.gecode.org/pipermail/users/2011-March/003177.html
   // * log, log10, exp, pow
   // * sqrt other than floor(sqrt())

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

   LinExpr VisitNumberOf(NumberOfExpr e);

   LinExpr VisitPow2(UnaryExpr e) {
     return sqr(Visit(e.arg()));
   }

   LinExpr VisitNumericConstant(NumericConstant c) {
     double value = c.value();
     if (static_cast<int>(value) != value)
       throw UnsupportedExprError("non-integer constant");
     return c.value();
   }

   LinExpr VisitVariable(Variable v) {
     return problem_.vars()[v.index()];
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

   BoolExpr VisitAtLeast(LogicalCountExpr e) {
     return Visit(e.value()) <= Visit(e.count());
   }

   BoolExpr VisitAtMost(LogicalCountExpr e) {
     return Visit(e.value()) >= Visit(e.count());
   }

   BoolExpr VisitExactly(LogicalCountExpr e) {
     return Visit(e.value()) == Visit(e.count());
   }

   BoolExpr VisitNotAtLeast(LogicalCountExpr e) {
     return Visit(e.value()) > Visit(e.count());
   }

   BoolExpr VisitNotAtMost(LogicalCountExpr e) {
     return Visit(e.value()) < Visit(e.count());
   }

   BoolExpr VisitNotExactly(LogicalCountExpr e) {
     return Visit(e.value()) != Visit(e.count());
   }

   BoolExpr VisitForAll(IteratedLogicalExpr e) {
     return Convert(Gecode::BOT_AND, e);
   }

   BoolExpr VisitExists(IteratedLogicalExpr e) {
     return Convert(Gecode::BOT_OR, e);
   }

   BoolExpr VisitImplication(ImplicationExpr e);

   BoolExpr VisitIff(BinaryLogicalExpr e) {
     return Visit(e.lhs()) == Visit(e.rhs());
   }

   BoolExpr VisitAllDiff(AllDiffExpr e);

   BoolExpr VisitLogicalConstant(LogicalConstant c) {
     bool value = c.value();
     return Gecode::BoolVar(problem_, value, value);
   }

 public:
  NLToGecodeConverter(int num_vars) : problem_(num_vars) {}

  LinExpr ConvertFullExpr(NumericExpr e) { return Visit(e); }
  BoolExpr ConvertFullExpr(LogicalExpr e, bool post = true);

  void Convert(const Problem &p);

  GecodeProblem &problem() { return problem_; }
};

// The Gecode driver for AMPL.
class GecodeDriver : public Driver {
 private:
  OptionInfo<GecodeDriver> options_;

 public:
  GecodeDriver();

  // Run the driver.
  int Run(char **argv);
};
}

#endif // AMPL_SOLVERS_GECODE_H

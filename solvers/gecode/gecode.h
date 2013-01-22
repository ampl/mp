/*
 AMPL solver interface to Gecode.

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

#ifndef AMPL_SOLVERS_GECODE_H
#define AMPL_SOLVERS_GECODE_H

#include <memory>
#include <gecode/minimodel.hh>
#include <gecode/search.hh>

#include "solvers/util/solver.h"

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

   BoolExpr Convert(Gecode::BoolOpType op, IteratedLogicalExpr e);

   static void RequireNonzeroConstRHS(
       BinaryExpr e, const std::string &func_name);

   template <typename Term>
   LinExpr ConvertExpr(LinearExpr<Term> linear, NumericExpr nonlinear);

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

template <typename Value>
struct OptionValue {
  const char *name;
  Value value;
};

template <typename T>
struct OptionInfo {
  const OptionValue<T> *values;
  T &value;

  OptionInfo(const OptionValue<T> *values, T &value)
  : values(values), value(value) {}
};

// Gecode solver.
class GecodeSolver : public Solver<GecodeSolver> {
 private:
  bool output_;
  Gecode::IntVarBranch var_branching_;
  Gecode::IntValBranch val_branching_;
  Gecode::Search::Options options_;
  double time_limit_; // Time limit in seconds.
  unsigned long node_limit_;
  unsigned long fail_limit_;
  std::size_t memory_limit_;

  void EnableOutput(const char *, bool enable) { output_ = enable; }

  template <typename T>
  void SetStrOption(const char *name, const char *value,
      const OptionInfo<T> &info);

  template <typename T>
  void SetIntOption(const char *name, int value, T *option);

  void SetDblOption(const char *, double value, double *option) {
    *option = value;
  }

  class Stop : public Gecode::Search::Stop {
   private:
    SignalHandler sh_;
    Gecode::Support::Timer timer_;
    const GecodeSolver &solver_;
    double time_limit_in_milliseconds_;
    bool has_limit_;

   public:
    explicit Stop(const GecodeSolver &s);

    bool stop(const Gecode::Search::Statistics &s,
              const Gecode::Search::Options &);
  };

 public:
  GecodeSolver();

  // Run the solver.
  int Run(char **argv);

  const Gecode::Search::Options &options() const { return options_; }
};
}

#endif // AMPL_SOLVERS_GECODE_H

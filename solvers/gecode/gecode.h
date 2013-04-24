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
#include <string>

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4200; disable: 4345; disable: 4800)
#endif
#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#ifdef _MSC_VER
# pragma warning(pop)
#endif

#include "solvers/util/solver.h"

namespace ampl {

typedef Gecode::LinIntExpr LinExpr;

class GecodeProblem: public Gecode::Space {
 private:
  Gecode::IntVarArray vars_;
  Gecode::IntVar obj_;
  Gecode::IntRelType obj_irt_; // IRT_NQ - no objective,
                               // IRT_LE - minimization, IRT_GR - maximization
  Gecode::IntConLevel icl_;

  Gecode::Space &space() { return *this; }

 public:
  GecodeProblem(int num_vars, Gecode::IntConLevel icl);
  GecodeProblem(bool share, GecodeProblem &s);

  Gecode::Space *copy(bool share);

  Gecode::IntVarArray &vars() { return vars_; }
  Gecode::IntVar &obj() { return obj_; }

  void SetObj(ObjType obj_type, const LinExpr &expr);

  virtual void constrain(const Gecode::Space &best);
};

// Converter of constraint programming problems from NL to Gecode format.
class NLToGecodeConverter :
  public ExprVisitor<NLToGecodeConverter, LinExpr, Gecode::BoolExpr> {
 private:
  GecodeProblem problem_;
  Gecode::IntConLevel icl_;

  typedef Gecode::BoolExpr BoolExpr;

  static int CastToInt(double value) {
    int int_value = static_cast<int>(value);
    if (int_value != value) {
      throw UnsupportedExprError::CreateFromMessage(
          fmt::Format("value {} can't be represented as int") << value);
    }
    return int_value;
  }

  BoolExpr Convert(Gecode::BoolOpType op, IteratedLogicalExpr e);

  typedef void (*VarArgFunc)(
      Gecode::Home, const Gecode::IntVarArgs &,
      Gecode::IntVar, Gecode::IntConLevel);

  LinExpr Convert(VarArgExpr e, VarArgFunc f);

  static void RequireNonzeroConstRHS(
      BinaryExpr e, const std::string &func_name);

  template<typename Term>
  LinExpr ConvertExpr(LinearExpr<Term> linear, NumericExpr nonlinear);

 public:
  NLToGecodeConverter(int num_vars, Gecode::IntConLevel icl)
  : problem_(num_vars, icl), icl_(icl) {}

  void Convert(const Problem &p);

  GecodeProblem &problem() { return problem_; }

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

  LinExpr VisitMin(VarArgExpr e) {
    return Convert(e, Gecode::min);
  }

  LinExpr VisitMax(VarArgExpr e) {
    return Convert(e, Gecode::max);
  }

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

  LinExpr VisitPowConstExp(BinaryExpr e) {
    return Gecode::pow(Visit(e.lhs()),
        CastToInt(Cast<NumericConstant>(e.rhs()).value()));
  }

  LinExpr VisitPow2(UnaryExpr e) {
    return sqr(Visit(e.arg()));
  }

  LinExpr VisitNumericConstant(NumericConstant c) {
    return CastToInt(c.value());
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
    return Visit(e.value()) <= VisitCount(e.count());
  }

  BoolExpr VisitAtMost(LogicalCountExpr e) {
    return Visit(e.value()) >= VisitCount(e.count());
  }

  BoolExpr VisitExactly(LogicalCountExpr e) {
    return Visit(e.value()) == VisitCount(e.count());
  }

  BoolExpr VisitNotAtLeast(LogicalCountExpr e) {
    return Visit(e.value()) > VisitCount(e.count());
  }

  BoolExpr VisitNotAtMost(LogicalCountExpr e) {
    return Visit(e.value()) < VisitCount(e.count());
  }

  BoolExpr VisitNotExactly(LogicalCountExpr e) {
    return Visit(e.value()) != VisitCount(e.count());
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
  double output_frequency_;
  unsigned output_count_;
  bool print_problem_;
  std::string header_;

  Gecode::IntConLevel icl_;
  Gecode::IntVarBranch::Select var_branching_;
  Gecode::IntValBranch val_branching_;
  double decay_;
  Gecode::Search::Options options_;
  double time_limit_;  // Time limit in seconds.
  unsigned long node_limit_;
  unsigned long fail_limit_;
  std::size_t memory_limit_;

  void SetBoolOption(const char *name, int value, bool *option);
  void SetOutputFrequency(const char *name, int value);

  template <typename T>
  void SetStrOption(const char *name, const char *value,
      const OptionInfo<T> &info);

  template <typename T, typename OptionT>
  void SetOption(const char *name, T value, OptionT *option);

  void SetDecay(const char *name, double value) {
    if (value > 0 && value <= 1)
      decay_ = value;
    else
      ReportInvalidOptionValue(name, value);
  }

  void SetDblOption(const char *, double value, double *option) {
    *option = value;
  }

  fmt::TempFormatter<fmt::Write> Output(fmt::StringRef format);

  class Stop : public Gecode::Search::Stop {
   private:
    SignalHandler sh_;
    Gecode::Support::Timer timer_;
    GecodeSolver &solver_;
    double time_limit_in_milliseconds_;
    double last_output_time_;
    bool output_or_limit_;

   public:
    explicit Stop(GecodeSolver &s);

    bool stop(const Gecode::Search::Statistics &s,
              const Gecode::Search::Options &);
  };

 protected:
  std::string GetOptionHeader();

 public:
  GecodeSolver();

  // Run the solver.
  int Run(char **argv);

  void Solve(Problem &p);

  Gecode::IntConLevel icl() const { return icl_; }
  Gecode::IntVarBranch::Select var_branching() const { return var_branching_; }
  Gecode::IntValBranch val_branching() const { return val_branching_; }
  const Gecode::Search::Options &options() const { return options_; }
  double decay() const { return decay_; }
};
}

#endif // AMPL_SOLVERS_GECODE_H

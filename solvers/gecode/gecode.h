/*
 AMPL solver interface to Gecode.

 Copyright (C) 2012 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#ifndef MP_SOLVERS_GECODE_H_
#define MP_SOLVERS_GECODE_H_

#include <memory>
#include <string>

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4200; disable: 4345; disable: 4800)
#endif
#include <gecode/driver.hh>
#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#ifdef _MSC_VER
# pragma warning(pop)
#endif

#include "mp/clock.h"
#include "mp/expr-visitor.h"
#include "mp/problem.h"
#include "mp/solver.h"

namespace mp {

typedef Gecode::LinIntExpr LinExpr;

class GecodeProblem: public Gecode::Space {
 private:
  Gecode::IntVarArray vars_;
  Gecode::IntVar obj_;
  Gecode::IntRelType obj_irt_; // IRT_NQ - no objective,
                               // IRT_LE - minimization, IRT_GR - maximization
  Gecode::IntPropLevel ipl_;

  Gecode::Space &space() { return *this; }

 public:
  GecodeProblem(int num_vars, Gecode::IntPropLevel ipl);
#if GECODE_VERSION_NUMBER > 600000
  GecodeProblem(GecodeProblem &s);
  Gecode::Space *copy();
#else
  GecodeProblem(bool share, GecodeProblem& s);
  Gecode::Space* copy(bool share);
#endif
  Gecode::IntVarArray &vars() { return vars_; }
  Gecode::IntVar &obj() { return obj_; }

  bool has_obj() const { return obj_irt_ != Gecode::IRT_NQ; }
  void SetObj(obj::Type obj_type, const LinExpr &expr);

  virtual void constrain(const Gecode::Space &best);
};

// Converter of constraint programming problems from MP to Gecode format.
class MPToGecodeConverter : public ExprVisitor<MPToGecodeConverter, LinExpr> {
 private:
   
  GecodeProblem problem_;
  Gecode::IntPropLevel ipl_;
  IntSuffix ipl_suffix_;
  std::vector<LinExpr> common_exprs_;

  typedef Gecode::BoolExpr BoolExpr;

  static int CastToInt(double value) {
    int int_value = static_cast<int>(value);
    if (int_value != value)
      throw Error("value {} can't be represented as int", value);
    return int_value;
  }

  BoolExpr Convert(Gecode::BoolOpType op, IteratedLogicalExpr e);

  typedef void (*VarArgFunc)(
      Gecode::Home, const Gecode::IntVarArgs &,
      Gecode::IntVar, Gecode::IntPropLevel);

  LinExpr Convert(IteratedExpr e, VarArgFunc f);

  static void RequireZeroRHS(BinaryExpr e, fmt::StringRef func_name);

  LinExpr ConvertExpr(const LinearExpr &linear, NumericExpr nonlinear);

  Gecode::IntPropLevel GetIPL(int con_index) const;

  class LogicalExprConverter :
      public ExprConverter<LogicalExprConverter, Gecode::BoolExpr> {
   private:
    MPToGecodeConverter &converter_;  // Main converter.

   public:
    explicit LogicalExprConverter(MPToGecodeConverter &c) : converter_(c) {}

    using ExprConverter<LogicalExprConverter, Gecode::BoolExpr>::Visit;

    LinExpr Visit(NumericExpr e) {
      return converter_.Visit(e);
    }

    BoolExpr VisitLogicalConstant(LogicalConstant c) {
      bool value = c.value();
      return Gecode::BoolVar(converter_.problem_, value, value);
    }

    BoolExpr VisitOr(BinaryLogicalExpr e) {
      return Visit(e.lhs()) || Visit(e.rhs());
    }

    BoolExpr VisitAnd(BinaryLogicalExpr e) {
      return Visit(e.lhs()) && Visit(e.rhs());
    }

    BoolExpr VisitLT(RelationalExpr e) {
      return Visit(e.lhs()) < Visit(e.rhs());
    }

    BoolExpr VisitLE(RelationalExpr e) {
      return Visit(e.lhs()) <= Visit(e.rhs());
    }

    BoolExpr VisitEQ(RelationalExpr e) {
      return Visit(e.lhs()) == Visit(e.rhs());
    }

    BoolExpr VisitGE(RelationalExpr e) {
      return Visit(e.lhs()) >= Visit(e.rhs());
    }

    BoolExpr VisitGT(RelationalExpr e) {
      return Visit(e.lhs()) > Visit(e.rhs());
    }

    BoolExpr VisitNE(RelationalExpr e) {
      return Visit(e.lhs()) != Visit(e.rhs());
    }

    BoolExpr VisitNot(NotExpr e) {
      return !Visit(e.arg());
    }

    BoolExpr VisitForAll(IteratedLogicalExpr e) {
      return converter_.Convert(Gecode::BOT_AND, e);
    }

    BoolExpr VisitExists(IteratedLogicalExpr e) {
      return converter_.Convert(Gecode::BOT_OR, e);
    }

    BoolExpr VisitImplication(ImplicationExpr e);

    BoolExpr VisitIff(BinaryLogicalExpr e) {
      return Visit(e.lhs()) == Visit(e.rhs());
    }

    BoolExpr VisitAllDiff(PairwiseExpr e);
    BoolExpr VisitNotAllDiff(PairwiseExpr e) { return VisitAllDiff(e); }
  };

  using ExprVisitor<MPToGecodeConverter, LinExpr>::Visit;

  BoolExpr Visit(LogicalExpr e) {
    return LogicalExprConverter(*this).Visit(e);
  }

 public:
  MPToGecodeConverter(int num_vars, Gecode::IntPropLevel ipl)
  : problem_(num_vars, ipl), ipl_(ipl) {}

  void Convert(const Problem &p);

  GecodeProblem &problem() { return problem_; }

  // The methods below perform conversion of AMPL NL expressions into
  // equivalent Gecode expressions. Gecode doesn't support the following
  // expressions/functions:
  // * division other than the integer one
  // * trigonometric & hyperbolic functions
  //   http://www.gecode.org/pipermail/users/2011-March/003177.html
  // * log, log10, exp, pow
  // * sqrt other than floor(sqrt())

  LinExpr VisitAdd(BinaryExpr e) {
    return Visit(e.lhs()) + Visit(e.rhs());
  }

  LinExpr VisitSub(BinaryExpr e) {
    return Visit(e.lhs()) - Visit(e.rhs());
  }

  LinExpr VisitMul(BinaryExpr e) {
    return Visit(e.lhs()) * Visit(e.rhs());
  }

  LinExpr VisitMod(BinaryExpr e) {
    return Visit(e.lhs()) % Visit(e.rhs());
  }

  LinExpr VisitLess(BinaryExpr e) {
    return max(Visit(e.lhs()) - Visit(e.rhs()), 0);
  }

  LinExpr VisitMin(IteratedExpr e) {
    return Convert(e, Gecode::min);
  }

  LinExpr VisitMax(IteratedExpr e) {
    return Convert(e, Gecode::max);
  }

  LinExpr VisitMinus(UnaryExpr e) {
    return -Visit(e.arg());
  }

  LinExpr VisitAbs(UnaryExpr e) {
    return abs(Visit(e.arg()));
  }

  LinExpr VisitFloor(UnaryExpr e);

  LinExpr VisitCeil(UnaryExpr e) {
    // ceil does nothing because Gecode supports only integer expressions.
    return Visit(e.arg());
  }

  LinExpr VisitIf(IfExpr e);

  LinExpr VisitSum(IteratedExpr e);

  LinExpr VisitTruncDiv(BinaryExpr e) {
    return Visit(e.lhs()) / Visit(e.rhs());
  }

  LinExpr VisitRound(BinaryExpr e) {
    // round does nothing because Gecode supports only integer expressions.
    RequireZeroRHS(e, "round");
    return Visit(e.lhs());
  }

  LinExpr VisitTrunc(BinaryExpr e) {
    // trunc does nothing because Gecode supports only integer expressions.
    RequireZeroRHS(e, "trunc");
    return Visit(e.lhs());
  }

  LinExpr VisitCount(CountExpr e);

  LinExpr VisitNumberOf(IteratedExpr e);

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

  LinExpr VisitVariable(Reference v) {
    return problem_.vars()[v.index()];
  }

  LinExpr VisitCommonExpr(Reference v) {
    return common_exprs_[v.index()];
  }
};

// Gecode solver.
class GecodeSolver : public SolverImpl<Problem> {
 private:
  bool output_;
  double output_frequency_;
  unsigned output_count_;
  std::string header_;
  int solve_code_;
  std::string status_;

  Gecode::IntPropLevel ipl_;
  Gecode::IntVarBranch::Select var_branching_;
  Gecode::IntValBranch::Select val_branching_;
  double decay_;
  Gecode::Search::Options options_;
  double time_limit_;  // Time limit in seconds.
  unsigned long node_limit_;
  unsigned long fail_limit_;
  unsigned solution_limit_;

  Gecode::RestartMode restart_;
  double restart_base_;
  unsigned long restart_scale_;

  void SetBoolOption(const SolverOption &opt, int value, bool *ptr);

  double GetOutputFrequency(const SolverOption &) const {
    return output_frequency_;
  }
  void SetOutputFrequency(const SolverOption &opt, double value);

  template <typename T>
  std::string GetEnumOption(const SolverOption &opt, T *ptr) const;

  template <typename T>
  void SetEnumOption(const SolverOption &opt, fmt::StringRef value, T *ptr);

  template <typename T, typename OptionT>
  T GetOption(const SolverOption &, OptionT *option) const {
    return static_cast<T>(*option);
  }

  template <typename T, typename OptionT>
  void SetNonnegativeOption(const SolverOption &opt, T value, OptionT *option);

  double GetDecay(const SolverOption &) const { return decay_; }
  void SetDecay(const SolverOption &opt, double value) {
    if (value <= 0 || value > 1)
      throw InvalidOptionValue(opt, value);
    decay_ = value;
  }

  void DoSetDblOption(const SolverOption &, double value, double *option) {
    *option = value;
  }

  void SetStatus(int solve_code, const char *status) {
    solve_code_ = solve_code;
    status_ = status;
  }

  void Output(fmt::CStringRef format, const fmt::ArgList &args);
  FMT_VARIADIC(void, Output, fmt::CStringRef)

  class Stop : public Gecode::Search::Stop {
   private:
    GecodeSolver &solver_;
    steady_clock::time_point end_time_;
    steady_clock::time_point next_output_time_;
    bool output_or_limit_;

    steady_clock::duration GetOutputInterval() const {
      return steady_clock::duration(
            static_cast<steady_clock::rep>(solver_.output_frequency_ *
                steady_clock::period::den / steady_clock::period::num));
    }

   public:
    explicit Stop(GecodeSolver &s);

    bool stop(const Gecode::Search::Statistics &s,
              const Gecode::Search::Options &);
  };

#ifdef MP_USE_UNIQUE_PTR
  typedef std::unique_ptr<GecodeProblem> ProblemPtr;
#else
  typedef std::auto_ptr<GecodeProblem> ProblemPtr;
#endif

  template<template<typename, template<typename> class> class Meta>
  ProblemPtr Search(Problem &p, GecodeProblem &gecode_problem,
                    Gecode::Search::Statistics &stats, SolutionHandler &sh);

 public:
  GecodeSolver();

  Gecode::IntPropLevel ipl() const { return ipl_; }
  Gecode::IntVarBranch::Select var_branching() const { return var_branching_; }
  Gecode::IntValBranch val_branching() const { return val_branching_; }
  const Gecode::Search::Options &options() const { return options_; }
  double decay() const { return decay_; }

  Gecode::RestartMode restart() const { return restart_; }
  double restart_base() const { return restart_base_; }
  unsigned long restart_scale() const { return restart_scale_; }

  void Solve(Problem &p, SolutionHandler &sh);
};
}

#endif  // MP_SOLVERS_GECODE_H_

/*
 IBM/ILOG CP solver for AMPL.

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

#ifndef SOLVERS_ILOGCP_ILOGCP_H_
#define SOLVERS_ILOGCP_ILOGCP_H_

#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#include <string.h> /* This and -fpermissive seem to be needed for MacOSX, */
                    /* at least with g++ 4.6.  Otherwise there are errors */
                    /* with iloconcert/iloenv.h . */
#include <limits.h> /* Needed for g++ -m32 on MacOSX. */
#include <memory>
#include <string>
#include <vector>

#include "solvers/util/solver.h"
#include "solvers/util/format.h"

namespace ampl {

class Optimizer : public Interruptible {
 private:
  IloRangeArray cons_;

 public:
  Optimizer(IloEnv env, const Problem &p);
  virtual ~Optimizer();

  IloRangeArray cons() const { return cons_; }

  virtual IloAlgorithm algorithm() const = 0;

  virtual void StartSearch() = 0;
  virtual void EndSearch() = 0;
  virtual bool FindNextSolution() = 0;

  virtual void GetSolutionInfo(fmt::Formatter &format_message,
      std::vector<double> &dual_values) const = 0;
};

class CPLEXOptimizer : public Optimizer {
 private:
  IloCplex cplex_;
  IloCplex::Aborter aborter_;
  bool started_;

 public:
  CPLEXOptimizer(IloEnv env, const Problem &p);

  IloCplex cplex() const { return cplex_; }
  IloAlgorithm algorithm() const { return cplex_; }

  void StartSearch() { started_ = true; }
  void EndSearch() { started_ = false; }

  bool FindNextSolution();

  void GetSolutionInfo(fmt::Formatter &format_message,
      std::vector<double> &dual_values) const;

  void Interrupt() { aborter_.abort(); }
};

class CPOptimizer : public Optimizer {
 private:
  IloCP solver_;

 public:
  CPOptimizer(IloEnv env, const Problem &p);

  IloCP solver() const { return solver_; }
  IloAlgorithm algorithm() const { return solver_; }

  void StartSearch() { solver_.startNewSearch(); }
  void EndSearch() { solver_.endSearch(); }
  bool FindNextSolution() { return solver_.next(); }

  void GetSolutionInfo(fmt::Formatter &format_message,
      std::vector<double> &dual_values) const;

  void Interrupt() { solver_.abortSearch(); }
};

class NLToConcertConverter;

typedef ExprVisitor<NLToConcertConverter, IloExpr, IloConstraint> Visitor;

// Converter of optimization problems from NL to Concert format.
class NLToConcertConverter : public Visitor {
 private:
  IloEnv env_;
  IloModel model_;
  IloNumVarArray vars_;
  bool use_numberof_;
  bool debug_;

  class CreateVar {
   private:
    IloEnv env_;

   public:
    CreateVar(IloEnv env) : env_(env) {}

    IloIntVar operator()() const {
      return IloIntVar(env_, IloIntMin, IloIntMax);
    }
  };

  typedef NumberOfMap<IloIntVar, CreateVar> IlogNumberOfMap;
  IlogNumberOfMap numberofs_;

  IloNumExprArray ConvertArgs(VarArgExpr e);

 public:
  NLToConcertConverter(
      IloEnv env, IloNumVarArray vars, bool use_numberof, bool debug);

  IloModel model() const { return model_; }

  IloExpr Visit(NumericExpr e) {
    if (debug_)
      printf("%s\n", e.opname());
    return Visitor::Visit(e);
  }

  IloConstraint Visit(LogicalExpr e) {
    if (debug_)
      printf("%s\n", e.opname());
    return Visitor::Visit(e);
  }

  IloExpr VisitPlus(BinaryExpr e) {
    return Visit(e.lhs()) + Visit(e.rhs());
  }

  IloExpr VisitMinus(BinaryExpr e) {
    return Visit(e.lhs()) - Visit(e.rhs());
  }

  IloExpr VisitMult(BinaryExpr e) {
    return Visit(e.lhs()) * Visit(e.rhs());
  }

  IloExpr VisitDiv(BinaryExpr e) {
    return Visit(e.lhs()) / Visit(e.rhs());
  }

  IloExpr VisitRem(BinaryExpr e) {
    IloNumExpr lhs(Visit(e.lhs())), rhs(Visit(e.rhs()));
    return lhs - IloTrunc(lhs / rhs) * rhs;
  }

  IloExpr VisitPow(BinaryExpr e) {
    return IloPower(Visit(e.lhs()), Visit(e.rhs()));
  }

  IloExpr VisitNumericLess(BinaryExpr e) {
    return IloMax(Visit(e.lhs()) - Visit(e.rhs()), 0.0);
  }

  IloExpr VisitMin(VarArgExpr e) {
    return IloMin(ConvertArgs(e));
  }

  IloExpr VisitMax(VarArgExpr e) {
    return IloMax(ConvertArgs(e));
  }

  IloExpr VisitFloor(UnaryExpr e) {
    return IloFloor(Visit(e.arg()));
  }

  IloExpr VisitCeil(UnaryExpr e) {
    return IloCeil(Visit(e.arg()));
  }

  IloExpr VisitAbs(UnaryExpr e) {
    return IloAbs(Visit(e.arg()));
  }

  IloExpr VisitUnaryMinus(UnaryExpr e) {
    return -Visit(e.arg());
  }

  IloExpr VisitIf(IfExpr e);

  IloExpr VisitTanh(UnaryExpr e) {
    IloNumExpr exp(IloExponent(2 * Visit(e.arg())));
    return (exp - 1) / (exp + 1);
  }

  IloExpr VisitTan(UnaryExpr e) {
    return IloTan(Visit(e.arg()));
  }

  IloExpr VisitSqrt(UnaryExpr e) {
    return IloPower(Visit(e.arg()), 0.5);
  }

  IloExpr VisitSinh(UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return (IloExponent(arg) - IloExponent(-arg)) * 0.5;
  }

  IloExpr VisitSin(UnaryExpr e) {
    return IloSin(Visit(e.arg()));
  }

  IloExpr VisitLog10(UnaryExpr e) {
    return IloLog10(Visit(e.arg()));
  }

  IloExpr VisitLog(UnaryExpr e) {
    return IloLog(Visit(e.arg()));
  }

  IloExpr VisitExp(UnaryExpr e) {
    return IloExponent(Visit(e.arg()));
  }

  IloExpr VisitCosh(UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return (IloExponent(arg) + IloExponent(-arg)) * 0.5;
  }

  IloExpr VisitCos(UnaryExpr e) {
    return IloCos(Visit(e.arg()));
  }

  IloExpr VisitAtanh(UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return IloLog((1 + arg) / (1 - arg)) * 0.5;
  }

  IloExpr VisitAtan2(BinaryExpr e);

  IloExpr VisitAtan(UnaryExpr e) {
    return IloArcTan(Visit(e.arg()));
  }

  IloExpr VisitAsinh(UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return IloLog(arg + IloPower(IloSquare(arg) + 1, 0.5));
  }

  IloExpr VisitAsin(UnaryExpr e) {
    return IloArcSin(Visit(e.arg()));
  }

  IloExpr VisitAcosh(UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return IloLog(arg + IloPower(arg + 1, 0.5) * IloPower(arg - 1, 0.5));
  }

  IloExpr VisitAcos(UnaryExpr e) {
    return IloArcCos(Visit(e.arg()));
  }

  IloExpr VisitSum(SumExpr e);

  IloExpr VisitIntDiv(BinaryExpr e) {
    return IloTrunc(Visit(e.lhs()) / Visit(e.rhs()));
  }

  IloExpr VisitRound(BinaryExpr e);

  IloExpr VisitTrunc(BinaryExpr e);

  IloExpr VisitCount(CountExpr e);

  IloExpr VisitNumberOf(NumberOfExpr e);

  IloExpr VisitPowConstExp(BinaryExpr e) {
    return IloPower(Visit(e.lhs()), Cast<NumericConstant>(e.rhs()).value());
  }

  IloExpr VisitPow2(UnaryExpr e) {
    return IloSquare(Visit(e.arg()));
  }

  IloExpr VisitPowConstBase(BinaryExpr e) {
    return IloPower(Cast<NumericConstant>(e.lhs()).value(), Visit(e.rhs()));
  }

  IloExpr VisitPLTerm(PiecewiseLinearTerm t);

  IloExpr VisitNumericConstant(NumericConstant n) {
    return IloExpr(env_, n.value());
  }

  IloExpr VisitVariable(Variable v) {
    return vars_[v.index()];
  }

  IloConstraint VisitLogicalConstant(LogicalConstant c) {
    return IloNumVar(env_, 1, 1) == c.value();
  }

  IloConstraint VisitLess(RelationalExpr e) {
    return Visit(e.lhs()) < Visit(e.rhs());
  }

  IloConstraint VisitLessEqual(RelationalExpr e) {
    return Visit(e.lhs()) <= Visit(e.rhs());
  }

  IloConstraint VisitEqual(RelationalExpr e) {
    return Visit(e.lhs()) == Visit(e.rhs());
  }

  IloConstraint VisitGreaterEqual(RelationalExpr e) {
    return Visit(e.lhs()) >= Visit(e.rhs());
  }

  IloConstraint VisitGreater(RelationalExpr e) {
    return Visit(e.lhs()) > Visit(e.rhs());
  }

  IloConstraint VisitNotEqual(RelationalExpr e) {
    return Visit(e.lhs()) != Visit(e.rhs());
  }

  IloConstraint VisitAtMost(LogicalCountExpr e) {
    return Visit(e.value()) >= VisitCount(e.count());
  }

  IloConstraint VisitNotAtMost(LogicalCountExpr e) {
    return !(Visit(e.value()) >= VisitCount(e.count()));
  }

  IloConstraint VisitAtLeast(LogicalCountExpr e) {
    return Visit(e.value()) <= VisitCount(e.count());
  }

  IloConstraint VisitNotAtLeast(LogicalCountExpr e) {
    return !(Visit(e.value()) <= VisitCount(e.count()));
  }

  IloConstraint VisitExactly(LogicalCountExpr e) {
    return Visit(e.value()) == VisitCount(e.count());
  }

  IloConstraint VisitNotExactly(LogicalCountExpr e) {
    return Visit(e.value()) != VisitCount(e.count());
  }

  IloConstraint VisitOr(BinaryLogicalExpr e) {
    return IloIfThen(env_, !Visit(e.lhs()), Visit(e.rhs()));
  }

  IloConstraint VisitExists(IteratedLogicalExpr e);

  IloConstraint VisitAnd(BinaryLogicalExpr e) {
    return Visit(e.lhs()) && Visit(e.rhs());
  }

  IloConstraint VisitForAll(IteratedLogicalExpr e);

  IloConstraint VisitNot(NotExpr e) {
    return !Visit(e.arg());
  }

  IloConstraint VisitIff(BinaryLogicalExpr e) {
    return Visit(e.lhs()) == Visit(e.rhs());
  }

  IloConstraint VisitImplication(ImplicationExpr e);

  IloConstraint VisitAllDiff(AllDiffExpr e);

  // Combines 'numberof' operators into IloDistribute constraints
  // which are much more useful to the solution procedure.
  void FinishBuildingNumberOf();
};

// IlogCP solver.
class IlogCPSolver : public Solver<IlogCPSolver> {
 private:
  IloEnv env_;
  std::auto_ptr<Optimizer> optimizer_;
  bool gotopttype_;

  // Do not implement.
  IlogCPSolver(const IlogCPSolver&);
  IlogCPSolver &operator=(const IlogCPSolver&);

 public:
  // Options accessible from AMPL.
  enum Option {
    DEBUGEXPR,
    OPTIMIZER,
    TIMING,
    USENUMBEROF,
    NUM_OPTIONS
  };

  // Values for the OPTIMIZER option.
  enum {
    AUTO  = -1,
    CP    =  0,
    CPLEX =  1
  };

 private:
  int options_[NUM_OPTIONS];

  std::string GetOptimizer(const char *name);
  void SetOptimizer(const char *name, const char *value);

  void SetBoolOption(const char *name, int value, Option opt);

  // Information about a constraint programming solver option.
  struct CPOptionInfo {
    IloCP::IntParam param;
    int start;           // start value for the enumerated options
    const char **values; // string values for enum options
    bool accepts_auto;   // true if the option accepts IloCP::Auto value

    CPOptionInfo(IloCP::IntParam p, int start = 0,
        const char **values = 0, bool accepts_auto = false)
    : param(p), start(start), values(values), accepts_auto(accepts_auto) {}
  };

  // Sets an integer option of the constraint programming optimizer.
  void SetCPOption(
      const char *name, const char *value, const CPOptionInfo &info);

  // Sets a double option of the constraint programming optimizer.
  void SetCPDblOption(const char *name, double value, IloCP::NumParam param);

  // Sets an integer option of the CPLEX optimizer.
  void SetCPLEXIntOption(const char *name, int value, int param);

  void CreateOptimizer(const Problem &p);

 protected:

  std::string GetOptionHeader();

 public:
  IlogCPSolver();
  virtual ~IlogCPSolver();

  IloEnv env() const { return env_; }

  IloAlgorithm alg() const {
    return optimizer_.get() ? optimizer_->algorithm() : IloAlgorithm();
  }

  Optimizer *optimizer() const { return optimizer_.get(); }

  bool ParseOptions(char **argv, unsigned flags = 0);

  int GetOption(Option opt) const {
    assert(opt >= 0 && opt < NUM_OPTIONS);
    return options_[opt];
  }

  void use_numberof(bool use = true) { options_[USENUMBEROF] = use; }

  void Solve(Problem &p);
};
}

#endif  // SOLVERS_ILOGCP_ILOGCP_H_

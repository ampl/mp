#ifndef AMPL_SOLVERS_CONCERT_H
#define AMPL_SOLVERS_CONCERT_H

#include <string.h>	/* This and -fpermissive seem to be needed for MacOSX, */
                    /* at least with g++ 4.6.  Otherwise there are errors */
                    /* with iloconcert/iloenv.h . */
#include <limits.h>	/* Needed for g++ -m32 on MacOSX. */
#include <memory>
#include <vector>
#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplex.h>
#include <ilcp/cp.h>

#include "solvers/util/expr.h"

struct keyword;
struct Option_Info;
struct ASL_fg;

namespace ampl {

class NumberOf {
 private:
  IloIntVarArray cards_;
  IloIntArray values_;
  IloIntVarArray vars_;
  NumberOfExpr expr_;

 public:
  NumberOf(IloIntVarArray cards, IloIntArray values,
      IloIntVarArray vars, NumberOfExpr e) :
    cards_(cards), values_(values), vars_(vars), expr_(e) {}

  IloInt num_vars() const { return vars_.getSize(); }
  NumberOfExpr expr() const { return expr_; }

  IloDistribute Convert(IloEnv env) const {
    return IloDistribute(env, cards_, values_, vars_);
  }

  IloIntVar Add(double value, IloEnv env);
};

class Optimizer {
 private:
  IloObjective obj_;
  IloNumVarArray vars_;
  IloRangeArray cons_;

 public:
  Optimizer(IloEnv env, ASL_fg *asl);
  virtual ~Optimizer();

  IloObjective obj() const { return obj_; }
  void set_obj(IloObjective obj) { obj_ = obj; }

  IloNumVarArray vars() const { return vars_; }
  IloRangeArray cons() const { return cons_; }

  virtual IloAlgorithm algorithm() const = 0;

  virtual void set_option(const void *key, int value) = 0;
  virtual void set_option(const void *key, double value) = 0;

  virtual void get_solution(ASL_fg *asl, char *message,
      std::vector<double> &primal, std::vector<double> &dual) const = 0;
};

class CPLEXOptimizer : public Optimizer {
 private:
  IloCplex cplex_;

 public:
  CPLEXOptimizer(IloEnv env, ASL_fg *asl) : Optimizer(env, asl), cplex_(env) {
    cplex_.setParam(IloCplex::MIPDisplay, 0);
  }

  IloCplex cplex() const { return cplex_; }
  IloAlgorithm algorithm() const { return cplex_; }

  void set_option(const void *key, int value);
  void set_option(const void *key, double value);

  void get_solution(ASL_fg *asl, char *message,
      std::vector<double> &primal, std::vector<double> &dual) const;
};

class CPOptimizer : public Optimizer {
 private:
  IloSolver solver_;

 public:
  CPOptimizer(IloEnv env, ASL_fg *asl) : Optimizer(env, asl), solver_(env) {
    solver_.setIntParameter(IloCP::LogVerbosity, IloCP::Quiet);
  }

  IloSolver solver() const { return solver_; }
  IloAlgorithm algorithm() const { return solver_; }

  void set_option(const void *key, int value);
  void set_option(const void *key, double value);

  void get_solution(ASL_fg *asl, char *message,
      std::vector<double> &primal, std::vector<double> &dual) const;
};

// The Ilogcp driver for AMPL.
class Driver : public ExprVisitor<Driver, IloExpr, IloConstraint> {
 private:
  IloEnv env_;
  IloModel mod_;
  IloNumVarArray vars_;
  std::auto_ptr<Optimizer> optimizer_;
  std::vector<NumberOf> numberofs_;
  ASL_fg *asl;
  std::vector<char> version_;
  std::auto_ptr<Option_Info> oinfo_;
  bool gotopttype;
  int n_badvals;
  static keyword keywords_[];

  // Do not implement.
  Driver(const Driver&);
  Driver &operator=(const Driver&);

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

  static char *set_optimizer(Option_Info *oi, keyword *kw, char *value);
  static char *set_int_option(Option_Info *oi, keyword *kw, char *value);
  static char *set_bool_option(Option_Info *oi, keyword *kw, char *value);

  void set_option(keyword *kw, int value);

  // Sets an integer option of the constraint programming optimizer.
  static char *set_cp_int_option(Option_Info *oi, keyword *kw, char *value);

  // Sets a double option of the constraint programming optimizer.
  static char *set_cp_dbl_option(Option_Info *oi, keyword *kw, char *value);

  // Sets an integer option of the CPLEX optimizer.
  static char *set_cplex_int_option(Option_Info *oi, keyword *kw, char *value);

  IloNumExprArray ConvertArgs(VarArgExpr e);
  
 public:
  Driver();
  virtual ~Driver();

  IloEnv env() const { return env_; }
  IloModel mod() const { return mod_; }
  ASL_fg *get_asl() const { return asl; }

  IloAlgorithm alg() const {
    return optimizer_.get() ? optimizer_->algorithm() : IloAlgorithm();
  }

  Optimizer *optimizer() const { return optimizer_.get(); }

  IloNumVarArray vars() const { return vars_; }
  void set_vars(IloNumVarArray vars) { vars_ = vars; }

  // Get and process ILOG Concert and driver options.
  bool parse_options(char **argv);

  int get_option(Option opt) const { return options_[opt]; }
  void use_numberof(bool use = true) { options_[USENUMBEROF] = use; }
  bool show_version() const;
  int wantsol() const;

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

  IloExpr VisitLess(BinaryExpr e) {
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

  IloExpr VisitMinus(UnaryExpr e) {
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
    return (IloLog(1 + arg) - IloLog(1 - arg)) * 0.5;
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

  IloExpr VisitConstExpPow(BinaryExpr e) {
    return IloPower(Visit(e.lhs()), Cast<NumericConstant>(e.rhs()).value());
  }

  IloExpr VisitPow2(UnaryExpr e) {
    return IloSquare(Visit(e.arg()));
  }

  IloExpr VisitConstBasePow(BinaryExpr e) {
    return IloPower(Cast<NumericConstant>(e.lhs()).value(), Visit(e.rhs()));
  }

  IloExpr VisitPLTerm(PiecewiseLinearTerm t);

  IloExpr VisitNumber(NumericConstant n) {
    return IloExpr(env_, n.value());
  }

  IloExpr VisitVariable(Variable v) {
    return vars_[v.index()];
  }

  IloConstraint VisitConstant(LogicalConstant c) {
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

  IloConstraint VisitAtMost(RelationalExpr e) {
    return Visit(e.lhs()) >= Visit(e.rhs());
  }

  IloConstraint VisitNotAtMost(RelationalExpr e) {
    return !(Visit(e.lhs()) >= Visit(e.rhs()));
  }

  IloConstraint VisitAtLeast(RelationalExpr e) {
    return Visit(e.lhs()) <= Visit(e.rhs());
  }

  IloConstraint VisitNotAtLeast(RelationalExpr e) {
    return !(Visit(e.lhs()) <= Visit(e.rhs()));
  }

  IloConstraint VisitExactly(RelationalExpr e) {
    return Visit(e.lhs()) == Visit(e.rhs());
  }

  IloConstraint VisitNotExactly(RelationalExpr e) {
    return Visit(e.lhs()) != Visit(e.rhs());
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

  // Runs the driver.
  int run(char **argv);
};
}

#endif // AMPL_SOLVERS_CONCERT_H

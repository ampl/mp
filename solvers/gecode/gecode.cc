#include "gecode.h"

#include <gecode/int.hh>
#include <gecode/minimodel.hh>
#include <gecode/search.hh>
#include <gecode/gist.hh>

#include <iostream>
#include <memory>
#include <ctime>

#include "solvers/util/expr.h"
#include "solvers/getstub.h"
#include "solvers/nlp.h"
#include "solvers/opcode.hd"

using std::cerr;
using std::endl;
using std::vector;

using ampl::Expr;
using ampl::NumericExpr;
using ampl::LogicalExpr;
using ampl::UnaryExpr;
using ampl::BinaryExpr;
using ampl::VarArgExpr;
using ampl::SumExpr;
using ampl::CountExpr;
using ampl::IfExpr;
using ampl::PiecewiseLinearTerm;
using ampl::NumericConstant;
using ampl::Variable;
using ampl::NumberOfExpr;
using ampl::LogicalConstant;
using ampl::RelationalExpr;
using ampl::NotExpr;
using ampl::BinaryLogicalExpr;
using ampl::ImplicationExpr;
using ampl::IteratedLogicalExpr;
using ampl::AllDiffExpr;
using ampl::Driver;

using Gecode::BoolExpr;
using Gecode::BAB;
using Gecode::DFS;
using Gecode::IntArgs;
using Gecode::IntVarArgs;
using Gecode::IntVar;
using Gecode::IntVarArray;
using Gecode::IntRelType;
using Gecode::IRT_NQ;
using Gecode::LinExpr;
using Gecode::linear;
using Gecode::Space;

namespace {

class GecodeProblem;
typedef ampl::ExprVisitor<GecodeProblem, LinExpr, BoolExpr> Visitor;

class GecodeProblem: public Space, public Visitor {
 private:
  IntVarArray vars_;
  IntVar obj_;
  IntRelType obj_irt_; // IRT_NQ - no objective,
                       // IRT_LE - minimization, IRT_GR - maximization
  static const BoolExpr DUMMY_EXPR;

 public:
  GecodeProblem(int num_vars) : vars_(*this, num_vars), obj_irt_(IRT_NQ) {}

  GecodeProblem(bool share, GecodeProblem &s) :
    Space(share, s), obj_irt_(s.obj_irt_) {
    vars_.update(*this, share, s.vars_);
    if (obj_irt_ != IRT_NQ)
      obj_.update(*this, share, s.obj_);
  }

  Space *copy(bool share);

  IntVarArray &vars() { return vars_; }
  IntVar &obj() { return obj_; }

  void SetObj(Driver::ObjType obj_type, const Gecode::LinExpr &expr) {
    obj_irt_ = obj_type == Driver::MAX ? Gecode::IRT_GR : Gecode::IRT_LE;
    obj_ = IntVar(*this, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
    rel(*this, obj_ == expr);
  }

  virtual void constrain(const Space &best);

  // The methods below perform conversion of AMPL expressions into
  // equivalent Gecode expressions. Gecode doesn't support the following
  // expressions/functions:
  // * division
  // * trigonometric functions
  //   http://www.gecode.org/pipermail/users/2011-March/003177.html
  // * log, log10, exp, pow

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
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitMax(VarArgExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitFloor(UnaryExpr e) {
    // floor does nothing because Gecode supports only integer expressions
    // currently.
    return Visit(e.arg());
  }

  LinExpr VisitCeil(UnaryExpr e) {
    // ceil does nothing because Gecode supports only integer expressions
    // currently.
    return Visit(e.arg());
  }

  LinExpr VisitAbs(UnaryExpr e) {
    return abs(Visit(e.arg()));
  }

  LinExpr VisitUnaryMinus(UnaryExpr e) {
    return -Visit(e.arg());
  }

  LinExpr VisitIf(IfExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitSqrt(UnaryExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitSum(SumExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitIntDiv(BinaryExpr e) {
    return Visit(e.lhs()) / Visit(e.rhs());
  }

  LinExpr VisitPrecision(BinaryExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitRound(BinaryExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitTrunc(BinaryExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

  LinExpr VisitCount(CountExpr e) {
    return VisitUnhandledNumericExpr(e); // TODO
  }

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
    return VisitUnhandledLogicalExpr(e); // TODO
  }

  BoolExpr VisitAtMost(RelationalExpr e) {
    return VisitUnhandledLogicalExpr(e); // TODO
  }

  BoolExpr VisitExactly(RelationalExpr e) {
    return VisitUnhandledLogicalExpr(e); // TODO
  }

  BoolExpr VisitNotAtLeast(RelationalExpr e) {
    return VisitUnhandledLogicalExpr(e); // TODO
  }

  BoolExpr VisitNotAtMost(RelationalExpr e) {
    return VisitUnhandledLogicalExpr(e); // TODO
  }

  BoolExpr VisitNotExactly(RelationalExpr e) {
    return VisitUnhandledLogicalExpr(e); // TODO
  }

  BoolExpr VisitForAll(IteratedLogicalExpr e) {
    return VisitUnhandledLogicalExpr(e); // TODO
  }

  BoolExpr VisitExists(IteratedLogicalExpr e) {
    return VisitUnhandledLogicalExpr(e); // TODO
  }

  BoolExpr VisitImplication(ImplicationExpr e) {
    return VisitUnhandledLogicalExpr(e); // TODO
  }

  BoolExpr VisitIff(BinaryLogicalExpr e) {
    return VisitUnhandledLogicalExpr(e); // TODO
  }

  BoolExpr VisitAllDiff(AllDiffExpr e) {
    int num_args = e.num_args();
    IntVarArgs x(num_args);
    for (int i = 0; i < num_args; ++i) {
      NumericExpr arg(e[i]);
      Variable var(ampl::Cast<Variable>(arg));
      if (var) {
        x[i] = vars_[var.index()];
      } else {
        IntVar gecode_var(*this,
            Gecode::Int::Limits::min, Gecode::Int::Limits::max);
        rel(*this, gecode_var == Visit(arg));
        x[i] = gecode_var;
      }
    }
    distinct(*this, x);
    return DUMMY_EXPR;
  }

  BoolExpr VisitLogicalConstant(LogicalConstant c) {
    return VisitUnhandledLogicalExpr(c); // TODO
  }
};

const BoolExpr GecodeProblem::DUMMY_EXPR((Gecode::BoolVar()));

Space *GecodeProblem::copy(bool share) {
  return new GecodeProblem(share, *this);
}

void GecodeProblem::constrain(const Space &best) {
  if (obj_irt_ != IRT_NQ)
    rel(*this, obj_, obj_irt_, static_cast<const GecodeProblem&>(best).obj_);
}
}

namespace ampl {

GecodeDriver::GecodeDriver() : oinfo_(new Option_Info()) {}

int GecodeDriver::run(char **argv) {
  if (!Read(argv, oinfo_.get()))
    return 1;

  // TODO: parse options
  /*if (!parse_options(argv))
   return 1;*/

  // Set up an optimization problem in Gecode.
  if (num_continuous_vars() != 0) {
    cerr << "Gecode doesn't support continuous variables" << endl;
    return 1;
  }
  std::auto_ptr<GecodeProblem> problem(new GecodeProblem(num_vars()));
  IntVarArray &vars = problem->vars();
  for (int j = 0, n = num_vars(); j < n; ++j) {
    double lb = GetVarLB(j), ub = GetVarUB(j);
    vars[j] = IntVar(*problem,
        lb <= negInfinity ? Gecode::Int::Limits::min : lb,
        ub >= Infinity ? Gecode::Int::Limits::max : ub);
  }

  // Post branching.
  branch(*problem, vars, Gecode::INT_VAR_SIZE_MIN, Gecode::INT_VAL_MIN);

  bool has_obj = num_objs() != 0;
  if (has_obj) {
    Gecode::LinExpr obj_expr(0);
    for (ograd *cg = GetObjGradient(0); cg; cg = cg->next)
      obj_expr = obj_expr + cg->coef * vars[cg->varno];
    if (NumericExpr expr = GetNonlinearObjExpr(0))
      obj_expr = obj_expr + problem->Visit(expr);
    problem->SetObj(GetObjType(0), obj_expr);
  }

  // Convert constraints.
  for (int i = 0, n = num_cons(); i < n; ++i) {
    Gecode::LinExpr con_expr(0);
    for (cgrad *cg = GetConGradient(i); cg; cg = cg->next)
      con_expr = con_expr + cg->coef * vars[cg->varno];
    double lb = GetConLB(i);
    double ub = GetConUB(i);
    if (i < num_nonlinear_cons())
      con_expr = con_expr + problem->Visit(GetNonlinearConExpr(i));
    if (lb <= negInfinity) {
      rel(*problem, con_expr <= ub);
    } else if (ub >= Infinity) {
      rel(*problem, con_expr >= lb);
    } else if (lb == ub) {
      rel(*problem, con_expr == lb);
    } else {
      rel(*problem, con_expr >= lb);
      rel(*problem, con_expr <= ub);
    }
  }

  // Convert logical constraints.
  for (int i = 0, n = num_logical_cons(); i < n; ++i) {
    LogicalExpr expr(GetLogicalConExpr(i));
    BoolExpr gecode_expr(problem->Visit(expr));
    if (!Cast<AllDiffExpr>(expr))
      rel(*problem, gecode_expr);
  }

  // TODO
  // finish_building_numberof();

  // Solve the problem.
  std::auto_ptr<GecodeProblem> solution;
  if (has_obj) {
    BAB<GecodeProblem> e(problem.get());
    problem.reset();
    while (GecodeProblem *next = e.next())
      solution.reset(next);
  } else {
    DFS<GecodeProblem> e(problem.get());
    problem.reset();
    solution.reset(e.next());
  }

  // Convert solution status.
  const char *status = 0;
  vector<real> primal;
  int solve_code = 0;
  if (solution.get()) {
    if (has_obj) {
      solve_code = 0;
      status = "optimal solution";
    } else {
      solve_code = 100;
      status = "feasible solution";
    }
    IntVarArray &vars = solution->vars();
    primal.resize(num_vars());
    for (int j = 0, n = num_vars(); j < n; ++j)
      primal[j] = vars[j].val();
  } else {
    solve_code = 200;
    status = "infeasible problem";
  }
  SetSolveCode(solve_code);

  char message[256];
  Sprintf(message, "%s: %s\n", oinfo_->bsname, status);
  WriteSolution(message, primal.empty() ? 0 : &primal[0], 0, oinfo_.get());
  return 0;
}
}

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
using ampl::Variable;
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
using Gecode::linear;
using Gecode::Space;

namespace {

class GecodeProblem: public Space,
  public ampl::ExprVisitor<GecodeProblem, IntVar, BoolExpr> {
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

  void SetObjType(Driver::ObjType obj_type,
      const IntArgs &c, const IntVarArgs &x) {
    obj_irt_ = obj_type == Driver::MAX ? Gecode::IRT_GR : Gecode::IRT_LE;
    obj_ = IntVar(*this, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
    linear(*this, c, x, Gecode::IRT_EQ, obj_);
  }

  virtual void constrain(const Space &best);

  IntVar VisitNumber(ampl::NumericConstant n) {
    double value = n.value();
    IntVar var(*this, value, value);
    rel(*this, var == value);
    return var;
  }

  IntVar VisitVariable(Variable v) {
    return vars_[v.index()];
  }

  BoolExpr VisitNotEqual(ampl::RelationalExpr e) {
      return Visit(e.lhs()) != Visit(e.rhs());
  }

  BoolExpr VisitAllDiff(ampl::AllDiffExpr e) {
    int num_args = e.num_args();
    IntVarArgs x(num_args);
    for (int i = 0; i < num_args; ++i) {
      NumericExpr arg(e[i]);
      Variable var(ampl::Cast<Variable>(arg));
      x[i] = var ? vars_[var.index()] : Visit(arg);
    }
    distinct(*this, x);
    return DUMMY_EXPR;
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
  for (int j = 0, n = num_vars(); j < n; ++j)
    vars[j] = IntVar(*problem, GetVarLB(j), GetVarUB(j));

  // Post branching.
  branch(*problem, vars, Gecode::INT_VAR_SIZE_MIN, Gecode::INT_VAL_MIN);

  bool has_obj = num_objs() != 0;
  if (has_obj) {
    NumericExpr expr(GetNonlinearObjExpr(0));
    NumericConstant constant(Cast<NumericConstant>(expr));
    int num_terms = 0;
    for (ograd *cg = GetObjGradient(0); cg; cg = cg->next)
      ++num_terms;
    IntArgs c(num_terms);
    IntVarArgs x(num_terms);
    int index = 0;
    for (ograd *cg = GetObjGradient(0); cg; cg = cg->next) {
      c[index] = cg->coef;
      x[index] = vars[cg->varno];
      ++index;
    }
    // TODO: constant and nonlinear part
    /*IloExpr ilo_expr(env_, constant ? constant.value() : 0);
    if (num_nonlinear_objs() > 0)
      ilo_expr += Visit(expr);*/
    problem->SetObjType(GetObjType(0), c, x);
  }

  // Convert constraints.
  for (int i = 0, n = num_cons(); i < n; ++i) {
    int num_terms = 0;
    for (cgrad *cg = GetConGradient(i); cg; cg = cg->next)
      ++num_terms;
    IntArgs c(num_terms);
    IntVarArgs x(num_terms);
    int index = 0;
    for (cgrad *cg = GetConGradient(i); cg; cg = cg->next) {
      c[index] = cg->coef;
      x[index] = vars[cg->varno];
      ++index;
    }
    double lb = GetConLB(i);
    double ub = GetConUB(i);
    if (lb <= negInfinity) {
      linear(*problem, c, x, Gecode::IRT_LE, ub);
    } else if (ub >= Infinity) {
      linear(*problem, c, x, Gecode::IRT_GQ, lb);
    } else if (lb == ub) {
      linear(*problem, c, x, Gecode::IRT_EQ, lb);
    } else {
      linear(*problem, c, x, Gecode::IRT_GQ, lb);
      linear(*problem, c, x, Gecode::IRT_LE, ub);
    }
    // TODO: nonlinear part
    //if (i < nlc)
    //   conExpr += model->ConvertExpr(con_de[i].e);
    //Con[i] = (LUrhs[i] <= conExpr <= Urhsx[i]);
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

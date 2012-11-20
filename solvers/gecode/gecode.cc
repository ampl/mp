#include "gecode.h"

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

namespace ampl {

const BoolExpr GecodeProblem::DUMMY_EXPR((Gecode::BoolVar()));

GecodeProblem::GecodeProblem(bool share, GecodeProblem &s) :
  Space(share, s), obj_irt_(s.obj_irt_) {
  vars_.update(*this, share, s.vars_);
  if (obj_irt_ != IRT_NQ)
    obj_.update(*this, share, s.obj_);
}

Space *GecodeProblem::copy(bool share) {
  return new GecodeProblem(share, *this);
}

void GecodeProblem::SetObj(
    Driver::ObjType obj_type, const Gecode::LinExpr &expr) {
  obj_irt_ = obj_type == Driver::MAX ? Gecode::IRT_GR : Gecode::IRT_LE;
  obj_ = IntVar(*this, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  rel(*this, obj_ == expr);
}

void GecodeProblem::constrain(const Space &best) {
  if (obj_irt_ != IRT_NQ)
    rel(*this, obj_, obj_irt_, static_cast<const GecodeProblem&>(best).obj_);
}

template <typename Grad>
Gecode::LinExpr GecodeProblem::ConvertExpr(Grad *grad, NumericExpr nonlinear) {
  Gecode::LinExpr expr;
  bool has_linear_part = grad != 0;
  if (has_linear_part)
    expr = grad->coef * vars_[grad->varno];
  for (grad = grad->next; grad; grad = grad->next)
    expr = expr + grad->coef * vars_[grad->varno];
  if (!nonlinear)
    return expr;
  if (has_linear_part)
    expr = expr + Visit(nonlinear);
  else
    expr = Visit(nonlinear);
  return expr;
}

LinExpr GecodeProblem::VisitMin(VarArgExpr e) {
  VarArgExpr::iterator i = e.begin();
  if (!*i)
    throw UnsupportedExprError("min with empty argument list");
  LinExpr result = Visit(*i);
  for (++i; *i; ++i)
    result = min(result, Visit(*i));
  return result;
}

LinExpr GecodeProblem::VisitMax(VarArgExpr e) {
  VarArgExpr::iterator i = e.begin();
  if (!*i)
    throw UnsupportedExprError("max with empty argument list");
  LinExpr result = Visit(*i);
  for (++i; *i; ++i)
    result = max(result, Visit(*i));
  return result;
}

LinExpr GecodeProblem::VisitFloor(UnaryExpr e) {
  // floor does nothing because Gecode supports only integer expressions
  // currently.
  NumericExpr arg = e.arg();
  if (arg.opcode() == OP_sqrt)
    return sqrt(Visit(Cast<UnaryExpr>(arg).arg()));
  return Visit(arg);
}

LinExpr GecodeProblem::VisitIf(IfExpr e) {
  Gecode::IntVar result(*this,
      Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  Gecode::BoolExpr condition = Visit(e.condition());
  rel(*this, result, Gecode::IRT_EQ,
      CreateVar(Visit(e.true_expr())), CreateVar(condition));
  rel(*this, result, Gecode::IRT_EQ,
      CreateVar(Visit(e.false_expr())), CreateVar(!condition));
  return result;
}

BoolExpr GecodeProblem::VisitAllDiff(AllDiffExpr e) {
  int num_args = e.num_args();
  IntVarArgs x(num_args);
  for (int i = 0; i < num_args; ++i) {
    NumericExpr arg(e[i]);
    if (Variable var = ampl::Cast<Variable>(arg)) {
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
    problem->SetObj(GetObjType(0),
        problem->ConvertExpr(GetObjGradient(0), GetNonlinearObjExpr(0)));
  }

  // Convert constraints.
  for (int i = 0, n = num_cons(); i < n; ++i) {
    Gecode::LinExpr con_expr(
        problem->ConvertExpr(GetConGradient(i), GetNonlinearConExpr(i)));
    double lb = GetConLB(i);
    double ub = GetConUB(i);
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
    BAB<GecodeProblem> engine(problem.get());
    problem.reset();
    while (GecodeProblem *next = engine.next())
      solution.reset(next);
  } else {
    DFS<GecodeProblem> engine(problem.get());
    problem.reset();
    solution.reset(engine.next());
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

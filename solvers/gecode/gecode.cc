#include "gecode.h"

#include <gecode/search.hh>
#include <gecode/gist.hh>

#include <iostream>
#include <limits>
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
    Problem::ObjType obj_type, const Gecode::LinExpr &expr) {
  obj_irt_ = obj_type == Problem::MAX ? Gecode::IRT_GR : Gecode::IRT_LE;
  obj_ = IntVar(*this, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  rel(*this, obj_ == expr);
}

void GecodeProblem::constrain(const Space &best) {
  if (obj_irt_ != IRT_NQ)
    rel(*this, obj_, obj_irt_, static_cast<const GecodeProblem&>(best).obj_);
}

BoolExpr NLToGecodeConverter::Convert(
    Gecode::BoolOpType op, IteratedLogicalExpr e) {
  Gecode::BoolVarArgs args(e.num_args());
  int index = 0;
  for (IteratedLogicalExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i, ++index) {
    args[index] = CreateVar(Visit(*i));
  }
  Gecode::BoolVar var(problem_, 0, 1);
  rel(problem_, op, args, var);
  return var;
}

void NLToGecodeConverter::RequireNonzeroConstRHS(
    BinaryExpr e, const std::string &func_name) {
  NumericConstant num = Cast<NumericConstant>(e.rhs());
  if (!num || num.value() != 0)
    throw UnsupportedExprError(func_name + " with nonzero second parameter");
}

template <typename Grad>
Gecode::LinExpr NLToGecodeConverter::ConvertExpr(
    Grad *grad, NumericExpr nonlinear) {
  IntVarArray &vars = problem_.vars();
  Gecode::LinExpr expr;
  bool has_linear_part = grad != 0;
  if (has_linear_part)
    expr = grad->coef * vars[grad->varno];
  for (grad = grad->next; grad; grad = grad->next)
    expr = expr + grad->coef * vars[grad->varno];
  if (!nonlinear)
    return expr;
  if (has_linear_part)
    expr = expr + Visit(nonlinear);
  else
    expr = Visit(nonlinear);
  return expr;
}

BoolExpr NLToGecodeConverter::ConvertFullExpr(LogicalExpr e, bool post) {
  AllDiffExpr alldiff = Cast<AllDiffExpr>(e);
  if (!alldiff) {
    BoolExpr result = ExprVisitor::Visit(e);
    if (post)
      rel(problem_, result);
    return result;
  }
  IntVarArray &vars = problem_.vars();
  int num_args = alldiff.num_args();
  IntVarArgs args(num_args);
  for (int i = 0; i < num_args; ++i) {
    NumericExpr arg(alldiff[i]);
    if (Variable var = ampl::Cast<Variable>(arg)) {
      args[i] = vars[var.index()];
    } else {
      IntVar gecode_var(problem_,
          Gecode::Int::Limits::min, Gecode::Int::Limits::max);
      rel(problem_, gecode_var == Visit(arg));
      args[i] = gecode_var;
    }
  }
  distinct(problem_, args);
  return Gecode::BoolVar();
}

void NLToGecodeConverter::Convert(const Problem &p) {
  if (p.num_continuous_vars() != 0)
    throw std::runtime_error("Gecode doesn't support continuous variables");

  IntVarArray &vars = problem_.vars();
  for (int j = 0, n = p.num_vars(); j < n; ++j) {
    double lb = p.GetVarLB(j), ub = p.GetVarUB(j);
    vars[j] = IntVar(problem_,
        lb <= negInfinity ? Gecode::Int::Limits::min : lb,
        ub >= Infinity ? Gecode::Int::Limits::max : ub);
  }

  // Post branching.
  branch(problem_, vars, Gecode::INT_VAR_SIZE_MIN, Gecode::INT_VAL_MIN);

  if (p.num_objs() != 0) {
    problem_.SetObj(p.GetObjType(0),
        ConvertExpr(p.GetLinearObjExpr(0), p.GetNonlinearObjExpr(0)));
  }

  // Convert constraints.
  for (int i = 0, n = p.num_cons(); i < n; ++i) {
    Gecode::LinExpr con_expr(
        ConvertExpr(p.GetLinearConExpr(i), p.GetNonlinearConExpr(i)));
    double lb = p.GetConLB(i);
    double ub = p.GetConUB(i);
    if (lb <= negInfinity) {
      rel(problem_, con_expr <= ub);
    } else if (ub >= Infinity) {
      rel(problem_, con_expr >= lb);
    } else if (lb == ub) {
      rel(problem_, con_expr == lb);
    } else {
      rel(problem_, con_expr >= lb);
      rel(problem_, con_expr <= ub);
    }
  }

  // Convert logical constraints.
  for (int i = 0, n = p.num_logical_cons(); i < n; ++i)
    ConvertFullExpr(p.GetLogicalConExpr(i));
}

LinExpr NLToGecodeConverter::VisitMin(VarArgExpr e) {
  VarArgExpr::iterator i = e.begin();
  if (!*i)
    throw UnsupportedExprError("min with empty argument list");
  IntVarArgs args;
  for (; *i; ++i)
    args << CreateVar(Visit(*i));
  Gecode::IntVar result(problem_,
      Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  min(problem_, args, result);
  return result;
}

LinExpr NLToGecodeConverter::VisitMax(VarArgExpr e) {
  VarArgExpr::iterator i = e.begin();
  if (!*i)
    throw UnsupportedExprError("max with empty argument list");
  IntVarArgs args;
  for (; *i; ++i)
    args << CreateVar(Visit(*i));
  Gecode::IntVar result(problem_,
      Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  max(problem_, args, result);
  return result;
}

LinExpr NLToGecodeConverter::VisitFloor(UnaryExpr e) {
  // floor does nothing because Gecode supports only integer expressions
  // currently.
  NumericExpr arg = e.arg();
  if (arg.opcode() == OP_sqrt)
    return sqrt(Visit(Cast<UnaryExpr>(arg).arg()));
  return Visit(arg);
}

LinExpr NLToGecodeConverter::VisitIf(IfExpr e) {
  Gecode::IntVar result(problem_,
      Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  Gecode::BoolExpr condition = Visit(e.condition());
  rel(problem_, result, Gecode::IRT_EQ,
      CreateVar(Visit(e.true_expr())), CreateVar(condition));
  rel(problem_, result, Gecode::IRT_EQ,
      CreateVar(Visit(e.false_expr())), CreateVar(!condition));
  return result;
}

LinExpr NLToGecodeConverter::VisitSum(SumExpr e) {
  SumExpr::iterator i = e.begin(), end = e.end();
  if (i == end)
    return 0;
  LinExpr sum = Visit(*i++);
  for (; i != end; ++i)
    sum = sum + Visit(*i);
  return sum;
}

LinExpr NLToGecodeConverter::VisitCount(CountExpr e) {
  Gecode::BoolVarArgs args(e.num_args());
  int index = 0;
  for (CountExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i, ++index) {
    args[index] = CreateVar(Visit(*i));
  }
  Gecode::IntVar result(problem_, 0, e.num_args());
  linear(problem_, args, Gecode::IRT_EQ, result);
  return result;
}

LinExpr NLToGecodeConverter::VisitNumberOf(NumberOfExpr e) {
  // Gecode only supports global cardinality (count) constraint where no other
  // values except those specified may occur, so we use only local count
  // constraints.
  Gecode::IntVar result(problem_,
      Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  int index = 0;
  Gecode::IntVarArgs args(e.num_args());
  for (NumberOfExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    args[index++] = CreateVar(Visit(*i));
  count(problem_, args, CreateVar(Visit(e.value())), Gecode::IRT_EQ, result);
  return result;
}

BoolExpr NLToGecodeConverter::VisitImplication(ImplicationExpr e) {
  BoolExpr condition = Visit(e.condition());
  LogicalConstant c = Cast<LogicalConstant>(e.false_expr());
  if (c && !c.value())
    return condition >> Visit(e.true_expr());
  return (condition && Visit(e.true_expr())) ||
        (!condition && Visit(e.false_expr()));
}

BoolExpr NLToGecodeConverter::VisitAllDiff(AllDiffExpr) {
  throw UnsupportedExprError("nested 'alldiff'");
  return BoolExpr();
}

GecodeDriver::GecodeDriver() : Driver(options_), options_(*this) {
  std::string name = str(fmt::Format("Gecode {0}") << GECODE_VERSION);
  options_.set_long_solver_name(name);
  options_.set_options_var_name("gecode_options");
  options_.set_version(name);
}

int GecodeDriver::Run(char **argv) {
  Problem &problem = Driver::problem();
  if (!problem.Read(argv, options_))
    return 1;

  if (!ParseOptions(argv))
    return 1;

  // Set up an optimization problem in Gecode.
  std::auto_ptr<NLToGecodeConverter>
    converter(new NLToGecodeConverter(problem.num_vars()));
  converter->Convert(problem);

  // Solve the problem.
  double obj_val = std::numeric_limits<double>::quiet_NaN();
  std::auto_ptr<GecodeProblem> solution;
  bool has_obj = problem.num_objs() != 0;
  Gecode::Search::Statistics stats;
  if (has_obj) {
    BAB<GecodeProblem> engine(&converter->problem());
    converter.reset();
    while (GecodeProblem *next = engine.next())
      solution.reset(next);
    if (solution.get())
      obj_val = solution->obj().val();
    stats = engine.statistics();
  } else {
    DFS<GecodeProblem> engine(&converter->problem());
    converter.reset();
    solution.reset(engine.next());
    stats = engine.statistics();
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
    int num_vars = problem.num_vars();
    primal.resize(num_vars);
    for (int j = 0; j < num_vars; ++j)
      primal[j] = vars[j].val();
  } else {
    solve_code = 200;
    status = "infeasible problem";
  }
  problem.SetSolveCode(solve_code);

  fmt::Formatter format;
  format("{0}: {1}\n") << options_.long_solver_name() << status;
  format("{0} nodes, {1} fails") << stats.node << stats.fail;
  if (has_obj)
    format(", objective {0}") << ObjPrec(obj_val);
  HandleSolution(format.c_str(), primal.empty() ? 0 : &primal[0], 0, obj_val);
  return 0;
}
}

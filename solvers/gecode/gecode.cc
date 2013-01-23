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

#include "gecode.h"

#include <limits>
#include <memory>
#include <string>
#include <vector>

using Gecode::BoolExpr;
using Gecode::IntVarArgs;
using Gecode::IntVar;
using Gecode::IntVarArray;
using Gecode::LinExpr;
namespace Search = Gecode::Search;

namespace {

const ampl::OptionValue<Gecode::IntVarBranch> VAR_BRANCHINGS[] = {
    {"none",            Gecode::INT_VAR_NONE},
    {"rnd",             Gecode::INT_VAR_RND},
    {"degree_min",      Gecode::INT_VAR_DEGREE_MIN},
    {"degree_max",      Gecode::INT_VAR_DEGREE_MAX},
    {"afc_min",         Gecode::INT_VAR_AFC_MIN},
    {"afc_max",         Gecode::INT_VAR_AFC_MAX},
    {"min_min",         Gecode::INT_VAR_MIN_MIN},
    {"min_max",         Gecode::INT_VAR_MIN_MAX},
    {"max_min",         Gecode::INT_VAR_MAX_MIN},
    {"max_max",         Gecode::INT_VAR_MAX_MAX},
    {"size_min",        Gecode::INT_VAR_SIZE_MIN},
    {"size_max",        Gecode::INT_VAR_SIZE_MAX},
    {"size_degree_min", Gecode::INT_VAR_SIZE_DEGREE_MIN},
    {"size_degree_max", Gecode::INT_VAR_SIZE_DEGREE_MAX},
    {"size_afc_min",    Gecode::INT_VAR_SIZE_AFC_MIN},
    {"size_afc_max",    Gecode::INT_VAR_SIZE_AFC_MAX},
    {"regret_min_min",  Gecode::INT_VAR_REGRET_MIN_MIN},
    {"regret_min_max",  Gecode::INT_VAR_REGRET_MIN_MAX},
    {"regret_max_min",  Gecode::INT_VAR_REGRET_MAX_MIN},
    {"regret_max_max",  Gecode::INT_VAR_REGRET_MAX_MAX},
    {}
};

const ampl::OptionValue<Gecode::IntValBranch> VAL_BRANCHINGS[] = {
    {"min",        Gecode::INT_VAL_MIN},
    {"med",        Gecode::INT_VAL_MED},
    {"max",        Gecode::INT_VAL_MAX},
    {"rnd",        Gecode::INT_VAL_RND},
    {"split_min",  Gecode::INT_VAL_SPLIT_MIN},
    {"split_max",  Gecode::INT_VAL_SPLIT_MAX},
    {"range_min",  Gecode::INT_VAL_RANGE_MIN},
    {"range_max",  Gecode::INT_VAL_RANGE_MAX},
    {"values_min", Gecode::INT_VALUES_MIN},
    {"values_max", Gecode::INT_VALUES_MAX},
    {}
};
}

namespace ampl {

GecodeProblem::GecodeProblem(bool share, GecodeProblem &s) :
  Gecode::Space(share, s), obj_irt_(s.obj_irt_) {
  vars_.update(*this, share, s.vars_);
  if (obj_irt_ != Gecode::IRT_NQ)
    obj_.update(*this, share, s.obj_);
}

Gecode::Space *GecodeProblem::copy(bool share) {
  return new GecodeProblem(share, *this);
}

void GecodeProblem::SetObj(Problem::ObjType obj_type, const LinExpr &expr) {
  obj_irt_ = obj_type == Problem::MAX ? Gecode::IRT_GR : Gecode::IRT_LE;
  obj_ = Gecode::expr(*this, expr);
}

void GecodeProblem::constrain(const Gecode::Space &best) {
  if (obj_irt_ != Gecode::IRT_NQ)
    rel(*this, obj_, obj_irt_, static_cast<const GecodeProblem&>(best).obj_);
}

BoolExpr NLToGecodeConverter::Convert(
    Gecode::BoolOpType op, IteratedLogicalExpr e) {
  Gecode::BoolVarArgs args(e.num_args());
  int index = 0;
  for (IteratedLogicalExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i, ++index) {
    args[index] = Gecode::expr(problem_, Visit(*i));
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

template <typename Term>
LinExpr NLToGecodeConverter::ConvertExpr(
    LinearExpr<Term> linear, NumericExpr nonlinear) {
  IntVarArray &vars = problem_.vars();
  LinExpr expr;
  typename LinearExpr<Term>::iterator i = linear.begin(), end = linear.end();
  bool has_linear_part = i != end;
  if (has_linear_part)
    expr = i->coef() * vars[i->var_index()];
  for (++i; i != end; ++i)
    expr = expr + i->coef() * vars[i->var_index()];
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
    if (Variable var = ampl::Cast<Variable>(arg))
      args[i] = vars[var.index()];
    else
      args[i] = Gecode::expr(problem_, Visit(arg));
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

  if (p.num_objs() != 0) {
    problem_.SetObj(p.GetObjType(0),
        ConvertExpr(p.GetLinearObjExpr(0), p.GetNonlinearObjExpr(0)));
  }

  // Convert constraints.
  for (int i = 0, n = p.num_cons(); i < n; ++i) {
    LinExpr con_expr(
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
    args << Gecode::expr(problem_, Visit(*i));
  IntVar result(problem_, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  min(problem_, args, result);
  return result;
}

LinExpr NLToGecodeConverter::VisitMax(VarArgExpr e) {
  VarArgExpr::iterator i = e.begin();
  if (!*i)
    throw UnsupportedExprError("max with empty argument list");
  IntVarArgs args;
  for (; *i; ++i)
    args << Gecode::expr(problem_, Visit(*i));
  IntVar result(problem_, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
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
  IntVar result(problem_, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  BoolExpr condition = Visit(e.condition());
  rel(problem_, result, Gecode::IRT_EQ,
      Gecode::expr(problem_, Visit(e.true_expr())),
      Gecode::expr(problem_, condition));
  rel(problem_, result, Gecode::IRT_EQ,
      Gecode::expr(problem_, Visit(e.false_expr())),
      Gecode::expr(problem_, !condition));
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
    args[index] = Gecode::expr(problem_, Visit(*i));
  }
  IntVar result(problem_, 0, e.num_args());
  Gecode::linear(problem_, args, Gecode::IRT_EQ, result);
  return result;
}

LinExpr NLToGecodeConverter::VisitNumberOf(NumberOfExpr e) {
  // Gecode only supports global cardinality (count) constraint where no other
  // values except those specified may occur, so we use only local count
  // constraints.
  IntVar result(problem_, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  int index = 0;
  IntVarArgs args(e.num_args());
  for (NumberOfExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    args[index++] = Gecode::expr(problem_, Visit(*i));
  count(problem_, args, Gecode::expr(problem_, Visit(e.value())),
      Gecode::IRT_EQ, result);
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

GecodeSolver::Stop::Stop(const GecodeSolver &s)
: sh_(s), solver_(s), time_limit_in_milliseconds_(s.time_limit_ * 1000) {
  has_limit_ = time_limit_in_milliseconds_ < DBL_MAX ||
      s.node_limit_ != ULONG_MAX || s.fail_limit_ != ULONG_MAX ||
      s.memory_limit_ != std::numeric_limits<std::size_t>::max();
  timer_.start();
}

bool GecodeSolver::Stop::stop(
    const Search::Statistics &s, const Search::Options &) {
  if (SignalHandler::stop()) return true;
  return has_limit_ && (timer_.stop() > time_limit_in_milliseconds_ ||
      s.node > solver_.node_limit_ || s.fail > solver_.fail_limit_ ||
      s.memory > solver_.memory_limit_);
}

template <typename T>
void GecodeSolver::SetStrOption(const char *name, const char *value,
    const OptionInfo<T> &info) {
  for (const OptionValue<T> *p = info.values; p->name; ++p) {
    if (std::strcmp(value, p->name) == 0) {
      info.value = p->value;
      return;
    }
  }
  ReportError("Invalid value {} for option {}") << value << name;
}

template <typename T, typename OptionT>
void GecodeSolver::SetOption(const char *name, T value, OptionT *option) {
  if (value < 0)
    ReportError("Invalid value {} for option {}") << value << name;
  else
    *option = value;
}

GecodeSolver::GecodeSolver()
: Solver<GecodeSolver>("gecode", "gecode " GECODE_VERSION), output_(false),
  var_branching_(Gecode::INT_VAR_SIZE_MIN),
  val_branching_(Gecode::INT_VAL_MIN),
  time_limit_(DBL_MAX), node_limit_(ULONG_MAX), fail_limit_(ULONG_MAX),
  memory_limit_(std::numeric_limits<std::size_t>::max()) {

  set_version("Gecode " GECODE_VERSION);

  AddIntOption("outlev",
      "0 or 1 (default 0):  Whether to print solution log.",
      &GecodeSolver::EnableOutput);

  AddStrOption("var_branching",
      "Variable branching.  Possible values:\n"
      "      none            - first unassigned\n"
      "      rnd             - random\n"
      "      degree_min      - smallest degree\n"
      "      degree_max      - largest degree\n"
      "      afc_min         - smallest accumulated failure count (AFC)\n"
      "      afc_max         - largest accumulated failure count (AFC)\n"
      "      min_min         - smallest minimum value\n"
      "      min_max         - largest minimum value\n"
      "      max_min         - smallest maximum value\n"
      "      max_max         - largest maximum value\n"
      "      size_min        - smallest domain size (default)\n"
      "      size_max        - largest domain size\n"
      "      size_degree_min - smallest domain size divided by degree\n"
      "      size_degree_max - largest domain size divided by degree\n"
      "      size_afc_min    - smallest domain size divided by AFC\n"
      "      size_afc_max    - largest domain size divided by AFC\n"
      "      regret_min_min  - smallest minimum-regret\n"
      "      regret_min_max  - largest minimum-regret\n"
      "      regret_max_min  - smallest maximum-regret\n"
      "      regret_max_max  - largest maximum-regret\n",
      &GecodeSolver::SetStrOption<Gecode::IntVarBranch>,
      OptionInfo<Gecode::IntVarBranch>(VAR_BRANCHINGS, var_branching_));

  AddStrOption("val_branching",
      "Value branching.  Possible values:\n"
      "      min        - smallest value (default)\n"
      "      med        - greatest value not greater than the median\n"
      "      max        - largest value\n"
      "      rnd        - random value\n"
      "      split_min  - values not greater than mean of smallest and\n"
      "                   largest value\n"
      "      split_max  - values greater than mean of smallest and largest\n"
      "                   value\n"
      "      range_min  - values from smallest range, if domain has several\n"
      "                   ranges; otherwise, values not greater than mean of\n"
      "                   smallest and largest value\n"
      "      range_max  - values from largest range, if domain has several\n"
      "                   ranges; otherwise, values greater than mean of\n"
      "                   smallest and largest value\n"
      "      values_min - all values starting from smallest\n"
      "      values_min - all values starting from largest\n",
      &GecodeSolver::SetStrOption<Gecode::IntValBranch>,
      OptionInfo<Gecode::IntValBranch>(VAL_BRANCHINGS, val_branching_));

  AddDblOption("threads",
      "The number of parallel threads to use.  Assume that your computer\n"
      "has m processing units and that the value for threads is n.\n"
      "    * If n = 0, then m threads are used (as many as available\n"
      "      processing units).\n"
      "    * If n >= 1, then n threads are used (absolute number of\n"
      "      threads to be used).\n"
      "    * If n <= −1, then m + n threads are used (absolute number\n"
      "      of processing units not to be used). For example, when\n"
      "      n = −6 and m = 8, then 2 threads are used.\n"
      "    * If 0 < n < 1, then n * m threads are used (relative number\n"
      "      of processing units to be used). For example, when n = 0.5\n"
      "      and m = 8, then 4 threads are used.\n"
      "    * If −1 < n < 0, then (1 + n) * m threads are used (relative\n"
      "      number of processing units not to be used). For example,\n"
      "      when n = −0.25 and m = 8, then 6 threads are used.\n"
      "All values are rounded and at least one thread is used.\n",
      &GecodeSolver::SetDblOption, &options_.threads);

  AddIntOption("c_d", "Commit recomputation distance.",
      &GecodeSolver::SetOption<int, unsigned>, &options_.c_d);
  AddIntOption("a_d", "Adaptive recomputation distance.",
      &GecodeSolver::SetOption<int, unsigned>, &options_.a_d);

  AddDblOption("timelimit", "Time limit.",
      &GecodeSolver::SetOption<double, double>, &time_limit_);
  AddIntOption("nodelimit", "Node limit.",
      &GecodeSolver::SetOption<int, unsigned long>, &node_limit_);
  AddIntOption("faillimit", "Fail limit.",
      &GecodeSolver::SetOption<int, unsigned long>, &fail_limit_);
  AddIntOption("memorylimit", "Memory limit.",
      &GecodeSolver::SetOption<int, std::size_t>, &memory_limit_);
}

int GecodeSolver::Run(char **argv) {
  if (!ProcessArgs(argv, *this))
    return 1;

  // Set up an optimization problem in Gecode.
  Problem &problem = Solver::problem();
  std::auto_ptr<NLToGecodeConverter>
    converter(new NLToGecodeConverter(problem.num_vars()));
  converter->Convert(problem);

  // Post branching.
  GecodeProblem &gecode_problem = converter->problem();
  branch(gecode_problem, gecode_problem.vars(),
      var_branching_, val_branching_);

  Stop stop(*this);
  options_.stop = &stop;

  // Solve the problem.
  double obj_val = std::numeric_limits<double>::quiet_NaN();
  std::auto_ptr<GecodeProblem> solution;
  bool has_obj = problem.num_objs() != 0;
  Search::Statistics stats;
  bool stopped = false;
  if (has_obj) {
    Gecode::BAB<GecodeProblem> engine(&gecode_problem, options_);
    converter.reset();
    while (GecodeProblem *next = engine.next()) {
      if (output_)
        fmt::Print("Best objective: {}\n") << next->obj().val();
      solution.reset(next);
    }
    if (solution.get())
      obj_val = solution->obj().val();
    stopped = engine.stopped();
    stats = engine.statistics();
  } else {
    Gecode::DFS<GecodeProblem> engine(&gecode_problem, options_);
    converter.reset();
    solution.reset(engine.next());
    stopped = engine.stopped();
    stats = engine.statistics();
  }

  // Convert solution status.
  const char *status = 0;
  int solve_code = 0;
  if (stopped) {
    solve_code = 600;
    status = "interrupted";
  } else if (!solution.get()) {
    solve_code = 200;
    status = "infeasible problem";
  } else if (has_obj) {
    solve_code = 0;
    status = "optimal solution";
  } else {
    solve_code = 100;
    status = "feasible solution";
  }
  problem.set_solve_code(solve_code);

  std::vector<real> primal;
  if (solution.get()) {
    IntVarArray &vars = solution->vars();
    int num_vars = problem.num_vars();
    primal.resize(num_vars);
    for (int j = 0; j < num_vars; ++j)
      primal[j] = vars[j].val();
  }

  fmt::Formatter format;
  format("{}: {}\n") << long_name() << status;
  format("{} nodes, {} fails") << stats.node << stats.fail;
  if (has_obj && solution.get())
    format(", objective {}") << ObjPrec(obj_val);
  HandleSolution(format.c_str(), primal.empty() ? 0 : &primal[0], 0, obj_val);
  return 0;
}
}

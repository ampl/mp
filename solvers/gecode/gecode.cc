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

#include "gecode/gecode.h"

#include <limits>
#include <memory>
#include <vector>

using Gecode::BoolExpr;
using Gecode::IntValBranch;
using Gecode::IntVarArgs;
using Gecode::IntVar;
using Gecode::IntVarArray;
using Gecode::IntVarBranch;
using Gecode::Reify;
namespace Search = Gecode::Search;

namespace {

const mp::OptionValueInfo INT_PROP_LEVELS[] = {
  {"val", "value propagation or consistency (naive)", Gecode::IPL_VAL},
  {"bnd", "bounds propagation or consistency",        Gecode::IPL_BND},
  {"dom", "domain propagation or consistency",        Gecode::IPL_DOM},
  {"def", "the default propagation for a constraint", Gecode::IPL_DEF}
};

const mp::OptionValueInfo VAR_BRANCHINGS[] = {
  {
    "none",
    "first unassigned",
    IntVarBranch::SEL_NONE
  },
  {
    "rnd",
    "random",
    IntVarBranch::SEL_RND
  },
  {
    "degree_min",
    "smallest degree",
    IntVarBranch::SEL_DEGREE_MIN
  },
  {
    "degree_max",
    "largest degree",
    IntVarBranch::SEL_DEGREE_MAX
  },
  {
    "afc_min",
    "smallest accumulated failure count (AFC)",
    IntVarBranch::SEL_AFC_MIN
  },
  {
    "afc_max",
    "largest accumulated failure count (AFC)",
    IntVarBranch::SEL_AFC_MAX
  },
  {
    "action_min",
    "lowest action",
    IntVarBranch::SEL_ACTION_MIN
  },
  {
    "action_max",
    "highest action",
    IntVarBranch::SEL_ACTION_MAX
  },
  {
    "min_min",
    "smallest minimum value",
    IntVarBranch::SEL_MIN_MIN
  },
  {
    "min_max",
    "largest minimum value",
    IntVarBranch::SEL_MIN_MAX
  },
  {
    "max_min",
    "smallest maximum value",
    IntVarBranch::SEL_MAX_MIN
  },
  {
    "max_max",
    "largest maximum value",
    IntVarBranch::SEL_MAX_MAX
  },
  {
    "size_min",
    "smallest domain size (default)",
    IntVarBranch::SEL_SIZE_MIN
  },
  {
    "size_max",
    "largest domain size",
    IntVarBranch::SEL_SIZE_MAX
  },
  {
    "degree_size_min",
    "smallest domain size divided by degree",
    IntVarBranch::SEL_DEGREE_SIZE_MIN
  },
  {
    "degree_size_max",
    "largest domain size divided by degree",
    IntVarBranch::SEL_DEGREE_SIZE_MAX
  },
  {
    "afc_size_min",
    "smallest domain size divided by AFC",
    IntVarBranch::SEL_AFC_SIZE_MIN
  },
  {
    "afc_size_max",
    "largest domain size divided by AFC",
    IntVarBranch::SEL_AFC_SIZE_MAX
  },
  {
    "action_size_min",
    "smallest action divided by domain size",
    IntVarBranch::SEL_ACTION_SIZE_MIN},
  {
    "action_size_max",
    "largest action divided by domain size",
    IntVarBranch::SEL_ACTION_SIZE_MAX
  },
  {
    "regret_min_min",
    "smallest minimum-regret",
    IntVarBranch::SEL_REGRET_MIN_MIN
  },
  {
    "regret_min_max",
    "largest minimum-regret",
    IntVarBranch::SEL_REGRET_MIN_MAX
  },
  {
    "regret_max_min",
    "smallest maximum-regret",
    IntVarBranch::SEL_REGRET_MAX_MIN
  },
  {
    "regret_max_max",
    "largest maximum-regret",
    IntVarBranch::SEL_REGRET_MAX_MAX
  }
};

const mp::OptionValueInfo VAL_BRANCHINGS[] = {
  {
    "min",
    "smallest value (default)",
    IntValBranch::SEL_MIN
  },
  {
    "med",
    "greatest value not greater than the median",
    IntValBranch::SEL_MED
  },
  {
    "max",
    "largest value",
    IntValBranch::SEL_MAX
  },
  {
    "rnd",
    "random value",
    IntValBranch::SEL_RND
  },
  {
    "split_min",
    "values not greater than mean of smallest and largest value",
    IntValBranch::SEL_SPLIT_MIN
  },
  {
    "split_max",
    "values greater than mean of smallest and largest value",
    IntValBranch::SEL_SPLIT_MAX
  },
  {
    "range_min",
    "values from smallest range, if domain has several ranges; "
    "otherwise, values not greater than mean of smallest and largest value",
    IntValBranch::SEL_RANGE_MIN
  },
  {
    "range_max",
    "values from largest range, if domain has several ranges; "
    "otherwise, values greater than mean of smallest and largest value",
    IntValBranch::SEL_RANGE_MAX
  },
  {
    "values_min",
    "all values starting from smallest",
    IntValBranch::SEL_VALUES_MIN
  },
  {
    "values_max",
    "all values starting from largest",
    IntValBranch::SEL_VALUES_MAX
  }
};

const mp::OptionValueInfo RESTART_MODES[] = {
  {"none",      "no restarts",                     Gecode::RM_NONE},
  {"constant",  "restart with constant sequence",  Gecode::RM_CONSTANT},
  {"linear",    "restart with linear sequence",    Gecode::RM_LINEAR},
  {"luby",      "restart with Luby sequence",      Gecode::RM_LUBY},
  {"geometric", "restart with geometric sequence", Gecode::RM_GEOMETRIC}
};
}

namespace mp {

GecodeProblem::GecodeProblem(int num_vars, Gecode::IntPropLevel ipl) :
  vars_(space(), num_vars), obj_irt_(Gecode::IRT_NQ), ipl_(ipl) {
}

#if GECODE_VERSION_NUMBER > 600000
GecodeProblem::GecodeProblem(GecodeProblem &s) :
  Gecode::Space(s), obj_irt_(s.obj_irt_), ipl_(s.ipl_) {
  vars_.update(*this, s.vars_);
  if (obj_irt_ != Gecode::IRT_NQ)
    obj_.update(*this, s.obj_);
}
Gecode::Space* GecodeProblem::copy() {
  return new GecodeProblem(*this);
}
#else
GecodeProblem::GecodeProblem(bool share, GecodeProblem& s) :
  Gecode::Space(share, s), obj_irt_(s.obj_irt_), ipl_(s.ipl_) {
  vars_.update(*this, share, s.vars_);
  if (obj_irt_ != Gecode::IRT_NQ)
    obj_.update(*this, share, s.obj_);
}
Gecode::Space* GecodeProblem::copy(bool share) {
  return new GecodeProblem(share, *this);
}
#endif



void GecodeProblem::SetObj(obj::Type obj_type, const LinExpr &expr) {
  obj_irt_ = obj_type == obj::MAX ? Gecode::IRT_GR : Gecode::IRT_LE;
  obj_ = Gecode::expr(*this, expr);
}

void GecodeProblem::constrain(const Gecode::Space &best) {
  if (obj_irt_ != Gecode::IRT_NQ) {
    rel(*this, obj_, obj_irt_,
        static_cast<const GecodeProblem&>(best).obj_, ipl_);
  }
}

BoolExpr MPToGecodeConverter::Convert(
    Gecode::BoolOpType op, IteratedLogicalExpr e) {
  Gecode::BoolVarArgs args(e.num_args());
  int index = 0;
  for (IteratedLogicalExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i, ++index) {
    args[index] = Gecode::expr(problem_, Visit(*i), ipl_);
  }
  Gecode::BoolVar var(problem_, 0, 1);
  rel(problem_, op, args, var, ipl_);
  return var;
}

LinExpr MPToGecodeConverter::Convert(IteratedExpr e, VarArgFunc f) {
  IntVarArgs args;
  for (VarArgExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    args << Gecode::expr(problem_, Visit(*i), ipl_);
  IntVar result(problem_, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  f(problem_, args, result, ipl_);
  return result;
}

void MPToGecodeConverter::RequireZeroRHS(
    BinaryExpr e, fmt::StringRef func_name) {
  if (!IsZero(e.rhs()))
    throw MakeUnsupportedError("{} with nonzero second parameter", func_name);
}

LinExpr MPToGecodeConverter::ConvertExpr(
    const LinearExpr &linear, NumericExpr nonlinear) {
  IntVarArray &vars = problem_.vars();
  LinExpr expr;
  LinearExpr::const_iterator i = linear.begin(), end = linear.end();
  bool has_linear_part = i != end;
  if (has_linear_part) {
    expr = CastToInt(i->coef()) * vars[i->var_index()];
    for (++i; i != end; ++i)
      expr = expr + CastToInt(i->coef()) * vars[i->var_index()];
  }
  if (!nonlinear)
    return expr;
  if (has_linear_part)
    expr = expr + Visit(nonlinear);
  else
    expr = Visit(nonlinear);
  return expr;
}

Gecode::IntPropLevel MPToGecodeConverter::GetIPL(int con_index) const {
  if (!ipl_suffix_)
    return ipl_;
  int value = ipl_suffix_.value(con_index);
  assert(value == Gecode::IPL_VAL || value == Gecode::IPL_BND ||
         value == Gecode::IPL_DOM || value == Gecode::IPL_DEF);
  if (value < 0 || value > Gecode::IPL_DEF)
    throw Error("Invalid value \"{}\" for suffix \"ipl\"", value);
  return static_cast<Gecode::IntPropLevel>(value);
}

void MPToGecodeConverter::Convert(const Problem &p) {
  double inf = std::numeric_limits<double>::infinity();
  IntVarArray &vars = problem_.vars();
  for (int j = 0, n = p.num_vars(); j < n; ++j) {
    Problem::Variable var = p.var(j);
    if (var.type() == mp::var::CONTINUOUS)
      throw Error("Gecode doesn't support continuous variables");
    double lb = var.lb(), ub = var.ub();
    vars[j] = IntVar(problem_,
        lb <= -inf ? Gecode::Int::Limits::min : CastToInt(lb),
        ub >=  inf ? Gecode::Int::Limits::max : CastToInt(ub));
  }
  int num_common_exprs = p.num_common_exprs();
  common_exprs_.resize(num_common_exprs);
  for (int i = 0; i < num_common_exprs; ++i) {
    Problem::CommonExpr expr = p.common_expr(i);
    common_exprs_[i] = ConvertExpr(expr.linear_expr(), expr.nonlinear_expr());
  }

  if (p.num_objs() > 0) {
    Problem::Objective obj = p.obj(0);
    problem_.SetObj(obj.type(),
                    ConvertExpr(obj.linear_expr(), obj.nonlinear_expr()));
  }

  ipl_suffix_ = p.suffixes(suf::CON).Find<int>("ipl");

  class IPLSetter {
   private:
    Gecode::IntPropLevel &ipl_;
    Gecode::IntPropLevel saved_value_;

   public:
    IPLSetter(Gecode::IntPropLevel &ipl, Gecode::IntPropLevel new_value) :
      ipl_(ipl), saved_value_(ipl) {
      ipl = new_value;
    }
    ~IPLSetter() { ipl_ = saved_value_; }
  };

  // Convert algebraic constraints.
  for (int i = 0, n = p.num_algebraic_cons(); i < n; ++i) {
    Problem::AlgebraicCon con = p.algebraic_con(i);
    LinExpr con_expr(
        ConvertExpr(con.linear_expr(), con.nonlinear_expr()));
    double lb = con.lb(), ub = con.ub();
    IPLSetter ipl_setter(ipl_, GetIPL(i));
    if (lb <= -inf) {
      rel(problem_, con_expr <= CastToInt(ub), ipl_);
      continue;
    }
    if (ub >= inf) {
      rel(problem_, con_expr >= CastToInt(lb), ipl_);
      continue;
    }
    int int_lb = CastToInt(lb), int_ub = CastToInt(ub);
    if (int_lb == int_ub) {
      rel(problem_, con_expr == int_lb, ipl_);
    } else {
      rel(problem_, con_expr >= int_lb, ipl_);
      rel(problem_, con_expr <= int_ub, ipl_);
    }
  }

  // Convert logical constraints.
  int num_logical_cons = p.num_logical_cons();
  for (int i = 0; i < num_logical_cons; ++i) {
    LogicalExpr e = p.logical_con(i).expr();
    IPLSetter ipl_setter(ipl_, GetIPL(p.num_algebraic_cons() + i));
    if (e.kind() != expr::ALLDIFF) {
      rel(problem_, Visit(e), ipl_);
      continue;
    }
    PairwiseExpr alldiff = Cast<PairwiseExpr>(e);
    int num_args = alldiff.num_args();
    IntVarArgs args(num_args);
    for (int i = 0; i < num_args; ++i) {
      NumericExpr arg = alldiff.arg(i);
      if (arg.kind() == expr::VARIABLE)
        args[i] = vars[Cast<Variable>(arg).index()];
      else
        args[i] = Gecode::expr(problem_, Visit(arg), ipl_);
    }
    distinct(problem_, args, ipl_);
  }
}

LinExpr MPToGecodeConverter::VisitFloor(UnaryExpr e) {
  // floor does nothing because Gecode supports only integer expressions
  // currently.
  NumericExpr arg = e.arg();
  if (arg.kind() == expr::SQRT)
    return sqrt(Visit(Cast<UnaryExpr>(arg).arg()));
  return Visit(arg);
}

LinExpr MPToGecodeConverter::VisitIf(IfExpr e) {
  IntVar result(problem_, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  BoolExpr condition = Visit(e.condition());
  NumericExpr then_expr = e.then_expr(), else_expr = e.else_expr();
  NumericConstant false_const = Cast<NumericConstant>(else_expr);
  if (false_const && false_const.value() == 0) {
    NumericConstant true_const = Cast<NumericConstant>(then_expr);
    if (true_const && true_const.value() == 1) {
      Gecode::channel(problem_, Gecode::expr(problem_, condition, ipl_), result);
      return result;
    }
  }
  rel(problem_, result, Gecode::IRT_EQ,
      Gecode::expr(problem_, Visit(then_expr), ipl_),
      Reify(Gecode::expr(problem_, condition, ipl_), Gecode::RM_IMP), ipl_);
  rel(problem_, result, Gecode::IRT_EQ,
      Gecode::expr(problem_, Visit(else_expr), ipl_),
      Reify(Gecode::expr(problem_, !condition, ipl_), Gecode::RM_IMP), ipl_);
  return result;
}

LinExpr MPToGecodeConverter::VisitSum(SumExpr e) {
  SumExpr::iterator i = e.begin(), end = e.end();
  if (i == end)
    return 0;
  LinExpr sum = Visit(*i++);
  for (; i != end; ++i)
    sum = sum + Visit(*i);
  return sum;
}

LinExpr MPToGecodeConverter::VisitCount(CountExpr e) {
  Gecode::BoolVarArgs args(e.num_args());
  int index = 0;
  for (CountExpr::iterator
      i = e.begin(), end = e.end(); i != end; ++i, ++index) {
    args[index] = Gecode::expr(problem_, Visit(*i), ipl_);
  }
  IntVar result(problem_, 0, e.num_args());
  Gecode::linear(problem_, args, Gecode::IRT_EQ, result, ipl_);
  return result;
}

LinExpr MPToGecodeConverter::VisitNumberOf(IteratedExpr e) {
  // Gecode only supports global cardinality (count) constraint where no other
  // values except those specified may occur, so we use only local count
  // constraints.
  IntVar result(problem_, Gecode::Int::Limits::min, Gecode::Int::Limits::max);
  int num_args = e.num_args();
  IntVarArgs args(num_args - 1);
  for (int i = 1; i < num_args; ++i)
    args[i - 1] = Gecode::expr(problem_, Visit(e.arg(i)), ipl_);
  count(problem_, args, Gecode::expr(problem_, Visit(e.arg(0)), ipl_),
      Gecode::IRT_EQ, result);
  return result;
}

BoolExpr MPToGecodeConverter::LogicalExprConverter::VisitImplication(
    ImplicationExpr e) {
  BoolExpr condition = Visit(e.condition());
  return (condition && Visit(e.then_expr())) ||
        (!condition && Visit(e.else_expr()));
}

BoolExpr MPToGecodeConverter::LogicalExprConverter::VisitAllDiff(
    PairwiseExpr e) {
  bool negate = e.kind() == expr::NOT_ALLDIFF;
  int n = e.num_args();
  std::vector<LinExpr> args(n);
  int index = 0;
  for (PairwiseExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
    args[index++] = Visit(*i);
  Gecode::BoolVarArgs logical_args(n * (n - 1) / 2);
  index = 0;
  GecodeProblem &problem = converter_.problem_;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      Gecode::BoolExpr expr = negate ? args[i] == args[j] : args[i] != args[j];
      logical_args[index++] = Gecode::expr(problem, expr, converter_.ipl_);
    }
  }
  Gecode::BoolVar var(problem, 0, 1);
  rel(problem, negate ? Gecode::BOT_OR : Gecode::BOT_AND,
      logical_args, var, converter_.ipl_);
  return var;
}

GecodeSolver::Stop::Stop(GecodeSolver &solver)
: solver_(solver) {
  output_or_limit_ = solver.output_ || solver.time_limit_ < DBL_MAX ||
      solver.node_limit_ != ULONG_MAX || solver.fail_limit_ != ULONG_MAX;
  steady_clock::time_point start = steady_clock::now();
  double end_time_in_ticks = start.time_since_epoch().count() +
      solver.time_limit_ * steady_clock::period::den /
      steady_clock::period::num;
  end_time_ = steady_clock::time_point(steady_clock::duration(
      end_time_in_ticks >= std::numeric_limits<steady_clock::rep>::max() ?
          std::numeric_limits<steady_clock::rep>::max() :
          static_cast<steady_clock::rep>(end_time_in_ticks)));
  next_output_time_ = start + GetOutputInterval();
}

bool GecodeSolver::Stop::stop(
    const Search::Statistics &s, const Search::Options &) {
  if (solver_.interrupter()->Stop()) {
    solver_.SetStatus(600, "interrupted");
    return true;
  }
  if (!output_or_limit_) return false;
  steady_clock::time_point time = steady_clock::now();
  if (solver_.output_ && time >= next_output_time_) {
    solver_.Output("{:10} {:10} {:10}\n", s.depth, s.node, s.fail);
    next_output_time_ += GetOutputInterval();
  }
  if (time > end_time_)
    solver_.SetStatus(400, "time limit");
  else if (s.node > solver_.node_limit_)
    solver_.SetStatus(401, "node limit");
  else if (s.fail > solver_.fail_limit_)
    solver_.SetStatus(402, "fail limit");
  else
    return false;
  return true;
}

void GecodeSolver::SetBoolOption(
    const SolverOption &opt, int value, bool *option) {
  if (value != 0 && value != 1)
    throw InvalidOptionValue(opt, value);
  *option = value != 0;
}

void GecodeSolver::SetOutputFrequency(const SolverOption &opt, double value) {
  if (value <= 0)
    throw InvalidOptionValue(opt, value);
  output_frequency_ = value;
}

template <typename T>
std::string GecodeSolver::GetEnumOption(const SolverOption &opt, T *ptr) const {
  for (mp::ValueArrayRef::iterator
      i = opt.values().begin(), e = opt.values().end(); i != e; ++i) {
    if (*ptr == i->data)
      return i->value;
  }
  return fmt::format("{}", *ptr);
}

template <typename T>
void GecodeSolver::SetEnumOption(
    const SolverOption &opt, fmt::StringRef value, T *ptr) {
  for (mp::ValueArrayRef::iterator
      i = opt.values().begin(), e = opt.values().end(); i != e; ++i) {
    if (value == i->value) {
      *ptr = static_cast<T>(i->data);
      return;
    }
  }
  throw InvalidOptionValue(opt, value);
}

template <typename T, typename OptionT>
void GecodeSolver::SetNonnegativeOption(
    const SolverOption &opt, T value, OptionT *option) {
  if (value < 0)
    throw InvalidOptionValue(opt, value);
  *option = value;
}

void GecodeSolver::Output(fmt::CStringRef format, const fmt::ArgList &args) {
  if (output_count_ == 0)
    Print("{}", header_);
  output_count_ = (output_count_ + 1) % 20;
  Print(format, args);
}

GecodeSolver::GecodeSolver()
: SolverImpl<Problem>(
    "gecode", "gecode " GECODE_VERSION, 20160205, MULTIPLE_SOL),
  output_(false), output_frequency_(1), output_count_(0), solve_code_(-1),
  ipl_(Gecode::IPL_DEF),
  var_branching_(IntVarBranch::SEL_SIZE_MIN),
  val_branching_(IntValBranch::SEL_MIN),
  decay_(1),
  time_limit_(DBL_MAX), node_limit_(ULONG_MAX), fail_limit_(ULONG_MAX),
  solution_limit_(UINT_MAX),
  restart_(Gecode::RM_NONE), restart_base_(1.5), restart_scale_(250) {

  set_version("Gecode " GECODE_VERSION);

  AddSuffix("ipl", 0, suf::CON);

  set_option_header(
      "Gecode Options for AMPL\n"
      "-----------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to "
      "the AMPL option ``gecode_options``. For example::\n"
      "\n"
      "  ampl: option gecode_options 'version nodelimit=30000 "
      "val_branching=min';\n");

  AddIntOption("log:lev outlev",
      "0 or 1 (default 0): Whether to print solution log.",
      &GecodeSolver::GetOption<int, bool>,
      &GecodeSolver::SetBoolOption, &output_);

  AddDblOption("log:freq logfreq outfreq",
      "Output frequency in seconds. The value should be a positive number.",
      &GecodeSolver::GetOutputFrequency, &GecodeSolver::SetOutputFrequency);

  AddStrOption("cp:ipl ipl",
      "Propagation level for integer propagators. Possible values:\n"
      "\n"
      ".. value-table::\n",
      &GecodeSolver::GetEnumOption<Gecode::IntPropLevel>,
      &GecodeSolver::SetEnumOption<Gecode::IntPropLevel>,
      &ipl_, INT_PROP_LEVELS);

  AddStrOption("cp:var_branching var_branching",
      "Variable branching. Possible values:\n"
      "\n"
      ".. value-table::\n",
      &GecodeSolver::GetEnumOption<IntVarBranch::Select>,
      &GecodeSolver::SetEnumOption<IntVarBranch::Select>,
      &var_branching_, VAR_BRANCHINGS);

  AddStrOption("cp:val_branching val_branching",
      "Value branching. Possible values:\n"
      "\n"
      ".. value-table::\n",
      &GecodeSolver::GetEnumOption<IntValBranch::Select>,
      &GecodeSolver::SetEnumOption<IntValBranch::Select>,
      &val_branching_, VAL_BRANCHINGS);

  AddDblOption("cp:decay decay",
      "Decay factor for AFC and activity branchings. Default = 1.",
      &GecodeSolver::GetDecay, &GecodeSolver::SetDecay);

  AddDblOption("tech:threads threads",
      "The number of parallel threads to use. Assume that your computer "
      "has m processing units and that the value for threads is n.\n"
      "\n"
      "* If n = 0, then m threads are used (as many as available "
      "  processing units).\n"
      "* If n >= 1, then n threads are used (absolute number of "
      "  threads to be used).\n"
      "* If n <= -1, then m + n threads are used (absolute number "
      "  of processing units not to be used). For example, when "
      "  n = -6 and m = 8, then 2 threads are used.\n"
      "* If 0 < n < 1, then n * m threads are used (relative number "
      "  of processing units to be used). For example, when n = 0.5 "
      "  and m = 8, then 4 threads are used.\n"
      "* If -1 < n < 0, then (1 + n) * m threads are used (relative "
      "  number of processing units not to be used). For example, "
      "  when n = -0.25 and m = 8, then 6 threads are used.\n"
      "\n"
      "All values are rounded and at least one thread is used.\n",
      &GecodeSolver::GetOption<double, double>,
      &GecodeSolver::DoSetDblOption, &options_.threads);

  AddIntOption("cp:c_d c_d", "Commit recomputation distance.",
      &GecodeSolver::GetOption<int, unsigned>,
      &GecodeSolver::SetNonnegativeOption<int, unsigned>, &options_.c_d);
  AddIntOption("cp:a_d a_d", "Adaptive recomputation distance.",
      &GecodeSolver::GetOption<int, unsigned>,
      &GecodeSolver::SetNonnegativeOption<int, unsigned>, &options_.a_d);

  AddDblOption("lim:time timelim timelimit", "Time limit in seconds.",
      &GecodeSolver::GetOption<double, double>,
      &GecodeSolver::SetNonnegativeOption<double, double>, &time_limit_);
  AddIntOption("lim:nodes nodelim nodelimit", "Node limit.",
      &GecodeSolver::GetOption<int, unsigned long>,
      &GecodeSolver::SetNonnegativeOption<int, unsigned long>, &node_limit_);
  AddIntOption("lim:fail faillim faillimit", "Fail limit.",
      &GecodeSolver::GetOption<int, unsigned long>,
      &GecodeSolver::SetNonnegativeOption<int, unsigned long>, &fail_limit_);

  AddStrOption("cp:restart restart",
      "Restart sequence type. Possible values:\n"
      "\n"
      ".. value-table::\n",
      &GecodeSolver::GetEnumOption<Gecode::RestartMode>,
      &GecodeSolver::SetEnumOption<Gecode::RestartMode>,
      &restart_, RESTART_MODES);

  AddDblOption("cp:restart_base restart_base",
      "Base for geometric restart sequence. Default = 1.5.",
      &GecodeSolver::GetOption<double, double>,
      &GecodeSolver::DoSetDblOption, &restart_base_);

  AddIntOption("cp:restart_scale restart_scale",
      "Scale factor for restart sequence. Default = 250.",
      &GecodeSolver::GetOption<int, unsigned long>,
      &GecodeSolver::SetNonnegativeOption<int, unsigned long>, &restart_scale_);

  AddIntOption("lim:sol solutionlimit",
      "Limit on the number of feasible solutions found before terminating "
      "a search. Leaving the solution limit unspecified will make the "
      "optimizer search for an optimal solution if there is an objective "
      "function or for a feasible solution otherwise.",
      &GecodeSolver::GetOption<int, unsigned>,
      &GecodeSolver::SetNonnegativeOption<int, unsigned>, &solution_limit_);
}

void GetSolution(GecodeProblem &gecode_problem, std::vector<double> &solution) {
  IntVarArray &vars = gecode_problem.vars();
  int num_vars = static_cast<int>(solution.size());
  for (int j = 0; j < num_vars; ++j)
    solution[j] = vars[j].val();
}

template<template<typename, template<typename> class> class Meta>
GecodeSolver::ProblemPtr GecodeSolver::Search(
    Problem &p, GecodeProblem &problem,
    Search::Statistics &stats, SolutionHandler &sh) {
  ProblemPtr final_problem;
  unsigned solution_limit = solution_limit_;
  unsigned num_solutions = 0;
  if (problem.has_obj()) {
    Meta<GecodeProblem, Gecode::BAB> engine(&problem, options_);
    while (GecodeProblem *next = engine.next()) {
      if (output_)
        Output("{:46}\n", next->obj().val());
      final_problem.reset(next);
      if (++num_solutions >= solution_limit_) {
        SetStatus(403, "solution limit");
        break;
      }
    }
    stats = engine.statistics();
  } else {
    if (solution_limit == UINT_MAX)
      solution_limit = 1;
    Meta<GecodeProblem, Gecode::DFS> engine(&problem, options_);
    std::vector<double> solution;
    bool multiple_sol = need_multiple_solutions();
    if (multiple_sol)
      solution.resize(p.num_vars());
    std::string feasible_sol_message =
        fmt::format("{}: feasible solution", long_name());
    while (GecodeProblem *next = engine.next()) {
      final_problem.reset(next);
      if (multiple_sol) {
        GetSolution(*final_problem, solution);
        sh.HandleFeasibleSolution(feasible_sol_message, solution.data(), 0, 0);
      }
      if (++num_solutions >= solution_limit_)
        break;
    }
    stats = engine.statistics();
  }
  return final_problem;
}

void GecodeSolver::Solve(Problem &p, SolutionHandler &sh) {
  steady_clock::time_point time = steady_clock::now();

  SetStatus(-1, "");

  // Set up an optimization problem in Gecode.
  MPToGecodeConverter converter(p.num_vars(), ipl_);
  converter.Convert(p);

  // Post branching.
  GecodeProblem &gecode_problem = converter.problem();
  IntVarBranch var_branch;
  switch (var_branching_) {
  case IntVarBranch::SEL_RND:
    var_branch = IntVarBranch(Gecode::Rnd(0));
    break;
  case IntVarBranch::SEL_AFC_MIN:
  case IntVarBranch::SEL_AFC_MAX:
  case IntVarBranch::SEL_ACTION_MIN:
  case IntVarBranch::SEL_ACTION_MAX:
  case IntVarBranch::SEL_AFC_SIZE_MIN:
  case IntVarBranch::SEL_AFC_SIZE_MAX:
  case IntVarBranch::SEL_ACTION_SIZE_MIN:
  case IntVarBranch::SEL_ACTION_SIZE_MAX:
    var_branch = IntVarBranch(var_branching_, decay_, 0);
    break;
  default:
    var_branch = IntVarBranch(var_branching_, 0);
    break;
  }
  IntValBranch val_branch = val_branching_ == IntValBranch::SEL_RND ?
      IntValBranch(Gecode::Rnd(0)) : IntValBranch(val_branching_);
  branch(gecode_problem, gecode_problem.vars(), var_branch, val_branch);

  Stop stop(*this);
  options_.stop = &stop;

  double setup_time = GetTimeAndReset(time);

  // Solve the problem.
  double obj_val = std::numeric_limits<double>::quiet_NaN();
  bool has_obj = p.num_objs() > 0;
  Search::Statistics stats;
  header_ = fmt::format("{:>10} {:>10} {:>10} {:>13}\n",
    "Max Depth", "Nodes", "Fails", (has_obj ? "Best Obj" : ""));
  output_count_ = 0;
  GecodeSolver::ProblemPtr solution;
  if (restart_ != Gecode::RM_NONE) {
    options_.cutoff = Gecode::Driver::createCutoff(*this);
    solution = Search<Gecode::RBS>(p, gecode_problem, stats, sh);
  } else {
    solution =
        Search<Gecode::Driver::EngineToMeta>(p, gecode_problem, stats, sh);
  }

  // Convert solution status.
  if (solution.get()) {
    if (!has_obj) {
      // If the problem has no objectives and there is a solution, report
      // it as solved even if some limit was reached.
      solve_code_ = sol::SOLVED;
      status_ = "feasible solution";
    } else if (solve_code_ == -1) {
      solve_code_ = sol::SOLVED;
      obj_val = solution->obj().val();
      status_ = "optimal solution";
    }
  } else if (solve_code_ == -1) {
    solve_code_ = sol::INFEASIBLE;
    status_ = "infeasible problem";
  }

  std::vector<double> final_solution;
  if (solution.get()) {
    final_solution.resize(p.num_vars());
    GetSolution(*solution, final_solution);
  }

  double solution_time = GetTimeAndReset(time);

  fmt::MemoryWriter w;
  w.write("{}: {}\n", long_name(), status_);
  w.write("{} nodes, {} fails", stats.node, stats.fail);
  if (has_obj && solution.get())
    w.write(", objective {}", FormatObjValue(obj_val));
  sh.HandleSolution(solve_code_, w.c_str(),
      final_solution.empty() ? 0 : final_solution.data(), 0, obj_val);

  double output_time = GetTimeAndReset(time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n",
          setup_time, solution_time, output_time);
  }
}

SolverPtr create_gecode(const char *) { return SolverPtr(new GecodeSolver()); }
}  // namespace mp

/*
 AMPL solver interface to LocalSolver.

 Copyright (C) 2014 AMPL Optimization Inc

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

#include "localsolver/localsolver.h"

#include <cmath>
#include <limits>
#include "mp/clock.h"
#include "mp/safeint.h"

#ifndef MP_DRIVER_DATE
# include "locsol_date.h"
# define MP_DRIVER_DATE YYYYMMDD
#endif

namespace {

enum { TERSE_VERBOSITY = -1 };

const mp::OptionValueInfo VERBOSITIES[] = {
  {"terse",  "Terse verbosity.", TERSE_VERBOSITY},
  {"quiet",  "All the traces are disabled (default).", 0},
  {"normal", "Normal verbosity.", 1},
  {"detailed",
   "Detailed verbosity. Displays extended statistics on the model.", 2}
};

// Returns the value of an expression.
inline double GetValue(localsolver::LSExpression e) {
  return e.isDouble() ? e.getDoubleValue() : e.getValue();
}

bool StopSolver(void *data) {
  localsolver::LocalSolver* solver =
      static_cast<localsolver::LocalSolver*>(data);
  if (solver->getState() != localsolver::S_Running)
    return false;
  solver->stop();
  return true;
}
}  // namespace

namespace mp {

// LocalSolver doesn't allow infinite bounds for variables, so use the max
// double value instead.
const double LSProblemBuilder::LS_INF = std::numeric_limits<double>::max();

LSProblemBuilder::HyperbolicTerms
    LSProblemBuilder::MakeHyperbolicTerms(ls::LSExpression arg) {
  HyperbolicTerms terms;
  terms.exp_x = model_.createExpression(ls::O_Exp, arg);
  terms.exp_minus_x = model_.createExpression(ls::O_Exp, Negate(arg));
  return terms;
}

LSProblemBuilder::LSProblemBuilder(LocalSolver &s)
  : model_(solver_.getModel()), num_cons_(0), pl_bigm_(s.pl_bigm()) {
  solver_.getParam().setVerbosity(0);
}

void LSProblemBuilder::SetInfo(const ProblemInfo &info) {
  vars_.reserve(info.num_vars);
  objs_.reserve(info.num_objs);
  cons_.reserve(info.num_algebraic_cons);
}

void LSProblemBuilder::AddVar(double lb, double ub, var::Type type) {
  vars_.push_back(ls::LSExpression());
  ls::LSExpression &var = vars_.back();
  double inf = std::numeric_limits<double>::infinity();
  if (type == var::CONTINUOUS) {
    var = model_.createExpression(
          ls::O_Float, lb == -inf ? -LS_INF : lb, ub == inf ? LS_INF : ub);
  } else if (lb == 0 && ub == 1) {
    var = model_.createExpression(ls::O_Bool);
  } else {
    ls::lsint int_lb = lb == -inf ?
          std::numeric_limits<int>::min() : ConvertToInt(lb);
    ls::lsint int_ub = ub == inf ?
          std::numeric_limits<int>::max() : ConvertToInt(ub);
    var = model_.createExpression(ls::O_Int, int_lb, int_ub);
  }
}

ls::LSExpression LSProblemBuilder::MakeUnary(
    expr::Kind kind, ls::LSExpression arg) {
  ls::LSOperator op = ls::O_Bool;
  switch (kind) {
  case expr::FLOOR:
    op = ls::O_Floor;
    break;
  case expr::CEIL:
    op = ls::O_Ceil;
    break;
  case expr::ABS:
    op = ls::O_Abs;
    break;
  case expr::MINUS:
    return Negate(arg);
  case expr::TANH: {
    HyperbolicTerms terms = MakeHyperbolicTerms(arg);
    return MakeBinary(ls::O_Div,
                      MakeBinary(ls::O_Sub, terms.exp_x, terms.exp_minus_x),
                      MakeBinary(ls::O_Sum, terms.exp_x, terms.exp_minus_x));
  }
  case expr::TAN:
    op = ls::O_Tan;
    break;
  case expr::SQRT:
    op = ls::O_Sqrt;
    break;
  case expr::SINH: {
    HyperbolicTerms terms = MakeHyperbolicTerms(arg);
    return Half(MakeBinary(ls::O_Sub, terms.exp_x, terms.exp_minus_x));
  }
  case expr::SIN:
    op = ls::O_Sin;
    break;
  case expr::LOG10:
    return MakeBinary(ls::O_Div,
                      model_.createExpression(ls::O_Log, arg), std::log(10.0));
  case expr::LOG:
    op = ls::O_Log;
    break;
  case expr::EXP:
    op = ls::O_Exp;
    break;
  case expr::COSH:{
    HyperbolicTerms terms = MakeHyperbolicTerms(arg);
    return Half(MakeBinary(ls::O_Sum, terms.exp_x, terms.exp_minus_x));
  }
  case expr::COS:
    return model_.createExpression(ls::O_Cos, arg);
  case expr::ATANH:
    arg = MakeBinary(ls::O_Div, Plus1(arg),                    // (1 + x) /
                     MakeBinary(ls::O_Sub, AsLSInt(1), arg));  // (1 - x)
    return Half(model_.createExpression(ls::O_Log, arg));
  case expr::ASINH: {
    ls::LSExpression arg2 = model_.createExpression(
          ls::O_Sqrt, Plus1(MakeBinary(ls::O_Pow, arg, AsLSInt(2))));
    return model_.createExpression(ls::O_Log, MakeBinary(ls::O_Sum, arg, arg2));
  }
  case expr::ACOSH: {
    ls::LSExpression x_minus_1 = MakeBinary(ls::O_Sub, arg, AsLSInt(1));
    ls::LSExpression arg2 = MakeBinary(
          ls::O_Prod, model_.createExpression(ls::O_Sqrt, Plus1(arg)),
          model_.createExpression(ls::O_Sqrt, x_minus_1));
    return model_.createExpression(ls::O_Log, MakeBinary(ls::O_Sum, arg, arg2));
  }
  case expr::POW2:
    return MakeBinary(ls::O_Pow, arg, AsLSInt(2));
  case expr::ATAN: case expr::ASIN: case expr::ACOS:
    // LocalSolver doesn't support these expressions.
    // Fall through.
  default:
    return Base::MakeUnary(kind, arg);
  }
  return model_.createExpression(op, arg);
}

ls::LSExpression LSProblemBuilder::MakeBinary(
    expr::Kind kind, ls::LSExpression lhs, ls::LSExpression rhs) {
  ls::LSOperator op = ls::O_Bool;
  switch (kind) {
  case expr::ADD:
    op = ls::O_Sum;
    break;
  case expr::SUB:
    op = ls::O_Sub;
    break;
  case expr::MUL:
    op = ls::O_Prod;
    break;
  case expr::DIV:
    op = ls::O_Div;
    break;
  case expr::TRUNC_DIV:
    return TruncDiv(lhs, rhs);
  case expr::MOD:
    op = ls::O_Mod;
    break;
  case expr::POW:
  case expr::POW_CONST_BASE:
  case expr::POW_CONST_EXP:
    op = ls::O_Pow;
    break;
  case expr::LESS:
    return MakeBinary(ls::O_Max, MakeBinary(ls::O_Sub, lhs, rhs), AsLSInt(0));
  case expr::ROUND:
    RequireZero(rhs, "round");
    return model_.createExpression(ls::O_Round, lhs);
  case expr::TRUNC:
    RequireZero(rhs, "trunc");
    return MakeIf(MakeBinary(ls::O_Geq, lhs, AsLSInt(0)),
                  model_.createExpression(ls::O_Floor, lhs),
                  model_.createExpression(ls::O_Ceil, lhs));
  case expr::PRECISION:
  case expr::ATAN2:
    // LocalSolver doesn't support these functions.
    // Fall through.
  default:
    return Base::MakeBinary(kind, lhs, rhs);
  }
  return MakeBinary(op, lhs, rhs);
}

ls::LSExpression LSProblemBuilder::EndPLTerm(
    const PLTermBuilder &builder, ls::LSExpression arg) {
  std::size_t num_breakpoints = static_cast<int>(builder.breakpoints.size());
  assert(builder.slopes.size() == num_breakpoints + 1);
  // Find the starting index of nonnegative breakpoints.
  std::size_t nonnegative_start = 0;
  for (; nonnegative_start < num_breakpoints; ++nonnegative_start) {
    if (builder.breakpoints[nonnegative_start] >= 0)
      break;
  }
  // Process nonnegative breakpoints.
  struct Converter {
    double lb, ub;
    double prev_x, prev_y;
    ls::LSExpression xs, ys;

    // Swaps elements i and j in array.
    static void swap(ls::LSExpression array, int i, int j) {
      ls::LSExpression temp = array.getOperand(i);
      array.setOperand(i, array.getOperand(j));
      array.setOperand(j, temp);
    }

    explicit Converter(ls::LSModel model, double lb, double ub)
      : lb(lb), ub(ub), prev_x(0), prev_y(0),
        xs(model.createExpression(ls::O_Array)),
        ys(model.createExpression(ls::O_Array)) {}

    void Convert(double x, double slope) {
      double y = slope * (x - prev_x) + prev_y;
      if (x >= lb && x <= ub) {
        xs.addOperand(x);
        ys.addOperand(y);
      }
      prev_x = x;
      prev_y = y;
    }

    // Reverse the order of points as xs must be nondecreasing.
    void Reverse() {
      for (int i = 0, j = xs.getNbOperands() - 1; i < j; ++i, --j) {
        swap(xs, i, j);
        swap(ys, i, j);
      }
    }
  };
  double lb = arg.getOperand(0).getDoubleValue();
  double ub = arg.getOperand(1).getDoubleValue();
  if (lb == -LS_INF || ub == LS_INF) {
    // If the argument is unbounded, impose bounds ourselves because
    // LocalSolver doesn't support unbounded piecewise-linear terms.
    double max_abs_bp = 0;
    for (std::size_t i = 0; i < num_breakpoints; ++i)
      max_abs_bp = std::max(max_abs_bp, std::abs(builder.breakpoints[i]));
    double bound = std::max(pl_bigm_, 2 * max_abs_bp);
    if (lb > -LS_INF) bound = std::max(bound, std::abs(lb));
    if (ub <  LS_INF) bound = std::max(bound, std::abs(ub));
    if (lb == -LS_INF) lb = -bound;
    if (ub ==  LS_INF) ub =  bound;
  }
  Converter converter(model_, lb, ub);
  // Convert negative breakpoints.
  double slope = builder.slopes[nonnegative_start];
  for (int i = val(SafeInt<int>(nonnegative_start)) - 1; i >= 0; --i) {
    double x = builder.breakpoints[i];
    if (x < lb) break;
    converter.Convert(x, slope);
    slope = builder.slopes[i];
  }
  if (ub < converter.prev_x)
    converter.Convert(ub, slope);
  if (lb < converter.prev_x)
    converter.Convert(lb, slope);
  converter.Reverse();
  // Convert nonnegative breakpoints.
  converter.prev_x = converter.prev_y = 0;
  slope = builder.slopes[nonnegative_start];
  for (std::size_t i = nonnegative_start; i < num_breakpoints; ++i) {
    double x = builder.breakpoints[i];
    if (x > ub) break;
    converter.Convert(x, slope);
    slope = builder.slopes[i + 1];
  }
  if (lb > converter.prev_x)
    converter.Convert(lb, slope);
  if (ub > converter.prev_x)
    converter.Convert(ub, slope);
  if (converter.xs.getNbOperands() == 1) {
    // LocalSolver requires at least two points.
    converter.xs.addOperand(converter.xs.getOperand(0));
    converter.ys.addOperand(converter.ys.getOperand(0));
  }
  // For compatibility with pre-5.0 versions:
  ls::LSOperator piecewise = static_cast<ls::LSOperator>(ls::O_Int + 1);
  return model_.createExpression(
        piecewise, converter.xs, converter.ys, arg);
}

LSProblemBuilder::ExprBuilder LSProblemBuilder::BeginIterated(
    expr::Kind kind, int num_args) {
  ls::LSOperator op = ls::O_Min;
  if (kind == expr::MAX)
    op = ls::O_Max;
  else if (kind != expr::MIN)
    Base::BeginIterated(kind, num_args);
  return ExprBuilder(model_.createExpression(op));
}

ls::LSExpression LSProblemBuilder::MakeBinaryLogical(
    expr::Kind kind, ls::LSExpression lhs, ls::LSExpression rhs) {
  ls::LSOperator op = ls::O_Bool;
  switch (kind) {
  case expr::OR:
    op = ls::O_Or;
    break;
  case expr::AND:
    op = ls::O_And;
    break;
  case expr::IFF:
    op = ls::O_Eq;
    break;
  default:
    return Base::MakeBinaryLogical(kind, lhs, rhs);
  }
  return MakeBinary(op, lhs, rhs);
}

ls::LSExpression LSProblemBuilder::MakeRelational(
    expr::Kind kind, ls::LSExpression lhs, ls::LSExpression rhs) {
  ls::LSOperator op = ls::O_Bool;
  switch (kind) {
  case expr::LT:
    op = ls::O_Lt;
    break;
  case expr::LE:
    op = ls::O_Leq;
    break;
  case expr::EQ:
    op = ls::O_Eq;
    break;
  case expr::GE:
    op = ls::O_Geq;
    break;
  case expr::GT:
    op = ls::O_Gt;
    break;
  case expr::NE:
    if (lhs.getOperator() == ls::O_Bool && IsConst(rhs, 0))
      return lhs;
    op = ls::O_Neq;
    break;
  default:
    return Base::MakeRelational(kind, lhs, rhs);
  }
  return MakeBinary(op, lhs, rhs);
}

ls::LSExpression LSProblemBuilder::MakeLogicalCount(
    expr::Kind kind, ls::LSExpression lhs, ls::LSExpression rhs) {
  ls::LSOperator op = ls::O_Bool;
  switch (kind) {
  case expr::ATLEAST:
    op = ls::O_Leq;
    break;
  case expr::ATMOST:
    op = ls::O_Geq;
    break;
  case expr::EXACTLY:
    op = ls::O_Eq;
    break;
  case expr::NOT_ATLEAST:
    op = ls::O_Gt;
    break;
  case expr::NOT_ATMOST:
    op = ls::O_Lt;
    break;
  case expr::NOT_EXACTLY:
    op = ls::O_Neq;
    break;
  default:
    return Base::MakeLogicalCount(kind, lhs, rhs);
  }
  return MakeBinary(op, lhs, rhs);
}

LSProblemBuilder::ExprBuilder LSProblemBuilder::BeginIteratedLogical(
    expr::Kind kind, int num_args) {
  ls::LSOperator op = ls::O_Or;
  if (kind == expr::FORALL)
    op = ls::O_And;
  else if (kind != expr::EXISTS)
    Base::BeginIteratedLogical(kind, num_args);
  return ExprBuilder(model_.createExpression(op));
}

ls::LSExpression LSProblemBuilder::EndPairwise(PairwiseExprBuilder builder) {
  std::vector<ls::LSExpression> &args = builder.args;
  ls::LSOperator logical_op = ls::O_And, comparison_op = ls::O_Neq;
  if (builder.kind == expr::NOT_ALLDIFF) {
    logical_op = ls::O_Or;
    comparison_op = ls::O_Eq;
  }
  ls::LSExpression alldiff = model_.createExpression(logical_op);
  for (std::size_t i = 0, n = args.size(); i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j)
      alldiff.addOperand(MakeBinary(comparison_op, args[i], args[j]));
  }
  return alldiff;
}

std::string LocalSolver::GetVerbosity(const SolverOption &opt) const {
  int value = options_[VERBOSITY];
  for (mp::ValueArrayRef::iterator
      i = opt.values().begin(), e = opt.values().end(); i != e; ++i) {
    if (i->data == value)
      return i->value;
  }
  return fmt::format("{}", value);
}

void LocalSolver::SetVerbosity(const SolverOption &opt, fmt::StringRef value) {
  char *end = 0;
  // Copy the value adding a terminating null for strtol.
  fmt::MemoryWriter writer;
  writer << value;
  const char *str = writer.c_str();
  long intval = std::strtol(str, &end, 0);
  if (!*end) {
    std::size_t i = 0, size = sizeof(VERBOSITIES) / sizeof(*VERBOSITIES);
    for (; i < size; ++i) {
      if (intval == verbosities_[i].data)
        break;
    }
    if (i == size)
      throw InvalidOptionValue(opt, str);
    options_[VERBOSITY] = intval;
    return;
  }
  for (mp::ValueArrayRef::iterator
      i = opt.values().begin(), e = opt.values().end(); i != e; ++i) {
    if (std::strcmp(i->value, str) == 0) {
      options_[VERBOSITY] = static_cast<int>(i->data);
      return;
    }
  }
  throw InvalidOptionValue(opt, str);
}

LocalSolver::LocalSolver()
  : SolverImpl<LSProblemBuilder>("locsol", 0, MP_DRIVER_DATE, MULTIPLE_OBJ) {
  options_[SEED] = 0;
  options_[THREADS] = 2;
  options_[ANNEALING_LEVEL] = 1;
  options_[VERBOSITY] = 0;
  options_[TIME_BETWEEN_DISPLAYS] = 1;
  options_[TIMELIMIT] = 10;
  iterlimit_ = std::numeric_limits<fmt::LongLong>::max();
  pl_bigm_ = 1e6;

  int ls_version = localsolver::LSVersion::getMajorVersionNumber();
  std::string version = fmt::format("{}.{}",
      ls_version, localsolver::LSVersion::getMinorVersionNumber());
  set_long_name("LocalSolver " + version);
  set_version("LocalSolver " + version);

  set_option_header(
      "LocalSolver Options for AMPL\n"
      "----------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to "
      "the AMPL option ``locsol_options``. For example::\n"
      "\n"
      "  ampl: option locsol_options 'version timelimit=30';\n");

  AddIntOption("seed",
      "Seed of the pseudo-random number generator used by the solver. "
      "Default = 0.",
      &LocalSolver::DoGetIntOption, &LocalSolver::DoSetIntOption<INT_MAX>,
      OptionInfo(SEED, 0));

  AddIntOption("threads",
      "Number of threads used to parallelize the search. Default = 2.",
      &LocalSolver::DoGetIntOption, &LocalSolver::DoSetIntOption<1024>,
      OptionInfo(THREADS, 1));

  AddIntOption("annealing_level",
      "Simulated annealing level. Default = 1.",
      &LocalSolver::DoGetIntOption, &LocalSolver::DoSetIntOption<9>,
      OptionInfo(ANNEALING_LEVEL, 0));

  FMT_STATIC_ASSERT(sizeof(VERBOSITIES) == sizeof(verbosities_),
                    "size mismatch");
  std::size_t num_verbosities = sizeof(VERBOSITIES) / sizeof(*VERBOSITIES);
  std::copy(VERBOSITIES, VERBOSITIES + num_verbosities, verbosities_);
  if (ls_version < 5)
    verbosities_[num_verbosities - 1].data = 10;
  AddStrOption("verbosity",
      "Verbosity level of the solver. Possible values:\n"
      "\n"
      ".. value-table::\n"
      "\n",
      &LocalSolver::GetVerbosity, &LocalSolver::SetVerbosity, verbosities_);

  AddIntOption("time_between_displays",
      "Time in seconds between two consecutive displays in console while "
      "the solver is running. Default = 1.",
      &LocalSolver::DoGetIntOption, &LocalSolver::DoSetIntOption<65535>,
      OptionInfo(TIME_BETWEEN_DISPLAYS, 1));

  AddStrOption("logfile",
      "Path of the LocalSolver log file. Default = \"\" (no log file).",
      &LocalSolver::GetLogFile, &LocalSolver::SetLogFile);

  AddIntOption("timelimit",
      "Time limit in seconds (nonnegative integer). Default = 10.",
      &LocalSolver::DoGetIntOption, &LocalSolver::DoSetIntOption<INT_MAX>,
      OptionInfo(TIMELIMIT, ls_version < 5 ? 1 : 0));

  AddIntOption("iterlimit",
      "Iteration limit (nonnegative integer). "
      "Default = largest positive integer.",
      &LocalSolver::GetIterLimit, &LocalSolver::SetIterLimit);

  AddDblOption("pl_bigm",
    "The artificial bound used for unbounded variables in piecewise-linear "
    "terms. Default = 1e6.", &LocalSolver::GetPLBigM, &LocalSolver::SetPLBigM);

  AddStrOption("envfile",
    "Path to the file where to export LocalSolver environment. "
    "Default = \"\" (no file)",
    &LocalSolver::GetEnvFile, &LocalSolver::SetEnvFile);
}

void LocalSolver::Solve(ProblemBuilder &builder, SolutionHandler &sh) {
  steady_clock::time_point time = steady_clock::now();

  ls::LocalSolver &solver = builder.solver();
  ls::LSModel model = solver.getModel();

  // LocalSolver requires at least one objective - create a dummy one.
  if (model.getNbObjectives() == 0)
    model.addObjective(model.createConstant(AsLSInt(0)), ls::OD_Minimize);
  model.close();

  // Set initial values.
  ls::LSExpression *vars = builder.vars();
  if (const double *initial_values = builder.initial_values()) {
    for (int i = 0; i < builder.num_vars(); ++i) {
      double value = initial_values[i];
      auto var = vars[i];
      if (var.getOperator() == ls::O_Float)
        var.setValue(value);
      else
        var.setValue(static_cast<ls::lsint>(value));
    }
  }

  int verbosity = options_[VERBOSITY];
  bool custom_output = verbosity == TERSE_VERBOSITY;

  // Set options. LS requires this to be done after the model is closed.
  ls::LSParam param = solver.getParam();
  param.setSeed(options_[SEED]);
  param.setNbThreads(options_[THREADS]);
  param.setAnnealingLevel(options_[ANNEALING_LEVEL]);
  param.setVerbosity(custom_output ? 0 : verbosity);
  param.setTimeBetweenDisplays(options_[TIME_BETWEEN_DISPLAYS]);
  if (!logfile_.empty())
    param.setLogFile(logfile_);
  ls::LSPhase phase = solver.createPhase();
  phase.setTimeLimit(options_[TIMELIMIT]);
  phase.setIterationLimit(iterlimit_);
  interrupter()->SetHandler(StopSolver, &solver);

  typedef LSProblemBuilder::Bound Bound;
  const std::vector<Bound> &obj_bounds = builder.obj_bounds();
  for (std::vector<Bound>::const_iterator
       i = obj_bounds.begin(), e = obj_bounds.end(); i != e; ++i) {
    int index = i->index();
    if (model.getObjective(index).isDouble())
      param.setDoubleObjectiveBound(index, i->dbl_value());
    else
      param.setObjectiveBound(index, i->int_value());
  }

  struct Callback {
    LocalSolver &solver;
    ls::LocalSolver &ls_solver;
    ls::LSExpression obj;
    bool custom_output;
    bool print_header;

    // The best objective value found so far, multiplied by obj_sign.
    double adjusted_obj_value;
    int obj_sign;
    int seconds_passed;
    int seconds_to_best_obj;
    fmt::LongLong iters_to_best_obj;

    Callback(LocalSolver &s, ls::LocalSolver &ls_solver, bool custom_output)
      : solver(s), ls_solver(ls_solver),
        custom_output(custom_output), print_header(true),
        adjusted_obj_value(std::numeric_limits<double>::infinity()),
        seconds_passed(0), seconds_to_best_obj(0), iters_to_best_obj(0) {
      ls::LSModel model = ls_solver.getModel();
      obj_sign = model.getObjectiveDirection(0) == ls::OD_Minimize ? 1 : -1;
      obj = model.getObjective(0);
    }

    static void Call(ls::LSCallbackType, void *data) {
      static_cast<Callback*>(data)->Call();
    }

    void Call() {
      ++seconds_passed;
      double obj_value = obj.isDouble() ? obj.getDoubleValue() : obj.getValue();
      ls::LSStatistics stats = ls_solver.getStatistics();
      if (custom_output) {
        if (print_header) {
          print_header = false;
          solver.Print(
                "\n"
                "                    |                Moves               |\n"
                "    Time       Iter |    Total  Infeas  Accepted  Improv |"
                "    Obj\n");
        }
        solver.Print("{:7}s {:10} {:10}  {:5.1f}%  {:5.1f}%{:10}   {:6}\n",
                     seconds_passed, stats.getNbIterations(),
                     stats.getNbMoves(), stats.getPercentInfeasibleMoves(),
                     stats.getPercentAcceptedMoves(),
                     stats.getNbImprovingMoves(), obj_value);
      }
      if (obj_sign * obj_value < adjusted_obj_value) {
        seconds_to_best_obj = seconds_passed;
        iters_to_best_obj = stats.getNbIterations();
        adjusted_obj_value = obj_sign * obj_value;
      }
    }
  } callback(*this, solver, custom_output);
  solver.addCallback(ls::CT_Ticked, &Callback::Call, &callback);

  if (!envfile_.empty())
    solver.saveEnvironment(envfile_);

  double setup_time = GetTimeAndReset(time);

  // Solve the problem.
  DoSolve(solver);

  // Convert solution status.
  int solve_code = sol::UNKNOWN;
  ls::LSSolutionStatus ls_status = solver.getSolution().getStatus();
  const char *status = "unknown";
  if (interrupter()->Stop()) {
    solve_code = sol::INTERRUPTED;
    status = "interrupted";
  } else {
    switch (ls_status) {
    case ls::SS_Inconsistent:
      solve_code = sol::INFEASIBLE;
      status = "infeasible problem";
      break;
    case ls::SS_Infeasible:
      // Solution is infeasible, but problem may be feasible.
      // This can only happen if stopped by a limit.
      solve_code = sol::LIMIT;
      status = "limit";
      break;
    case ls::SS_Feasible:
      solve_code = sol::UNCERTAIN;
      status = "feasible solution";
      break;
    case ls::SS_Optimal:
      solve_code = sol::SOLVED;
      status = builder.num_objs() > 0 ?
            "optimal solution" : "feasible solution";
      break;
    default:
      solve_code = sol::FAILURE;
      status = "unknown solution status";
      break;
    }
  }

  fmt::MemoryWriter w;
  w.write("{}: {}", long_name(), status);
  double obj_val = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> solution;
  if (ls_status == ls::SS_Optimal || ls_status == ls::SS_Feasible) {
    if (builder.num_objs() != 0) {
      obj_val = GetValue(model.getObjective(0));
      w.write("; objective {}", FormatObjValue(obj_val));
    }
    int num_vars = builder.num_vars();
    solution.resize(num_vars);
    for (int i = 0; i < num_vars; ++i)
      solution[i] = GetValue(vars[i]);
  }
  w << "\n";
  double solution_time = GetTimeAndReset(time);

  w.write("Best solution found at {} second(s) and {} iteration(s)\n",
          callback.seconds_to_best_obj, callback.iters_to_best_obj);

  ls::LSStatistics stats = solver.getStatistics();
  if (stats.getRunningTime() >= options_[TIMELIMIT])
    w.write("Stopped at time limit of {} second(s)\n", options_[TIMELIMIT]);
  w.write("{}", stats.toString());
  sh.HandleSolution(solve_code, w.c_str(),
                    solution.empty() ? 0 : solution.data(), 0, obj_val);
  double output_time = GetTimeAndReset(time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n",
          setup_time, solution_time, output_time);
  }
}

SolverPtr create_localsolver(const char *) {
  return SolverPtr(new LocalSolver());
}
}  // namespace mp

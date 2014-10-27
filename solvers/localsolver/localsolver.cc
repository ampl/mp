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

namespace {
// Returns the value of an expression.
inline double GetValue(localsolver::LSExpression e) {
  return e.isDouble() ? e.getDoubleValue() : e.getValue();
}

inline localsolver::lsint ConvertToInt(double value) {
  localsolver::lsint int_value = value;
  if (int_value != value)
    throw mp::Error("value {} can't be represented as int", value);
  return int_value;
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

LSProblemBuilder::HyperbolicTerms
    LSProblemBuilder::MakeHyperbolicTerms(ls::LSExpression arg) {
  HyperbolicTerms terms;
  terms.exp_x = model_.createExpression(ls::O_Exp, arg);
  terms.exp_minus_x = model_.createExpression(ls::O_Exp, Negate(arg));
  return terms;
}

LSProblemBuilder::LSProblemBuilder(LocalSolver &)
  : model_(solver_.getModel()), num_objs_(0), num_cons_(0) {
  solver_.getParam().setVerbosity(0);
}

void LSProblemBuilder::SetInfo(const ProblemInfo &info) {
  vars_.reserve(info.num_vars);
}

void LSProblemBuilder::AddVar(double lb, double ub, var::Type type) {
  vars_.push_back(ls::LSExpression());
  ls::LSExpression &var = vars_.back();
  double inf = std::numeric_limits<double>::infinity();
  if (type == var::CONTINUOUS) {
    // LocalSolver doesn't allow infinite bounds, so use min an max double
    // values instead.
    if (lb == -inf)
      lb = std::numeric_limits<double>::min();
    if (ub == inf)
      ub = std::numeric_limits<double>::max();
    var = model_.createExpression(ls::O_Float, lb, ub);
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

LSProblemBuilder::LinearObjBuilder
    LSProblemBuilder::AddObj(obj::Type type, ls::LSExpression expr, int) {
  ++num_objs_;
  ls::LSObjectiveDirection dir =
      type == obj::MIN ? ls::OD_Minimize : ls::OD_Maximize;
  ls::LSExpression sum = model_.createExpression(ls::O_Sum);
  // TODO: multiobjective optimization support
  model_.addObjective(sum, dir);
  return LinearObjBuilder(*this, expr, sum);
}

LSProblemBuilder::LinearConBuilder
    LSProblemBuilder::AddCon(ls::LSExpression expr, double lb, double ub, int) {
  ++num_cons_;
  double inf = std::numeric_limits<double>::infinity();
  ls::LSExpression sum = model_.createExpression(ls::O_Sum);
  if (lb <= -inf) {
    model_.addConstraint(model_.createExpression(ls::O_Leq, sum, ub));
  } else if (ub >= inf) {
    model_.addConstraint(model_.createExpression(ls::O_Geq, sum, lb));
  } else if (lb == ub) {
    model_.addConstraint(model_.createExpression(ls::O_Eq, sum, lb));
  } else {
    model_.addConstraint(model_.createExpression(ls::O_Geq, sum, lb));
    model_.addConstraint(model_.createExpression(ls::O_Leq, sum, ub));
  }
  return LinearConBuilder(*this, expr, sum);
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
  case expr::INT_DIV:
    return IntDiv(lhs, rhs);
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

LSProblemBuilder::ArgHandler LSProblemBuilder::BeginVarArg(
    expr::Kind kind, int num_args) {
  ls::LSOperator op = ls::O_Min;
  if (kind == expr::MAX)
    op = ls::O_Max;
  else if (kind != expr::MIN)
    Base::BeginVarArg(kind, num_args);
  return ArgHandler(model_.createExpression(op));
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

LSProblemBuilder::ArgHandler LSProblemBuilder::BeginIteratedLogical(
    expr::Kind kind, int num_args) {
  ls::LSOperator op = ls::O_Or;
  if (kind == expr::FORALL)
    op = ls::O_And;
  else if (kind != expr::EXISTS)
    Base::BeginIteratedLogical(kind, num_args);
  return ArgHandler(model_.createExpression(op));
}

ls::LSExpression LSProblemBuilder::EndAllDiff(AllDiffArgHandler handler) {
  std::vector<ls::LSExpression> &args = handler.args;
  ls::LSExpression alldiff = model_.createExpression(ls::O_And);
  for (std::size_t i = 0, n = args.size(); i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j)
      alldiff.addOperand(MakeBinary(ls::O_Neq, args[i], args[j]));
  }
  return alldiff;
}

LocalSolver::LocalSolver()
  : SolverImpl<LSProblemBuilder>("localsolver", 0, 20140710, MULTIPLE_OBJ),
    timelimit_(0) {
  std::string version = fmt::format("{}.{}",
      localsolver::LSVersion::getMajorVersionNumber(),
      localsolver::LSVersion::getMinorVersionNumber());
  set_long_name("localsolver " + version);
  set_version("LocalSolver " + version);

  set_option_header(
      "LocalSolver Options for AMPL\n"
      "----------------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to "
      "the AMPL option ``localsolver_options``. For example::\n"
      "\n"
      "  ampl: option localsolver_options 'version timelimit=30;\n");

  AddIntOption("timelimit",
      "Time limit in seconds (positive integer). Default = no limit.",
      &LocalSolver::GetTimeLimit, &LocalSolver::SetTimeLimit);
}

void LocalSolver::Solve(ProblemBuilder &builder, SolutionHandler &sh) {
  steady_clock::time_point time = steady_clock::now();

  ls::LocalSolver &solver = builder.solver();
  ls::LSModel model = solver.getModel();

  // LocalSolver requires at least one objective - create a dummy one.
  if (model.getNbObjectives() == 0)
    model.addObjective(model.createConstant(AsLSInt(0)), ls::OD_Minimize);
  model.close();

  // Set options. LS requires this to be done after the model is closed.
  ls::LSPhase phase = solver.createPhase();
  if (timelimit_ != 0)
    phase.setTimeLimit(timelimit_);
  interrupter()->SetHandler(StopSolver, &solver);

  // TODO: more options

  double setup_time = GetTimeAndReset(time);

  // Solve the problem.
  solver.solve();

  // Convert solution status.
  int solve_code = sol::UNKNOWN;
  ls::LSSolution sol = solver.getSolution();
  const char *status = "unknown";
  if (interrupter()->Stop()) {
    solve_code = sol::INTERRUPTED;
    status = "interrupted";
  } else {
    switch (sol.getStatus()) {
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
      solve_code = sol::UNSOLVED;
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

  std::vector<double> solution;
  if (solve_code < sol::INFEASIBLE) {
    int num_vars = builder.num_vars();
    const ls::LSExpression *vars = builder.vars();
    solution.resize(num_vars);
    for (int i = 0; i < num_vars; ++i)
      solution[i] = GetValue(vars[i]);
  }
  double solution_time = GetTimeAndReset(time);

  fmt::MemoryWriter w;
  w.write("{}: {}\n", long_name(), status);
  w.write("{}", solver.getStatistics().toString());
  double obj_val = std::numeric_limits<double>::quiet_NaN();
  if (builder.num_objs() != 0) {
    obj_val = GetValue(model.getObjective(0));
    w.write("objective {}", FormatObjValue(obj_val));
  }
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

SolverPtr CreateSolver(const char *) { return SolverPtr(new LocalSolver()); }
}  // namespace mp

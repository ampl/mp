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

#include "localsolver-direct/localsolver-direct.h"

#include <cmath>
#include "mp/clock.h"
#include "asl/problem.h"

namespace {
// Returns the value of an expression.
inline double GetValue(localsolver::LSExpression e) {
  return e.isDouble() ? e.getDoubleValue() : e.getValue();
}
}

namespace mp {

LSProblemBuilder::LSProblemBuilder(LocalSolver &solver)
  : model_(solver.model()), num_continuous_vars_(0) {
}

void LSProblemBuilder::BeginBuild(const NLHeader &header) {
  vars_.resize(header.num_vars);
  objs_.resize(header.num_objs);
  cons_.resize(header.num_algebraic_cons + header.num_logical_cons);
  num_continuous_vars_ = header.num_continuous_vars();
  for (int i = 0; i < num_continuous_vars_; ++i)
    vars_[i] = model_.createExpression(ls::O_Float);
  for (int i = num_continuous_vars_; i < header.num_vars; ++i)
    vars_[i] = model_.createExpression(ls::O_Int);
}

void LSProblemBuilder::EndBuild() {
  // Add objectives.
  for (std::size_t i = 0, n = objs_.size(); i < n; ++i) {
    const ObjInfo &obj = objs_[i];
    model_.addObjective(obj.expr, obj.direction);
  }
  // Add constraints.
  double inf = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0, n = cons_.size(); i < n; ++i) {
    const ConInfo &con = cons_[i];
    ls::LSExpression expr;
    if (con.lb <= -inf) {
      expr = model_.createExpression(ls::O_Leq, con.expr, con.ub);
    } else if (con.ub >= inf) {
      expr = model_.createExpression(ls::O_Geq, con.expr, con.lb);
    } else if (con.lb == con.ub) {
      expr = model_.createExpression(ls::O_Eq, con.expr, con.lb);
    } else {
      expr = model_.createExpression(ls::O_Geq, con.expr, con.lb);
      expr = model_.createExpression(ls::O_Leq, con.expr, con.ub);
    }
    model_.addConstraint(expr);
  }
  model_.close();
}

ls::LSExpression LSProblemBuilder::MakeUnary(
    expr::Kind kind, ls::LSExpression arg) {
  switch (kind) {
    case expr::FLOOR: return model_.createExpression(ls::O_Floor, arg);
    case expr::CEIL:  return model_.createExpression(ls::O_Ceil, arg);
    case expr::ABS:   return model_.createExpression(ls::O_Abs, arg);
    case expr::MINUS: return Negate(arg);
    case expr::TANH: {
      ls::LSExpression exp_x = model_.createExpression(ls::O_Exp, arg);
      ls::LSExpression exp_minus_x =
          model_.createExpression(ls::O_Exp, Negate(arg));
      return model_.createExpression(
            ls::O_Div, model_.createExpression(ls::O_Sub, exp_x, exp_minus_x),
            model_.createExpression(ls::O_Sum, exp_x, exp_minus_x));
    }
    case expr::TAN:   return model_.createExpression(ls::O_Tan, arg);
    case expr::SQRT:  return model_.createExpression(ls::O_Sqrt, arg);
    case expr::SINH:  break; // TODO
    case expr::SIN:   return model_.createExpression(ls::O_Sin, arg);
    case expr::LOG10: break; // TODO
    case expr::LOG:   return model_.createExpression(ls::O_Log, arg);
    case expr::EXP:   return model_.createExpression(ls::O_Exp, arg);
    case expr::COSH:  break; // TODO
    case expr::COS:   return model_.createExpression(ls::O_Cos, arg);
    case expr::ATANH:
    case expr::ASINH:
    case expr::ACOSH:
      break; // TODO
    case expr::POW2:
      return model_.createExpression(ls::O_Pow, arg, MakeInt(2));
    case expr::ATAN: case expr::ASIN: case expr::ACOS:
      // LocalSolver doesn't support atan, asin and acos.
      // Fall through.
    default:
      break;
  }
  return Base::MakeUnary(kind, arg);
}

ls::LSExpression LSProblemBuilder::MakeBinary(
    expr::Kind kind, ls::LSExpression lhs, ls::LSExpression rhs) {
  switch (kind) {
    case expr::ADD:
      return model_.createExpression(ls::O_Sum, lhs, rhs);
    case expr::SUB:
      return model_.createExpression(ls::O_Sub, lhs, rhs);
    case expr::MUL:
      return model_.createExpression(ls::O_Prod, lhs, rhs);
    case expr::DIV:
      return model_.createExpression(ls::O_Div, lhs, rhs);
    case expr::INT_DIV:
    case expr::MOD:
      break; // TODO
    case expr::POW:
    case expr::POW_CONST_BASE:
    case expr::POW_CONST_EXP:
      return model_.createExpression(ls::O_Pow, lhs, rhs);
    case expr::LESS:
    case expr::PRECISION:
    case expr::ROUND:
    case expr::TRUNC:
      break; // TODO
    case expr::ATAN2:
      // LocalSolver doesn't support atan2.
      // Fall through.
    default:
      break;
  }
  return Base::MakeBinary(kind, lhs, rhs);
}

/*
void NLToLocalSolverConverter::Convert(const Problem &p) {
  // Convert logical constraints.
/*  int num_logical_cons = p.num_logical_cons();
  for (int i = 0; i < num_logical_cons; ++i) {
    LogicalExpr e = p.logical_con_expr(i);
    AllDiffExpr alldiff = Cast<AllDiffExpr>(e);
    ICLSetter icl_setter(icl_, GetICL(p.num_cons() + i));
    if (!alldiff) {
      rel(problem_, Visit(e), icl_);
      continue;
    }
    int num_args = alldiff.num_args();
    IntVarArgs args(num_args);
    for (int i = 0; i < num_args; ++i) {
      NumericExpr arg(alldiff[i]);
      if (Variable var = ampl::Cast<Variable>(arg))
        args[i] = vars[var.index()];
      else
        args[i] = Gecode::expr(problem_, Visit(arg), icl_);
    }
    distinct(problem_, args, icl_);
  }*/
/*}

ls::LSExpression NLToLocalSolverConverter::VisitLog10(UnaryExpr e) {
  return model_.createExpression(
      ls::O_Div, ConvertUnary(ls::O_Log, e), std::log(10.0));
}*/

LocalSolver::LocalSolver()
  : SolverImpl<LSProblemBuilder>("localsolver", 0, 20140710), timelimit_(0) {
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

/*ls::LSExpression NLToLocalSolverConverter::VisitAllDiff(AllDiffExpr e) {
  ls::LSExpression result = model_.createExpression(ls::O_And);
  int num_args = e.num_args();
  std::vector<ls::LSExpression > args(num_args);
  for (int i = 0; i < num_args; ++i)
    args[i] = Visit(e[i]);
  for (int i = 0; i < num_args; ++i) {
    for (int j = i + 1; j < num_args; ++j)
      result.addOperand(model_.createExpression(ls::O_Neq, args[i], args[j]));
  }
  return result;
}*/

void LocalSolver::Solve(LSProblemBuilder &pb) {
  steady_clock::time_point time = steady_clock::now();

  pb.EndBuild();

  // Set options. LS requires this to be done after the model is closed.
  ls::LSPhase phase = solver_.createPhase();
  if (timelimit_ != 0)
    phase.setTimeLimit(timelimit_);

  double setup_time = GetTimeAndReset(time);

  // Solve the problem.
  solver_.solve();

  // Convert solution status.
  int solve_code = 0;
  ls::LSSolution sol = solver_.getSolution();
  const char *status = "unknown";
  switch (sol.getStatus()) {
  case ls::SS_Inconsistent:
    solve_code = INFEASIBLE;
    status = "infeasible problem";
    break;
  case ls::SS_Infeasible:
    // Solution is infeasible, but problem may be feasible.
    // This can only happen if stopped by a limit.
    solve_code = LIMIT;
    status = "infeasible solution";
    break;
  case ls::SS_Feasible:
    solve_code = SOLVED_MAYBE;
    status = "feasible solution";
    break;
  case ls::SS_Optimal:
    solve_code = SOLVED;
    status = "optimal solution";
    break;
  default:
    solve_code = FAILURE;
    status = "unknown solution status";
    break;
  }
  // TODO
  //p.set_solve_code(solve_code);

  int num_vars = pb.num_vars();
  ls::LSExpression const *vars = 0; //converter.vars();
  std::vector<double> solution(num_vars);
  int num_continuous_vars = pb.num_continuous_vars();
  for (int i = 0; i < num_continuous_vars; ++i)
    solution[i] = vars[i].getDoubleValue();
  for (int i = num_continuous_vars; i < num_vars; ++i)
    solution[i] = vars[i].getValue();
  double solution_time = GetTimeAndReset(time);

  fmt::Writer w;
  w.write("{}: {}\n", long_name(), status);
  w.write("{}", solver_.getStatistics().toString());
  double obj_val = std::numeric_limits<double>::quiet_NaN();
  if (pb.num_objs() != 0) {
    obj_val = GetValue(solver_.getModel().getObjective(0));
    w.write("objective {}", FormatObjValue(obj_val));
  }
  // TODO: write solution
  //HandleSolution(pb, w.c_str(),
  //    solution.empty() ? 0 : solution.data(), 0, obj_val);
  double output_time = GetTimeAndReset(time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n",
          setup_time, solution_time, output_time);
  }
}

// TODO
SolverPtr CreateSolver(const char *) { return SolverPtr(new LocalSolver()); }
}  // namespace mp

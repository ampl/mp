/*
 AMPL solver interface to Sulum.

 Copyright (C) 2013 AMPL Optimization Inc

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

#include "solvers/sulum/sulum.h"
#include "solvers/util/clock.h"
#include <sulumcpp.h>

namespace {
inline void Check(SlmReturn ret) {
  if (ret != SlmRetOk)
    throw Slm::SlmException(ret);
}
}

namespace ampl {

std::string SulumSolver::GetOptionHeader() {
  // TODO
  return
      "Sulum Directives for AMPL\n"
      "--------------------------\n"
      "\n"
      "To set these directives, assign a string specifying their values to "
      "the AMPL option sulum_options.  For example:\n"
      "\n"
      "  ampl: option sulum_options 'version nodelimit=30000 "
      "val_branching=min';\n";
}

SulumSolver::SulumSolver() : Solver("sulum", "", 20130820) {
  int major = 0, minor = 0, interim = 0;
  SlmGetSulumVersion(&major, &minor, &interim);
  std::string version = str(
      fmt::Format("sulum {}.{}.{}") << major << minor << interim);
  set_long_name(version);
  version[0] = 'S';
  set_version(version);
  
  // TODO: register options and suffixes
}

void SulumSolver::DoSolve(Problem &p) {
  steady_clock::time_point time = steady_clock::now();

  if (p.num_nonlinear_objs() != 0 || p.num_nonlinear_cons() != 0)
    throw Error("Sulum doesn't support nonlinear problems");

  class SulumModel : Noncopyable {
   private:
    SlmEnv_t env_;
    SlmModel_t model_;

   public:
    SulumModel() : env_(), model_() {
      Check(SlmMakeEnv(&env_));
      SlmReturn ret = SlmMakeModel(env_, &model_);
      if (ret != SlmRetOk)
        SlmFreeEnv(&env_);
      Check(ret);
    }
    ~SulumModel() {
      SlmFreeModel(env_, &model_);
      SlmFreeEnv(&env_);
    }

    operator SlmModel_t() const { return model_; }
  };
  SulumModel model;

  // Convert variables.
  int num_vars = p.num_vars();
  std::vector<SlmBoundKey> var_status(num_vars);
  std::vector<double> var_lb(num_vars);
  std::vector<double> var_ub(num_vars);
  for (int i = 0; i < num_vars; ++i) {
    double lb = p.var_lb(i), ub = p.var_ub(i);
    if (lb <= -Infinity)
      lb = -SlmInfinity;
    if (ub >= Infinity)
      ub = SlmInfinity;
    var_lb[i] = lb;
    var_ub[i] = ub;
  }
  // TODO: handle integer variables

  std::vector<double> obj_coefs(num_vars);
  if (p.num_objs() > 0) {
    // Convert objective.
    Check(SlmSetIntParam(model, SlmPrmIntObjSense,
        p.obj_type(0) == MIN ? SlmObjSenseMin : SlmObjSenseMax));
    LinearObjExpr expr = p.linear_obj_expr(0);
    for (LinearObjExpr::iterator
        i = expr.begin(), end = expr.end(); i != end; ++i) {
      obj_coefs[i->var_index()] = i->coef();
    }
  }

  // Convert constraints.
  int num_cons = p.num_cons();
  std::vector<SlmBoundKey> con_status(num_cons);
  std::vector<double> con_lb(num_cons);
  std::vector<double> con_ub(num_cons);
  for (int i = 0; i < num_cons; ++i) {
    double lb = p.con_lb(i), ub = p.con_ub(i);
    if (lb <= -Infinity)
      lb = -SlmInfinity;
    if (ub >= Infinity)
      ub = SlmInfinity;
    con_lb[i] = lb;
    con_ub[i] = ub;
  }

  Problem::ColMatrix matrix = p.col_matrix();
  Check(SlmSetAllData(model, num_cons, num_vars, 0,
      ptr(con_status), ptr(con_lb), ptr(con_ub),
      ptr(var_status), ptr(obj_coefs), ptr(var_lb), ptr(var_ub),
      1, matrix.col_starts(), matrix.row_indices(), matrix.values()));

  double setup_time = GetTimeAndReset(time);

  // Solve the problem.
  Check(SlmOptimize(model));

  // TODO: convert the solution status
  int solve_code = 0;
  const char *status = "";
  p.set_solve_code(solve_code);

  double solution_time = GetTimeAndReset(time);

  fmt::Writer w;
  w.Format("{}: {}\n") << long_name() << status;
  double obj_val = 0; // TODO
  std::vector<double> solution(num_vars);
  Check(SlmGetSolPrimVars(model, &solution[0]));
  // TODO
  //w.Format("{} nodes, {} fails") << stats.node << stats.fail;
  //if (has_obj && solution.get())
  //  w.Format(", objective {}") << ObjPrec(obj_val);
  HandleSolution(p, w.c_str(), ptr(solution), 0, obj_val);

  double output_time = GetTimeAndReset(time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n")
              << setup_time << solution_time << output_time;
  }
}
}

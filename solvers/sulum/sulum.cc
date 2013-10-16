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
template <typename Param>
struct OptionInfo {
  Param param;
  const char *name;
  const char *description;
};

const OptionInfo<SlmParamInt> INT_OPTION_INFO[] = {
#define INT_OPTION(param, name, description) \
  {param, name, description},
#define DBL_OPTION(param, name, description)
#include "optsulum.ampl"
#undef INT_OPTION
#undef DBL_OPTION
};

const OptionInfo<SlmParamDb> DBL_OPTION_INFO[] = {
#define INT_OPTION(param, name, description)
#define DBL_OPTION(param, name, description) \
    {param, name, description},
#include "optsulum.ampl"
};

inline void Check(SlmReturn ret) {
  if (ret != SlmRetOk)
    throw Slm::SlmException(ret);
}

inline SlmBoundKey GetBoundKey(double lb, double ub) {
  if (lb <= -Infinity)
    return ub >= Infinity ? SlmBndFr : SlmBndUp;
  if (ub >= Infinity)
    return SlmBndLo;
  return lb == ub ? SlmBndFx : SlmBndRa;
}
}

namespace ampl {

class SulumSolver::IntSulumOption : public TypedSolverOption<int> {
 private:
  SulumSolver *solver_;
  SlmParamInt param_;

 public:
  IntSulumOption(const OptionInfo<SlmParamInt> &info, SulumSolver *s)
  : TypedSolverOption<int>(info.name, info.description),
    solver_(s), param_(info.param) {}

  int GetValue() const {
    int value = 0;
    Check(SlmGetIntParam(solver_->model_, param_, &value));
    return value;
  }
  void SetValue(int value) {
    Check(SlmSetIntParam(solver_->model_, param_, value));
  }
};

class SulumSolver::DblSulumOption : public TypedSolverOption<double> {
 private:
  SulumSolver *solver_;
  SlmParamDb param_;

 public:
  DblSulumOption(const OptionInfo<SlmParamDb> &info, SulumSolver *s)
  : TypedSolverOption<double>(info.name, info.description),
    solver_(s), param_(info.param) {}

  double GetValue() const {
    double value = 0;
    Check(SlmGetDbParam(solver_->model_, param_, &value));
    return value;
  }
  void SetValue(double value) {
    Check(SlmSetDbParam(solver_->model_, param_, value));
  }
};

std::string SulumSolver::GetOptionHeader() {
  return
      "Sulum Directives for AMPL\n"
      "--------------------------\n"
      "\n"
      "To set these directives, assign a string specifying their values to "
      "the AMPL option sulum_options.  For example:\n"
      "\n"
      "  ampl: option sulum_options 'version loglevel=10 "
      "simmaxiter=100';\n";
}

SulumSolver::SulumSolver() : Solver("sulum", "", 20130908), env_(), model_() {
  int major = 0, minor = 0, interim = 0;
  SlmGetSulumVersion(&major, &minor, &interim);
  std::string version = str(
      fmt::Format("sulum {}.{}.{}") << major << minor << interim);
  set_long_name(version);
  version[0] = 'S';
  set_version(version);
  set_read_flags(Problem::READ_COLUMNWISE);

  Check(SlmMakeEnv(&env_));
  SlmReturn ret = SlmMakeModel(env_, &model_);
  if (ret != SlmRetOk)
    SlmFreeEnv(&env_);
  Check(ret);

  Check(SlmSetIntParam(model_, SlmPrmIntLogPrefix, SlmOff));

  size_t num_int_options = sizeof(INT_OPTION_INFO) / sizeof(*INT_OPTION_INFO);
  for (size_t i = 0; i < num_int_options; ++i)
    AddOption(OptionPtr(new IntSulumOption(INT_OPTION_INFO[i], this)));

  size_t num_dbl_options = sizeof(DBL_OPTION_INFO) / sizeof(*DBL_OPTION_INFO);
  for (size_t i = 0; i < num_dbl_options; ++i)
    AddOption(OptionPtr(new DblSulumOption(DBL_OPTION_INFO[i], this)));

  // TODO: register suffixes
}

SulumSolver::~SulumSolver() {
  SlmFreeModel(env_, &model_);
  SlmFreeEnv(&env_);
}

void SulumSolver::DoSolve(Problem &p) {
  steady_clock::time_point time = steady_clock::now();

  if (p.num_nonlinear_objs() != 0 || p.num_nonlinear_cons() != 0)
    throw Error("Sulum doesn't support nonlinear problems");

  // Convert variables.
  int num_vars = p.num_vars();
  std::vector<SlmBoundKey> var_bound_keys(num_vars);
  for (int i = 0; i < num_vars; ++i)
    var_bound_keys[i] = GetBoundKey(p.var_lb(i), p.var_ub(i));

  std::vector<double> obj_coefs(num_vars);
  if (p.num_objs() > 0) {
    // Convert objective.
    Check(SlmSetIntParam(model_, SlmPrmIntObjSense,
        p.obj_type(0) == MIN ? SlmObjSenseMin : SlmObjSenseMax));
    LinearObjExpr expr = p.linear_obj_expr(0);
    for (LinearObjExpr::iterator
        i = expr.begin(), end = expr.end(); i != end; ++i) {
      obj_coefs[i->var_index()] = i->coef();
    }
  }

  // Convert constraints.
  int num_cons = p.num_cons();
  std::vector<SlmBoundKey> con_bound_keys(num_cons);
  for (int i = 0; i < num_cons; ++i)
    con_bound_keys[i] = GetBoundKey(p.con_lb(i), p.con_ub(i));

  // Convert constraint matrix.
  Problem::ColMatrix matrix = p.col_matrix();
  std::vector<int> col_starts;
  col_starts.reserve(num_vars + 1);
  col_starts.assign(matrix.col_starts(), matrix.col_starts() + num_vars);
  col_starts.push_back(p.num_con_nonzeros());
  Check(SlmSetAllData(model_, num_cons, num_vars, 0,
      con_bound_keys.data(), p.con_lb(), p.con_ub(),
      var_bound_keys.data(), obj_coefs.data(), p.var_lb(), p.var_ub(),
      1, &col_starts[0], matrix.row_indices(), matrix.values()));

  // Set up integer variables.
  if (int num_int_vars = p.num_integer_vars()) {
    std::vector<SlmVarType> var_types(num_int_vars, SlmVarTypeInt);
    for (int i = p.num_continuous_vars(), j = 0; i < num_vars; ++i, ++j) {
      if (p.var_lb(i) == 0 && p.var_ub(i) == 1)
        var_types[j] = SlmVarTypeBin;
    }
    SlmSetTypeVarsFromTo(model_,
        p.num_continuous_vars(), num_vars, &var_types[0]);
  }

  double setup_time = GetTimeAndReset(time);

  // Make Sulum update info items.
  Check(SlmSetIntParam(model_, SlmPrmIntUpdateSolQuality, SlmOn));

  // Solve the problem.
  // No need to handle SIGINT because Sulum does it.
  SlmReturn ret = SlmOptimize(model_);
  int solve_code = NOT_SOLVED;
  const char *status = "";
  if (ret == SlmRetUserTerm) {
    solve_code = 600;
    status = "interrupted";
    // TODO: make sure solve code is not overwritten later
  } else {
    Check(ret);
  }

  SlmSolStatus sulum_status = SlmSolStatUnk;
  Check(SlmGetSolStatus(model_, &sulum_status));
  switch (sulum_status) {
  default:
    assert(0 && "unknown solution status");
    // Fall through.
  case SlmSolStatUnk:
    solve_code = 500;
    status = "unknown";
    break;
  case SlmSolStatOpt:
    solve_code = SOLVED;
    status = "optimal solution";
    break;
  case SlmSolStatPrimFeas:
    solve_code = SOLVED_MAYBE;
    status = "feasible solution";
    break;
  case SlmSolStatDualFeas:
    solve_code = SOLVED_MAYBE + 1;
    status = "dual feasible solution";
    break;
  case SlmSolStatPrimInf:
    solve_code = INFEASIBLE;
    status = "infeasible problem";
    break;
  case SlmSolStatDualInf:
    solve_code = INFEASIBLE + 1;
    status = "infeasible or unbounded";
    break;
  case SlmSolStatIntFeas:
    solve_code = SOLVED_MAYBE;
    status = "integer feasible solution";
    break;
  case SlmSolStatIntInf:
    solve_code = INFEASIBLE + 1;
    status = "integer infeasible";
    break;
  case SlmSolStatIntUnBndInf:
    solve_code = INFEASIBLE + 2;
    status = "integer infeasible or unbounded";
    break;
  }
  p.set_solve_code(solve_code);

  double solution_time = GetTimeAndReset(time);

  fmt::Writer w;
  w.Format("{}: {}\n") << long_name() << status;
  double obj_val = 0;
  Check(SlmGetDbInfo(model_, SlmInfoDbPrimObj, &obj_val));
  std::vector<double> solution(num_vars);
  Check(SlmGetSolPrimVars(model_, solution.data()));
  std::vector<double> dual_solution(num_cons);
  Check(SlmGetSolDualCons(model_, dual_solution.data()));
  w << status;
  if (p.num_objs() > 0)
    w.Format("; objective {}") << ObjPrec(obj_val);
  HandleSolution(p, w.c_str(), solution.data(), dual_solution.data(), obj_val);

  double output_time = GetTimeAndReset(time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n")
              << setup_time << solution_time << output_time;
  }
}
}

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

#include "sulum/sulum.h"
#include "mp/utils-clock.h"

#include <limits>

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
  if (ret != SlmRetOk) {
    // TODO: get error message
    throw mp::Error("Sulum error {}", ret);
  }
}

inline SlmBoundKey GetBoundKey(double lb, double ub) {
  double inf = std::numeric_limits<double>::infinity();
  if (lb <= -inf)
    return ub >= inf ? SlmBndFr : SlmBndUp;
  if (ub >= inf)
    return SlmBndLo;
  return lb == ub ? SlmBndFx : SlmBndRa;
}
}

namespace mp {

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

SulumSolver::SulumSolver()
  : SolverImpl<ColProblem>("sulum", "", 20160218), env_(), model_() {
  int major = 0, minor = 0, interim = 0;
  SlmGetSulumVersion(&major, &minor, &interim);
  std::string version = fmt::format("sulum {}.{}.{}", major, minor, interim);
  set_long_name(version);
  version[0] = 'S';
  set_version(version);

  Check(SlmMakeEnv(&env_));
  SlmReturn ret = SlmMakeModel(env_, &model_);
  if (ret != SlmRetOk)
    SlmFreeEnv(&env_);
  Check(ret);

  Check(SlmSetIntParam(model_, SlmPrmIntLogPrefix, SlmOff));

  set_option_header(
      "Sulum Options for AMPL\n"
      "----------------------\n"
      "\n"
      "To set these options, assign a string specifying their values to "
      "the AMPL option ``sulum_options``. For example::\n"
      "\n"
      "  ampl: option sulum_options 'version loglevel=10 simmaxiter=100';\n");

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

void SulumSolver::Solve(ColProblem &p, SolutionHandler &sh) {
  std::chrono::steady_clock::time_point time = std::chrono::steady_clock::now();

  if (p.has_nonlinear_cons() || p.num_logical_cons() != 0)
    throw Error("Sulum doesn't support nonlinear problems");

  // Convert variables.
  int num_vars = p.num_vars();
  std::vector<SlmVarType> var_types;
  std::vector<SlmBoundKey> var_bound_keys(num_vars);
  std::vector<double> var_lb(num_vars), var_ub(num_vars);
  for (int i = 0; i < num_vars; ++i) {
    Problem::Variable var = p.var(i);
    if (var.type() == mp::var::INTEGER) {
      if (var_types.empty())
        var_types.resize(num_vars);
      var_types[i] = var.lb() == 0 && var.ub() == 1 ?
            SlmVarTypeBin : SlmVarTypeInt;
    }
    var_bound_keys[i] = GetBoundKey(var.lb(), var.ub());
    var_lb[i] = var.lb();
    var_ub[i] = var.ub();
  }

  std::vector<double> obj_coefs(num_vars);
  if (p.num_objs() > 0) {
    // Convert objective.
    Problem::Objective obj = p.obj(0);
    if (obj.nonlinear_expr())
      throw Error("Sulum doesn't support nonlinear problems");
    Check(SlmSetIntParam(model_, SlmPrmIntObjSense,
        obj.type() == obj::MIN ? SlmObjSenseMin : SlmObjSenseMax));
    const LinearExpr &expr = obj.linear_expr();
    for (LinearExpr::iterator
        i = expr.begin(), end = expr.end(); i != end; ++i) {
      obj_coefs[i->var_index()] = i->coef();
    }
  }

  // Convert constraints.
  int num_cons = p.num_algebraic_cons();
  std::vector<SlmBoundKey> con_bound_keys(num_cons);
  std::vector<double> con_lb, con_ub;
  for (int i = 0; i < num_cons; ++i) {
    Problem::AlgebraicCon con = p.algebraic_con(i);
    con_bound_keys[i] = GetBoundKey(con.lb(), con.ub());
  }

  // Create Sulum problem.
  Check(SlmSetAllData(model_, num_cons, num_vars, 0,
      con_bound_keys.data(), con_lb.data(), con_ub.data(),
      var_bound_keys.data(), obj_coefs.data(), var_lb.data(), var_ub.data(),
      1, p.col_starts(), p.row_indices(), p.values()));

  // Set up integer variables.
  if (!var_types.empty())
    SlmSetTypeVars(model_, &var_types[0]);

  double setup_time = GetTimeAndReset(time);

  // Make Sulum update info items.
  Check(SlmSetIntParam(model_, SlmPrmIntUpdateSolQuality, SlmOn));

  // Solve the problem.
  // No need to handle SIGINT because Sulum does it.
  SlmReturn ret = SlmOptimize(model_);
  int solve_code = sol::UNKNOWN;
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
    solve_code = sol::SOLVED;
    status = "optimal solution";
    break;
  case SlmSolStatPrimFeas:
    solve_code = sol::UNCERTAIN;
    status = "feasible solution";
    break;
  case SlmSolStatDualFeas:
    solve_code = sol::UNCERTAIN + 1;
    status = "dual feasible solution";
    break;
  case SlmSolStatPrimInf:
    solve_code = sol::INFEASIBLE;
    status = "infeasible problem";
    break;
  case SlmSolStatDualInf:
    solve_code = sol::INFEASIBLE + 1;
    status = "infeasible or unbounded";
    break;
  case SlmSolStatIntFeas:
    solve_code = sol::UNCERTAIN + 2;
    status = "integer feasible solution";
    break;
  case SlmSolStatIntInf:
    solve_code = sol::INFEASIBLE + 2;
    status = "integer infeasible";
    break;
  case SlmSolStatIntUnBndInf:
    solve_code = sol::INFEASIBLE + 3;
    status = "integer infeasible or unbounded";
    break;
  }

  double solution_time = GetTimeAndReset(time);

  fmt::MemoryWriter w;
  w.write("{}: {}\n", long_name(), status);
  double obj_val = 0;
  Check(SlmGetDbInfo(model_, SlmInfoDbPrimObj, &obj_val));
  std::vector<double> solution(num_vars);
  Check(SlmGetSolPrimVars(model_, solution.data()));
  std::vector<double> dual_solution(num_cons);
  Check(SlmGetSolDualCons(model_, dual_solution.data()));
  w << status;
  if (p.num_objs() > 0)
    w.write("; objective {}", FormatObjValue(obj_val));
  sh.HandleSolution(solve_code, w.c_str(), solution.data(),
                    dual_solution.data(), obj_val);

  double output_time = GetTimeAndReset(time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n",
          setup_time, solution_time, output_time);
  }
}

SolverPtr create_sulum(const char *) { return SolverPtr(new SulumSolver()); }
}

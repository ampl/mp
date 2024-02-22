/**
 * C "easy API" example using NLW2_NLModel_C
 * to write NL file, execute solver, and read SOL file
 * for a MILP model.
 *
 * For full API, see documentation of NLW2
 * and related tests/examples.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "api/c/nl-solver-c.h"

/// Solution checking data
typedef struct NLW2_SolCheckData_C {
  double *x_ref_;
  double obj_val_ref_;
} NLW2_SolCheckData_C;

/// Build model
NLW2_NLModel_C BuildModel(NLW2_SolCheckData_C* p_solchkdata) {
  /*
var x1 >=0, <=14, integer;
var x2 >=-3, <=15;
var x3 >=-1, <=16;

minimize O: x1+x2+x3;

C1: 15*x1 + 17*x3 >= 19;
C2: 33*x1 + 28*x2 >= 37;
   */
  static const char* prob_name_  = "NLSOL_C_EasyAPI_MILP";
  static double var_lb_[] = {0, -3, -1};
  static double var_ub_[] = {14, 15, 16};
  static int var_type_[]  = {1,  0,  0};
  static const char* var_names_[] = {"x1_3", "x2_1", "x3_2"};
  static int A_format_ = NLW2_MatrixFormatRowwise;
  static size_t A_start_[] = {0, 2};
  static int A_index_[] = {0,2,0,1};
  static double A_value_[] = {15,17,33,28};
  static double row_lb_[] = {19, 37};
  static double row_ub_[] = {INFINITY, INFINITY};
  static const char* row_names_[] = {"C1", "C2"};
  static NLW2_ObjSense obj_sense_ = NLW2_ObjSenseMinimize;
  static double obj_c0_ = 0.0;
  static double obj_c_[] = {1,1,1};
  static const char* obj_name_ = NULL;

  int num_col = sizeof(var_lb_) / sizeof(var_lb_[0]);
  int num_row = sizeof(row_lb_) / sizeof(row_lb_[0]);

  /// Solution
  static double x_ref[] = {3, -2.21429, -1};
  static double obj_val_ref  = -2.142857142857e-01;

  NLW2_NLModel_C nlme
      = NLW2_MakeNLModel_C(prob_name_);

  NLW2_SetCols_C(&nlme, num_col,
                var_lb_, var_ub_, var_type_);
  NLW2_SetColNames_C(&nlme, var_names_);

  NLW2_SetRows_C(&nlme, num_row, row_lb_, row_ub_,
                 A_format_, sizeof(A_index_) / sizeof(A_index_[0]),
                 A_start_, A_index_, A_value_);
  NLW2_SetRowNames_C(&nlme, row_names_);

  NLW2_SetLinearObjective_C(&nlme, obj_sense_, obj_c0_, obj_c_);
  NLW2_SetObjName_C(&nlme, obj_name_);

  p_solchkdata->obj_val_ref_ = obj_val_ref;
  p_solchkdata->x_ref_       = x_ref;

  return nlme;
}

/// Approx equal?
static int ApproxEqual(double n, double m) {
  return fabs(n-m)
      <= 1e-5 * fmin(1.0, fabs(n)+fabs(m));
}

/// Check solution
int NLW2_CheckSolution_C(
    NLW2_NLSolution_C* p_sol, NLW2_SolCheckData_C* p_sol_ref) {
  if (!ApproxEqual(p_sol->obj_val_, p_sol_ref->obj_val_ref_)) {
    printf("MILP 1: wrong obj val (%.17g !~ %.17g)\n",
           p_sol->obj_val_, p_sol_ref->obj_val_ref_);
    return 0;
  }
  for (int i=p_sol->n_primal_values_; i--; )
    if (!ApproxEqual(p_sol_ref->x_ref_[i], p_sol->x_[i])) {
      printf("MILP 1: wrong x[%d] (%.17g !~ %.17g)\n",
             i+1, p_sol_ref->x_ref_[i], p_sol->x_[i]);
      return 0;
    }
  printf("MILP 1: solution check ok.\n");
  return 1;
}

/// Solver with given parameters
int SolveAndCheck(const char* solver, const char* sopts,
                  int binary, const char* stub) {
  int result = 1;
  NLW2_SolCheckData_C solchkdata;
  NLW2_NLModel_C nlme = BuildModel(&solchkdata);
  NLW2_NLOptionsBasic_C nlopts = NLW2_MakeNLOptionsBasic_C_Default();
  nlopts.n_text_mode_ = !binary;
  nlopts.want_nl_comments_ = 1;
  NLW2_NLSolver_C nlse = NLW2_MakeNLSolver_C(NULL);
  NLW2_SetNLOptions_C(&nlse, nlopts);
  NLW2_SetFileStub_C(&nlse, stub);
  NLW2_NLSolution_C sol = NLW2_SolveNLModel_C(&nlse, &nlme, solver, sopts);
  if (sol.solve_result_ >= -1) {
    if (!NLW2_CheckSolution_C(&sol, &solchkdata)) {
      printf("Solution check failed.\n");
      result = 0;
    }
  } else {
    printf("%s\n", NLW2_GetErrorMessage_C(&nlse));
    result = 0;
  }

  // Destroy NLME... etc
  NLW2_DestroyNLSolver_C(&nlse);
  NLW2_DestroyNLModel_C(&nlme);

  return result;
}


/// main()
int main(int argc, const char* const* argv) {
  if (argc<2) {
    printf("%s\n",
           "AMPL NL writer C API example.\n"
           "Usage:\n"
           "  <this_exe> <solver> [\"<solver_options>\" [binary/text [<stub>]]],\n\n"
           "where <solver> is ipopt, gurobi, minos, ...;\n"
           "binary/text is the NL format (default: binary.)\n"
           "Examples:\n"
           "  <this_exe> highs \"writeprob=/tmp/stub.lp\" text /tmp/stub\n"
           "  <this_exe> gurobi \"mip:return_gap=1\"");
    exit(0);
  }

  const char* solver = (argc>1) ? argv[1] : "highs";
  const char* sopts = (argc>2) ? argv[2] : "";
  int binary = (argc<=3) || strcmp("text", argv[3]);
  const char* stub = (argc>4) ? argv[4] : "";

  if (!SolveAndCheck(solver, sopts, binary, stub))
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

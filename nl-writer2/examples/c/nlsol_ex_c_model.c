#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "nlsol_ex_c_model.h"

CAPIExample MakeCAPIExample_Linear_01(void) {
  /// Variables.
  /// Put y first because continuous variables
  /// come before integer ones.
  static const double
      var_lb[] = {-17, 0};
  static const double
      var_ub[] = {504, INFINITY};
  static const char* const
      var_name[] = {"y", "x"};

  /// Constraints
  static const double
      con_lb[] =   {-INFINITY, -INFINITY};
  static const double
      con_ub[] =   {      3e4,       3e5};
  static const char* const
      con_name[] = {     "C2",      "C3"};

  /// Constraint matrix
  static const SparseEntry row_01[]
      = { {0, 0.6}, {1, 3700} };
  static const SparseEntry row_02[]
      = { {0, 14536}, {1, 22} };
  static const SparseEntry* const rows[] = {row_01, row_02};

  static int row_nnz[] = {2, 2};
  static int col_sizes[] = {2};

  /// Objective
  static const SparseEntry obj_linpart[]
      = { {0, 13}, {1, 1.0} };

  static const double ini_x[] = {30.15, 15.11};
  static const double ini_y[] = {-10, -120};

  static const SparseEntry suf_val[] = {{0, 5.3}};

  static int n_suf = 1;
  static const Suffix suf[] = {
    {"zork", 1+4,   // constraints, float-valued
     1, suf_val}
  };

  /// Solution
  static double sol_dual[] = {NAN, NAN};
  static double sol_primal[] = {NAN, NAN};

  CAPIExample result = {
    .n_var = 2,
    .n_var_int = 1,
    .n_con = 2,
    .n_obj = 1,

    .var_lb = var_lb,
    .var_ub = var_ub,
    .var_name = var_name,

    .con_lb = con_lb,
    .con_ub = con_ub,
    .con_name = con_name,

    .con_linpart = rows,
    .row_nnz = row_nnz,
    .col_sizes = col_sizes,
    .n_con_nz = 4,

    .obj_sense = 1,
    .obj_linpart = obj_linpart,
    .n_obj_nz = 2,
    .obj_name = "TotalSum",

    .ini_x = ini_x,
    .ini_y = ini_y,

    .n_suf = n_suf,
    .suf = suf,

    .binary_nl = 1,

    .sol_dual_ = sol_dual,
    .sol_primal_ = sol_primal,
    .objno_ = -2,
    .solve_code_ = -100
  };
  return result;
}

void DestroyCAPIExample_Linear_01(CAPIExample* pEx) {
  pEx->var_lb = NULL;
  pEx->var_ub = NULL;
  pEx->var_name = NULL;

  // ...
}

void PrintSolution_C(CAPIExample* pex, const char* stub) {
  printf(
        "\n     ********** SOLUTION (%s.sol) ***********\n",
        stub);
  printf("%s\n", "DUALS.");
  for (int i=0; i<pex->n_con; ++i)
    printf("   %10s = %.17g\n",
           pex->con_name[i], pex->sol_dual_[i]);
  printf("%s\n", "PRIMALS.");
  for (int i=0; i<pex->n_con; ++i)
    printf("      %7s = %.17g\n",
           pex->var_name[i], pex->sol_primal_[i]);

  printf("\nObjno used: %d, solve_result_num: %d\n",
         pex->objno_+1, pex->solve_code_);
}

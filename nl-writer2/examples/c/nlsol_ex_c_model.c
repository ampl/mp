#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "nlsol_ex_c_model.h"

CAPIExample MakeCAPIExample_Linear_01() {
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
    .obj_name = "TotalSum"
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
  assert(0);
}

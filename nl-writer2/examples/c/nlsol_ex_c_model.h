#ifndef NLSOL_EX_C_MODEL_H
#define NLSOL_EX_C_MODEL_H

// #include "api/c/nl-header.h"
//#include "mp/nl-opcodes.h"        // When non-linearities

/**
 * A linear model for NL Writer C API example.
 * Illustrates NL variable order (continuous -> integer),
 * linear constraints, var/con/obj names.
 *

## To write NL file and name files in AMPL, use commands
##    ampl: option presolve 0;      # otherwise it's simplified
##    ampl: option nl_comments 1;
##    ampl: option auxfiles rc;
##    ampl: write gmodel;

var x >=0, integer;
var y >=-17, <=504;
maximize TotalSum:
    x + 13*y;
subj to C2:
       3700*x + 0.6*y <= 3e4;
subj to C3:
       22*x + 14536*y <= 3e5;

## Initial guess
let x := 1.5;
let y := 0.11;

## A suffix
suffix zork;
let C2.zork := 5;

 *
 */


/// Typedef SparseEntry
typedef struct SparseEntry {
  int index_;
  double value_;
} SparseEntry;


/// C API example data
typedef struct CAPIExample {
  const int n_var;
  const int n_var_int;
  const int n_con;
  const int n_obj;

  /////////////////////////////////////////////////
  //////////// Bounds and linear parts ////////////
  /////////////////////////////////////////////////

  /// Variables.
  /// Put y first because continuous variables
  /// come before integer ones.
  const double* var_lb;
  const double* var_ub;
  const char* const* var_name;

  /// Algebraic constraint bounds.
  /// If we had nonlinear constraints,
  /// constraints would be reordered: nonlinear first.
  /// If we had logical constraints, they'd go last
  /// (represented by expressions only.)
  const double* con_lb;
  const double* con_ub;
  const char* const* con_name;

  /// Linear constraints.
  ///
  /// If we had non-linear constraints, here would be
  /// their linear parts.
  /// They are presented sparse, but should include entries
  /// which can become nonzero in the Jacobian.


  /// Sparse vector per constraint
  const int* row_sizes;
  const SparseEntry* const* con_linpart;

  /// Sizes of all Jacobian columns except the last
  const int* col_sizes;
  const int n_con_nz;

  /// Obj sense: min/max
  const int obj_sense;
  /// Objective: linear part.
  ///
  /// Similar to the Jacobian,
  /// need all elements which can become nonzero.
  const SparseEntry* obj_linpart;
  int n_obj_nz;
  const char* obj_name;

  /// Some technical stuff
  int binary_nl;

} CAPIExample;

/// Create linear example data
CAPIExample MakeCAPIExample_Linear_01();

/// Destroy linear example data
void DestroyCAPIExample_Linear_01(CAPIExample* );

/// Print solution
void PrintSolution_C(CAPIExample* pex, const char* stub);

#endif // NLSOL_EX_C_MODEL_H

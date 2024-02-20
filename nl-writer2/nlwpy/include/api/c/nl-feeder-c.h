/*
 C API: extern "C" wrappers for the mp::NLFeeder interface,
 as well as for NLWriter2 calls.

 Copyright (C) 2024 AMPL Optimization Inc.

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

 Author: Gleb Belov */

#ifndef NLFeeder_C_H
#define NLFeeder_C_H

#include "mp/nl-header-c.h"


#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif

/// Declare callbacks

/// Write sparse vector(int) entry
void NLW2_WriteSparseIntEntry(
    void* p_api_data_, int index, int value);
/// Write sparse vector(double) entry
void NLW2_WriteSparseDblEntry(
    void* p_api_data_, int index, double value);
/// Write next variable's Lb, Ub
void NLW2_WriteVarLbUb(void* p_api_data, double lb, double ub);
/// Write next Jacobian column size
void NLW2_WriteColSize(void* p_api_data, int sz);

/// \rst
/// Algebraic constraint bounds (for a single constraint):
/// either range (lb, ub),
/// or complementarity info (k, cvar), when k>0.
///
/// For a complementarity constraint to hold, if cvar is at
///	its lower bound, then body >= 0; if cvar is at its upper
/// bound, then body <= 0;
///	and if cvar is strictly between its bounds, then body = 0.
/// The integer k in a complementarity constraint line indicates
/// which bounds on cvar are finite: 1 and 3 imply a finite
/// lower bound; 2 and 3 imply a finite upper bound; 0 (which
///	should not occur) would imply no finite bounds, i.e.,
/// body = 0 must always hold.
///
/// Example:
///
/// .. code-block:: ampl
///
///    ampl: var x; var y; var z;
///	   ampl: s.t. Compl1: x+y >= 3 complements x-z <= 15;
///	   ampl: s.t. Compl2: -2 <= 2*y+3*z <= 13 complements 6*z-2*x;
///	   ampl: expand;
///	   subject to Compl1:
///					3 <= x + y
///			 complements
///					x - z <= 15;
///
///	   subject to Compl2:
///					-2 <= 2*y + 3*z <= 13
///			 complements
///					-2*x + 6*z;
///
///	   ampl: solexpand;
///	   Nonsquare complementarity system:
///					4 complementarities including 2 equations
///					5 variables
///	   subject to Compl1.L:
///					x + y + Compl1$cvar = 0;
///
///	   subject to Compl1.R:
///					-15 + x - z <= 0
///			 complements
///					Compl1$cvar <= -3;
///
///	   subject to Compl2.L:
///					2*y + 3*z - Compl2$cvar = 0;
///
///	   subject to Compl2.R:
///					-2*x + 6*z
///			 complements
///					-2 <= Compl2$cvar <= 13;
///
/// \endrst
typedef struct NLW2_AlgConRange_C {
  double L, U;
  int k, cvar;    // k>0 means complementarity to cvar
} NLW2_AlgConRange_C;

/// Callback: write next constraint's range
void NLW2_WriteAlgConRange(void* , NLW2_AlgConRange_C*);

/// Callback: start int suffix
void* NLW2_StartIntSuffix(void* p_api_1,
    const char* suf_name, int kind, int nnz);
/// Callback: start dbl suffix
void* NLW2_StartDblSuffix(void* p_api_1,
    const char* suf_name, int kind, int nnz);

/// Callback: write model item name
void NLW2_WriteName(void* p_api_data, const char* name);
/// Callback: write fixed var name and comment
void NLW2_WriteNameAndComment(
    void* p_api_data, const char* name, const char* comment);
/// Callback: write obj name and offset
void NLW2_WriteNameAndNumber(
    void* p_api_data, const char* name, double val);


/** Wrap mp::NLFeeder for C API.

  NLW2_NLFeeder_C: writes model details on request
  via provided callback objects.
  See the examples folder.

  To fill some **default values and methods**,
  e.g., options and some methods like name feeders,
  call NLW2_MakeNLFeeder_C_Default() / NLW2_Destroy...().

  2023-11: CURRENT IMPLEMENTATION SUPPORTS LINEAR MODELS.

  For the NL format, variables and constraints must have certain order.

  **Variable ordering:**
    first continuous, then integer.
  Some solvers might require more elaborate ordering, see NLHeader_C.

  **Constraint ordering:**
    first algebraic (including complementarity), then logical.
  Some solvers might require nonlinear constraints first.
 */
typedef struct NLW2_NLFeeder_C {
  /// User data, provided as the 1st argument to the methods
  void* p_user_data_;

  ///////////////////////////////////////////////////////////////
  /// Set the below function pointers
  ///////////////////////////////////////////////////////////////

  ///////////////////// 1. NL HEADER AND OPTIONS /////////////////
  /** Provide NLHeader.
   *
   *	This method is called first.
   *
   *  NLHeader summarizes the model and provides some
   *  technical parameters,
   *  such as text/binary NL format. */
  NLHeader_C (*Header)(void* p_user_data);

  /// Options. Safe to leave as by the default generator.

  /// NL comments?
  int want_nl_comments_;

  /// The maximum number of significant digits written.
  /// The default value 0 requests full precision, which
  /// might be the shortest representation that, when
  /// converted to binary and properly rounded, will
  /// give exactly the binary value stored in the computer.
  int output_precision_;

  /// Write bounds first?
  /// The default is 1 (yes) in AMPL, controlled by
  /// (the value of option nl_permute) & 32
  /// (the bit is 0 for yes).
  /// Changing this option is deprecated, see
  /// https://netlib.org/ampl/changes.
  int want_bounds_first_;

  /// Want Jacobian column sizes?
  /// Required by some nonlinear solvers.
  /// Options: 0 - none, 1 - cumulative (default),
  /// 2 - non-cumulative.
  /// This option controls how ColSizeWriter
  /// writes the provided sizes (which should be
  /// non-cumulative).
  int want_column_sizes_;


  ///////////////////// 2. OBJECTIVES /////////////////////
  /** Description for objective function \a i
   *    (\a i in 0..num_objs-1).
   *  With WantNLComments()==true, this is
   *  written to text-format NL as a comment. */
  const char* (*ObjDescription)(void* p_user_data, int i);

  /** Provide type of objective \a i.
   *  0 - minimization;
   *  1 - maximization. */
  int (*ObjType)(void* p_user_data, int i);

  /** Number of nonzeros in the gradient for objective \a i.
   *  Should include entries for all potentially
   *  nonzero elements (sparsity pattern). */
  int (*ObjGradientNNZ)(void* p_user_data, int i);

  /** Feed gradient for objective \a i.
   *  Should include entries for all potentially
   *  nonzero elements (sparsity pattern).
   *
   *  Implementation skeleton:
   *      for (size_t j=0; j<obj_grad_size[i]; ++j)
   *        NLW2_WriteSparseDblEntry(p_api_data,
   *            obj_grad_index[i][j], obj_grad_value[i][j]);
   */
  void (*FeedObjGradient)(
      void* p_user_data, int i, void* p_api_data);

  /** Feed nonlinear expression of objective \a i.
   *
   *  The default puts a constant 0.
   *
   *  Implementation example:
   *      ew.EPut(obj_root_expr[i]); */
//  template <class ObjExprWriter>
//  void FeedObjExpression(int i, ObjExprWriter& ) { }


  ///////////////////// 3. DEFINED VARIABLES /////////////////////
  /** Defined variables.
   *
   *  Classical NL writes first the defined variables
   *  which are used in several places (constraints and/or
   *  objectives). Defined variables used in a single place
   *  (1 constraint, or 1 objective), are written
   *  just before the expression tree of their usage.
   *
   *  For most solvers, this requirement can be ignored
   *  and this method can return all defined variables
     *  in the first group (for \a i=0).
   *
   *	The method is guaranteed to be called in the following order:
   *		1. For \a i=0;
   *		2. For \a i>0, increasing, before constraint \a (i-1)'s expression;
   *		3. For \a i<0, decreasing, before objective \a (-i-1)'s expression.
   *
   *  @param i:
   *		- For \a i=0, feed a sequence of defined variables
   *			used in several constraints and/or objectives.
   *		- For \a i>0, feed the defined variables used solely
   *			in constraint \a i-1.
   *		- For \a i<0, feed the defined variables used solely
   *			in objective \a -i-1.
   *
   *  Implementation skeleton:
   *      // dvar_index in num_vars..num_vars+num_defvars-1.
   *      for (int dvar_index: dvar_indexes[i]) {
   *        auto dv = dvw.StartDefVar(dvar_index, lin_nnz, name_or_comment);
   *        /////////// Write the linear part:
   *        auto linw = dv.GetLinExprWriter();
   *        for (int i=0; i<lin_nnz; ++i)
   *          linw.Write(linexp_var[i], linexp_coef[i]);
   *        /////////// Write the expression tree:
   *        auto ew = dv.GetExprWriter();
   *        ew.EPut(root_expr);
   *      }
   */
//  template <class DefVarWriterFactory>
//  void FeedDefinedVariables(int i, DefVarWriterFactory& ) { }


  ///////////////////// 4. VARIABLE BOUNDS /////////////////////
  /** Bounds for variables (except defined variables).
   *  Use +-inf for missing lower and/or upper bounds.
   *  Note that variable type is given by variable ordering,
   *  see NLHeader.
   *
   *  Implementation skeleton:
   *      for (int i = 0; i < hdr.num_vars; i++)
   *        NLW2_WriteVarLbUb(p_api_data, lb[i], ub[i]);
   */
  void (*FeedVarBounds)(void* p_user_data, void* p_api_data);


  ///////////////// 5. CONSTRAINT BOUNDS & COMPLEMENTARITY ///////

  /** Bounds/complementarity for all algebraic constraints
   *  (\a num_algebraic_cons).
   *
   *  Implementation skeleton:
   *      for (int j=0; j<hdr.num_algebraic_cons; j++) {
   *        NLW2_AlgConRange_C bnd;
   *        if (compl_var && compl_var[j]) {
   *          j = compl_var[j]-1;
   *          bnd.k = 0;
   *          if (vlb[j] > negInfinity)
   *            bnd.k = 1;
   *          if (vub[j] < Infinity)
   *            bnd.k |= 2;
   *          assert(bnd.k);
   *          bnd.cvar = j;
   *        } else {
   *          bnd.k = 0;
   *          bnd.L = clb[j];
   *          bnd.U = cub[j];
   *        }
   *        NLW2_WriteAlgConRange(p_api_data, &bnd);
   *      }
   */
  void (*FeedConBounds)(void* p_user_data, void* p_api_data);


  ///////////////////// 6. CONSTRAINTS /////////////////////
  /** Description of constraint \a i
   *    (\a i in 0..num_algebraic_cons+num_logical_cons-1).
   *  With WantNLComments()==true, this is
   *  written to text-format NL as a comment. */
  const char* (*ConDescription)(void *p_user_data, int );

  /** Number of nonzeros in the linear part of constraint \a i.
   *  Should include entries for all potentially
   *  nonzero elements (sparsity pattern). */
  int (*LinearConExprNNZ)(void* p_user_data, int i);

  /** Feed the linear part of algebraic constraint \a i.
    * For smooth solvers, should contain entries for all
    * potential nonzeros (Jacobian sparsity pattern).
    *
    *  Implementation skeleton:
    *    for (size_t j=0; j<con_grad.size(); ++j)
    *      NLW2_WriteSparseDblEntry(p_api_data,
    *        con_grad[j].var_index, con_grad[j].coef);
    *    }
    */
  void (*FeedLinearConExpr)(
      void* p_user_data, int i, void* p_api_data);

  /** Feed nonlinear expression of constraint \a i.
   *  Algebraic constraints (num_algebraic_cons)
   *  come before logical (num_logical_cons).
   *  For linear constraints, the expression should be
   *  constant 0.
   */
//  void FeedConExpression(int i, ConExprWriter& ) { }


  ///////////////////// 7. EXPRESSIONS /////////////////////
  /** Feed native expression.
     *  This method is recursively called from NLWriter,
     *  when Feeder uses ExprWriter::EPut().
     *  Feeder should not call this method itself.
     *
     *  Details of ExprWriter: see NLWriter2.
   */
//  void FeedExpr(Expr e, ExprWriter& ) { }


  ///////////////////// 8. PL-SOS CONSTRAINTS ////////////
  /**
   *  The below feature is for AMPL's internal
   *  linearization of piecewise-linear functions.
   *  For user-definable SOS constraints, use suffixes
   *  .sosno/.ref.
   *
   *  The below is a feeder interface
   *  for .sos/.sosref suffixes.
   *  The feeder can provide 3 sparse vectors:
   *  - .sos for variables:
   *    Each nonzero value defines SOS group number.
   *    Negative means SOS Type 2, positive - SOS Type 1.
   *  - .sos for constraints:
   *    Each nonzero value denotes a constraint used in a
   *    linearization of an SOS. The constraint can be deleted
   *    by the solver driver if using solver's SOS.
   *  - .sosref for variables:
   *    SOS weights. Variables participating in an SOS having
   *    zero weights are involved in linearization and can be
   *    deleted if the solver accepts SOS natively.
   *
   *  Implementation:
   *      auto sosv = plsos.StartSOSVars(nvsos);
   *      for (int i=0; i<nvsos; ++i)
   *        sosv.Write(i, vsos[i]);
   *      if (ncsos) {
   *        auto sosc = plsos.StartSOSCons(ncsos);
   *        for ....
   *      }
   *      auto sosrefv = plsos.StartSOSREFVars(ac->nsosref);
   *      ....
  */
//  void FeedPLSOS(PLSOSWriter& ) { }


  ///////////////////// 9. FUNCTIONS /////////////////////
  /** Function definition. */
//  struct FuncDef {
//    const char* Name() { return ""; }
//    int NumArgs() { return 0; }
//    /** Function type.
//     *  0 - numeric;
//     *  1 - symbolic. */
//    int Type() { return 0; }
//  };

  /** Provide definition
   *  of function \a i, i=0..num_funcs-1. */
//  FuncDef Function(int i) { return {}; }


  ///////////////////// 10. RANDOM VARIABLES /////////////////////
  /// Random variables.
  /// Undocumented feature. SNL2006.
  /// Example:
  /// var z >= 0;
  ///	let z.stage := 1;
  ///	var x{0..1, 0..1} random := Uniform(0,2);
  ///	for {i in 0..1, j in 0..1} {let x[i,j].stage := 1;};
  ///	display z.stage, x.stage;
  ///	c: z * sum{i in 0..1, j in 0..1} x[i,j] <= 3 + Sample(Uniform(0,2));
  ///
  /// Feed random variables.
  /// Indexes: num_vars+num_common_exprs
  ///   .. num_vars+num_common_exprs+num_rand_vars-1.
  ///
  /// Implementation skeleton:
  ///     for(j = num_vars+num_common_exprs;
  ///         j < num_vars+num_common_exprs+num_rand_vars; j++) {
  ///       auto ew = rvw.StartRandVar(j, rand_var_comment(j));
  ///       ew.EPut(rand_var_root_expr(j));
  ///     }
//  void FeedRandomVariables(RandVarWriterFactory& ) { }


  ///////////////////// 11. COLUMN SIZES /////////////////////

  /** Jacobian column sizes (with potential nonzeros).
   *  Should feed column sizes
   *  for all but the last variable.
   *
   *  Implementation skeleton:
   *      if (WantColumnSizes())
   *        for (int i=0; i < num_vars+num_rand_vars-1; ++i)
   *          NLW2_WriteColSize(col_size[i]);
   */
  void (*FeedColumnSizes)(void* p_user_data, void* p_api_data);


  ///////////////////// 12. INITIAL GUESSES /////////////////////
  /** Initial primal guesses.
   *
   *  Implementation: write all meaningful entries (incl. zeros.)
   *      for (size_t i=0; i<n_ini_guess; ++i)
   *        NLW2_WriteSparseDblEntry(
   *          p_api_data, ini_index[i], ini_value[i]);
   */
  int (*InitialGuessesNNZ)(void* p_user_data);
  void (*FeedInitialGuesses)(void* p_user_data, void* p_api_data);

  /// Initial dual guesses
  int (*InitialDualGuessesNNZ)(void* p_user_data);
  void (*FeedInitialDualGuesses)(
      void* p_user_data, void* p_api_data);


  ///////////////////// 13. SUFFIXES /////////////////////
  /** Feed suffixes.
   *
   *  For constraints, assume ordering:
   *  first algebraic, then logical.
   *
   *  Implementation: write all non-0 entries (0 is the default.)
   *      while (....) {
   *        void* p_api_2 = NLW2_StartIntSuffix(  // or ...DblSuffix
   *          p_api_data, suf_name, kind, n_nonzeros);
   *        for (int i=0; i<n_nonzeros; ++i)
   *          NLW2_WriteSparseIntEntry(p_api_2,   // or ...DblEntry
   *            index[i], value[i]);              // ^<- p_api_2 here
   *      }
   */
  void (*FeedSuffixes)(void* p_user_data, void* p_api_data);


  //////////////////// 14. ROW/COLUMN NAMES ETC /////////////////////

  /// Want row/obj names?
  int want_row_and_obj_names_;
  /// Want del row names?
  int want_del_row_names_;
  /// Want col names?
  int want_col_names_;
  /// Want unused var names?
  int want_unused_var_names_;
  /// Want fixed var names?
  int want_fixed_var_names_;
  /// Want objective offsets?
  int want_obj_adj_;

  /** FeedRowAndObjNames:
   *  Provide constraint, then objective names.
   *
   *  Implementation:
   *    for (i: {algcons, logcons, objs})
   *      NLW2_WriteName( p_api_data, con_obj_name[i] );
   */
  void (*FeedRowAndObjNames)(void* p_user_data, void* p_api_data);

  /** Provide deleted row names.*/
  void (*FeedDelRowNames)(void* p_user_data, void* p_api_data);

  /** Provide variable names. */
  void (*FeedColNames)(void* p_user_data, void* p_api_data);

  /** Provide unused variable names. */
  void (*FeedUnusedVarNames)(void* p_user_data, void* p_api_data);

  /** Provide {fixed variable, extra info} pairs.
   *  This includes defined eliminated variables.
   *
   *  Implementation:
   *        for (....)
   *          NLW2_WriteNameAndComment(
   *            p_api_data, name[i], comment[i] );
   */
  void (*FeedFixedVarNames)(void* p_user_data, void* p_api_data);

  /** Provide {obj name, constant term} pairs.
   *
   *  Implementation:
   *        for (....)
   *          NLW2_WriteNameAndNumber(
   *            p_api_data, name[i], obj_offset[i] );
   */
  void (*FeedObjAdj)(void* p_user_data, void* p_api_data);

} NLW2_NLFeeder_C;


/// Return NLW2_NLFeeder_C with default options / methods
NLW2_NLFeeder_C NLW2_MakeNLFeeder_C_Default(void);

/// Destroy NLW2_NLFeeder_C created by NLW2_MakeNLFeeder_C_default()
void NLW2_DestroyNLFeeder_C_Default(NLW2_NLFeeder_C* );

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLFeeder_C_H

#include <stdlib.h>
#include <assert.h>

#include "nlsol_ex_c_nl.h"

NLHeader_C CAPI_ex_Header(void* pex_void) {
  CAPIExample* pex = (CAPIExample*)pex_void;
  // Get default header
  NLHeader_C hdr = MakeNLHeader_C_Default();

  hdr.pi.num_vars = pex->n_var;
  hdr.pi.num_algebraic_cons = pex->n_con;
  hdr.pi.num_objs = pex->n_obj;
  hdr.pi.num_ranges = 0;
  hdr.pi.num_eqns = 0;
  hdr.pi.num_logical_cons = 0;

  // Setting some other common data (more may be needed)

  /** Total number of nonlinear constraints. */
  hdr.pi.num_nl_cons = 0;
  hdr.pi.num_nl_objs = 0;
  hdr.pi.num_compl_conds = 0;
  hdr.pi.num_nl_compl_conds = 0;
  hdr.pi.num_compl_dbl_ineqs = 0;
  hdr.pi.num_compl_vars_with_nz_lb = 0;

  /** Number of nonlinear network constraints. */
  hdr.pi.num_nl_net_cons = 0;
  hdr.pi.num_linear_net_cons = 0;

  /**
    Number of nonlinear variables in constraints including nonlinear
    variables in both constraints and objectives.
   */
  hdr.pi.num_nl_vars_in_cons = 0;

  /**
    Number of nonlinear variables in objectives including nonlinear
    variables in both constraints and objectives.
   */
  hdr.pi.num_nl_vars_in_objs = 0;

  /** Number of nonlinear variables in both constraints and objectives. */
  hdr.pi.num_nl_vars_in_both = 0;

  // Miscellaneous
  // -------------

  /** Number of linear network variables (arcs). */
  hdr.pi.num_linear_net_vars = 0;

  /** Number of functions. */
  hdr.pi.num_funcs = 0;

  // Information about discrete variables
  // ------------------------------------

  /** Number of linear binary variables. */
  hdr.pi.num_linear_binary_vars = 0;

  /** Number of linear non-binary integer variables. */
  hdr.pi.num_linear_integer_vars = pex->n_var_int;

  /**
    Number of integer nonlinear variables in both constraints and objectives.
   */
  hdr.pi.num_nl_integer_vars_in_both = 0;

  /** Number of integer nonlinear variables just in constraints. */
  hdr.pi.num_nl_integer_vars_in_cons = 0;

  /** Number of integer nonlinear variables just in objectives. */
  hdr.pi.num_nl_integer_vars_in_objs = 0;

  // Information about nonzeros
  // --------------------------

  /** Number of nonzeros in constraints' Jacobian. */
  hdr.pi.num_con_nonzeros = pex->n_con_nz;

  /** Number of nonzeros in all objective gradients. */
  hdr.pi.num_obj_nonzeros = pex->n_obj_nz;

  // Information about names
  // -----------------------

  /** Length of longest con/obj name if names are available. */
  hdr.pi.max_con_name_len = 0;    // no need to set

  /** Length of longest variable name if names are available. */
  hdr.pi.max_var_name_len = 0;    // no need to set

  // Information about common expressions
  // ------------------------------------

  /**
    Number of common expressions that appear both in constraints
    and objectives.
   */
  hdr.pi.num_common_exprs_in_both = 0;

  /**
    Number of common expressions that appear in multiple constraints
    and don't appear in objectives.
   */
  hdr.pi.num_common_exprs_in_cons = 0;

  /**
    Number of common expressions that appear in multiple objectives
    and don't appear in constraints.
   */
  hdr.pi.num_common_exprs_in_objs = 0;

  /**
    Number of common expressions that only appear in a single constraint
    and don't appear in objectives.
   */
  hdr.pi.num_common_exprs_in_single_cons = 0;

  /**
    Number of common expressions that only appear in a single objective
    and don't appear in constraints.
   */
  hdr.pi.num_common_exprs_in_single_objs = 0;


  // Technical
  hdr.nli.format
      = pex->binary_nl ? NL_FORMAT_BINARY : NL_FORMAT_TEXT;

  hdr.nli.prob_name = "c_api_example_model";

  return hdr;
}

const char* ObjDescription(void* p_user_data, int i) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  return pex->obj_name;
}

int ObjType(void* p_user_data, int i) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  return pex->obj_sense;
}

int ObjGradientNNZ(void* p_user_data, int i) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  return pex->n_obj_nz;
}

void FeedObjGradient(
    void* p_user_data, int i, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  assert(0==i);
  for (int j=0; j<pex->n_obj_nz; ++j)
    NLW2_WriteSparseDblEntry(p_api_data,
                             pex->obj_linpart[j].index_,
                             pex->obj_linpart[j].value_);
}

void FeedVarBounds(void* p_user_data, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  for (int i = 0; i < pex->n_var; i++)
    NLW2_WriteVarLbUb(p_api_data,
                      pex->var_lb[i], pex->var_ub[i]);
}

void FeedConBounds(void* p_user_data, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  NLW2_AlgConRange_C bnd;
  for (int j=0; j<pex->n_con; j++) {
    bnd.k = 0;
    bnd.L = pex->con_lb[j];
    bnd.U = pex->con_ub[j];
    NLW2_WriteAlgConRange(p_api_data, &bnd);
  }
}

const char* ConDescription(void* p_user_data, int i) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  return pex->con_name[i];
}

int LinearConExprNNZ(void* p_user_data, int i) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  return pex->row_nnz[i];
}

void FeedLinearConExpr(
    void* p_user_data, int i, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  assert(i<pex->n_con);
  for (int j=0; j<pex->row_nnz[i]; ++j)
    NLW2_WriteSparseDblEntry(p_api_data,
                             pex->con_linpart[i][j].index_,
                             pex->con_linpart[i][j].value_);
}

//  void FeedConExpression(int i, ConExprWriter& ) { }


  ///////////////////// 7. EXPRESSIONS /////////////////////
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
void FeedColumnSizes(void* p_user_data, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  if (1 /*WantColumnSizes()*/)
    for (int i=0; i < pex->n_var-1; ++i)
      NLW2_WriteColSize(p_api_data, pex->col_sizes[i]);
}


  ///////////////////// 12. INITIAL GUESSES /////////////////////

int CountNNZ(const double* parray, int len) {
  int result = 0;
  for (int i=0; i<len; ++i)
    if (parray[i])
      ++result;
  return result;
}
void WriteDenseAsSparse(
    const double* parray, int len, void* p_api_data) {
  for (int i=0; i<len; ++i)
    if (parray[i])
      NLW2_WriteSparseDblEntry(p_api_data, i, parray[i]);
}

  /**
   *  Implementation:
   *      for (size_t i=0; i<n_ini_guess; ++i)
   *        NLW2_WriteSparseDblEntry(p_api_data, i, ini_guess[i]);
   */
int InitialGuessesNNZ(void* p_user_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  return CountNNZ(pex->ini_x, pex->n_var);
}
void FeedInitialGuesses(void* p_user_data, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  WriteDenseAsSparse(pex->ini_x, pex->n_var, p_api_data);
}

/// Initial dual guesses
int InitialDualGuessesNNZ(void* p_user_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  return CountNNZ(pex->ini_y, pex->n_con);
}
void FeedInitialDualGuesses(
    void* p_user_data, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  WriteDenseAsSparse(pex->ini_y, pex->n_con, p_api_data);
}


  ///////////////////// 13. SUFFIXES /////////////////////
  /** Feed suffixes.
   *
   *  For constraints, assume ordering:
   *  first algebraic, then logical.
   *
   *  Implementation:
   *      while (....) {
   *        void* p_api_2 = NLW2_StartIntSuffix(  // or ...DblSuffix
   *          p_api_data, suf_name, kind, n_nonzeros);
   *        for (int i=0; i<n_nonzeros; ++i)
   *          NLW2_WriteSparseIntEntry(p_api_2,   // <- new API pointer
   *            index[i], value[i]);              // or ...DblEntry
   *      }
   */
void FeedSuffixes(void* p_user_data, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  for (int i_suf=0; i_suf<pex->n_suf; ++i_suf) {
    const Suffix* ps = &pex->suf[i_suf];
    void* p_api_2 = NLW2_StartDblSuffix(
          p_api_data, ps->name_, ps->kind_, ps->n_val_);
    for (int i=0; i<ps->n_val_; ++i)
      NLW2_WriteSparseDblEntry(p_api_2,
                               ps->values_[i].index_,
                               ps->values_[i].value_);
  }
}


  //////////////////// 14. ROW/COLUMN NAMES ETC /////////////////////
void FeedRowAndObjNames(void* p_user_data, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  for (int i=0; i<pex->n_con; ++i)
    NLW2_WriteName(p_api_data, pex->con_name[i]);
  NLW2_WriteName(p_api_data, pex->obj_name);
}

  /** Provide deleted row names.*/
//  void FeedDelRowNames(DelRowNameWriter& ) { }

  /** Provide variable names. */
void FeedColNames(void* p_user_data, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  for (int i=0; i<pex->n_var; ++i)
    NLW2_WriteName(p_api_data, pex->var_name[i]);
}

  /** Provide unused variable names. */
//  void FeedUnusedVarNames(UnusedVarNameWriter& ) { }

  /** Provide {fixed variable, extra info} pairs.
   *  This includes defined eliminated variables.
   *
   *  Implementation:
   *      if ((output_desired) && wrt)
   *        for (....)
   *          wrt << typename Writer::StrStrValue
   *          { name[i].c_str(), comment[i].c_str() };
   */
//  void FeedFixedVarNames(FixedVarNameWriter& ) { }

  /** Provide {obj name, constant term} pairs.
   *
   *  Implementation:
   *      if (wrt)
   *        for (....)
   *          wrt << typename Writer::StrDblValue
   *          { name[i].c_str(), (double)obj_offset[i] };
   */
//  void FeedObjAdj(ObjOffsetWriter& ) { }


NLW2_NLFeeder_C MakeNLFeeder_C(
    CAPIExample* pex, int binary) {
  NLW2_NLFeeder_C result        // Fill with default values
      = NLW2_MakeNLFeeder_C_Default();

  result.p_user_data_ = pex;
  pex->binary_nl = binary;

  // Header feeder
  result.Header = CAPI_ex_Header;

  // Change some options
  result.want_nl_comments_ = 1;

  // Objective
  result.ObjDescription = ObjDescription;
  result.ObjType = ObjType;
  result.ObjGradientNNZ = ObjGradientNNZ;
  result.FeedObjGradient = FeedObjGradient;

  // Reuse default-generated FeedObjExpr, FeedDefVars

  result.FeedVarBounds = FeedVarBounds;
  result.FeedConBounds = FeedConBounds;

  // Constraints
  result.ConDescription = ConDescription;
  result.LinearConExprNNZ = LinearConExprNNZ;
  result.FeedLinearConExpr = FeedLinearConExpr;

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
  result.FeedColumnSizes = FeedColumnSizes;


    ///////////////////// 12. INITIAL GUESSES /////////////////////
    /** Initial primal guesses.
     *
     *  Implementation:
     *      if (ini_guess.size()) {
     *        auto ig = igw.MakeVectorWriter(ini_guess.size());
     *        for (size_t i=0; i<ini_guess.size(); ++i)
     *          ig.Write(i, ini_guess[i]);
     *      }
     */
  result.InitialGuessesNNZ = InitialGuessesNNZ;
  result.InitialDualGuessesNNZ = InitialDualGuessesNNZ;
  result.FeedInitialGuesses = FeedInitialGuesses;
  result.FeedInitialDualGuesses = FeedInitialDualGuesses;

    ///////////////////// 13. SUFFIXES /////////////////////
    /** Feed suffixes.
     *
     *  For constraints, assume ordering:
     *  first algebraic, then logical.
     *
     *  Implementation:
     *      while (....) {
     *        auto sw = swf.StartIntSuffix(  // or ...DblSuffix
     *          suf_name, kind, n_nonzeros);
     *        for (int i=0; i<n_nonzeros; ++i)
     *          sw.Write(index[i], value[i]);
     *      }
     */
  result.FeedSuffixes = FeedSuffixes;

    //////////////////// 14. ROW/COLUMN NAMES ETC /////////////////////
  result.want_row_and_obj_names_ = 1;
  result.want_col_names_ = 1;

  result.FeedRowAndObjNames = FeedRowAndObjNames;
  result.FeedColNames = FeedColNames;


  return result;
}

void DestroyNLFeeder_C(NLW2_NLFeeder_C* pf) {
  pf->p_user_data_ = NULL;

  NLW2_DestroyNLFeeder_C_Default(pf);
}

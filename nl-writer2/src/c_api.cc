/**
 * NL Writer C API implementation
 *
 */

#include <cstdlib>
#include <functional>
#include <cassert>

#include "api/c/nl-feeder2-c.h"
#include "api/c/sol-handler2-c.h"
#include "api/c/nl-writer2-misc-c.h"
#include "api/c/nlsol-c.h"

#include "api/c/nl-feeder2-c-impl.h"
#include "api/c/sol-handler2-c-impl.h"
#include "api/c/nl-writer2-misc-c-impl.h"
#include "mp/nlsol.h"

#include "mp/nl-writer2.h"
#include "mp/nl-writer2.hpp"
#include "mp/sol-reader2.h"
#include "mp/sol-reader2.hpp"

#ifdef __cplusplus  // Implementing C API from C++
extern "C" {
#endif

///////////////////////// NLFeeder_C ///////////////////////////

/// Sparse vector writers
void NLW2_WriteSparseDblEntry(
    void* svw, int index, double value) {
  auto& f = *(std::function<void(int, double)>*)(svw);
  f(index, value);
}

/// Var bound writer
void NLW2_WriteVarLbUb(void* vbw, double lb, double ub) {
  auto& f = *(std::function<void(double, double)>*)(vbw);
  f(lb, ub);
}

/// Algebraic constraint writer
void NLW2_WriteAlgConRange(void* crw, NLW2_AlgConRange_C* pbnd) {
  auto& f = *(std::function<void(NLW2_AlgConRange_C* )>*)(crw);
  f(pbnd);
}

void NLW2_WriteColSize(void* csw, int sz) {
  auto& f = *(std::function<void(int s)>*)(csw);
  f(sz);
}


/// Default implementations
static const char* NLW2_ObjDescription_C_Default(void* , int )
{ return ""; }
static int NLW2_ObjType_C_Default(void* , int ) { return 0; }
static int NLW2_ObjGradientNNZ_C_Default(void* , int ) { return 0; }
static void NLW2_FeedObjGradient_C_Default(void* , int , void* ) { }

// ObjExpr

// DefVars

static void NLW2_FeedVarBounds_C_Default(void* , void* ) { }
static void NLW2_FeedConBounds_C_Default(void* , void* ) { }

static const char* NLW2_ConDescription_C_Default(void *, int )
{ return ""; }
static int NLW2_LinearConExprNNZ_C_Default(void* , int ) { return 0; }
static void NLW2_FeedLinearConExpr_C_Default(void* , int , void* ) { }

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

  /** Jacobian column sizes.
   *  Should feed LP column sizes
   *  for all but the last variable.
   *
   *  Implementation skeleton:
   *      if (WantColumnSizes())
   *        for (int i=0; i < num_vars+num_rand_vars-1; ++i)
   *          csw.Write(col_size[i]);
   */
static void NLW2_FeedColumnSizes_C_Default(void* , void* )
{ assert(0 && "this needs implementation"); }


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
//  void FeedInitialGuesses(IGWriter& ) { }

  /** Initial dual guesses. */
//  void FeedInitialDualGuesses(IDGWriter& ) { }


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
//  void FeedSuffixes(SuffixWriterFactory& ) { }


  //////////////////// 14. ROW/COLUMN NAMES ETC /////////////////////
  /** FeedRowAndObjNames:
   *  Provide constraint, then objective names.
   *  Name information is optional.
   *
   *  Implementation:
   *      if ((output_desired) && wrt)
   *        for (i: ....)
   *          wrt << name[i].c_str();
   */
//  void FeedRowAndObjNames(RowObjNameWriter& wrt) { }

  /** Provide deleted row names.*/
//  void FeedDelRowNames(DelRowNameWriter& ) { }

  /** Provide variable names. */
//  void FeedColNames(ColNameWriter& ) { }

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


NLW2_NLFeeder2_C NLW2_MakeNLFeeder2_C_Default() {
  NLW2_NLFeeder2_C result;

  std::memset(&result, 0, sizeof(result));       // all 0

  result.p_user_data_ = NULL;

  result.Header = NULL;           // User should provide this

  // Default options
  result.want_nl_comments_ = 0;
  result.output_precision_ = 0;
  result.want_bounds_first_ =1;
  result.want_column_sizes_ =1;

  // Objectives
  result.ObjDescription = NLW2_ObjDescription_C_Default;
  result.ObjType = NLW2_ObjType_C_Default;
  result.ObjGradientNNZ = NLW2_ObjGradientNNZ_C_Default;
  result.FeedObjGradient = NLW2_FeedObjGradient_C_Default;
  // ObjExpr... relying on NLFeeder2's default (0)

  // DefVars...

  result.FeedVarBounds = NLW2_FeedVarBounds_C_Default;
  result.FeedConBounds = NLW2_FeedConBounds_C_Default;

  // Constraints
  result.ConDescription = NLW2_ConDescription_C_Default;
  result.LinearConExprNNZ = NLW2_LinearConExprNNZ_C_Default;
  result.FeedLinearConExpr = NLW2_FeedLinearConExpr_C_Default;
  // ConExpr...

  //  void FeedConExpression(int i, ConExprWriter& ) { }


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
  result.FeedColumnSizes = NLW2_FeedColumnSizes_C_Default;


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
  //  void FeedInitialGuesses(IGWriter& ) { }

    /** Initial dual guesses. */
  //  void FeedInitialDualGuesses(IDGWriter& ) { }


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
  //  void FeedSuffixes(SuffixWriterFactory& ) { }


    //////////////////// 14. ROW/COLUMN NAMES ETC /////////////////////
    /** FeedRowAndObjNames:
     *  Provide constraint, then objective names.
     *  Name information is optional.
     *
     *  Implementation:
     *      if ((output_desired) && wrt)
     *        for (i: ....)
     *          wrt << name[i].c_str();
     */
  //  void FeedRowAndObjNames(RowObjNameWriter& wrt) { }

    /** Provide deleted row names.*/
  //  void FeedDelRowNames(DelRowNameWriter& ) { }

    /** Provide variable names. */
  //  void FeedColNames(ColNameWriter& ) { }

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

  return result;
}

void NLW2_DestroyNLFeeder2_C_Default(NLW2_NLFeeder2_C* )
{ }

///////////////////////// NLUtils_C ///////////////////////////
/// log message
static void NLW2_log_message_C_Default(
    void* p_api_data, const char* format, ...) {
  va_list args;
  va_start (args, format);
  std::vprintf (format, args);
  va_end (args);
}
/// log warning
static void NLW2_log_warning_C_Default(
    void* p_api_data, const char* format, ...) {
  std::fprintf(stderr, "WARNING: ");
  va_list args;
  va_start (args, format);
  std::vfprintf (stderr, format, args);
  va_end (args);
  std::fprintf(stderr, "\n");
}
/// Override this to your error handler.
/// Not using exceptions by default.
/// Only called with wrong output format string
/// (internal error.)
static void NLW2_myexit_C_Default(
    void* p_api_data, const char* msg) {
  using namespace std;
  fprintf(stderr, "%s\n", msg);
  exit(1);
}


NLW2_NLUtils_C NLW2_MakeNLUtils_C_Default() {
  NLW2_NLUtils_C result;

  result.p_user_data_ = NULL;

  result.log_message = NLW2_log_message_C_Default;
  result.log_warning = NLW2_log_warning_C_Default;
  result.myexit = NLW2_myexit_C_Default;

  return result;
}

void NLW2_DestroyNLUtils_C_Default(NLW2_NLUtils_C* )
{ }


/////////////////// SOLHandler2_C /////////////////////

/// Callbacks

double NLW2_ReadSolVal(void* p_api_data) {
  // In contrast to NLWriter2_C,
  // here we know the type
  return ((mp::VecReader<double>*)p_api_data)
      ->ReadNext();
}


/// Default methods
static void NLW2_OnSolveMessage_C_Default(
    void* p_user_data, const char* s, int nbs) {
  if (nbs < (int)strlen(s)) {
    printf("%s\n", s+nbs);
    fflush(stdout);
  }
}
static int NLW2_OnAMPLOptions_C_Default(
    void* p_user_data, AMPLOptions_C ) { return 0; }
static void NLW2_OnDualSolution_C_Default(
    void* p_user_data, int nvals, void* p_api_data) { }
static void NLW2_OnPrimalSolution_C_Default(
    void* p_user_data, int nvals, void* p_api_data) { }
static void NLW2_OnObjno_C_Default(void* p_user_data, int ) { }
static void NLW2_OnSolveCode_C_Default(void* p_user_data, int ) { }
//  void OnIntSuffix(SuffixReader& ) { }
//  void OnDblSuffix(SuffixReader& ) { }


NLW2_SOLHandler2_C NLW2_MakeSOLHandler2_C_Default() {
  NLW2_SOLHandler2_C result;

  std::memset(&result, 0, sizeof(result));       // all 0

  result.p_user_data_ = NULL;

  result.Header = NULL;           // User should provide this

  result.OnSolveMessage = NLW2_OnSolveMessage_C_Default;
  result.OnAMPLOptions = NLW2_OnAMPLOptions_C_Default;
  result.OnDualSolution = NLW2_OnDualSolution_C_Default;
  result.OnPrimalSolution = NLW2_OnPrimalSolution_C_Default;
  result.OnObjno = NLW2_OnObjno_C_Default;
  result.OnSolveCode = NLW2_OnSolveCode_C_Default;
  //  void OnIntSuffix(SuffixReader& ) { }
  //  void OnDblSuffix(SuffixReader& ) { }

  return result;
}

void NLW2_DestroySOLHandler2_C_Default(NLW2_SOLHandler2_C* )
{ }


//////////// NLSOL_C API //////////////

namespace mp {

/// Typedef our specialization of NLSOL
using NLSOL_C_Impl
  = mp::NLSOL<mp::NLFeeder2_C_Impl,
      mp::SOLHandler2_C_Impl>;

}  // namespace mp

/// Construct.
///
/// Note that the argument objects are stored by value.
NLW2_NLSOL_C NLW2_MakeNLSOL_C(
    NLW2_NLFeeder2_C* pnlf, NLW2_SOLHandler2_C* psolh, NLW2_NLUtils_C* putl) {
  NLW2_NLSOL_C result;

  result.p_nlf_ = new mp::NLFeeder2_C_Impl(pnlf);
  result.p_solh_ = new mp::SOLHandler2_C_Impl(psolh);
  result.p_utl_ = new mp::NLUtils_C_Impl(putl);
  result.p_nlsol_
      = new mp::NLSOL_C_Impl(
        *(mp::NLFeeder2_C_Impl*)result.p_nlf_,
        *(mp::SOLHandler2_C_Impl*)result.p_solh_,
        *(mp::NLUtils_C_Impl*)result.p_utl_);

  return result;
}

/// Destroy
void NLW2_DestroyNLSOL_C(NLW2_NLSOL_C* pnls) {
  delete (mp::NLUtils_C_Impl*)(pnls->p_utl_);
  delete (mp::SOLHandler2_C_Impl*)(pnls->p_solh_);
  delete (mp::NLFeeder2_C_Impl*)(pnls->p_nlf_);
  delete (mp::NLSOL_C_Impl*)(pnls->p_nlsol_);
}

/// Set solver, such as "gurobi", "highs", "ipopt"
void NLW2_SetSolver(NLW2_NLSOL_C* pnls, const char* solver) {
  ((mp::NLSOL_C_Impl*)(pnls->p_nlsol_))
      ->SetSolver(solver);
}

/// Set solver options, such as "outlev=1 lim:time=500"
void NLW2_SetSolverOptions(NLW2_NLSOL_C* pnls, const char* sopts) {
  ((mp::NLSOL_C_Impl*)(pnls->p_nlsol_))
      ->SetSolverOptions(sopts);
}

/// Solve.
/// @param filestub: filename stub to be used
/// for input files (.nl, .col., .row, etc.),
/// and output files (.sol).
/// @return true if all ok.
int NLW2_Solve(NLW2_NLSOL_C* pnls, const char* filestub) {
  return ((mp::NLSOL_C_Impl*)(pnls->p_nlsol_))
      ->Solve(filestub);
}

/// Get error message.
const char* NLW2_GetErrorMessage(NLW2_NLSOL_C* pnls) {
  return ((mp::NLSOL_C_Impl*)(pnls->p_nlsol_))
      ->GetErrorMessage();
}

/// Substep: write NL and any accompanying files.
int NLW2_WriteNLFile(NLW2_NLSOL_C* pnls, const char* filestub) {
  return ((mp::NLSOL_C_Impl*)(pnls->p_nlsol_))
      ->WriteNLFile(filestub);
}

/// Substep: invoke chosen solver for \a filestub.
int NLW2_InvokeSolver(NLW2_NLSOL_C* pnls, const char* filestub) {
  return ((mp::NLSOL_C_Impl*)(pnls->p_nlsol_))
      ->InvokeSolver(filestub);
}

/// Substep: read solution.
/// @param filename: complete file name,
/// normally (stub).sol.
int NLW2_ReadSolution(NLW2_NLSOL_C* pnls, const char* filename) {
  return ((mp::NLSOL_C_Impl*)(pnls->p_nlsol_))
      ->ReadSolution(filename);
}


#ifdef __cplusplus
}  // extern "C"
#endif
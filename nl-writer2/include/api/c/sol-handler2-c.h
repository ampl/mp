/**
 * C API: extern "C" wrappers for the SOLHandler2 interface,
 * as well as SOLReader2 calls.
 *
 */

#ifndef SOLHANDLER2C_H
#define SOLHANDLER2C_H

#include "mp/nl-header-c.h"

#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif

/// Callback: read next dual/variable value
double NLW2_ReadSolVal(void* p_api_data);

/**
 * AMPL internal options.
 */
typedef struct AMPLOptions_C {
  /// Actual number of options
  int n_options_;
  /// Option values
  long options_[MAX_AMPL_OPTIONS];
  /// Whether vbtol specified
  int has_vbtol_;
  /// vbtol value
  double vbtol_;
} AMPLOptions_C;


/// Wrap mp::SOLHandler2 for C API.
///
/// NLW2_SOLHandler2_C: reads solution details on request
/// via provided callback objects.
/// See the examples folder.
///
/// To fill some **default methods**,
/// call NLW2_MakeSOLHandler2_C_Default() / NLW2_Destroy...().
typedef struct NLW2_SOLHandler2_C {
  /// User data, provided to the methods as the 1st arg
  void* p_user_data_;

  /** The NLHeader used to write the NL file. */
  NLHeader_C (*Header)(void* p_user_data);

  /** Receive solve message.
   *  The message always ends with '\n'.
   *
   *  @param nbs: number of backspaces
   *  in the original solve message.
   *  So many characters should be skipped
   *  from the message if printed straightaway.
   *  AMPL solver drivers can supply the message
   *  with initial backspaces to indicate
   *  that so many characters should be skipped
   *  when printing. For example, if the driver prints
   *  MINOS 5.51:
   *  and exits, and the message starts with that again,
   *  this part should be skipped.
   */
  void (*OnSolveMessage)(
      void* p_user_data, const char* s, int nbs);

  /**
   * Can be ignored by external systems.
   * @return non-zero to stop solution input.
   */
  int (*OnAMPLOptions)(
      void* p_user_data, AMPLOptions_C );

  /**
   * Dual values for algebraic constraints,
   * if provided in the solution.
   * Number of values <= NumAlgCons().
   * Implementation:
   *
   *   duals.reserve(nvals);
   *   while (nvals--)
   *     duals.push_back(NLW2_ReadSolVal(p_api_data));
   */
  void (*OnDualSolution)(
      void* p_user_data, int nvals, void* p_api_data);

  /**
   * Variable values, if provided.
   * Number of values <= NumVars().
   */
  void (*OnPrimalSolution)(
      void* p_user_data, int nvals, void* p_api_data);

  /**
   * Receive notification of the objective index
   * used by the driver (solver option 'objno').
   */
  void (*OnObjno)(void* p_user_data, int );

  /**
   * Receive notification of the solve code.
   */
  void (*OnSolveCode)(void* p_user_data, int );

  /**
   * OnIntSuffix().
   *
   * For constraints, can include values for
   * logical constraints (after algebraic.)
   * Sparse representation - can be empty
   * (i.e., all values zero.)
   *
   * const auto& si = sr.SufInfo();
   * int kind = si.Kind();
   * int nmax = nitems_max[kind & 3];
   * const std::string& name = si.Name();
   * const std::string& table = si.Table();
   * while (sr.Size()) {
   *   std::pair<int, int> val = sr.ReadNext();
   *   if (val.first<0 || val.first>=nmax) {
   *     sr.SetError(mp::SOL_Read_Bad_Suffix,
   *       "bad suffix element index");
   *     return;
   *   }
   *   suf[val.first] = val.second;
   * }
   * if (mp::SOL_Read_OK == sr.ReadResult())    // Can check
   *   RegisterSuffix(kind, name, table, suf);
   */
//  template <class SuffixReader>
//  void OnIntSuffix(SuffixReader& ) { }

  /**
   * Same as OnIntSuffix(), but
   * sr.ReadNext() returns pair<int, double>
   */
//  template <class SuffixReader>
//  void OnDblSuffix(SuffixReader& ) { }

} NLW2_SOLHandler2_C;

/// Create an SOLHandler2_C with default methods
NLW2_SOLHandler2_C NLW2_MakeSOLHandler2_C_Default();

/// Destroy an SOLHandler2_C
/// created with NLW2_Make...Default()
void NLW2_DestroySOLHandler2_C_Default(NLW2_SOLHandler2_C* );

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // SOLHANDLER2C_H

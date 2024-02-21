/*
 C API: extern "C" wrappers for the SOLHandler interface,
 as well as SOLReader2 calls.

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

 Author: Gleb Belov
 */

#ifndef SOLHandlerC_H
#define SOLHandlerC_H

#include "mp/nl-header-c.h"

#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif

/// Callback: read next dual/variable value
double NLW2_ReadSolVal(void* p_api_data);

/// Suffix information
typedef struct NLW2_SuffixInfo_C {
  /// Suffix kind
  int kind_;
  /// Name
  const char* name_;
  /// Suffix table, if provided
  const char* table_;
} NLW2_SuffixInfo_C;

/// Number of suffix non-zero elements
int NLW2_IntSuffixNNZ(void* p_api_data);
/// Number of suffix non-zero elements
int NLW2_DblSuffixNNZ(void* p_api_data);
/// Read suffix entry
void NLW2_ReadIntSuffixEntry(
    void* p_api_data, int* , int*);
/// Read suffix entry
void NLW2_ReadDblSuffixEntry(
    void* p_api_data, int* , double*);
/// Report suffix error.
/// This causes NLW2_DblSuffixNNZ() to return 0.
void NLW2_ReportDblSuffixError(
    void* p_api_data, const char* msg);
/// Report suffix error.
/// This causes NLW2_IntSuffixNNZ() to return 0.
void NLW2_ReportIntSuffixError(
    void* p_api_data, const char* msg);
/// Check suffix read result
int NLW2_IntSuffixReadOK(void* p_api_data);
/// Check suffix read result
int NLW2_DblSuffixReadOK(void* p_api_data);


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


/// Wrap mp::SOLHandler for C API.
///
/// NLW2_SOLHandler_C: reads solution details on request
/// via provided callback objects.
/// See the examples folder.
///
/// To fill some **default methods**,
/// call NLW2_MakeSOLHandler_C_Default() / NLW2_Destroy...().
typedef struct NLW2_SOLHandler_C {
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
   * used by the driver (solver option 'objno'-1).
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
   * int kind = si.kind_;
   * int nmax = nitems_max[kind & 3];  // {vars, cons, objs, 1}
   * const char* name = si.name_;
   * const char* table = si.table_;
   * int i;
   * int v;
   * while (NLW2_IntSuffixNNZ(p_api_data)) {
   *   NLW2_ReadIntSuffixEntry(p_api_data, &i, &v);
   *   if (i<0 || i>=nmax) {
   *     NLW2_ReportIntSuffixError(
   *       p_api_data, "bad suffix element index");
   *     return;
   *   }
   *   suf[i] = v;
   * }
   * if (NLW2_IntSuffixReadOK(p_api_data))    // Can check
   *   RegisterSuffix(kind, name, table, suf);
   */
  void (*OnIntSuffix)(
      void* p_user_data, NLW2_SuffixInfo_C si, void* p_api_data);

  /**
   * Same as OnIntSuffix(), but
   * use NLW2_ReadDblSuffixEntry() etc.
   */
  void (*OnDblSuffix)(
      void* p_user_data, NLW2_SuffixInfo_C si, void* p_api_data);

} NLW2_SOLHandler_C;

/// Create an NLW2_SOLHandler_C with default methods
NLW2_SOLHandler_C NLW2_MakeSOLHandler_C_Default(void);

/// Destroy an NLW2_SOLHandler_C
/// created with NLW2_Make...Default()
void NLW2_DestroySOLHandler_C_Default(NLW2_SOLHandler_C* );

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // SOLHandlerC_H

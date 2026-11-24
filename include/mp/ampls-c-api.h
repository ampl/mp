#ifndef AMPLS_C_INTERFACE_H
#define AMPLS_C_INTERFACE_H

#include <stddef.h>
#include "mp/ampls-ccallbacks.h"

/*
 * A generic C API for MP
 */
#ifdef _WIN32
#define AMPLS_C_EXPORT __declspec(dllexport)
#else
#define AMPLS_C_EXPORT __attribute__((visibility("default")))
#endif


// Macro to be used in each solver to automatically declare the API functions
#define DECLARE_SOLVER_API_FUNCTIONS(name) \
    AMPLS_MP_Solver* Open_##name(CCallbacks cb); \
    AMPLS_C_EXPORT void AMPLSClose_##name(AMPLS_MP_Solver* slv); \
    AMPLS_C_EXPORT void* AMPLSGetModel_##name(AMPLS_MP_Solver* slv);

/// An AMPLS solver instance
typedef struct AMPLS_MP_Solver_T {
  /// AMPLS internal info
  void *internal_info_;
  /// Extra info, managed by the specific solver
  void *solver_info_;
  /// User info, free to assign
  void *user_info_;
} AMPLS_MP_Solver;

  /*
   * Return values of the calls that can fail (AMPLSLoadNLModel(),
   * AMPLSSolve(), AMPLSReportResults(), AMPLSReportError()).
   *
   * They tell only whether the call succeeded, not how the model
   * was solved:
   * - 0: success. This is not a solve result: a solve that ended,
   *   e.g., infeasible (200) or at a limit (400) also returns 0.
   *   The solve result is in the .sol file that AMPLSReportResults()
   *   writes, and AMPLSGetSolveResult() gives it.
   * - > 0: failure, and the error carries a solve result code
   *   (AMPL's solve_result_num), which the caller can report with
   *   AMPLSReportError(); e.g., 200 when the model conversion proved
   *   the model infeasible, 150 for a solution check failed with
   *   chk:fail. Such a code describes the failed call; no solution
   *   comes with it.
   * - -1: failure without a solve result code (report it as 500).
   * The error message is among AMPLSGetMessages().
   */

  /// Load model incl. suffixes.
  /// The method assumes that slv has been initialized by
  /// a solver-specific API using functions of ampls-cpp-api.h.
  /// @return 0 on success; on failure, the error's solve result
  ///   code if it has one (e.g., 200: the converter proved the model
  ///   infeasible), otherwise -1. See "Return values" above.
  AMPLS_C_EXPORT int AMPLSLoadNLModel(AMPLS_MP_Solver* slv,
    const char* nl_filename, char** options);

  /// Load suffixes and other option-dependent info
  AMPLS_C_EXPORT void AMPLSReadExtras(AMPLS_MP_Solver* slv);

  /// Solve the model, as the driver does (e.g., with emulated
  /// multiple objectives). Can be called again, e.g., after changing
  /// the model through the solver's own API.
  /// @return 0 if the solve ran, whatever its result (see
  ///   AMPLSGetSolveResult()); on failure, the error's solve result
  ///   code if it has one, otherwise -1. See "Return values" above.
  ///   After a failure, do not call
  ///   AMPLSReportResults(), which would report the previous solve;
  ///   report the error with AMPLSReportError() if needed.
  AMPLS_C_EXPORT int AMPLSSolve(AMPLS_MP_Solver* slv);

  /// Report results.
  /// The kind of results reported is influenced by solver option 
  /// `wantsol`.
  /// @return 0 on success (not the solve result, which is in the
  ///   .sol file); on failure, the error's solve result code if it has
  ///   one (e.g., 150: a solution check failed with chk:fail),
  ///   otherwise -1. See "Return values" above.
  AMPLS_C_EXPORT int AMPLSReportResults(AMPLS_MP_Solver* slv, const char* solFileName);

  /// Report an error as a driver does: write the .sol file
  /// with \a solve_result (e.g., 500) and \a msg, without exiting.
  /// Needs the model's file name, i.e., a call of AMPLSLoadNLModel()
  /// that got past reading the NL header (e.g., failed on an option).
  /// @return 0 on success; -1 or a code on failure, as above.
  ///   Does not change AMPLSGetSolveResult().
  AMPLS_C_EXPORT int AMPLSReportError(AMPLS_MP_Solver* slv,
                                      int solve_result, const char* msg);

  /// Solve result code (AMPL's solve_result_num) of the current
  /// state, unlike the return values above also after a success:
  /// if the last call of AMPLSLoadNLModel(), AMPLSSolve() or
  /// AMPLSReportResults() failed, that error's code, or 500 (failure)
  /// if it has none; otherwise the result of the last solve
  /// (e.g., 0 solved, 200 infeasible, 400 limit), or -200 before
  /// any solve.
  AMPLS_C_EXPORT int AMPLSGetSolveResult(AMPLS_MP_Solver* slv);

  /// Add message
  AMPLS_C_EXPORT void AMPLSAddMessage(AMPLS_MP_Solver* slv, const char* msg);

  /// Retrieve messages, 0-terminated array
  AMPLS_C_EXPORT const char* const* AMPLSGetMessages(AMPLS_MP_Solver* slv);


  typedef struct AMPLSCOption_T {
    const char* name;
    const char* description;
    int type; // 0=int, 1=bool, 2=double, 3=string, 4=undefined - aliases
  } AMPLS_C_Option;
  
  AMPLS_C_EXPORT AMPLS_C_Option* AMPLSGetOptions(AMPLS_MP_Solver* slv);
  /// Return 0 on success
  AMPLS_C_EXPORT int AMPLSSetIntOption(AMPLS_MP_Solver*, const char* name, int value);
  AMPLS_C_EXPORT int AMPLSGetIntOption(AMPLS_MP_Solver*, const char* name, int *value);
  AMPLS_C_EXPORT int AMPLSSetDblOption(AMPLS_MP_Solver* slv, const char* name, double value);
  AMPLS_C_EXPORT int AMPLSGetDblOption(AMPLS_MP_Solver* slv, const char* name, double* value);
  AMPLS_C_EXPORT int AMPLSSetStrOption(AMPLS_MP_Solver* slv, const char* name, const char* value);
  AMPLS_C_EXPORT int AMPLSGetStrOption(AMPLS_MP_Solver* slv, const char* name, const char** value);

#endif // AMPLS_C_INTERFACE_H

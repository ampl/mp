#ifndef AMPLS_C_INTERFACE_H
#define AMPLS_C_INTERFACE_H

#include <stddef.h>

/*
 * A generic C API for MP
 */

#include <stddef.h>    // for size_t

#ifdef _WIN32
#define AMPLS_C_EXPORT __declspec(dllexport)
#else
#define AMPLS_C_EXPORT
#endif

/// An AMPLS solver instance
typedef struct AMPLS_MP_Solver_T {
  /// AMPLS internal info
  void *internal_info_;
  /// Extra info, managed by the specific solver
  void *solver_info_;
  /// User info, free to assign
  void *user_info_;
} AMPLS_MP_Solver;


/// Load model incl. suffixes.
/// The method assumes the \a slv has been initialized by
/// a solver-specific API using functions of ampls-cpp-api.h.
/// @return 0 on success
AMPLS_C_EXPORT int AMPLSLoadNLModel(AMPLS_MP_Solver* slv,
                     const char* nl_filename);

/// No AMPLSOptimize(), user does it with an internal handler
/// such as that obtained by GetGRBmodel()
/// int AMPLSOptimize(AMPLS_MP_Solver* slv);

/// Report results.
/// The kind of results reported is influenced by solver option
/// `wantsol`.
/// @return 0 on success
AMPLS_C_EXPORT int AMPLSReportResults(AMPLS_MP_Solver* slv);

/// Add message
AMPLS_C_EXPORT void AMPLSAddMessage(AMPLS_MP_Solver* slv, const char* msg);

/// Retrieve messages, 0-terminated array
const char* const * AMPLSGetMessages(AMPLS_MP_Solver* slv);


/// Set of callbacks provided to a driver for licensing issues
typedef struct CCallbacks_T {
  /// If provided, replacement function to create solver environment
  /// for a given solver.
  void* (*init)(); 
  /// Check if model sizes are allowed; throw if not.
  /// Parameters: n_vars, n_alg_con, n_log_con
  void (*check)(size_t, size_t, size_t);
  /// Return additional text to be displayed on '-v' output
  const char* (*additionalText)();
  /// Function called after failures to provide additional diagnostics
  void (*diagnostics)();
} CCallbacks;

#endif // AMPLS_C_INTERFACE_H

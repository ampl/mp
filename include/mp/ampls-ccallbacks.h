#ifndef AMPLSCCALLCACKS_H
#define AMPLSCCALLCACKS_H

#ifdef __cplusplus
extern "C" {
#endif

/// Model traits for license check
typedef struct AMPLS_ModelTraits_T {
	long long
	n_vars,
	n_alg_con,
	n_log_con,
	n_quad_con,
	n_conic_con;
} AMPLS_ModelTraits;


typedef void (*Checker_AMPLS_ModeltTraits)(const AMPLS_ModelTraits*);


/// Set of callbacks provided to a driver for licensing issues
typedef struct CCallbacks_T {
  /// If provided, replacement function to create solver environment
  /// for a given solver.
  void* (*init)();
  /// Check if model sizes are allowed; throw if not.
  /// Parameters: n_vars, n_alg_con, n_log_con
	Checker_AMPLS_ModeltTraits check;
  /// Return additional text to be displayed on '-v' output
  const char* (*additionalText)();
  /// Function called after failures to provide additional diagnostics
  void (*diagnostics)();
} CCallbacks;

#ifdef __cplusplus
}
#endif

#endif

#ifndef AMPLSCCALLCACKS_H
#define AMPLSCCALLCACKS_H

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

#endif
/**
 * C API: extern "C" wrappers for
 * NLUtils etc.
 *
 */

#ifndef NLWRITER2MISCC_H
#define NLWRITER2MISCC_H


#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif

/// Wrap mp::NLUtils for C API.
///
/// NL writer and SOL reader utilities.
/// It provides facilities for logging
/// and error handling.
///
/// To fill **default methods**,
/// call NLW2_MakeNLUtils_C_Default() / NLW2_Destroy...().
///
/// The default error handler exit()s.
typedef struct NLUtils_C {

  /// Use this pointer if you need to store your data.
  /// It is provided as the 1st argument to the methods.
  void* p_user_data_;

  /// log message
  void (*log_message)(
      void* p_user_data, const char* format, ...);
  /// log warning
  void (*log_warning)(
      void* p_user_data, const char* format, ...);
  /// Override this to your error handler.
  /// Not using exceptions by default.
  /// Only called on a wrong output format string
  /// (internal error.)
  void (*myexit)(void* p_user_data, const char* msg);

} NLUtils_C;

/// Create a default NLUtils_C wrapper object.
/// User application might change some methods
/// and use the p_user_data_ pointer.
NLUtils_C NLW2_MakeNLUtils_C_Default();

/// Destroy the NLUtils_C object created by the API
void NLW2_DestroyNLUtils_C_Default(NLUtils_C*);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLWRITER2MISCC_H

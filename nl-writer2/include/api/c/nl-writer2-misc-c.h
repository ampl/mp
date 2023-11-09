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

/// Wrap mp::NLUtils.
///
/// NL writer and SOL reader utilities.
/// It provides default facilities for logging
/// and error handling.
/// The default error handler exit()s.
typedef struct NLUtils_C {
  void* p_;

} NLUtils_C;

/// Create a default NLUtils_C wrapper object.
/// User application might want to customize
/// using the below functions.
NLUtils_C NLW2_MakeNLUtils_C_Default();


#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLWRITER2MISCC_H

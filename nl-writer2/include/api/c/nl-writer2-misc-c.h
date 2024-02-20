/*
 C API: extern "C" wrappers for
 mp::NLUtils etc.

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
typedef struct NLW2_NLUtils_C {

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

} NLW2_NLUtils_C;

/// Create a default NLW2_NLUtils_C wrapper object.
/// User application might change some methods
/// and use the p_user_data_ pointer.
NLW2_NLUtils_C NLW2_MakeNLUtils_C_Default(void);

/// Destroy the NLW2_NLUtils_C object created by the API
void NLW2_DestroyNLUtils_C_Default(NLW2_NLUtils_C*);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLWRITER2MISCC_H

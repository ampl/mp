/**
 * C API: extern "C" wrappers for the NLFeeder2 interface,
 * as well as NLWriter2 calls.
 *
 */

#ifndef NLFEEDER2_C_H
#define NLFEEDER2_C_H


#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif

/// Wrap mp::NLFeeder2
typedef struct NLFeeder2_C {
  void* p_;

} NLFeeder2_C;

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLFEEDER2_C_H

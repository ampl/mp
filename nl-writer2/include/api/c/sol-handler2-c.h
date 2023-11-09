/**
 * C API: extern "C" wrappers for the SOLHandler2 interface,
 * as well as SOLReader2 calls.
 *
 */

#ifndef SOLHANDLER2C_H
#define SOLHANDLER2C_H


#ifdef __cplusplus  // Can be used from C++
extern "C" {
#endif

/// Wrap mp::SOLHandler2
typedef struct SOLHandler2_C {
  void* p_;

} SOLHandler2_C;

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // SOLHANDLER2C_H

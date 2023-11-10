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

/** Wrap mp::NLFeeder2 for C API.

  NLFeeder2_C: writes model details on request
  via provided callback objects.
  See the examples folder.

  For the NL format, variables and constraints must have certain order.

  **Variable ordering:**
    first continuous, then integer.
  Some solvers might require more elaborate ordering, see NLHeader_C.

  **Constraint ordering:**
    first algebraic (including complementarity), then logical.
  Some solvers might require nonlinear constraints first.
 */
typedef struct NLFeeder2_C {
  /// User data
  void* p_user_data_;

} NLFeeder2_C;

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLFEEDER2_C_H

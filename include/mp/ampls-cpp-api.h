#ifndef AMPLSCPPAPI_H
#define AMPLSCPPAPI_H

/*
 * --------------------- Internal AMPLS API -----------------------
 * --------- Only to be used by solver-specific AMPLS API ---------
 */

#include <memory>

#include "mp/backend-base.h"

extern "C" {
  #include "mp/ampls-c-api.h"
}

/// AMPLS internal initialize.
/// To be used by Init<Gurobi>(AMPLS_MP_Solver*, char*) etc.
/// This function is not for the user.
///
/// @param p_be: GurobiBackend etc.
/// @return slv, also on failure: the error message is then among
///   AMPLSGetMessages() (callers such as ampls-api read it there).
///   The instance should then not be used further, only closed.
AMPLS_MP_Solver* AMPLS__internal__Open(std::unique_ptr<mp::BasicBackend> p_be,
  CCallbacks cb);

/// Shut down solver instance. Internal only.
/// Deallocate the mp-related info and the structure itself.
/// Before calling this function, a solver-specific implementation should
/// delete solver_info_ and user_info_ if allocated.
void AMPLS__internal__Close(AMPLS_MP_Solver* slv);

/// Retrieve the Backend
mp::BasicBackend* AMPLSGetBackend(AMPLS_MP_Solver* slv);

/// Record the outcome of the last AMPLS call for AMPLSGetSolveResult():
/// 0 on success, otherwise the solve result code to report (> 0)
void AMPLS__internal__SetCallResult(AMPLS_MP_Solver* slv, int code);

/// Wrapper for try/catch.
/// Can be used from C++ functions returning \a int for C API.
/// With \a record, the outcome is kept for AMPLSGetSolveResult().
/// @return 0 if success; on failure, the solve result code of an
///   mp::Error if it has a nonzero one (e.g., 200 when the converter
///   proved the model infeasible), so that the caller can report it
///   as the solve result; otherwise -1 (an error without a code,
///   another exception, or code 0, which would mean success)
template <class Fn> inline
int AMPLS__internal__TryCatchWrapper(
    AMPLS_MP_Solver* slv, Fn fn, bool record=true) {
  int result = 0;
  try {
    fn();
  } catch (const mp::Error& err) {
    AMPLSAddMessage(slv, err.what());
    const int code = err.exit_code();
    result = code > 0 ? code : -1;
  } catch (const std::exception& exc) {
    AMPLSAddMessage(slv, exc.what());
    result = -1;
  } catch (...) {
    AMPLSAddMessage(slv, "Unknown exception");
    result = -1;
  }
  if (record)
    AMPLS__internal__SetCallResult(slv, result < 0 ? 500 : result);
  return result;
}

#endif // AMPLSCPPAPI_H

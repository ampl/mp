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
/// @param slv: pointer to struct AMPLS_MP_Solver to be populated.
/// @param p_be: GurobiBackend etc.
/// @param slv_opt: a string of solver options
/// (normally provided in the <solver>_options variable).
/// Can be \a nullptr.
/// @return 0 on success, otherwise use AMPLSGetMessages()
/// (well they can contain warnings in any case).
int AMPLS__internal__Open(AMPLS_MP_Solver* slv,
                          std::unique_ptr<mp::BasicBackend> p_be,
                          const char* slv_opt);

/// Shut down solver instance. Internal only.
void AMPLS__internal__Close(AMPLS_MP_Solver* slv);

/// Retrieve the Backend
mp::BasicBackend* AMPLSGetBackend(AMPLS_MP_Solver* slv);

/// Wrapper for try/catch.
/// Can be used from C++ functions returning \a int for C API.
/// @return 0 if success
template <class Fn> inline
int AMPLS__internal__TryCatchWrapper(AMPLS_MP_Solver* slv, Fn fn) {
  try {
    fn();
  } catch (const std::exception& exc) {
    AMPLSAddMessage(slv, exc.what());
    return -1;
  } catch (...) {
    AMPLSAddMessage(slv, "Unknown exception");
    return -1;
  }
  return 0;
}

#endif // AMPLSCPPAPI_H

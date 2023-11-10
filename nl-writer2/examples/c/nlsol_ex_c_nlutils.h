/**
 * Custom NLUtils_C for the C API example
 *
 */
#ifndef NLSOL_EX_C_NLUTILS_H
#define NLSOL_EX_C_NLUTILS_H

#include "api/c/nl-writer2-misc-c.h"

/// Fill NLUtils_C for the C API example
NLUtils_C MakeNLUtils_C();

/// Destroy custom NLUtils_C
void DestroyNLUtils_C(NLUtils_C* );

#endif // NLSOL_EX_C_NLUTILS_H

/**
 * Custom NLFeeder2_C for the C API example
 *
 */

#ifndef NLSOL_EX_C_NL_H
#define NLSOL_EX_C_NL_H

#include "api/c/nl-feeder2-c.h"

#include "nlsol_ex_c_model.h"

/// Fill an NLFeeder2_C for the C API example
NLFeeder2_C MakeNLFeeder2_C(CAPIExample* , int binary);

/// Destroy custom NLFeeder2_C
void DestroyNLFeeder2_C(NLFeeder2_C* );

#endif // NLSOL_EX_C_NL_H

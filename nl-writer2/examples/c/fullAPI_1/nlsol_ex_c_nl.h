/**
 * Custom NLFeeder_C for the C API example
 *
 */

#ifndef NLSOL_EX_C_NL_H
#define NLSOL_EX_C_NL_H

#include "api/c/nl-feeder-c.h"

#include "nlsol_ex_c_model.h"

/// Fill an NLFeeder_C for the C API example
NLW2_NLFeeder_C MakeNLFeeder_C(CAPIExample* , int binary);

/// Destroy custom NLFeeder_C
void DestroyNLFeeder_C(NLW2_NLFeeder_C* );

#endif // NLSOL_EX_C_NL_H

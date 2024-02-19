/**
 * Custom SOLHandler_C for the C API example
 *
 */
#ifndef NLSOL_EX_C_SOL_H
#define NLSOL_EX_C_SOL_H

#include "api/c/sol-handler-c.h"

#include "nlsol_ex_c_model.h"

/// Fill a SOLHandler_C for the C API example
NLW2_SOLHandler_C MakeSOLHandler_C(CAPIExample* );

/// Destroy custom SOLHandler_C
void DestroySOLHandler_C(NLW2_SOLHandler_C* );

#endif // NLSOL_EX_C_SOL_H

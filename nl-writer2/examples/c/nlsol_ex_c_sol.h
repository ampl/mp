/**
 * Custom SOLHandler2_C for the C API example
 *
 */
#ifndef NLSOL_EX_C_SOL_H
#define NLSOL_EX_C_SOL_H

#include "api/c/sol-handler2-c.h"

#include "nlsol_ex_c_model.h"

/// Fill a SOLHandler2_C for the C API example
NLW2_SOLHandler2_C MakeSOLHandler2_C(CAPIExample* );

/// Destroy custom SOLHandler2_C
void DestroySOLHandler2_C(NLW2_SOLHandler2_C* );

#endif // NLSOL_EX_C_SOL_H

#include <stdlib.h>

#include "nlsol_ex_c_sol.h"

/// Declare, implementation in nlsol_ex_c_nl.c
NLHeader_C CAPI_ex_Header(void* pex_void);

void OnDualSolution(
    void* p_user_data, int nvals, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  for (int i=0; i<nvals; ++i)
    pex->sol_dual_[i] = NLW2_ReadSolVal(p_api_data);
}

void OnPrimalSolution(
    void* p_user_data, int nvals, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  for (int i=0; i<nvals; ++i)
    pex->sol_primal_[i] = NLW2_ReadSolVal(p_api_data);
}

void OnObjno(void* p_user_data, int on) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  pex->objno_  = on;
}

void OnSolveCode(void* p_user_data, int sc) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  pex->solve_code_ = sc;
}

NLW2_SOLHandler2_C MakeSOLHandler2_C(CAPIExample* pex) {
  NLW2_SOLHandler2_C result
      = NLW2_MakeSOLHandler2_C_Default();

  result.p_user_data_ = pex;

  result.Header = CAPI_ex_Header;
  // solve message: use default (print it)
  // AMPL options: use default
  result.OnDualSolution = OnDualSolution;
  result.OnPrimalSolution = OnPrimalSolution;
  result.OnObjno = OnObjno;
  result.OnSolveCode = OnSolveCode;

  return result;
}

void DestroySOLHandler2_C(NLW2_SOLHandler2_C* p) {
  p->p_user_data_ = NULL;
  NLW2_DestroySOLHandler2_C_Default(p);
}

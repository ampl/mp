#include <stdlib.h>

#include "nlsol_ex_c_sol.h"

SOLHandler2_C MakeSOLHandler2_C(CAPIExample* pex) {
  SOLHandler2_C result;

  result.p_user_data_ = pex;

  return result;
}

void DestroySOLHandler2_C(SOLHandler2_C* p) {
  p->p_user_data_ = NULL;
}

#include <stdlib.h>

#include "nlsol_ex_c_nl.h"

NLFeeder2_C MakeNLFeeder2_C(
    CAPIExample* pex, int binary) {
  NLFeeder2_C result;

  result.p_user_data_ = pex;
  pex->binary_nl = binary;

  return result;
}

void DestroyNLFeeder2_C(NLFeeder2_C* pf) {
  pf->p_user_data_ = NULL;
}

#include "nlsol_ex_c_nlutils.h"

NLUtils_C MakeNLUtils_C() {
  // Just return the API's default config:
  return NLW2_MakeNLUtils_C_Default();
}

void DestroyNLUtils_C(NLUtils_C* p) {
  NLW2_DestroyNLUtils_C_Default(p);
}

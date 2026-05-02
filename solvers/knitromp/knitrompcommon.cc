#include "mp/format.h"
#include "knitrompcommon.h"

namespace mp {


double KnitrompCommon::getDblAttr(const char* name) const  {
  double value = 0;
  /* TODO Utility function to get the value of an integer attribute
 * from the solver API
  KNITROMP_CCALL(KNITROMP_GetDblAttr(lp_, name, &value)); */
  return value;
}

int KnitrompCommon::NumCons() const {
    int num_cons;
	KN_get_number_cons(lp(), &num_cons);
    return num_cons;
}
int KnitrompCommon::NumVars() const {
    int num_vars;
    KN_get_number_vars(lp(), &num_vars);
    return num_vars;
}
int KnitrompCommon::NumObjs() const {
    return 1;
}

void KnitrompCommon::GetSolverOption(int key, int &value) const {
  KNITROMP_CCALL( KN_get_int_param(lp(), key, &value));
}

void KnitrompCommon::SetSolverOption(int key, int value) {
  KNITROMP_CCALL(KN_set_int_param(lp(), key, value));
}

void KnitrompCommon::GetSolverOption(int key, double &value) const {
    KNITROMP_CCALL(KN_get_double_param(lp(), key, &value));
}

void KnitrompCommon::SetSolverOption(int key, double value) {
    KNITROMP_CCALL(KN_set_double_param(lp(), key, value));
}

void KnitrompCommon::GetSolverOption(int key, std::string &value) const {
}

void KnitrompCommon::SetSolverOption(int key, const std::string& value) {
    KNITROMP_CCALL(KN_set_char_param(lp(), key, value.data()));
}


} // namespace mp

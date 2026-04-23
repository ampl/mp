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

int KnitrompCommon::NumLinCons() const {
    return 0;
  //return getIntAttr(Solver::NCONS_TYPE, Solver::ConsType::CONS_LIN);
  // TODO Get number of linear constraints using solver API
  // return getIntAttr(KNITROMP_INTATTR_ROWS);
}

int KnitrompCommon::NumVars() const {
    int num_vars;
    KN_get_number_vars(lp(), &num_vars);
    return num_vars;
}
int KnitrompCommon::NumObjs() const {
    return 1;
}

int KnitrompCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(KNITROMP_INTATTR_QCONSTRS);
    return 0;// getIntAttr(Solver::NCONS_TYPE, Solver::ConsType::CONS_QUAD);
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

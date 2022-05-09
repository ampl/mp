#include "mp/format.h"
#include "highscommon.h"

namespace mp {

void HighsCommon::OpenSolver() {
  int status = 0;
  void* prob = Highs_create();
  set_lp(prob); // Assign it
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
  // HIGHS_CCALL(HIGHS_SetIntParam(prob, "Logging", 0));

}

void HighsCommon::CloseSolver() {
  Highs_destroy(lp_);
}

void HighsCommon::copy_handlers_from_other_highs() {
  assert(other_highs());
  lp_ = other_highs()->lp(); 
}

void HighsCommon::copy_handlers_to_other_highs() {
  assert(other_highs());
  other_highs()->set_lp(lp_); 
}
int64_t HighsCommon::getInt64Attr(const char* name)  const {
  int64_t value = 0;
  HIGHS_CCALL(Highs_getInt64InfoValue(lp_, name, &value));
  return value;
}
int HighsCommon::getIntAttr(const char* name)  const {
  int value = 0;
  HIGHS_CCALL(Highs_getIntInfoValue(lp_, name, &value)); 
  return value;
}
double HighsCommon::getDblAttr(const char* name) const  {
  double value = 0;
  HIGHS_CCALL(Highs_getDoubleInfoValue(lp_, name, &value));
  return value;
}

int HighsCommon::NumLinCons() const {
  return Highs_getNumRows(lp());
}

int HighsCommon::NumVars() const {
  return Highs_getNumCols(lp());
}

int HighsCommon::NumObjs() const {
  // TODO Get number of objectives using solver API

  return 1;
}

int HighsCommon::NumQPCons() const {
  return 0;
}

int HighsCommon::NumSOSCons() const {
  return 0;
}

int HighsCommon::NumIndicatorCons() const {
  return 0;
}

void HighsCommon::GetSolverOption(const char* key, int &value) const {
  HIGHS_CCALL(Highs_getIntOptionValue(lp_, key, &value) );
}

void HighsCommon::SetSolverOption(const char* key, int value) {
  HIGHS_CCALL(Highs_setIntOptionValue(lp_, key, value));
}

void HighsCommon::GetSolverOption(const char* key, double &value) const {
  HIGHS_CCALL(Highs_getDoubleOptionValue(lp_, key, &value) );
}

void HighsCommon::SetSolverOption(const char* key, double value) {
  HIGHS_CCALL(Highs_setDoubleOptionValue(lp_, key, value));
}

void HighsCommon::GetSolverOption(const char* key, std::string &value) const {
  char option[256];
  HIGHS_CCALL(Highs_getStringOptionValue(lp_, key, option));
  value = option;
}

void HighsCommon::SetSolverOption(const char* key, const std::string& value) {
  HIGHS_CCALL(Highs_setStringOptionValue(lp_, key, value.data()));
}


} // namespace mp

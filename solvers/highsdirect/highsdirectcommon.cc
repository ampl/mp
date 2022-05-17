#include "mp/format.h"
#include "highsdirectcommon.h"

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
  Highs_destroy(lp());
}

int64_t HighsCommon::getInt64Attr(const char* name)  const {
  int64_t value = 0;
  HIGHS_CCALL(Highs_getInt64InfoValue(lp(), name, &value));
  return value;
}
int HighsCommon::getIntAttr(const char* name)  const {
  int value = 0;
  HIGHS_CCALL(Highs_getIntInfoValue(lp(), name, &value));
  return value;
}
double HighsCommon::getDblAttr(const char* name) const  {
  double value = 0;
  HIGHS_CCALL(Highs_getDoubleInfoValue(lp(), name, &value));
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


void checkOption(int retvalue, const char* key) {
  if (retvalue != HIGHS_RETCODE_OK)
    throw std::runtime_error(fmt::format("while setting option '{}'", key));
  
}
void HighsCommon::GetSolverOption(const char* key, int& value) const {
  int type;
  Highs_getOptionType(lp(), key, &type);
  if (type == kHighsOptionTypeBool)
    HIGHS_CCALL(Highs_getBoolOptionValue(lp(), key, &value));
  else
    HIGHS_CCALL(Highs_getIntOptionValue(lp(), key, &value));
}

void HighsCommon::SetSolverOption(const char* key, int value) {
  int type;
  int ret;
  Highs_getOptionType(lp(), key, &type);
  if (type == kHighsOptionTypeBool)
    ret = Highs_setBoolOptionValue(lp(), key, value);
  else
    ret =Highs_setIntOptionValue(lp(), key, value);
  checkOption(ret, key);
}

void HighsCommon::GetSolverOption(const char* key, double &value) const {
  HIGHS_CCALL(Highs_getDoubleOptionValue(lp(), key, &value) );
}

void HighsCommon::SetSolverOption(const char* key, double value) {
  int ret = Highs_setDoubleOptionValue(lp(), key, value);
  checkOption(ret, key);
}

void HighsCommon::GetSolverOption(const char* key, std::string &value) const {
  char option[256];
  HIGHS_CCALL(Highs_getStringOptionValue(lp(), key, option));
  value = option;
}

void HighsCommon::SetSolverOption(const char* key, const std::string& value) {
  int ret = Highs_setStringOptionValue(lp(), key, value.data());
  checkOption(ret, key);
}


} // namespace mp

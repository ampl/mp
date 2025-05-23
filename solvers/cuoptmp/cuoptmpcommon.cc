#include "mp/format.h"
#include "cuoptmpcommon.h"

namespace mp {


int CuoptmpCommon::NumLinCons() const {
  ProblemData* problem = lp();
  return problem->num_constraints;
}

int CuoptmpCommon::NumVars() const {
  ProblemData* problem = lp();
  return problem->num_variables;
}

int CuoptmpCommon::NumObjs() const {
  return 1;
}

int CuoptmpCommon::NumQPCons() const {
  return 0;
}

void CuoptmpCommon::GetSolverOption(const char* key, int &value) const {
  int value_out;
  CUOPTLP_CCALL(cuOptGetIntegerParameter(lp_->settings, key, &value_out));
  value = value_out;
}

void CuoptmpCommon::SetSolverOption(const char* key, int value) {
  CUOPTLP_CCALL(cuOptSetIntegerParameter(lp_->settings, key, value));
}

void CuoptmpCommon::GetSolverOption(const char* key, double &value) const {
  double value_out;
  CUOPTLP_CCALL(cuOptGetFloatParameter(lp_->settings, key, &value_out));
  value = value_out;
}

void CuoptmpCommon::SetSolverOption(const char* key, double value) {
  CUOPTLP_CCALL(cuOptSetFloatParameter(lp_->settings, key, value) );
}

void CuoptmpCommon::GetSolverOption(const char* key, std::string &value) const {
  constexpr size_t BUFFER_SIZE = 1024;
  char buffer[BUFFER_SIZE];
  CUOPTLP_CCALL(cuOptGetParameter(lp_->settings, key, BUFFER_SIZE, buffer));
  value = std::string(buffer);
}

void CuoptmpCommon::SetSolverOption(const char* key, const std::string& value) {
  CUOPTLP_CCALL(cuOptSetParameter(lp_->settings, key, value.c_str()));
}


} // namespace mp

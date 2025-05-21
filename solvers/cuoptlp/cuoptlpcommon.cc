#include "mp/format.h"
#include "cuoptlpcommon.h"

namespace Solver {
  SolverModel* CreateSolverModel() {
    return new SolverModel();
  }
}
namespace mp {


int CuoptlpCommon::NumLinCons() const {
  ProblemData* problem = lp();
  return problem->num_constraints;
}

int CuoptlpCommon::NumVars() const {
  ProblemData* problem = lp();
  return problem->num_variables;
}

int CuoptlpCommon::NumObjs() const {
  return 1;
}

int CuoptlpCommon::NumQPCons() const {
  return 0;
}

void CuoptlpCommon::GetSolverOption(const char* key, int &value) const {
  printf("GetSolverOption Int %s", key);
  int value_out;
  CUOPTLP_CCALL(cuOptGetIntegerParameter(lp_->settings, key, &value_out));
  printf(" %d\n", value_out);
  value = value_out;
}

void CuoptlpCommon::SetSolverOption(const char* key, int value) {
  printf("SetSolverOption Int %s, %d\n", key, value);
  CUOPTLP_CCALL(cuOptSetIntegerParameter(lp_->settings, key, value));
}

void CuoptlpCommon::GetSolverOption(const char* key, double &value) const {
  printf("GetSolverOption double %s, %f\n", key, value);
  double value_out;
  CUOPTLP_CCALL(cuOptGetFloatParameter(lp_->settings, key, &value_out));
  printf(" %f\n", value_out);
  value = value_out;
}

void CuoptlpCommon::SetSolverOption(const char* key, double value) {
  printf("SetSolverOption double %s, %f\n", key, value);
  CUOPTLP_CCALL(cuOptSetFloatParameter(lp_->settings, key, value) );
}

void CuoptlpCommon::GetSolverOption(const char* key, std::string &value) const {
  printf("GetSolverOption string %s\n", key);
#define BUFFER_SIZE 1024
  char buffer[BUFFER_SIZE];
  CUOPTLP_CCALL(cuOptGetParameter(lp_->settings, key, BUFFER_SIZE, buffer));
  printf(" %s\n", buffer);
  value = std::string(buffer);
}

void CuoptlpCommon::SetSolverOption(const char* key, const std::string& value) {
  printf("SetSolverOption string %s, %s\n", key, value.c_str());
  CUOPTLP_CCALL(cuOptSetParameter(lp_->settings, key, value.c_str()));
}


} // namespace mp

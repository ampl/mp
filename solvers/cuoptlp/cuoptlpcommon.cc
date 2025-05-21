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
  //CUOPTLP_CCALL( CUOPTLP_GetIntParam(lp_, key, &value) );
}

void CuoptlpCommon::SetSolverOption(const char* key, int value) {
  printf("\n SetSolverOption Int \n");
  printf("\n %s, %d \n", key, value);
  CUOPTLP_CCALL(cuOptSetIntegerParameter(lp_->settings, key, value));
}

void CuoptlpCommon::GetSolverOption(const char* key, double &value) const {
  //CUOPTLP_CCALL(CUOPTLP_GetDblParam(lp_, key, &value) );
}

void CuoptlpCommon::SetSolverOption(const char* key, double value) {
  printf("\n SetSolverOption double \n");
  printf("\n %s, %f \n", key, value);
  CUOPTLP_CCALL(cuOptSetFloatParameter(lp_->settings, key, value) );
}

void CuoptlpCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void CuoptlpCommon::SetSolverOption(const char* key, const std::string& value) {
  throw std::runtime_error("Not implemented"); // TODO
}


} // namespace mp

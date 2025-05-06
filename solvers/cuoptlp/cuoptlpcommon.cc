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
  //CUOPTLP_CCALL(CUOPTLP_SetIntParam(lp_, key, value));
}

void CuoptlpCommon::GetSolverOption(const char* key, double &value) const {
  //CUOPTLP_CCALL(CUOPTLP_GetDblParam(lp_, key, &value) );
}

void CuoptlpCommon::SetSolverOption(const char* key, double value) {
 // CUOPTLP_CCALL(CUOPTLP_SetDblParam(lp_, key, value) );
}

void CuoptlpCommon::GetSolverOption(const char* key, std::string &value) const {
  throw std::runtime_error("Not implemented"); // TODO
}

void CuoptlpCommon::SetSolverOption(const char* key, const std::string& value) {
  throw std::runtime_error("Not implemented"); // TODO
}


} // namespace mp

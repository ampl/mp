#include "mp/format.h"
#include "mp/error.h"
#include "cplexcommon.h"

namespace mp {


int CplexCommon::NumLinCons() const {
  return CPXgetnumrows (env(), lp());
}

int CplexCommon::NumVars() const {
  return CPXgetnumcols (env(), lp());
}

int CplexCommon::NumObjs() const {
  return CPXgetnumobjs (env(), lp());
}

int CplexCommon::NumQPCons() const {
  return CPXgetnumqconstrs(env(), lp());
}

int CplexCommon::NumIndicatorCons() const {
  return CPXgetnumindconstrs(env(), lp());
}

int CplexCommon::NumSOSCons() const {
  return CPXgetnumsos(env(), lp());
}
bool CplexCommon::HasQObj() const {
  return CPXgetnumquad(env(), lp()) > 0;
}
int CplexCommon::ModelSense() const {
  return CPXgetobjsen(env(), lp());
}

void CplexCommon::GetSolverOption(int key, int &value) const {
  value = 0;
  CPLEX_CALL( CPXgetintparam(env(), key, &value) );
}

void CplexCommon::SetSolverOption(int key, int value) {
  if (int e = CPXsetintparam(env(), key, value))
    RaiseRejectedParam(key, e, fmt::format("{}", value));
}

void CplexCommon::GetSolverOption(int key, double &value) const {
  value = 0.0;
  CPLEX_CALL( CPXgetdblparam(env(), key, &value) );
}

void CplexCommon::SetSolverOption(int key, double value) {
  if (int e = CPXsetdblparam(env(), key, value))
    RaiseRejectedParam(key, e, fmt::format("{}", value));
}

void CplexCommon::GetSolverOption(int key, std::string &value) const {
  char buffer[CPX_STR_PARAM_MAX] = "";
  CPLEX_CALL( CPXgetstrparam(env(), key, buffer) );
  value = buffer;
}

void CplexCommon::SetSolverOption(int key, const std::string& value) {
  if (int e = CPXsetstrparam(env(), key, value.c_str()))
    RaiseRejectedParam(key, e, value);
}

std::string CplexCommon::ParamName(int key) const {
  char name[CPX_STR_PARAM_MAX] = "";
  if (CPXgetparamhiername(env(), key, name) && CPXgetparamname(env(), key, name))
    return fmt::format("parameter {}", key);
  return name;
}

/// Raise an error for a rejected parameter value, in the style of
/// the ASL driver: "Rejecting <name>=<value>; must be between L and U".
void CplexCommon::RaiseRejectedParam(
    int key, int status, const std::string& value) const {
  std::string reason;
  int type = CPX_PARAMTYPE_NONE;
  CPXgetparamtype(env(), key, &type);
  if (CPXERR_PARAM_TOO_SMALL == status || CPXERR_PARAM_TOO_BIG == status) {
    switch (type) {
    case CPX_PARAMTYPE_INT: {
      CPXINT d = 0, lo = 0, hi = 0;
      if (!CPXinfointparam(env(), key, &d, &lo, &hi))
        reason = fmt::format("; must be between {} and {}", lo, hi);
      break;
    }
    case CPX_PARAMTYPE_LONG: {
      CPXLONG d = 0, lo = 0, hi = 0;
      if (!CPXinfolongparam(env(), key, &d, &lo, &hi))
        reason = fmt::format("; must be between {} and {}",
                             (long long)lo, (long long)hi);
      break;
    }
    case CPX_PARAMTYPE_DOUBLE: {
      double d = 0, lo = 0, hi = 0;
      if (!CPXinfodblparam(env(), key, &d, &lo, &hi))
        reason = fmt::format("; must be between {} and {}", lo, hi);
      break;
    }
    default: break;
    }
  }
  if (reason.empty()) {
    char buffer[CPXMESSAGEBUFSIZE] = "";
    const char* msg = CPXgeterrorstring(env(), status, buffer);
    std::string m = msg ? msg : "";
    while (!m.empty() && (m.back() == '\n' || m.back() == ' '))
      m.pop_back();
    reason = fmt::format(": {}", m.empty() ? fmt::format("CPLEX error {}", status) : m);
  }
  MP_RAISE(fmt::format("Rejecting {}={}{}", ParamName(key), value, reason));
}

double CplexCommon::GetCPLEXDblParam(int param) {
  double value;
  CPLEX_CALL(CPXgetdblparam(env(), param, &value));
  return value;
}
int CplexCommon::GetCPLEXIntParam(int param) {
  int value;
  CPLEX_CALL(CPXgetintparam(env(), param, &value));
  return value;
}
void CplexCommon::SetCPLEXParam(int param, int value) {
  CPLEX_CALL(CPXsetintparam(env(), param, value));
}
void CplexCommon::SetCPLEXParam(int param, double value) {
  CPLEX_CALL(CPXsetdblparam(env(), param, value));
}
} // namespace mp

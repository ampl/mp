#include <climits>

#include "mp/error.h"
#include "mp/format.h"
#include "scipmpcommon.h"

static
SCIP_DECL_PROBDELORIG(probdataDelOrigNl)
{
  int i;

  assert((*probdata)->vars != NULL || (*probdata)->nvars == 0);
  assert((*probdata)->linconss.size() || (*probdata)->nlinconss == 0);

  for( i = 0; i < (*probdata)->nlinconss; ++i )
  {
    SCIP_CCALL( SCIPreleaseCons(scip, &(*probdata)->linconss[i]) );
  }

  for (auto& pnlc: (*probdata)->nlconss)
    SCIP_CCALL( SCIPreleaseCons(scip, &pnlc) );

  for( i = 0; i < (*probdata)->nvars; ++i )
  {
    if ((*probdata)->var_exprs[i])
      SCIP_CCALL( SCIPreleaseExpr(scip, &(*probdata)->var_exprs[i]) );
    SCIP_CCALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
  }

  for (auto& e: (*probdata)->exprs)
    SCIP_CCALL( SCIPreleaseExpr(scip, &e) );

  if ((*probdata)->dummyexpr)
    SCIP_CCALL( SCIPreleaseExpr(scip, &(*probdata)->dummyexpr) );

  SCIPfreeBlockMemoryArray(scip, &(*probdata)->vars, (*probdata)->nvars);

  SCIPsetMessagehdlrQuiet(scip, true);    // Spurious warnings and failure in Debug build
                                          // https://github.com/scipopt/scip/issues/71
  SCIPfreeMemory(scip, probdata);

  return SCIP_OKAY;
}

namespace mp {

void ScipCommon::OpenSolver() {
  int status = 0;
  SCIP* scip = NULL;
  SCIP_PROBDATA* probdata = NULL;

  // initialize SCIP
  status = SCIPcreate(&scip);
  setSCIP(scip); // Assign it

  // include default SCIP plugins
  SCIP_CCALL( SCIPincludeDefaultPlugins(scip) );

  // initialize empty SCIP problem
  SCIP_CCALL( SCIPallocClearMemory(scip, &probdata) );
  SCIP_CCALL( SCIPcreateProb(scip, "", probdataDelOrigNl, NULL, NULL, NULL, NULL, NULL, probdata) );
  setPROBDATA(probdata);

  if (status != 1)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
  SetSolverOption("display/verblevel", 0);
}

void ScipCommon::CloseSolver() {
  SCIP* scip = getSCIP();

  // free SCIP
  SCIP_CCALL( SCIPfree(&scip) );
}

int ScipCommon::NumLinCons() const {
  return getPROBDATA()->nlinconss;
}

int ScipCommon::NumVars() const {
  return getPROBDATA()->nvars;
}

int ScipCommon::NumObjs() const {
  return 1;
}

int ScipCommon::NumQPCons() const {
  int count = 0;

  for (int i = 0; i < SCIPgetNOrigConss(getSCIP()); i++) {
    SCIP_Bool isquadratic;
    SCIP_CCALL( SCIPcheckQuadraticNonlinear(getSCIP(), SCIPgetOrigConss(getSCIP())[i], &isquadratic) );
    if (isquadratic == true)
      count++;
  }

  return count;
}

int ScipCommon::NumSOSCons() const {
  int count = 0;
  
  for (int i = 0; i < SCIPgetNOrigConss(getSCIP()); i++) {
    if (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(SCIPgetOrigConss(getSCIP())[i])), "SOS1") == 0 ||
      strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(SCIPgetOrigConss(getSCIP())[i])), "SOS2") == 0)
      count++;
  }

  return count;
}

int ScipCommon::NumIndicatorCons() const {
  int count = 0;
  
  for (int i = 0; i < SCIPgetNOrigConss(getSCIP()); i++) {
    if (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(SCIPgetOrigConss(getSCIP())[i])), "indicator") == 0)
      count++;
  }

  return count;
}

SCIP_PARAM* ScipCommon::Param(const char* key) const {
  auto* param = SCIPgetParam(getSCIP(), key);
  if (!param)
    throw std::runtime_error(
        fmt::format("Unknown SCIP parameter '{}'", key));
  return param;
}

SCIP_PARAMTYPE ScipCommon::ParamType(const char* key) const {
  return SCIPparamGetType(Param(key));
}

std::string ScipCommon::AmplOptionName(const char* key) const {
  auto it = optNames_.find(key);
  return it == optNames_.end() ? std::string(key) : it->second;
}

/// Valid values of a SCIP parameter, for the rejection message
static std::string ParamRange(SCIP_PARAM* param) {
  switch (SCIPparamGetType(param)) {
  case SCIP_PARAMTYPE_BOOL:
    return "0 or 1";
  case SCIP_PARAMTYPE_INT:
    return fmt::format("in [{}, {}]",
                       SCIPparamGetIntMin(param), SCIPparamGetIntMax(param));
  case SCIP_PARAMTYPE_LONGINT:
    return fmt::format("in [{}, {}]",
                       SCIPparamGetLongintMin(param), SCIPparamGetLongintMax(param));
  case SCIP_PARAMTYPE_REAL:
    return fmt::format("in [{}, {}]",
                       SCIPparamGetRealMin(param), SCIPparamGetRealMax(param));
  case SCIP_PARAMTYPE_CHAR:
    if (auto* allowed = SCIPparamGetCharAllowedValues(param))
      return fmt::format("one of '{}'", allowed);
    return "a single character";
  default:
    return "a string";
  }
}

void ScipCommon::CheckOptionSet(SCIP_RETCODE code, const char* key,
                                const std::string& value) const {
  if (SCIP_OKAY == code)
    return;
  throw std::runtime_error(fmt::format(
      "Rejecting {}={} (SCIP parameter '{}'); the value must be {}",
      AmplOptionName(key), value, key, ParamRange(Param(key))));
}

void ScipCommon::GetSolverOption(const char* key, int &value) const {
  switch (ParamType(key)) {
  case SCIP_PARAMTYPE_BOOL: {
    SCIP_Bool buffer;
    SCIP_CCALL( SCIPgetBoolParam(getSCIP(), key, &buffer) );
    value = (int)buffer;
    break;
  }
  case SCIP_PARAMTYPE_LONGINT: {
    SCIP_Longint buffer;
    SCIP_CCALL( SCIPgetLongintParam(getSCIP(), key, &buffer) );
    MP_ASSERT_ALWAYS(buffer >= INT_MIN && buffer <= INT_MAX,
                     fmt::format("Cannot represent {} as int", buffer));
    value = (int)buffer;
    break;
  }
  default:
    SCIP_CCALL( SCIPgetIntParam(getSCIP(), key, &value) );
  }
}

void ScipCommon::SetSolverOption(const char* key, int value) {
  SCIP_RETCODE code;
  switch (ParamType(key)) {
  case SCIP_PARAMTYPE_BOOL:
    code = SCIPsetBoolParam(getSCIP(), key, value);
    break;
  case SCIP_PARAMTYPE_LONGINT:
    code = SCIPsetLongintParam(getSCIP(), key, value);
    break;
  default:
    code = SCIPsetIntParam(getSCIP(), key, value);
  }
  CheckOptionSet(code, key, std::to_string(value));
}

void ScipCommon::GetSolverOption(const char* key, SCIP_Longint &value) const {
  SCIP_CCALL( SCIPgetLongintParam(getSCIP(), key, &value) );
}

void ScipCommon::SetSolverOption(const char* key, SCIP_Longint value) {
  CheckOptionSet(SCIPsetLongintParam(getSCIP(), key, value),
                 key, std::to_string(value));
}

void ScipCommon::GetSolverOption(const char* key, double &value) const {
  SCIP_CCALL( SCIPgetRealParam(getSCIP(), key, &value) );
}

void ScipCommon::SetSolverOption(const char* key, double value) {
  CheckOptionSet(SCIPsetRealParam(getSCIP(), key, value),
                 key, fmt::format("{}", value));
}

void ScipCommon::GetSolverOption(const char* key, std::string &value) const {
  if (SCIP_PARAMTYPE_CHAR == ParamType(key)) {
    char buffer = 0;
    SCIP_CCALL( SCIPgetCharParam(getSCIP(), key, &buffer) );
    value.assign(1, buffer);
  } else {
    char* internal = nullptr;        // SCIP's own storage, not to be freed
    SCIP_CCALL( SCIPgetStringParam(getSCIP(), key, &internal) );
    value = internal ? internal : "";
  }
}

void ScipCommon::SetSolverOption(const char* key, const std::string& value) {
  SCIP_RETCODE code;
  if (SCIP_PARAMTYPE_CHAR == ParamType(key))
    code = value.size() == 1
        ? SCIPsetCharParam(getSCIP(), key, value[0])
        : SCIP_PARAMETERWRONGVAL;
  else
    code = SCIPsetStringParam(getSCIP(), key, value.c_str());
  CheckOptionSet(code, key, value);
}

/// Set a parameter of any type from its textual value
/// (tech:optionnative, parameter files)
void ScipCommon::SetNativeParam(const std::string& name,
                                const std::string& value) {
  switch (ParamType(name.c_str())) {
  case SCIP_PARAMTYPE_BOOL:
  case SCIP_PARAMTYPE_INT:
    SetSolverOption(name.c_str(), std::stoi(value));
    break;
  case SCIP_PARAMTYPE_LONGINT:
    SetSolverOption(name.c_str(), (SCIP_Longint)std::stoll(value));
    break;
  case SCIP_PARAMTYPE_REAL:
    SetSolverOption(name.c_str(), std::stod(value));
    break;
  default:
    SetSolverOption(name.c_str(), value);
  }
}


double ScipCommon::Infinity() const { 
  double inf;
  GetSolverOption("numerics/infinity", inf);
  return inf;  
}

double ScipCommon::MinusInfinity() { return -Infinity(); }

bool ScipCommon::IsContinuous() {
  // True iff scip has only continuous variables
  for (int i = 0; i < getPROBDATA()->nvars; i++) {
    if (SCIPvarGetType(getPROBDATA()->vars[i]) != SCIP_VARTYPE_CONTINUOUS)
      return false;
  }
  return true;
}

} // namespace mp

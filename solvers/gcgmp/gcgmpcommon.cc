#include <climits>

#include "mp/format.h"
#include "mp/error.h"

#include "gcgmpcommon.h"
#include "gcg/gcgplugins.h"
#include "gcg/struct_gcg.h"      // GCGgetOrigprob() is a macro in release
                                 // builds only; the library exports no function

static
SCIP_DECL_PROBDELORIG(probdataDelOrigNl)
{
  int i;

  assert((*probdata)->vars != NULL || (*probdata)->nvars == 0);
  assert((*probdata)->linconss.size() || (*probdata)->nlinconss == 0);

  for( i = 0; i < (*probdata)->nlinconss; ++i )
  {
    SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->linconss[i]) );
  }

  for( i = 0; i < (*probdata)->nvars; ++i )
  {
    GCG_CCALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
  }
  SCIPfreeBlockMemoryArray(scip, &(*probdata)->vars, (*probdata)->nvars);

  SCIPfreeMemory(scip, probdata);

  return SCIP_OKAY;
}
namespace mp {

void GcgCommon::OpenSolver() {
  int status = 0;
  GCG* gcg = NULL;
  SCIP* scip = NULL;
  SCIP_PROBDATA* probdata = NULL;

  // initialize SCIP
  status = GCGcreate(&gcg);
  setGCG(gcg);
  scip = gcg->origprob;
  setSCIP(scip); // Assign it

  // include default GCG plugins
  //GCG_CCALL( GCGincludeGcgPlugins(gcg) );

  // initialize empty SCIP problem
  GCG_CCALL( SCIPallocClearMemory(scip, &probdata) );
  GCG_CCALL( SCIPcreateProb(scip, "", probdataDelOrigNl, NULL, NULL, NULL, NULL, NULL, probdata) );
  setPROBDATA(probdata);

  if (status != 1)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
  SetSolverOption("display/verblevel", 0);
}

void GcgCommon::CloseSolver() {
  GCG* gcg = getGCG();

  // free GCG
  GCG_CCALL( GCGfree(&gcg) );
}

int GcgCommon::NumLinCons() const {
  return getPROBDATA()->nlinconss;
}

int GcgCommon::NumVars() const {
  return getPROBDATA()->nvars;
}

int GcgCommon::NumObjs() const {
  return 1;
}

int GcgCommon::NumQPCons() const {
  return 0;
}

int GcgCommon::NumSOSCons() const {
  return 0;
}

int GcgCommon::NumIndicatorCons() const {
  return 0;
}


SCIP_PARAM* GcgCommon::Param(const char* key) const {
  auto* param = SCIPgetParam(getSCIP(), key);
  if (!param)
    throw std::runtime_error(
        fmt::format("Unknown GCG/SCIP parameter '{}'", key));
  return param;
}

SCIP_PARAMTYPE GcgCommon::ParamType(const char* key) const {
  return SCIPparamGetType(Param(key));
}

std::string GcgCommon::AmplOptionName(const char* key) const {
  auto it = optNames_.find(key);
  return it == optNames_.end() ? std::string(key) : it->second;
}

/// Valid values of a parameter, for the rejection message
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

void GcgCommon::CheckOptionSet(SCIP_RETCODE code, const char* key,
                               const std::string& value) const {
  if (SCIP_OKAY == code)
    return;
  throw std::runtime_error(fmt::format(
      "Rejecting {}={} (GCG parameter '{}'); the value must be {}",
      AmplOptionName(key), value, key, ParamRange(Param(key))));
}

void GcgCommon::GetSolverOption(const char* key, int &value) const {
  switch (ParamType(key)) {
  case SCIP_PARAMTYPE_BOOL: {
    SCIP_Bool buffer;
    GCG_CCALL( SCIPgetBoolParam(getSCIP(), key, &buffer) );
    value = (int)buffer;
    break;
  }
  case SCIP_PARAMTYPE_LONGINT: {
    SCIP_Longint buffer;
    GCG_CCALL( SCIPgetLongintParam(getSCIP(), key, &buffer) );
    MP_ASSERT_ALWAYS(buffer >= INT_MIN && buffer <= INT_MAX,
                     fmt::format("Cannot represent {} as int", buffer));
    value = (int)buffer;
    break;
  }
  default:
    GCG_CCALL( SCIPgetIntParam(getSCIP(), key, &value) );
  }
}

void GcgCommon::SetSolverOption(const char* key, int value) {
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

void GcgCommon::GetSolverOption(const char* key, SCIP_Longint &value) const {
  GCG_CCALL( SCIPgetLongintParam(getSCIP(), key, &value) );
}

void GcgCommon::SetSolverOption(const char* key, SCIP_Longint value) {
  CheckOptionSet(SCIPsetLongintParam(getSCIP(), key, value),
                 key, std::to_string(value));
}

void GcgCommon::GetSolverOption(const char* key, double &value) const {
  GCG_CCALL( SCIPgetRealParam(getSCIP(), key, &value) );
}

void GcgCommon::SetSolverOption(const char* key, double value) {
  CheckOptionSet(SCIPsetRealParam(getSCIP(), key, value),
                 key, fmt::format("{}", value));
}

void GcgCommon::GetSolverOption(const char* key, std::string &value) const {
  if (SCIP_PARAMTYPE_CHAR == ParamType(key)) {
    char buffer = 0;
    GCG_CCALL( SCIPgetCharParam(getSCIP(), key, &buffer) );
    value.assign(1, buffer);
  } else {
    char* internal = nullptr;        // SCIP's own storage, not to be freed
    GCG_CCALL( SCIPgetStringParam(getSCIP(), key, &internal) );
    value = internal ? internal : "";
  }
}

void GcgCommon::SetSolverOption(const char* key, const std::string& value) {
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
/// (tech:optionnative)
void GcgCommon::SetNativeParam(const std::string& name,
                               const std::string& value) {
  switch (ParamType(name.c_str())) {
  case SCIP_PARAMTYPE_BOOL:
  case SCIP_PARAMTYPE_INT:
    SetSolverOption(name.c_str(), std::stoi(value));
    break;
  case SCIP_PARAMTYPE_LONGINT:
    SetSolverOption(name.c_str(), std::stoi(value));
    break;
  case SCIP_PARAMTYPE_REAL:
    SetSolverOption(name.c_str(), std::stod(value));
    break;
  default:
    SetSolverOption(name.c_str(), value);
  }
}


double GcgCommon::Infinity() const { 
  double inf;
  GetSolverOption("numerics/infinity", inf);
  return inf;  
}

double GcgCommon::MinusInfinity() const { 
  double inf;
  GetSolverOption("numerics/infinity", inf);
  return -inf;  
}

bool GcgCommon::IsContinuous() {
  // True iff gcg has only continuous variables
  for (int i = 0; i < getPROBDATA()->nvars; i++) {
    if (SCIPvarGetType(getPROBDATA()->vars[i]) != SCIP_VARTYPE_CONTINUOUS)
      return false;
  }
  return true;
}

} // namespace mp

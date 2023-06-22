#include "mp/format.h"
#include "scipmpcommon.h"

static
SCIP_DECL_PROBDELORIG(probdataDelOrigNl)
{
  int i;

  assert((*probdata)->vars != NULL || (*probdata)->nvars == 0);
  assert((*probdata)->linconss != NULL || (*probdata)->nlinconss == 0);

  for( i = 0; i < (*probdata)->nlinconss; ++i )
  {
    SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->linconss[i]) );
  }
  SCIPfreeBlockMemoryArray(scip, &(*probdata)->linconss, (*probdata)->nlinconss);

  for( i = 0; i < (*probdata)->nvars; ++i )
  {
    SCIP_CCALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
  }
  SCIPfreeBlockMemoryArray(scip, &(*probdata)->vars, (*probdata)->nvars);

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

void ScipCommon::GetSolverOption(const char* key, int &value) const {
  if (SCIPparamGetType(SCIPgetParam(getSCIP(), key))==SCIP_PARAMTYPE_BOOL) {
    SCIP_Bool buffer;
    SCIP_CCALL( SCIPgetBoolParam(getSCIP(), key, &buffer) );
    value = (int)buffer;
  }
  else
    SCIP_CCALL( SCIPgetIntParam(getSCIP(), key, &value) );
}

void ScipCommon::SetSolverOption(const char* key, int value) {
  if (SCIPparamGetType(SCIPgetParam(getSCIP(), key))==SCIP_PARAMTYPE_BOOL)
    SCIP_CCALL( SCIPsetBoolParam(getSCIP(), key, value) );
  else
    SCIP_CCALL( SCIPsetIntParam(getSCIP(), key, value) );
}

void ScipCommon::GetSolverOption(const char* key, double &value) const {
  SCIP_CCALL( SCIPgetRealParam(getSCIP(), key, &value) );
}

void ScipCommon::SetSolverOption(const char* key, double value) {
  SCIP_CCALL( SCIPsetRealParam(getSCIP(), key, value) );
}

void ScipCommon::GetSolverOption(const char* key, std::string &value) const {
  char* buffer;
  SCIP_CCALL( SCIPallocBufferArray(getSCIP(), &buffer, SCIP_MAXSTRLEN) );
  if (SCIPparamGetType(SCIPgetParam(getSCIP(), key))==SCIP_PARAMTYPE_CHAR)
    SCIP_CCALL( SCIPgetCharParam(getSCIP(), key, buffer) );
  else
    SCIP_CCALL( SCIPgetStringParam(getSCIP(), key, &buffer) );
  value = buffer;
  SCIPfreeBufferArray(getSCIP(), &buffer);
}

void ScipCommon::SetSolverOption(const char* key, const std::string& value) {
  if (SCIPparamGetType(SCIPgetParam(getSCIP(), key))==SCIP_PARAMTYPE_CHAR)
    SCIP_CCALL( SCIPsetCharParam(getSCIP(), key, value.c_str()[0]) );
  else
    SCIP_CCALL( SCIPsetStringParam(getSCIP(), key, value.c_str()) );
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

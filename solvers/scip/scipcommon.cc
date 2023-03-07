#include "mp/format.h"
#include "scipcommon.h"

static
SCIP_DECL_PROBDELORIG(probdataDelOrigNl)
{
   int i;

   assert((*probdata)->vars != NULL || (*probdata)->nvars == 0);
   //assert((*probdata)->conss != NULL || (*probdata)->conss == 0);

   //for( i = 0; i < (*probdata)->nconss; ++i )
   //{
   //   SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->conss[i]) );
   //}
   //SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->conss, (*probdata)->nconss);

   for( i = 0; i < (*probdata)->nvars; ++i )
   {
      SCIP_CCALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->vars, (*probdata)->nvars);

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
  // Get number of linear constraints using solver API
  int nlinconss = 0;
  for (int i = 0; i < SCIPgetNOrigConss(getSCIP()); i++) {
    SCIP_CONS* cons = SCIPgetOrigConss(getSCIP())[i];
    if (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "linear") == 1)
      nlinconss++;
  }
  return nlinconss;
}

int ScipCommon::NumVars() const {
  // Get number of active problem variables
  return getPROBDATA()->nvars;
}

int ScipCommon::NumObjs() const {
  // Get number of objectives using solver API
  return 1;
}

int ScipCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(SCIP_INTATTR_QCONSTRS);
  return 0;
}

int ScipCommon::NumSOSCons() const {
  // TODO Get number of SOS constraints using solver API
  // return getIntAttr(SCIP_INTATTR_SOSS);
  return 0;
}

int ScipCommon::NumIndicatorCons() const {
  // TODO Get number of indicator constraints using solver API
  // return getIntAttr(SCIP_INTATTR_INDICATORS);
  return 0;
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
  SCIP_CCALL( SCIPallocBlockMemoryArray(getSCIP(), &buffer, SCIP_MAXSTRLEN) );
  if (SCIPparamGetType(SCIPgetParam(getSCIP(), key))==SCIP_PARAMTYPE_CHAR)
    SCIP_CCALL( SCIPgetCharParam(getSCIP(), key, buffer) );
  else
    SCIP_CCALL( SCIPgetStringParam(getSCIP(), key, &buffer) );
  value = buffer;
  SCIPfreeBlockMemoryArray(getSCIP(), &buffer, SCIP_MAXSTRLEN);
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

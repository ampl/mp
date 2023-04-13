#include "mp/format.h"
#include "gcgcommon.h"

#include "gcg/gcgplugins.h"

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
    GCG_CCALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
  }
  SCIPfreeBlockMemoryArray(scip, &(*probdata)->vars, (*probdata)->nvars);

  SCIPfreeMemory(scip, probdata);

  return SCIP_OKAY;
}
namespace mp {

void GcgCommon::OpenSolver() {
  int status = 0;
  SCIP* scip = NULL;
  SCIP_PROBDATA* probdata = NULL;

  // initialize SCIP
  status = SCIPcreate(&scip);
  setSCIP(scip); // Assign it

  // include default GCG plugins
  GCG_CCALL( SCIPincludeGcgPlugins(scip) );

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
  SCIP* scip = getSCIP();

  // free SCIP
  GCG_CCALL( SCIPfree(&scip) );
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


void GcgCommon::GetSolverOption(const char* key, int &value) const {
  if (SCIPparamGetType(SCIPgetParam(getSCIP(), key))==SCIP_PARAMTYPE_BOOL) {
    SCIP_Bool buffer;
    GCG_CCALL( SCIPgetBoolParam(getSCIP(), key, &buffer) );
    value = (int)buffer;
  }
  else
    GCG_CCALL( SCIPgetIntParam(getSCIP(), key, &value) );
}

void GcgCommon::SetSolverOption(const char* key, int value) {
  if (SCIPparamGetType(SCIPgetParam(getSCIP(), key))==SCIP_PARAMTYPE_BOOL)
    GCG_CCALL( SCIPsetBoolParam(getSCIP(), key, value) );
  else
    GCG_CCALL( SCIPsetIntParam(getSCIP(), key, value) );
}

void GcgCommon::GetSolverOption(const char* key, double &value) const {
  GCG_CCALL( SCIPgetRealParam(getSCIP(), key, &value) );
}

void GcgCommon::SetSolverOption(const char* key, double value) {
  GCG_CCALL( SCIPsetRealParam(getSCIP(), key, value) );
}

void GcgCommon::GetSolverOption(const char* key, std::string &value) const {
  char* buffer;
  GCG_CCALL( SCIPallocBufferArray(getSCIP(), &buffer, SCIP_MAXSTRLEN) );
  if (SCIPparamGetType(SCIPgetParam(getSCIP(), key))==SCIP_PARAMTYPE_CHAR)
    GCG_CCALL( SCIPgetCharParam(getSCIP(), key, buffer) );
  else
    GCG_CCALL( SCIPgetStringParam(getSCIP(), key, &buffer) );
  value = buffer;
  SCIPfreeBufferArray(getSCIP(), &buffer);
}

void GcgCommon::SetSolverOption(const char* key, const std::string& value) {
  if (SCIPparamGetType(SCIPgetParam(getSCIP(), key))==SCIP_PARAMTYPE_CHAR)
    GCG_CCALL( SCIPsetCharParam(getSCIP(), key, value.c_str()[0]) );
  else
    GCG_CCALL( SCIPsetStringParam(getSCIP(), key, value.c_str()) );
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

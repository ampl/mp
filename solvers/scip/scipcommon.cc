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

namespace Solver {

  SolverModel* CreateSolverModel() {
    return new SolverModel();
  }
  int addContVar(SolverModel& s) {
    return s.addEntity(Solver::VARS_INT);
  }
  int addIntVar(SolverModel& s) {
    return s.addEntity(Solver::VARS);
  }
  int addLinCon(SolverModel& s) {
    return s.addEntity(Solver::CONS_LIN);
  }
  int addQuadCon(SolverModel& s) {
    return s.addEntity(Solver::CONS_QUAD);
  }
  int addIndicCon(SolverModel& s) {
    return s.addEntity(Solver::CONS_INDIC);
  }
  int addSOSCon(SolverModel& s) {
    return s.addEntity(Solver::CONS_SOS);
  }
  int getSolverIntAttr(SolverModel* s, int attr) {
    return s->getN(attr);
  }
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

  // release all variables
  //for (int i = 0; i < getPROBDATA()->nvars; i++) {
  //  SCIP_CCALL( SCIPreleaseVar(getSCIP(), &getPROBDATA()->vars[i]) );
  //}

  // free SCIP
  SCIP_CCALL( SCIPfree(&scip) );
}

int ScipCommon::getIntAttr(int name)  const {
  int value = 0;
  /* TODO Utility function to get the value of an integer attribute 
  * from the solver API 
  SCIP_CCALL(SCIP_GetIntAttr(lp_, name, &value)); */
  return value;//getSolverIntAttr(getSCIP(), name);
}
double ScipCommon::getDblAttr(const char* name) const  {
  double value = 0;
  /* TODO Utility function to get the value of an integer attribute
 * from the solver API
  SCIP_CCALL(SCIP_GetDblAttr(lp_, name, &value)); */
  return value;
}

int ScipCommon::NumLinCons() const {
  // TODO Get number of linear constraints using solver API
  // return getIntAttr(SCIP_INTATTR_ROWS);
  return 0;
}

int ScipCommon::NumVars() const {
  // gets number of active problem variables
  return SCIPgetNVars(getSCIP());
}

int ScipCommon::NumObjs() const {
  // TODO Get number of objectives using solver API
  //return SCIPgetnumobjs (env_, lp_);
  return 0;
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
  SCIP_CCALL( SCIPgetIntParam(getSCIP(), key, &value) );
}

void ScipCommon::SetSolverOption(const char* key, int value) {
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
  SCIP_CCALL( SCIPgetStringParam(getSCIP(), key, &buffer) );
  value = buffer;
  SCIPfreeBlockMemoryArray(getSCIP(), &buffer, SCIP_MAXSTRLEN);
}

void ScipCommon::SetSolverOption(const char* key, const std::string& value) {
  SCIP_CCALL( SCIPsetStringParam(getSCIP(), key, value.c_str()) );
}


double ScipCommon::Infinity() const { 
  double inf;
  GetSolverOption("numerics/infinity", inf);
  return inf;  
}

double ScipCommon::MinusInfinity() { return -Infinity(); }


} // namespace mp

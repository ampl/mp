#ifndef SCIPCOMMON_H
#define SCIPCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"

//extern "C" {
  #include "scip/scip.h"
  #include "scip/scipdefplugins.h"
//}

/// problem data stored in SCIP
struct SCIP_ProbData
{
   SCIP_VAR**            vars;               /**< variables in the order given by AMPL */
   int                   nvars;              /**< number of variables */

   SCIP_CONS**           conss;              /**< constraints in the order given by AMPL */
   int                   nconss;             /**< number of constraints */
};

/// Instead, faking a typical solver namespace and defs:
namespace Solver {
  enum TYPE {
    VARS = 0,
    VARS_INT,
    CONS_LIN,
    CONS_QUAD,
    CONS_INDIC,
    CONS_SOS,

    OBJ
  };

  class SolverModel {
    int nEntities[7];
  public:
    int addEntity(int type) {
      return ++nEntities[type];
    }
    int getN(int type) {
      return nEntities[type];
    }
  };

  SolverModel* CreateSolverModel();
  int addContVar(SolverModel& s);
  int addIntVar(SolverModel& s);
  int addLinCon(SolverModel& s);
  int addQuadCon(SolverModel& s);
  int addIndicCon(SolverModel& s);
  int addSOSCon(SolverModel& s);
  int getSolverIntAttr(SolverModel* s, int attr);
}


/// The below would go into actual ...common.h:

#include "mp/format.h"

namespace mp {

/// Information shared by both
/// `ScipBackend` and `ScipModelAPI`
struct ScipCommonInfo {
  SCIP* getSCIP() const { return scip_; }
  void setSCIP(SCIP* scip) { scip_ = scip; }

  SCIP_PROBDATA* getPROBDATA() { return probdata_; }
  void setPROBDATA(SCIP_PROBDATA* probdata) { probdata_ = probdata; }

private:
  SCIP* scip_ = NULL;
  SCIP_PROBDATA* probdata_;
};


/// Common API for Scip classes
class ScipCommon :
    public Backend2ModelAPIConnector<ScipCommonInfo> {
public:
  /// These methods access Scip options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  /// TODO Typically solvers define their own infinity; use them here
  double Infinity() const;
  double MinusInfinity();

protected:
  void OpenSolver();
  void CloseSolver();

  int getIntAttr(int name) const;
  double getDblAttr(const char* name) const;

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;

protected:
  // TODO if desirable, provide function to create the solver's environment
  // with own license
  // int (*createEnv) (solver_env**) = nullptr;
  
};


/// Convenience macro
// TODO This macro is useful to automatically throw an error if a function in the 
// solver API does not return a valid errorcode. In this mock driver, we define it 
// ourselves, normally this constant would be defined in the solver's API.
#define SCIP_RETCODE_OK 1
#define SCIP_CCALL( call ) do { if (int e = (call) != SCIP_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // SCIPCOMMON_H

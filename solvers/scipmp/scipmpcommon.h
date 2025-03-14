#ifndef SCIPCOMMON_H
#define SCIPCOMMON_H

#include <vector>
#include <string>

#include "mp/backend-to-model-api.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "mp/format.h"

/// problem data stored in SCIP
struct SCIP_ProbData
{
  SCIP_VAR**            vars;               /**< variables in the order given by AMPL */
  std::vector<SCIP_EXPR*>  var_exprs;       /**< expressions for the variables, when non-zero */
  int                   nvars;              /**< number of variables */

  std::vector<SCIP_CONS*>  linconss;        /**< linear constraints in the order given by AMPL */
  int                   i = 0;              /**< shows free slot of linear constraints */
  int                   nlinconss = 0;      /**< number of linear constraints */

  std::vector<SCIP_CONS*>  nlconss;         /**< NL linear constraints in the order given by AMPL */

  std::vector<SCIP_EXPR*> exprs;            /**< expressions*/

  SCIP_EXPR* dummyexpr {nullptr};
};

namespace mp {

/// Information shared by both
/// `ScipBackend` and `ScipModelAPI`
struct ScipCommonInfo {
  SCIP* getSCIP() const { return scip_; }
  void setSCIP(SCIP* scip) { scip_ = scip; }

  SCIP_PROBDATA* getPROBDATA() const { return probdata_; }
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
	void GetSolverOption(const char* key, SCIP_Longint& value) const;
	void SetSolverOption(const char* key, SCIP_Longint value);
	void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  /// SCIP own infinity
  double Infinity() const;
  double MinusInfinity();

  bool IsContinuous();

protected:
  void OpenSolver();
  void CloseSolver();

  int NumLinCons() const;
  int NumVars() const;
  int NumObjs() const;
  int NumQPCons() const;
  int NumSOSCons() const;
  int NumIndicatorCons() const;
};


/// Convenience macro
// This macro is useful to automatically throw an error if a function in the 
// solver API does not return a valid errorcode. In this mock driver, we define it 
// ourselves, normally this constant would be defined in the solver's API.
#define SCIP_RETCODE_OK 1
#define SCIP_CCALL( scip__call_ ) do { if (int err__code_ = (scip__call_) != SCIP_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #scip__call_, err__code_ )); } while (0)

} // namespace mp

#endif // SCIPCOMMON_H

#ifndef GCGCOMMON_H
#define GCGCOMMON_H

#include <string>

#include "mp/backend-to-model-api.h"

#include "scip/scip.h"
#include "gcg/gcg.h"

#include "gcg/class_partialdecomp.h"

#include "mp/format.h"

/// problem data stored in SCIP
struct SCIP_ProbData
{
   SCIP_VAR**            vars;               /**< variables in the order given by AMPL */
   int                   nvars;              /**< number of variables */

   SCIP_CONS**           linconss;           /**< linear constraints in the order given by AMPL */
   int                   i;                  /**< shows free slot of linear constraints */
   int                   nlinconss;          /**< number of linear constraints */

   gcg::PARTIALDECOMP*   decomp;             /**< user partialdecomp */
};

namespace mp {

/// Information shared by both
/// `GcgBackend` and `GcgModelAPI`
struct GcgCommonInfo {
  SCIP* getSCIP() const { return scip_; }
  void setSCIP(SCIP* scip) { scip_ = scip; }

  SCIP_PROBDATA* getPROBDATA() const { return probdata_; }
  void setPROBDATA(SCIP_PROBDATA* probdata) { probdata_ = probdata; }

private:
  SCIP* scip_ = NULL;
  SCIP_PROBDATA* probdata_;
};


/// Common API for Gcg classes
class GcgCommon :
    public Backend2ModelAPIConnector<GcgCommonInfo> {
public:
  /// These methods access Gcg options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

  /// GCG own infinity
  double Infinity() const;
  double MinusInfinity() const;

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
#define GCG_RETCODE_OK 1
#define GCG_CCALL( call ) do { if (int e = (call) != GCG_RETCODE_OK) \
  throw std::runtime_error( \
    fmt::format("  Call failed: '{}' with code {}", #call, e )); } while (0)

} // namespace mp

#endif // GCGCOMMON_H

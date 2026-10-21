#ifndef MP_CPSAT_BACKEND_H_
#define MP_CPSAT_BACKEND_H_

#include <map>
#include <string>
#include <variant>
#include <vector>
#include <utility>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "cpsatcommon.h"

namespace mp {

/// AMPL driver for Google OR-Tools CP-SAT.
class CpsatBackend :
    public FlatBackend< MIPBackend<CpsatBackend> >,
    public CpsatCommon
{
  using BaseBackend = FlatBackend< MIPBackend<CpsatBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  /// Construct
  CpsatBackend();
  /// Destruct
  ~CpsatBackend();

  /// Prefix used for the <prefix>_options environment variable
  static const char* GetAMPLSolverName() { return "cpsat"; }
  /// AMPL driver name displayed in messages
  static const char* GetAMPLSolverLongName() { return "AMPL-CPSAT"; }
  /// Solver name displayed in messages
  static const char* GetSolverName() { return "CP-SAT"; }
  /// Version displayed with -v
  std::string GetSolverVersion();
  /// External libraries displayed with -v
  std::string set_external_libs() override { return ""; }

  /// Name for diagnostic messages
  static const char* GetBackendName();
  /// "long name", rarely used
  static const char* GetBackendLongName() { return nullptr; }

  /// Init custom driver options, such as outlev, writeprob
  void InitCustomOptions() override;
  /// Chance for the Backend to init solver environment, etc.
  void InitOptionParsing() override;
  /// Chance to consider options immediately
  void FinishOptionParsing() override;


  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  USING_STD_FEATURES;

  /// Report all improving solutions found during the search
  ALLOW_STD_FEATURE(MULTISOL, true)

  /// Get MIP gap (adds option mip:return_gap)
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, true)
  double MIPGap() override;
  double MIPGapAbs() override;

  /// Get best objective bound (adds option mip:bestbound)
  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, true)
  double BestDualBound() override;

  /// Model export (adds options tech:writeprob, tech:writeprobindex):
  /// the CP-SAT model in protobuf text format
  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override;

  /////////////////////////// Model attributes /////////////////////////
  bool IsMIP() const override { return true; }
  bool IsQCP() const override { return false; }

  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

  /// Solve, no model modification any more.
  void Solve() override;

  /// Objective values (single objective)
  ArrayRef<double> GetObjectiveValues() override
  { return std::vector<double>{ObjectiveValue()}; }


  //////////////////// [[ Implementation details ]] //////////////////////
protected:
  void OpenSolver();
  void CloseSolver();

  double ObjectiveValue() const;
  bool HasSolution() const;

  /// PrimalSolution() for flat backends.
  /// @return empty vector if no primal solution.
  ArrayRef<double> PrimalSolution() override;
  /// CP-SAT has no duals
  pre::ValueMapDbl DualSolution() override { return {}; }

  void ReportResults() override;
  void ReportCpsatResults();
  void ReportCpsatPool();

  std::pair<int, std::string> GetSolveResult() override;
  void AddCpsatMessages();

  /// Statistics for option tech:stats
  std::map<std::string, std::variant<int, double, std::string>>
      SolutionStats() override;

private:
  /// These options are stored in the class
  struct Options {
    std::string params_;         ///< SatParameters in text format
    int outlev_ = 0;
    double varbound_ = kCpsatDefaultVarBound;
  };
  Options storedOptions_;

  /// Improving solutions found, for MULTISOL
  std::vector< std::pair<std::vector<double>, double> > pool_;
};

}  // namespace mp

#endif  // MP_CPSAT_BACKEND_H_

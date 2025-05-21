#ifndef MP_CUOPTLP_BACKEND_H_
#define MP_CUOPTLP_BACKEND_H_

#include <vector>
#include <list>
#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "cuoptlpcommon.h"

namespace mp {

class CuoptlpBackend :
    public FlatBackend< MIPBackend<CuoptlpBackend> >,
    public CuoptlpCommon
{
  using BaseBackend = FlatBackend< MIPBackend<CuoptlpBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  /// Construct
  CuoptlpBackend();
  /// Destruct
  ~CuoptlpBackend();

  /// Prefix used for the <prefix>_options environment variable
  static const char* GetAMPLSolverName() { return "cuoptlp"; }

  /// AMPL driver name displayed in messages
  static const char* GetAMPLSolverLongName() { return "AMPL-CUOPT"; }
  /// Solver name displayed in messages
  static const char* GetSolverName() { return "CUOPT"; }
  /// Version displayed with -v
  std::string GetSolverVersion();
  /// External libraries displayed with -v
  std::string set_external_libs() override { return ""; };

  /// Name for diagnostic messages
  static const char* GetBackendName();
  /// "long name", rarely used
  static const char* GetBackendLongName() { return nullptr; }

  /// Init custom driver options, such as outlev, writeprob
  void InitCustomOptions() override;
  /// Chance for the Backend to init solver environment, etc.
  void InitOptionParsing() override { }
  /// Chance to consider options immediately (open cloud, etc)
  void FinishOptionParsing() override;

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  // Use this section to declare and implement some standard features
  // that may or may not need additional functions.
  // For a full list of features possible,
  // grep "STD_FEATURE".
  USING_STD_FEATURES;

  /**
  * EXPORT SOLUTION
  **/
  ALLOW_STD_FEATURE(WRITE_SOLUTION, true)
  void DoWriteSolution(const std::string& name) override { }

  /**
  * Get MIP Gap
  **/
  // (adds option mip:return_gap)
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, true)
  double MIPGap() override;
  double MIPGapAbs() override;

  /**
  * Get MIP dual bound
  **/
  // (adds option mip:bestbound)
  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, true)
  double BestDualBound() override;

  /////////////////////////// Model attributes /////////////////////////

  bool IsMIP() const override;

  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;


public:  // public for static polymorphism
  /// Solve, to be overloaded by the solver.
  /// No model modification any more.
  /// @note If using STD_FEATURE( MULTISOL ),
  /// can report intermediate results
  /// via ReportIntermediateSolution() during this
  /// (check if (need_multiple_solutions())),
  /// otherwise afterwards.
  void Solve() override;

  /// Default impl of GetObjValues()
  ArrayRef<double> GetObjectiveValues() override
  { return std::vector<double>{ObjectiveValue()}; }


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
protected:
  void OpenSolver();
  void CloseSolver();

  double ObjectiveValue() const;

	/// PrimalSolution() for flat backends.
	/// @return empty vector if no primal solution.
	ArrayRef<double> PrimalSolution() override;
	/// DualSolution() for flat backends.
	/// @return empty map if no dual solution.
	pre::ValueMapDbl DualSolution() override;
	/// Dual solution for the LP part only.
	/// @return empty vector if none.
  ArrayRef<double> DualSolution_LP();

  void WindupCUOPTLPSolve();

  void ReportResults() override;
  void ReportCUOPTLPResults();


  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> GetSolveResult() override;
  void AddCUOPTLPMessages();
};

}  // namespace mp

#endif  // MP_CUOPTLP_BACKEND_H_

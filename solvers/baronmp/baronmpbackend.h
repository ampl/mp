#ifndef MP_BARONMP_BACKEND_H_
#define MP_BARONMP_BACKEND_H_

#include <vector>
#include <list>
#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "baronmpcommon.h"



namespace mp {


  
class BaronmpBackend :
    public FlatBackend< MIPBackend<BaronmpBackend> >,
    public BaronmpCommon
{
  

  int errorLevel = 0;
  TimFileData timFileData_;
  ResFileData resFileData_;

  using BaseBackend = FlatBackend< MIPBackend<BaronmpBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  /// Construct
  BaronmpBackend();
  /// Destruct
  ~BaronmpBackend();

  /// Prefix used for the <prefix>_options environment variable
  static const char* GetAMPLSolverName() { return "baron"; }

  /// AMPL driver name displayed in messages
  static const char* GetAMPLSolverLongName() { return "AMPL-BaronMP"; }
  /// Solver name displayed in messages
  static const char* GetSolverName() { return "BaronMP"; }
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
  * EXPORT PROBLEM
  **/
  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override { }

  /**
  * EXPORT SOLUTION
  **/
  ALLOW_STD_FEATURE(WRITE_SOLUTION, true)
  void DoWriteSolution(const std::string& name) override { }

  /**
  * Compute the IIS and obtain relevant values
  **/
  ALLOW_STD_FEATURE(IIS, false) // TODO
  /// Compute IIS
  void ComputeIIS() override;
  /// Retrieve IIS. Elements correspond to IISStatus
  IIS GetIIS() override;

  /**
  * General warm start
  * set primal/dual initial guesses for continuous case
  **/
  ALLOW_STD_FEATURE( WARMSTART, true )
  void AddPrimalDualStart(Solution sol0) override;
  /**
  * Specifically, MIP warm start
  **/
  ALLOW_STD_FEATURE( MIPSTART, true )
  void AddMIPStart(
      ArrayRef<double> x0, ArrayRef<int> s0) override;

  /**
   * MULTISOL support.
   *  If (need_multiple_solutions()),
   *  call ReportIntermediateSolution() during solve or after.
   **/
  ALLOW_STD_FEATURE(MULTISOL, true)

  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, true)
  double BestDualBound() override;
  /**
   * Get MIP Gap
   **/
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, true)
  double MIPGap() override;
  double MIPGapAbs() override;


    /////////////////////////// Model attributes /////////////////////////

  /// Reimplement if the solver gives more information
  /// than just the number of non-fixed integer variables
  /// (e.g., the solver might consider if it has PL expressions.)
  bool IsMIP() const override;
  bool IsQCP() const override;

  std::map<std::string, std::variant<int, double, std::string>>
	  SolutionStats() override;
  
  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;


public: 
  /// Solve, to be overloaded by the solver.
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

  void WindupBARONMPSolve();

  void ReportResults() override;
  void ReportBARONMPResults();

  void ReportBARONMPPool();

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> GetSolveResult() override;
  void AddBARONMPMessages();

  ArrayRef<int> VarsIIS();
  pre::ValueMapInt ConsIIS();



};

}  // namespace mp

#endif  // MP_BARONMP_BACKEND_H_

#ifndef MP_MOSEK_BACKEND_H_
#define MP_MOSEK_BACKEND_H_

#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "mosekcommon.h"

namespace mp {

class MosekBackend :
    public FlatBackend< MIPBackend<MosekBackend> >,
    public MosekCommon
{
  using BaseBackend = FlatBackend< MIPBackend<MosekBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  MosekBackend();
  ~MosekBackend();

  /*----------------------------------------------------*\
  | Standard and optional methods to provide or retrieve |
  | information to/from or manipulate the solver. Most   |
  | of them override placeholders from base classes.     |
  \*----------------------------------------------------*/

  /// Name displayed in messages
  static const char* GetSolverName() { return "MOSEK"; }
  std::string GetSolverVersion();

  /// AMPL solver name is used to parse solver options
  /// for the [name]_options environment variable.
  /// This is only done if the [executable_name]_options
  /// variable is not provided.
  static const char* GetAMPLSolverName() { return "mosek"; }
  static const char* GetAMPLSolverLongName() { return "AMPL-MOSEK"; }

  static const char* GetBackendName();
  static const char* GetBackendLongName() { return nullptr; }

  /// Chance for the Backend to init solver environment, etc
  void InitOptionParsing() override;
  /// Chance to consider options immediately (open cloud, etc)
  void FinishOptionParsing() override;


  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  // Use this section to declare and implement some standard features
  // that may or may not need additional functions. 
  USING_STD_FEATURES;

  // Export problem
  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override;

  // LP basis info, status keys
  ALLOW_STD_FEATURE(BASIS, true)
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;

  // LP hotstart
  // Set primal/dual initial guesses for continuous case
  ALLOW_STD_FEATURE(WARMSTART, true)
  void AddPrimalDualStart(Solution) override;

  // MIP hotstart
  ALLOW_STD_FEATURE( MIPSTART, true )
  void AddMIPStart(ArrayRef<double> x0) override;

  // Obtain inf/unbounded rays
  ALLOW_STD_FEATURE(RAYS, true)
  ArrayRef<double> Ray() override;
  ArrayRef<double> DRay() override;

  // MIP gap
  // (adds option mip:return_gap)
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, true)
  double MIPGap() override;
  double MIPGapAbs() override;

  // MIP dual bound
  // (adds option mip:bestbound)
  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, true)
  double BestDualBound() override;

  // Sensitivity analysis
  // Report sensitivity analysis suffixes (postsolved)
  ALLOW_STD_FEATURE(SENSITIVITY_ANALYSIS, true)
  SensRangesPresolved GetSensRangesPresolved() override;

  /////////////////////////// Model attributes /////////////////////////
  bool IsMIP() const override;
  bool IsQCP() const override;
  
  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

  /// Solve, no model modification any more.
  /// Can report intermediate results via HandleFeasibleSolution() during this,
  /// otherwise in ReportResults()
  void Solve() override;

  ArrayRef<double> GetObjectiveValues() override
  {
    return std::vector<double>{ObjectiveValue()};
  }


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
public:  // public for static polymorphism
  void InitCustomOptions() override;

protected:
  void OpenSolver();
  void CloseSolver();

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
	ArrayRef<double> DualSolution_LP();
	ArrayRef<double> DualSolution_Cones();

  void WindupMOSEKSolve();

  void ReportResults() override;
  void ReportMOSEKResults();

  void ReportMOSEKPool();

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  /// Solution + termination status
  std::pair<int, std::string> ConvertMOSEKStatus();
  /// Text to add for termination status
  std::string ConvertMOSEKTermStatus();
  void AddMOSEKMessages();

  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  void VarStatii(ArrayRef<int>);
  void ConStatii(ArrayRef<int>);

  //ArrayRef<int> VarsIIS();
  //pre::ValueMapInt ConsIIS();

  // utility function used to decide which solution to fetch after optimization
  // sets member solToFetch_ as side effect
  MSKsoltypee GetSolutionTypeToFetch();

private:
  MSKsoltypee solToFetch_;
  MSKrescodee termCode_;
  MSKprostae proSta_;
  MSKsolstae solSta_;

  /// These options are stored in the class
  struct Options {
    // Whether to set MSK_IPAR_MIO_CONSTRUCT_SOL
    int MIPConstructSol_=0;
  };
  Options storedOptions_;

protected:
  /**** Option accessors ****/
  int Mosek_mip_construct_sol() const { return storedOptions_.MIPConstructSol_; }

};

}  // namespace mp

#endif  // MP_MOSEK_BACKEND_H_

#ifndef MP_HIGHS_BACKEND_H_
#define MP_HIGHS_BACKEND_H_

#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "highsmpcommon.h"

namespace mp {

class HighsBackend :
    public FlatBackend< MIPBackend<HighsBackend> >,
    public HighsCommon
{
  using BaseBackend = FlatBackend< MIPBackend<HighsBackend> >;

  std::vector<int> conStatiii_;
  //////////////////// [[ The public interface ]] //////////////////////
public:
  HighsBackend();
  ~HighsBackend();

  /// Name displayed in messages
  static const char* GetSolverName() { return "HiGHS"; }
  std::string GetSolverVersion();
  
  static const char* GetAMPLSolverName() { return "highs"; }
  static const char* GetAMPLSolverLongName() { return "AMPL-HiGHS"; }
  static const char* GetBackendName();
  static const char* GetBackendLongName() { return nullptr; }



  /// This can actually modify the model -- e.g., suffixes
  void InputExtras() override;


  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  // Use this section to declare and implement some standard features
  // that may or may not need additional functions. 
  USING_STD_FEATURES;

  /**
  * MULTIOBJ
  **/
  // Note: native MO is currently disabled as its implementation is 
  // somewhat inconsitent with other solvers
  ALLOW_STD_FEATURE( MULTIOBJ, false )
  ArrayRef<double> GetObjectiveValues() override;
  void ObjPriorities(ArrayRef<int>) override;
  void ObjWeights(ArrayRef<double>) override;
  void ObjAbsTol(ArrayRef<double>) override;
  void ObjRelTol(ArrayRef<double>) override;
  /**
  * EXPORT PROBLEM
  **/
  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override;

  /**
  * EXPORT SOLUTION
  **/
  ALLOW_STD_FEATURE(WRITE_SOLUTION, true)
  void DoWriteSolution(const std::string& name) override;

  /**
  * General warm start:
  * set primal/dual initial guesses for continuous case
  **/
  ALLOW_STD_FEATURE(WARMSTART, true)
    void AddPrimalDualStart(Solution sol0) override;
  /**
  * Get/Set AMPL var/con statii
  **/
  ALLOW_STD_FEATURE(BASIS, true)
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;

 /**
  * Get MIP Gap
  **/
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, true)
  double MIPGap() override;
  double MIPGapAbs() override;
  /**
  * Get MIP dual bound
  **/
  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, true)
  double BestDualBound() override;

  /**
  * Obtain inf/unbounded rays
  **/
  ALLOW_STD_FEATURE(RAYS, true)
    ArrayRef<double> Ray() override;
  ArrayRef<double> DRay() override;

  /**
  * Stop the MIP search on a plateau (mip:plateau* options)
  **/
  ALLOW_STD_FEATURE( PLATEAU_STOP, true )
  ALLOW_STD_FEATURE(PLATEAU_STOP_BOUND, true)
  void SetupPlateauCallbacks() override;

  /////////////////////////// Model attributes /////////////////////////
  bool IsQCP() const override;
  
  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

  /// Solve, no model modification any more.
  /// Can report intermediate results via ReportIntermediateSolution() during this,
  /// otherwise in ReportResults()
  void Solve() override;

  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
public:  
  void InitCustomOptions() override;
  void FinishOptionParsing() override;

protected:

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();

  void WindupHIGHSSolve();

  /// Native callback implementing PLATEAU_STOP / PLATEAU_STOP_BOUND:
  /// forwards new incumbents (kHighsCallbackMipImprovingSolution) to
  /// ReportIncumbentForPlateau(), and gap updates (from
  /// kHighsCallbackMipInterrupt and kHighsCallbackMipSolution --
  /// absgap/relgap derived from mip_primal_bound/mip_dual_bound/mip_gap)
  /// to ReportGapForPlateau(), requesting termination via
  /// data_in->user_interrupt when either signals a plateau.
  static void DoPlateauCallback(
      int callback_type, const char* message,
      const HighsCallbackDataOut* data_out, HighsCallbackDataIn* data_in,
      void* user_callback_data);

  void ReportResults() override;
  void ReportHIGHSResults();

  /// Solution attributes
  int NodeCount() const;
  int SimplexIterations() const;
  int BarrierIterations() const;
  int PdlpIterations() const;

  std::map<std::string, std::variant<int, double, std::string>>
      SolutionStats() override;

  std::pair<int, std::string> GetSolveResult() override;
  void AddHIGHSMessages();
  
  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  void VarConStatii(ArrayRef<int>, ArrayRef<int>);

private:
  /// These options are stored in the class
  struct Options {
    std::string lpmethod_ = "choose";
    bool onGPU() { return lpmethod_ == "pdlp-gpu" || lpmethod_ == "hipdlp-gpu"; }
  };
  Options storedOptions_;


};

}  // namespace mp

#endif  // MP_HIGHS_BACKEND_H_

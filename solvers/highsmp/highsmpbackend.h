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

  /// Chance for the Backend to init solver environment, etc
  void InitOptionParsing() override { }
  /// Chance to consider options immediately (open cloud, etc)
  void FinishOptionParsing() override;



  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  // Use this section to declare and implement some standard features
  // that may or may not need additional functions. 
  USING_STD_FEATURES;

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


  /////////////////////////// Model attributes /////////////////////////
  bool IsQCP() const override;
  
  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

  /// Solve, no model modification any more.
  /// Can report intermediate results via HandleFeasibleSolution() during this,
  /// otherwise in ReportResults()
  void Solve() override;

  ArrayRef<double> GetObjectiveValues() override
  { return std::vector<double>{ObjectiveValue()}; } 


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
public:  // public for static polymorphism
  void InitCustomOptions() override;

protected:

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();

  void WindupHIGHSSolve();

  void ReportResults() override;
  void ReportHIGHSResults();

  void ReportHIGHSPool();

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> ConvertHIGHSStatus();
  void AddHIGHSMessages();
  
  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  void VarConStatii(ArrayRef<int>, ArrayRef<int>);

private:
  /// These options are stored in the class
  struct Options {
  };
  Options storedOptions_;


};

}  // namespace mp

#endif  // MP_HIGHS_BACKEND_H_

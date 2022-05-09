#ifndef MP_HIGHS_BACKEND_H_
#define MP_HIGHS_BACKEND_H_

#if __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunused-parameter"
# pragma clang diagnostic ignored "-Wunused-private-field"
#elif _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4244)
#endif

#if __clang__
# pragma clang diagnostic pop
#elif _MSC_VER
# pragma warning(pop)
#endif

#include <string>

#include "mp/backend-mip.h"
#include "mp/backend-with-pre.h"
#include "highscommon.h"

namespace mp {

class HighsBackend :
    public MIPBackend<HighsBackend>,
    public BackendWithPresolver,
    public HighsCommon
{
  using BaseBackend = MIPBackend<HighsBackend>;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  HighsBackend();
  ~HighsBackend();

  /// Name displayed in messages
  static const char* GetSolverName() { return "x-HIGHS"; }
  std::string GetSolverVersion();
  
  static const char* GetAMPLSolverName() { return "highs"; }
  static const char* GetAMPLSolverLongName() { return "AMPL-HIGHS"; }
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
 * MULTISOL support
 * No API, see ReportIntermediateSolution()
**/
  ALLOW_STD_FEATURE(MULTISOL, false)

  /**
   /**
  * Get/Set AMPL var/con statii
  **/
  ALLOW_STD_FEATURE(BASIS, true)
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;

  /**
  * MIP warm start
  **/
  // TODO If MIP warm start is supported, implement the function below
  // to set a non-presolved starting solution
  ALLOW_STD_FEATURE(MIPSTART, false)
  void AddMIPStart(ArrayRef<double> x0) override;


 /**
  * Get MIP Gap
  **/
  // TODO Implement to return MIP gap
  // (adds option mip:return_gap)
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, false)
  double MIPGap() override;
  double MIPGapAbs() override;
  /**
  * Get MIP dual bound
  **/
  // TODO Implement to return the best dual bound value
  // (adds option mip:bestbound)
  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, false)
  double BestDualBound() override;


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

  Solution GetSolution() override;
  ArrayRef<double> GetObjectiveValues() override
  { return std::vector<double>{ObjectiveValue()}; } 


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
public:  // public for static polymorphism
  void InitCustomOptions() override;

protected:

  void ExportModel(const std::string& file);

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution();
  pre::ValueMapDbl DualSolution();
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
  void VarStatii(ArrayRef<int>);
  void ConStatii(ArrayRef<int>);

private:
  /// These options are stored in the class
  struct Options {
    std::string exportFile_;
  };
  Options storedOptions_;


};

}  // namespace mp

#endif  // MP_HIGHS_BACKEND_H_

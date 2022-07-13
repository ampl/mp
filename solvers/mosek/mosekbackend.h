#ifndef MP_MOSEK_BACKEND_H_
#define MP_MOSEK_BACKEND_H_

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
#include "mp/valcvt-base.h"
#include "mosekcommon.h"

namespace mp {

class MosekBackend :
    public MIPBackend<MosekBackend>,
    public BasicValuePresolverKeeper,
    public MosekCommon
{
  using BaseBackend = MIPBackend<MosekBackend>;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  MosekBackend();
  ~MosekBackend();

  /// Name displayed in messages
  static const char* GetSolverName() { return "MOSEK"; }
  std::string GetSolverVersion();
  
  static const char* GetAMPLSolverName() { return "mosek"; }
  static const char* GetAMPLSolverLongName() { return "AMPL-MOSEK"; }
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

  ALLOW_STD_FEATURE(BASIS, true)
  // TODO If getting/setting a basis is supported, implement the 
  // accessor and the setter below
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;

 /**
  * Get MIP Gap
  **/
  // TODO Implement to return MIP gap
  // (adds option mip:return_gap)
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, true)
  double MIPGap() override;
  double MIPGapAbs() override;
  /**
  * Get MIP dual bound
  **/
  // TODO Implement to return the best dual bound value
  // (adds option mip:bestbound)
  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, true)
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
  void OpenSolver();
  void CloseSolver();

  void ExportModel(const std::string& file);

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution();
  pre::ValueMapDbl DualSolution();
  ArrayRef<double> DualSolution_LP();

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

  std::pair<int, std::string> ConvertMOSEKStatus();
  void AddMOSEKMessages();

  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  void VarStatii(ArrayRef<int>);
  void ConStatii(ArrayRef<int>);

  ArrayRef<int> VarsIIS();
  pre::ValueMapInt ConsIIS();

  // utility function used to decide which solution to fetch after optimization
  // sets member solToFetch_ as side effect
  MSKsoltypee GetSolutionTypeToFetch();

private:
  MSKsoltypee solToFetch_;
  MSKrescodee termCode_;
  MSKsolstae solSta_;

  /// These options are stored in the class
  struct Options {
    std::string exportFile_;
  };
  Options storedOptions_;


};

}  // namespace mp

#endif  // MP_MOSEK_BACKEND_H_

#ifndef MP_ORTOOLS_BACKEND_H_
#define MP_ORTOOLS_BACKEND_H_

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
#include "ortoolsmpcommon.h"

namespace mp {

class OrtoolsBackend :
    public MIPBackend<OrtoolsBackend>,
    public BasicValuePresolverKeeper,
    public OrtoolsCommon
{
  using BaseBackend = MIPBackend<OrtoolsBackend>;


  //////////////////// [[ The public interface ]] //////////////////////
public:
  OrtoolsBackend();
  ~OrtoolsBackend();

  /// Name displayed in messages
  static const char* GetSolverName() { return "ortools"; }
  std::string GetSolverVersion();
  
  static const char* GetAMPLSolverName() { return "ortools"; }
  static const char* GetAMPLSolverLongName() { return "AMPL-ortools"; }
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

  ALLOW_STD_FEATURE(MULTISOL, true)

    /**
   * Get/Set AMPL var/con statii
   **/
  ALLOW_STD_FEATURE(BASIS, true)
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;

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

  void WindupORTOOLSSolve();

  void ReportResults() override;
  void ReportORTOOLSResults();
  void ReportORTOOLSPool();

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> ConvertORTOOLSStatus();
  void AddORTOOLSMessages();

  // For basis
  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();

private:
  /// These options are stored in the class
  struct Options {
    std::string exportFile_;
    std::string solver_ = "cbc";
    int outlev_ = 0;
    int threads_ = 0;
    double timelimit_ = 0;
  };
  Options storedOptions_;

  operations_research::MPSolver::ResultStatus status_;

};

}  // namespace mp

#endif  // MP_ORTOOLS_BACKEND_H_

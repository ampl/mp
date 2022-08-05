#ifndef MP_CPLEX_BACKEND_H_
#define MP_CPLEX_BACKEND_H_

#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "cplexmpcommon.h"

namespace mp {

class CplexBackend :
    public MIPBackend<CplexBackend>,
    public BasicFlatBackend,
    public CplexCommon
{
  using BaseBackend = MIPBackend<CplexBackend>;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  CplexBackend();
  ~CplexBackend();

  /// Name displayed in messages
  static const char* GetSolverName() { return "x-CPLEX"; }
  std::string GetSolverVersion();
  /// Reuse cplex_options:
  static const char* GetAMPLSolverName() { return "cplex"; }
  static const char* GetAMPLSolverLongName() { return "AMPL-CPLEX"; }
  static const char* GetBackendName();
  static const char* GetBackendLongName() { return nullptr; }

  /// Chance for the Backend to init solver environment, etc
  void InitOptionParsing() override { }
  /// Chance to consider options immediately (open cloud, etc)
  void FinishOptionParsing() override;


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
  { return std::vector<double>{ObjectiveValue()}; } // TODO


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

  void WindupCPLEXSolve();

  void ReportResults() override;
  void ReportCPLEXResults();

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> ConvertCPLEXStatus();
  void AddCPLEXMessages();

private:
  /// These options are stored in the class
  struct Options {
    std::string exportFile_;
    std::string logFile_;
  };
  Options storedOptions_;


};

}  // namespace mp

#endif  // MP_CPLEX_BACKEND_H_

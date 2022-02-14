#ifndef MP_CPLEX_BACKEND_H_
#define MP_CPLEX_BACKEND_H_

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
#include "cplexcommon.h"

namespace mp {

class CplexBackend :
    public MIPBackend<CplexBackend>,
    public BackendWithPresolver,
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
  void SetInterrupter(mp::Interrupter* inter) override;
  void SolveAndReportIntermediateResults() override;

  Solution GetSolution() override;
  ArrayRef<double> GetObjectiveValues() override
  { return std::vector<double>{ObjectiveValue()}; } // TODO


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

  void WindupCPLEXSolve();

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
  };
  Options storedOptions_;


};

}  // namespace mp

#endif  // MP_CPLEX_BACKEND_H_

#ifndef MP_KNITROMP_BACKEND_H_
#define MP_KNITROMP_BACKEND_H_

#include <vector>
#include <list>
#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "knitrompcommon.h"

namespace mp {

class KnitrompBackend :
    public FlatBackend< MIPBackend<KnitrompBackend> >,
    public KnitrompCommon
{
  using BaseBackend = FlatBackend< MIPBackend<KnitrompBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  /// Construct
  KnitrompBackend();
  /// Destruct
  ~KnitrompBackend();

  /// Prefix used for the <prefix>_options environment variable
  static const char* GetAMPLSolverName() { return "knitromp"; }

  /// AMPL driver name displayed in messages
  static const char* GetAMPLSolverLongName() { return "AMPL-KNITROMP"; }
  /// Solver name displayed in messages
  static const char* GetSolverName() { return "x-KNITROMP"; }
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
  void InitOptionParsing() override;
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

  ALLOW_STD_FEATURE(WARMSTART, true)
  void AddPrimalDualStart(Solution sol0) override;

  /////////////////////////// Model attributes /////////////////////////

  /// Reimplement if the solver gives more information
  /// than just the number of non-fixed integer variables
  /// (e.g., the solver might consider if it has PL expressions.)
  bool IsMIP() const override;
  bool IsQCP() const override;
  
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


  void WindupKNITROMPSolve();

  void ReportResults() override;
  void ReportKNITROMPResults();


  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> GetSolveResult() override;
  void AddKNITROMPMessages();



private:
  /// These options are stored in the class
  struct Options {
    
    int outlev = 1;
  };
  Options storedOptions_;

protected:
  int solstatus_;

};

}  // namespace mp

#endif  // MP_KNITROMP_BACKEND_H_

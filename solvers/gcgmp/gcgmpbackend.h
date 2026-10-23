#ifndef MP_GCG_BACKEND_H_
#define MP_GCG_BACKEND_H_

#include <string>
#include <utility>
#include <vector>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "gcgmpcommon.h"

namespace mp {

class GcgBackend :
    public FlatBackend< MIPBackend<GcgBackend> >,
    public GcgCommon
{
  using BaseBackend = FlatBackend< MIPBackend<GcgBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  GcgBackend();
  ~GcgBackend();

  /// Prefix used for the <prefix>_options environment variable
  static const char* GetAMPLSolverName() { return "gcg"; }

  /// AMPL driver name displayed in messages
  static const char* GetAMPLSolverLongName() { return "AMPL-GCG"; }
  /// Solver name displayed in messages
  static const char* GetSolverName() { return "GCG"; }
  /// Version displayed with -v
  std::string GetSolverVersion();
  /// External libraries displayed with -v
  std::string set_external_libs() override;
  
  /// Name for diagnostic messages
  static const char* GetBackendName();
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
  USING_STD_FEATURES;

  /**
 * MULTISOL support
 * No API, see ReportIntermediateSolution()
**/
  ALLOW_STD_FEATURE(MULTISOL, true)

  /**
  * EXPORT PROBLEM / SOLUTION
  **/
  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override;
  ALLOW_STD_FEATURE(WRITE_SOLUTION, true)
  void DoWriteSolution(const std::string& name) override;

  /// No basis statuses: GCG solves the master problem (BASIS not allowed)

 /**
  * Get MIP Gap
  **/
  // (adds option mip:return_gap)
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, true)
  double MIPGap() override;
  double MIPGapAbs() override;
  /**
  * Get MIP dual bound
  **/
  // (adds option mip:bestbound)
  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, true)
  double BestDualBound() override;

  /////////////////////////// Model attributes /////////////////////////
  bool IsMIP() const override;
  bool IsQCP() const override;
  
  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

public:  // public for static polymorphism
  /// Solve, no model modification any more (such as feasrelax).
  /// Can report intermediate results via ReportIntermediateSolution() during this,
  /// otherwise/finally via ReportResults()
  void Solve() override;

  /// Default impl of GetObjValues()
  ArrayRef<double> GetObjectiveValues() override
  { return std::vector<double>{ObjectiveValue()}; } 


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
protected:
  void ExportModel(const std::string& file);

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();

  void WindupGCGSolve();

  void ReportResults() override;
  void ReportGCGResults();

  void ReportGCGPool();
  /// Whether a pool solution passes sol:poolgap / sol:poolgapabs
  bool AcceptPoolSolution(double obj, double best) const;

  /// tech:optionnative: "name=value" or "name value" pairs
  void ApplyInlineParams();
  /// tech:outlev, tech:logfile: verbosity and log file
  void SetUpLogging();
  /// heu/cut/pre:settings, pre:solve
  void ApplyParamSettings();

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  int NodeCount() const;
  int SimplexIterations() const;
  int BarrierIterations() const;
  std::map<std::string, std::variant<int, double, std::string>>
      SolutionStats() override;

  std::pair<int, std::string> GetSolveResult() override;
  void AddGCGMessages();




private:
  /// Register a parameter as an AMPL option and remember
  /// the AMPL name for rejection messages.
  /// Integer bounds of a SCIP_Longint parameter register a 64-bit option.
  template <class... Args>
  void AddGcgOption(const char* names, const char* descr,
                    const char* key, Args&&... args) {
    RecordOptionName(key, PrimaryName(names));
    if constexpr (2 == sizeof...(Args)
                  && (std::is_same_v<std::decay_t<Args>, int> && ...)) {
      if (SCIP_PARAMTYPE_LONGINT == ParamType(key)) {
        AddSolverOption(names, descr, key, (SCIP_Longint)0, (SCIP_Longint)0);
        return;
      }
    }
    AddSolverOption(names, descr, key, std::forward<Args>(args)...);
  }
  static std::string PrimaryName(const char* names) {
    std::string s(names);
    return s.substr(0, s.find(' '));
  }

  void InputDecomposition();
  
  /// These options are stored in the class
  struct Options {
    std::string logFile_, paramRead_, paramWrite_;
    std::vector<std::string> inlineParams_;
    int heuristics_ = 0;
    int cuts_ = 0;
    int presolvings_ = 0;
    int presolve_ = 1;
    int outlev_ = 0;
    int poolLimit_ = 10;
    double poolGap_ = -1.0;     // < 0: no limit
    double poolGapAbs_ = -1.0;
  };
  Options storedOptions_;


};

}  // namespace mp

#endif  // MP_GCG_BACKEND_H_

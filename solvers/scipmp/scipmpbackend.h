#ifndef MP_SCIP_BACKEND_H_
#define MP_SCIP_BACKEND_H_

#include <string>
#include <utility>
#include <vector>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "scipmpcommon.h"

namespace mp {

class ScipBackend :
    public FlatBackend< MIPBackend<ScipBackend> >,
    public ScipCommon
{
  using BaseBackend = FlatBackend< MIPBackend<ScipBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  ScipBackend();
  ~ScipBackend();

  /// Prefix used for the <prefix>_options environment variable
  static const char* GetAMPLSolverName() { return "scip"; }

  /// AMPL driver name displayed in messages
  static const char* GetAMPLSolverLongName() { return "AMPL-SCIP"; }
  /// Solver name displayed in messages
  static const char* GetSolverName() { return "SCIP"; }
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
  ALLOW_STD_FEATURE(FIX_MODEL, true)
  bool IsFixedModelEmulated() const override { return true; }
  FixedLPResult SolveFixedLP(ArrayRef<double>) override;

  /**
 * MULTISOL support
 * No API, see ReportIntermediateSolution()
**/
  ALLOW_STD_FEATURE(MULTISOL, true)

  /**
  * Compute IIS (SCIP 10 IIS finder)
  **/
  ALLOW_STD_FEATURE(IIS, true)
  void ComputeIIS() override;
  IIS GetIIS() override;

  /**
  * Get/Set AMPL var/con statii
  **/
  ALLOW_STD_FEATURE(BASIS, false)
  // TODO If getting/setting a basis is supported, implement the 
  // accessor and the setter below
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;
  /**
  * General warm start, e.g.,
  * set primal/dual initial guesses for continuous case
  **/
  ALLOW_STD_FEATURE( WARMSTART, true )
  void AddPrimalDualStart(Solution sol) override;
  /**
  * MIP warm start
  **/
  // If MIP warm start is supported, implement the function below
  // to set a non-presolved starting solution
  ALLOW_STD_FEATURE(MIPSTART, true)
  void AddMIPStart(ArrayRef<double> x0,
									 ArrayRef<int> sparsity) override;

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
  * Get MIP Gap
  **/
  // return MIP gap
  // (adds option mip:return_gap)
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, true)
  double MIPGap() override;
  double MIPGapAbs() override;
  /**
  * Get MIP dual bound
  **/
  // return the best dual bound value
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
  { return std::vector<double>{HasFixedLP() ? FixedLP().objective : ObjectiveValue()}; }


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
protected:
  void ExportModel(const std::string& file);

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();

  void WindupSCIPSolve();

  void ReportResults() override;
  void ReportSCIPResults();

  void ReportSCIPPool();
  /// Whether a pool solution passes sol:poolgap / sol:poolgapabs
  bool AcceptPoolSolution(double obj, double best) const;
  /// sol:poolmode > 0: enumerate the feasible solutions of a
  /// pure-integer model with SCIP's solution counting, on a copy
  void EnumerateSolutions();
  bool AllVarsInteger() const;

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// tech:optionnative: "name=value" or "name value" pairs
  void ApplyInlineParams();
  /// tech:threads: concurrent solve and LP threads
  void ApplyThreads();
  /// tech:outlev, tech:logfile: verbosity and log file
  void SetUpLogging();
  /// alg:emphasis, heu/cut/pre:settings, pre:solve
  void ApplyParamSettings();
  bool UseConcurrentSolve() const;

  /// Solution attributes
  int NodeCount() const;
  int SimplexIterations() const;
  int BarrierIterations() const;
  std::map<std::string, std::variant<int, double, std::string>>
      SolutionStats() override;

  std::pair<int, std::string> GetSolveResult() override;
  void AddSCIPMessages();

  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  void VarStatii(ArrayRef<int>);
  void ConStatii(ArrayRef<int>);


private:
  /// Register a SCIP parameter as an AMPL option and remember
  /// the AMPL name for rejection messages.
  /// Integer bounds of a SCIP_Longint parameter register a 64-bit option.
  template <class... Args>
  void AddScipOption(const char* names, const char* descr,
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

  /// These options are stored in the class
  struct Options {
    std::string logFile_, paramRead_, paramWrite_;
    std::vector<std::string> inlineParams_;
    int threads_ = 0;
    int concurrent_ = -1;       // -1: with tech:threads > 1
    int heuristics_ = 0;
    int cuts_ = 0;
    int presolvings_ = 0;
    int presolve_ = 1;
    int emphasis_ = 0;
    int outlev_ = 0;
    int poolLimit_ = 10;
    int poolMode_ = 2;
    double poolGap_ = -1.0;     // < 0: no limit
    double poolGapAbs_ = -1.0;
  };
  Options storedOptions_;

  /// IIS statuses of the variables and linear constraints
  std::vector<int> iisVars_, iisCons_;

  /// Solutions enumerated by sol:poolmode > 0, best first
  struct CountedSolution {
    double obj;
    std::vector<double> x;
  };
  std::vector<CountedSolution> countedSols_;


};

}  // namespace mp

#endif  // MP_SCIP_BACKEND_H_

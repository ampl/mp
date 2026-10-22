 #ifndef MP_XPRESSMP_BACKEND_H_
#define MP_XPRESSMP_BACKEND_H_

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
#include <list>
#include <variant>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "xpresscommon.h"

#include "optprm.h"
#include "nlprm.h"

namespace mp {

class XpressmpBackend :
    public FlatBackend< MIPBackend<XpressmpBackend> >,
    public XpressmpCommon,
    public CompiledOptimizerOptions<XpressmpBackend>,
    public CompiledNonlinearOptions<XpressmpBackend>
{
  using BaseBackend = FlatBackend< MIPBackend<XpressmpBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  XpressmpBackend();
  ~XpressmpBackend();

  /// Name displayed in messages
  static const char* GetSolverName() { return "XPRESS"; }
  std::string GetSolverVersion();
  
  static const char* GetAMPLSolverName() { return "xpress"; }
  static const char* GetAMPLSolverLongName() { return "AMPL-XPRESSMP"; }
  static const char* GetBackendName();
  static const char* GetBackendLongName() { return nullptr; }

  void InputExtras() override;
  void InputXPRESSExtras();

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  // Use this section to declare and implement some standard features
  // that may or may not need additional functions. 
  USING_STD_FEATURES;

  /**
 * MULTIOBJ
**/
  ALLOW_STD_FEATURE(MULTIOBJ, true)
  ArrayRef<double> GetObjectiveValues() override;
  void ObjPriorities(ArrayRef<int>) override;
  void ObjWeights(ArrayRef<double>) override;
  void ObjAbsTol(ArrayRef<double>) override;
  void ObjRelTol(ArrayRef<double>) override;

  /**
 * MULTISOL support
 * No API, see ReportIntermediateSolution()
**/
  ALLOW_STD_FEATURE(MULTISOL, true)

  /**
  * Get/Set AMPL var/con statii
  **/
  ALLOW_STD_FEATURE(BASIS, true)
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;

  ALLOW_STD_FEATURE(WARMSTART, true)
  void AddPrimalDualStart(Solution sol0) override;
  /**
  * MIP warm start
  **/
  ALLOW_STD_FEATURE(MIPSTART, true)
	void AddMIPStart(ArrayRef<double> x0,
									 ArrayRef<int> sparsity) override;

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
  * Compute the IIS and obtain relevant values
  **/
  ALLOW_STD_FEATURE(IIS, true)
  /// Compute IIS
  void ComputeIIS() override;
  /// Retrieve IIS elements
  IIS GetIIS() override;

  /**
  * FixModel - duals, basis, and sensitivity for MIP.
  **/
  ALLOW_STD_FEATURE(FIX_MODEL, true)
  bool FixedModelReady() const override { return !fixed_model_failed_; }

  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override;

  /// Suported result format extensions
  ResultExtensionsSet NativeResultExtensions() const override
  { return {".sol", ".prt", ".iis", ".bas", ".bss"}; }

  ALLOW_STD_FEATURE(WRITE_SOLUTION, true)
  void DoWriteSolution(const std::string& name) override;

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
public:  // public for static polymorphism
  /// Solution pool and enumerator controls, see MspControl / MseControl
  using XpressmpCommon::GetSolverOption;
  using XpressmpCommon::SetSolverOption;
  void GetSolverOption(MspControl key, int& value) const;
  void SetSolverOption(MspControl key, int value);
  void GetSolverOption(MspControl key, double& value) const;
  void SetSolverOption(MspControl key, double value);
  void GetSolverOption(MseControl key, int& value) const;
  void SetSolverOption(MseControl key, int value);
  void GetSolverOption(MseControl key, double& value) const;
  void SetSolverOption(MseControl key, double value);


  // Overriding the initialization to make sure a problem is created
  // even for -v command line switch, as we need it to get the version number
  void Init(char** argv) override;

  /// See implementation on how to compile full params list
  void InitCustomOptions() override;
  bool IsMIP() const override {
    return BaseBackend::IsMIP()
        || NumSOSCons()|| NumIndicatorCons()
        || NumPWLs() || NumGenCons();
  }

  
protected:
  void OpenSolver();
  void CloseSolver();

  void FinishOptionParsing() override;

  int GetGlobalFlag(const SolverOption& ) const;
  void SetGlobalFlag(const SolverOption& , int );

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();

  /// Attach the solution pool to the problem and prepare the enumerator
  void CreateSolutionPoolEnvironment();
  /// Enumerate the sol:poollimit best solutions into the pool
  void EnumerateSolutions();
  /// After enumeration, the problem holds the last solution the
  /// enumerator visited. Keep the pool's best as the primary solution.
  void TakePoolBestAsPrimary();
  /// mip:basis in pool mode: fix the integer columns to the pool's best
  bool FixIntegersToPoolBest();
  /// Objective value for the MIP gap: of the last pass
  /// in a native multi-objective solve
  double MIPGapObjValue() const;
  /// Pool solution ids, best objective first
  std::vector<int> PoolSolutionIdsBestFirst();
  /// Is a pool solution within sol:poolgap / sol:poolgapabs of the best?
  bool IsWithinPoolGap(double obj, double best) const;
  /// Apply the tech:optionnative entries ("NAME=value", "NAME=?")
  void ApplyNativeOptions();
  void ApplyNativeOption(const std::string& name, const std::string& value);
  /// Load the basis file given by alg:readbasis
  void ReadBasisFile();

  std::string XPRESSSolveFlags();
  void WindupXPRESSMPSolve();

  void ReportResults() override;
  void ReportXPRESSMPResults();

  void ReportXPRESSMPPool();
  void DoXPRESSTune();

  void ConsiderXpressFixedModel();
  std::string DoXpressFixedModel();

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  int NodeCount() const;
  int SimplexIterations() const;
  int BarrierIterations() const;
  int CrossoverIterations() const;
  int SLPIterations() const;

  std::map<std::string, std::variant<int, double, std::string>>
      SolutionStats() override;

  std::pair<int, std::string> GetSolveResult() override;
  std::pair<int, std::string> GetSolveResult_LP();
  std::pair<int, std::string> GetSolveResult_MIP();
  std::pair<int, std::string> GetSolveResult_NLP();
  std::pair<int, std::string> GetSolveResult_General();

  void AddXPRESSMPMessages();

  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  std::vector<int> VarStatii(ArrayRef<int>);
  std::vector<int> ConStatii(ArrayRef<int>);

  ArrayRef<int> VarsIIS();
  pre::ValueMapInt ConsIIS();

  IISStatus IIS_VarToAMPL(char c);

  /* Xpress-MP callback in case the user wants some output from Optimizer */
  static int outlev_;
  
  static void xpdisplay(XPRSprob prob, void* data, const char* ch, int n, int msglvl);
  static int xp_mse_display(XPRSobject o, void* context, void* thread,
    const char* ch, int msglvl, int msgnumber);

  const std::string& tunebase() { return storedOptions_.tunebase_; }
  const std::string& tunermethodfileread() { return storedOptions_.tunermethodfileread_; }
  const std::string& tunermethodfilewrite() { return storedOptions_.tunermethodfilewrite_; }
  const std::string& tunename() { return storedOptions_.tunename_; }

private:
  int global_flag_ {1};
  XPRSprob model_fixed_ = nullptr;
  bool fixed_model_failed_ = false;
  std::vector<double> mip_primal_before_fixed_, mip_objectives_before_fixed_;
  /// Pool's best solution after enumeration (sol:poolmode)
  std::vector<double> pool_best_x_;
  double pool_best_obj_ = 0.0;

  /// These options are stored in the class
  struct Options {
    int nbest_ = 10;                // sol:poollimit
    int pooldualred_;
    int pooldupcol_;
    int nPoolMode_ = 2;             // sol:poolmode
    double poolGap_ = 1e100;        // sol:poolgap, relative
    double poolGapAbs_ = 1e100;     // sol:poolgapabs
    int nonconvex_ = -1;            // deprecated qp:nonconvex, -1: unset
    bool fBarrier_ = 0;
    bool fPrimal_ = 0;
    bool fDual_ = 0;
    bool fNetwork_ = 0;
    std::string tunebase_;
    std::string tunermethodfileread_;
    std::string tunermethodfilewrite_;
    std::string tunename_;
    std::string logFile_;
    std::string readBasis_;
    std::list<std::string> nativeOptions_;   // "NAME=value" entries
  };
  
  bool msgCallbackSet_ = false;
  Options storedOptions_;
  
  /// Solution pool and enumerator, created with the problem;
  /// the pool is attached to the problem only for multi-solution runs
  XPRSmipsolpool msp_;
  XPRSmipsolenum mse_;
};

}  // namespace mp

#endif  // MP_XPRESSMP_BACKEND_H_

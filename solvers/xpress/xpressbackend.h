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

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "xpresscommon.h"

namespace mp {

class XpressmpBackend :
    public FlatBackend< MIPBackend<XpressmpBackend> >,
    public XpressmpCommon
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

  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override;

  ALLOW_STD_FEATURE(WRITE_SOLUTION, true)
  void DoWriteSolution(const std::string& name) override;

  /////////////////////////// Model attributes /////////////////////////
  bool IsQCP() const override;
  
  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

  /// Solve, no model modification any more.
  /// Can report intermediate results via HandleFeasibleSolution() during this,
  /// otherwise in ReportResults()
  void Solve() override;



  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
public:  // public for static polymorphism

  // Overriding the initialization to make sure a problem is created
  // even for -v command line switch, as we need it to get the version number
  void Init(char** argv) override;

  void InitCustomOptions() override;
  bool IsMIP() const override {
    return BaseBackend::IsMIP() || NumSOSCons()|| NumIndicatorCons();
  }

  
protected:
  void OpenSolver();
  void CloseSolver();

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();

  void CreateSolutionPoolEnvironment();

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
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> ConvertXPRESSMPStatus();
  void AddXPRESSMPMessages();

  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  std::vector<int> VarStatii(ArrayRef<int>);
  std::vector<int> ConStatii(ArrayRef<int>);

  ArrayRef<int> VarsIIS();
  pre::ValueMapInt ConsIIS();

  /* Xpress-MP callback in case the user wants some output from Optimizer */
  static int outlev_;
  
  static void xpdisplay(XPRSprob prob, void* data, const char* ch, int n, int msglvl);
  static int xp_mse_display(XPRSobject o, void* context, void* thread,
    const char* ch, int msglvl, int msgnumber);

  const std::string& tunebase() { return storedOptions_.tunebase_; }
  const std::string& tunename() { return storedOptions_.tunename_; }

private:
  XPRSprob model_fixed_ = nullptr;

  /// These options are stored in the class
  struct Options {
    int nbest_ = 0;
    int pooldualred_;
    int pooldupcol_;
    int poollimit_ = 10;
    int nPoolMode_ = 0;
    int fBarrier_ = 0;
    int fPrimal_ = 0;
    int fDual_ = 0;
    int fNetwork_ = 0;
    std::string tunebase_;
    std::string tunename_;
    std::string logFile_;
  };
  
  bool msgCallbackSet_ = false;
  Options storedOptions_;
  
  // For multisol
  XPRSmipsolpool msp_;
  XPRSmipsolenum mse_;
};

}  // namespace mp

#endif  // MP_XPRESSMP_BACKEND_H_

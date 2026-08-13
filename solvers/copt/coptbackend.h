#ifndef MP_COPT_BACKEND_H_
#define MP_COPT_BACKEND_H_

#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "coptcommon.h"

namespace mp {

class CoptBackend :
    public FlatBackend< MIPBackend<CoptBackend> >,
    public CoptCommon
{
  using BaseBackend = FlatBackend< MIPBackend<CoptBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  CoptBackend();
  ~CoptBackend();

  /// Name displayed in messages
  static const char* GetSolverName() { return "COPT"; }
  std::string GetSolverVersion();
  
  static const char* GetAMPLSolverName() { return "copt"; }
  static const char* GetAMPLSolverLongName() { return "AMPL-COPT"; }
  static const char* GetBackendName();
  static const char* GetBackendLongName() { return nullptr; }

  /// Chance for the Backend to init solver environment, etc
  void InitOptionParsing() override;
  /// Chance to consider options immediately (open cloud, etc)
  void FinishOptionParsing() override;



  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  USING_STD_FEATURES;

  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override;

  /**
  * EXPORT SOLUTION
  **/
  ALLOW_STD_FEATURE(WRITE_SOLUTION, true)
  void DoWriteSolution(const std::string& name) override;

  /// Suported result format extensions
  ResultExtensionsSet NativeResultExtensions() const override
  { return {".sol", ".iis", ".bas"}; }

    /**
   * MULTIOBJ
  **/
  ALLOW_STD_FEATURE(MULTIOBJ, true)
  ArrayRef<double> GetObjectiveValues() override;
  void ObjPriorities(ArrayRef<int>) override;
  void ObjWeights(ArrayRef<double>) override;
  void ObjAbsTol(ArrayRef<double>) override;
  void ObjRelTol(ArrayRef<double>) override;
  void SetMultiobjOptions(BasicObjOptionSetter*) override;
  
  /**
  * MULTISOL support
  * No API, use ReportIntermediateSolution()
  **/
  ALLOW_STD_FEATURE(MULTISOL, true)

  /**
* Get/Set AMPL var/con statii
**/
  ALLOW_STD_FEATURE(BASIS, true)
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;

  /**
  * General warm start, e.g.,
  * set primal/dual initial guesses for continuous case
  **/
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
  /// Retrieve IIS. Elements correspond to IISStatus
  IIS GetIIS() override;

  /**
* Obtain inf/unbounded rays
**/
  ALLOW_STD_FEATURE(RAYS, true)
  ArrayRef<double> Ray() override;
  ArrayRef<double> DRay() override;

  /**
  * FeasRelax
  **/
  ALLOW_STD_FEATURE(FEAS_RELAX, true)

  /////////////////////////// Model attributes /////////////////////////
  bool IsMIP() const override;
  bool IsQCP() const override;
  

  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

  /// Solve, no model modification any more.
  /// Can report intermediate results via ReportIntermediateSolution() during this,
  /// otherwise in ReportResults()
  void Solve() override;

  /// This can modify the model
  void InputExtras() override;




  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
public:  // public for static polymorphism
  void InitCustomOptions() override;
  void OpenSolver();
  void CloseSolver();

protected:
  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();


  void InputCOPTExtras();

  void DoCOPTFeasRelax();
  void WindupCOPTSolve();

  void ReportResults() override;
  void ReportCOPTResults();

  void ReportCOPTPool();

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  int NodeCount() const;
  int SimplexIterations() const;
  int BarrierIterations() const;

  std::map<std::string, std::variant<int, double, std::string> >
	  SolutionStats() override;
  std::pair<int, std::string> GetSolveResult() override;
  void AddCOPTMessages();

  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  /// Should be in 1 call?
  void SetCoptBasis(ArrayRef<int> vs, ArrayRef<int> cs);

  ArrayRef<int> VarsIIS();
  pre::ValueMapInt ConsIIS();

  /// Native callback implementing PLATEAU_STOP / PLATEAU_STOP_BOUND:
  /// on COPT_CBCONTEXT_MIPSOL (new incumbent found) forwards the
  /// updated "BestObj" info to ReportIncumbentForPlateau(); on both
  /// COPT_CBCONTEXT_MIPSOL and COPT_CBCONTEXT_MIPNODE (periodic, keeps
  /// the plateau clock alive between incumbents) derives absgap/relgap
  /// from "BestObj"/"BestBnd" to report via ReportGapForPlateau().
  /// Terminates the solve via COPT_Interrupt() when either signals a
  /// plateau (COPT's callback has no per-objective partial-stop
  /// equivalent to Gurobi's GRBcbstoponemultiobj/Xpress's
  /// XPRS_STOP_NEXTOBJECTIVE, so this always stops the whole solve).
  static int COPT_CALL DoPlateauCallback(
      copt_prob* prob, void* cbdata, int cbctx, void* usrdata);


private:
  struct Options {
    std::string logFile_;
    int mempeak_ = 0;
    int numericalStats_ = 0;;

	int mempeak() const { return mempeak_; }
	int numericalStats() const { return numericalStats_; }
  };
  Options storedOptions_;


};

}  // namespace mp

#endif  // MP_COPT_BACKEND_H_

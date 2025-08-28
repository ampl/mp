#ifndef MP_CPLEX_BACKEND_H_
#define MP_CPLEX_BACKEND_H_

#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "cplexcommon.h"
#include <map>

namespace mp {

class CplexBackend :
    public FlatBackend< MIPBackend<CplexBackend> >,
    public CplexCommon
{
  using BaseBackend = FlatBackend< MIPBackend<CplexBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  CplexBackend();
  ~CplexBackend();

  /// Name displayed in messages
  static const char* GetSolverName() { return "CPLEX"; }
  std::string GetSolverVersion();
  static const char* GetAMPLSolverName() { return "cplex"; }
  static const char* GetAMPLSolverLongName() { return "AMPL-CPLEX"; }
  static const char* GetBackendName();
  static const char* GetBackendLongName() { return nullptr; }

  void InitOptionParsing() override;

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  USING_STD_FEATURES;

  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override;

  ALLOW_STD_FEATURE(WRITE_SOLUTION, true)
  void DoWriteSolution(const std::string& name) override;

  /**
 * MULTISOL support.
 * No API, use ReportIntermediateSolution()
**/
  ALLOW_STD_FEATURE(MULTISOL, true)
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
* Get/Set AMPL var/con statii
**/
  ALLOW_STD_FEATURE(BASIS, true)
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;
  ArrayRef<int> ConStatii();
  ArrayRef<int> VarStatii();
  void VarConStatii(ArrayRef<int>, ArrayRef<int>);
  /**
  * General warm start, e.g.,
  * set primal/dual initial guesses for continuous case
  **/
  ALLOW_STD_FEATURE(WARMSTART, true)
  void AddPrimalDualStart(Solution sol0) override;
  /**
  * Specifically, MIP warm start
  **/
  ALLOW_STD_FEATURE(MIPSTART, true)
  void AddMIPStart(ArrayRef<double> x0, ArrayRef<int> s0) override;

  /**
  * Get basis condition value (kappa)
  **/
  ALLOW_STD_FEATURE(KAPPA, true)
  double Kappa() override;

  /**
  * FixModel - duals, basis, and sensitivity for MIP.
  * No API to overload,
  * Impl should check need_fixed_MIP()
  **/
  ALLOW_STD_FEATURE(FIX_MODEL, true)
  void ConsiderCplexFixedModel();
  std::string DoCplexFixedModel();
  /**
  * Report sensitivity analysis suffixes (postsolved)
  **/
  ALLOW_STD_FEATURE(SENSITIVITY_ANALYSIS, true)
  SensRanges GetSensRanges() override;

  ALLOW_STD_FEATURE(IIS, true)
  void ComputeIIS() override;
  IIS GetIIS() override;
  ArrayRef<int> VarsIIS();
  pre::ValueMapInt ConsIIS();

  /**
  * Obtain inf/unbounded rays
  **/
  ALLOW_STD_FEATURE(RAYS, true)
  ArrayRef<double> Ray() override;
  ArrayRef<double> DRay() override;

  ALLOW_STD_FEATURE(FEAS_RELAX, true);
  void DoCplexFeasRelax();

  /////////////////////////// Model attributes /////////////////////////
  bool IsMIP() const override;
  bool IsQCP() const override;
  


  //////////////////////////// SOLVING ///////////////////////////////

    /// This can actually modify the model -- e.g., suffixes
  void InputExtras() override;
  void SetInterrupter(mp::Interrupter* inter) override;


  void InputCPLEXExtras();
  /// Solve, no model modification any more.
  /// Can report intermediate results via ReportIntermediateSolution() during this,
  /// otherwise in ReportResults()
  void Solve() override;

  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////


  /// For "obj:*:method" etc
public:
  using ObjNParamKey = std::pair< std::string, std::string >;
  template <class T>
  using ObjNParam = std::pair< ObjNParamKey, T >;
  enum CplexObjParams {
    OBJ_PRIORITY,
    OBJ_WEIGHT,
    OBJ_ABSTOL,
    OBJ_RELTOL,
    OBJ_NOTVALID
  };
  void CplexPlayObjNParams();
private:
  std::vector< ObjNParam<int> > objnparam_int_;
  std::vector< ObjNParam<double> > objnparam_dbl_;

  void CplexSetObjIntParam(const SolverOption& opt, int val);
  void CplexSetObjDblParam(const SolverOption& opt, double val);
  int CplexGetObjIntParam(const SolverOption& opt) const;
  double CplexGetObjDblParam(const SolverOption& opt) const;

 /// End "obj:*:method" etc

public:  // public for static polymorphism
  void InitCustomOptions() override;
  /// Chance to consider options immediately (for certain stored options)
  void FinishOptionParsing() override;
protected:
  void OpenSolver();
  void CloseSolver();

  

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();

  void WindupCPLEXSolve();

  void ReportResults() override;


  const SuffixDef<double> sufBestNodeObj = { "bestnode", suf::OBJ | suf::OUTONLY };
  const SuffixDef<double> sufBestNodeProb = { "bestbound", suf::PROBLEM | suf::OUTPUT };
  void ReportCPLEXResults();

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> GetSolveResult() override;
  void AddCPLEXMessages();
  void ReportCPLEXPool();

private:

  // returns {objlow, objhigh}
  std::pair<ArrayRef<double>, ArrayRef<double>> Sensobj() const;

  // returns { rhslow, rhshigh}
  std::pair<ArrayRef<double>, ArrayRef<double>> Sensrhs() const;

  // returns { {lb_low, lb_high}, {ub_low, ub_high} }
  std::pair<std::pair<ArrayRef<double>, ArrayRef<double>>, std::pair<ArrayRef<double>, ArrayRef<double>> > Sensbounds() const;


  bool silenceOutput_ = false;
  void setSilenceOutput(bool doSilence) {
    silenceOutput_ = doSilence;
  }
  
  void mymsgfunc(const char* msg, const char* )
  {
    if (!silenceOutput_)
      fmt::print("{}", msg);
  }
  void RedirectOutput();

  int hasSolution_ = -1;
  bool HasSolution();
  void ReadBendersSuffix();
  void setSolutionMethod();
  int original_model_type_ = -1;
  
  typedef struct
    CutInfo {
    const char* cutname;
    int cuttype;
  } CutInfo;
  CutInfo* CI;
  const static CutInfo Cut_Info[];

  /// These options are stored in the class
  struct Options {
    std::string paramRead_, paramWrite_;
    std::string exportFile_;
    std::string logFile_;
    std::string cpuMask_;
    std::string endBasis_;
    std::string mipStart_;

    std::string workDir_ = "";
    int outlev_ = 0;
    int noSolve_ = 0;
    int nPoolMode_=2;
    int populate_ = -1;
    int poolIntensity_ = -1;
    int algMethod_ = CPX_ALG_AUTOMATIC;
    int nodeMethod_ = CPX_ALG_AUTOMATIC;
    int cpxMethod_ = CPX_ALG_AUTOMATIC; // to store actual method
    int crossover_ = 0;
    int solutionType_ = 0;
    int bestnode_ = 0;
    int dropTol_ = 0;

    bool fBarrier_ = false;
    bool fPrimal_ = false;
    bool fDual_ = false;
    bool fNetwork_ = false;
    int netopt_ = 0;
    bool fSifting_ = false;
    bool fBenders_ = false;

    int cuts_ = -2;
    int cutstats_ = 0;

    int numcores_ = 0;
		bool dualprob_ = false;
		int predual_ = 0;
		bool primalprob_dummy_ = false;
	};
  Options storedOptions_;
  // to store IIS
  std::vector<int> iisColIndices, iisColValues,
                    iisRowIndices, iisRowValues;

};

}  // namespace mp

#endif  // MP_CPLEX_BACKEND_H_

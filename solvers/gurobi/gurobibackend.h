#ifndef MP_GUROBI_BACKEND_H_
#define MP_GUROBI_BACKEND_H_

#include <string>
#include <list>

#include "gurobicommon.h"
#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"

namespace mp {

class GurobiBackend :
    public FlatBackend< MIPBackend<GurobiBackend> >,
    public GurobiCommon
{
  using BaseBackend = FlatBackend< MIPBackend<GurobiBackend> >;


public:
  /// Constructor
  GurobiBackend();

  /// Destructor
  ~GurobiBackend();

  ////////////////////////////////////////////////////////////
  /// Standard and optional methods to provide or retrieve ///
  /// information to/from or manipulate the solver. Most  ////
  /// of them override placeholders from base classes.   /////
  ////////////////////////////////////////////////////////////

  /// Name displayed in messages
  static const char* GetSolverName() { return "Gurobi"; }
  /// Version
  static std::string GetSolverVersion();
  /// Use 'gurobi_options'
  static const char* GetAMPLSolverName() { return "gurobi"; }
  /// In long messages
  static const char* GetAMPLSolverLongName() { return "AMPL-Gurobi"; }

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  USING_STD_FEATURES;

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
   * MULTIOBJ
  **/
  ALLOW_STD_FEATURE( MULTIOBJ, true )
  ArrayRef<double> GetObjectiveValues() override;
  void ObjPriorities(ArrayRef<int>) override;
  void ObjWeights(ArrayRef<double>) override;
  void ObjAbsTol(ArrayRef<double>) override;
  void ObjRelTol(ArrayRef<double>) override;
  /**
   * MULTISOL support.
   * No API, use ReportIntermediateSolution()
  **/
  ALLOW_STD_FEATURE( MULTISOL, true )
  /**
  * Set lazy/user cut attributes.
  * Negative suffix values are "user cuts"
  * Check lazy_/user_cuts() to see which kinds are allowed
  **/
  ALLOW_STD_FEATURE( LAZY_USER_CUTS, true )
  void MarkLazyOrUserCuts(ArrayRef<int> ) override;
  /**
  * Get/Set AMPL var/con statii
  **/
  ALLOW_STD_FEATURE( BASIS, true )
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis ) override;
  /**
  * General warm start, e.g.,
  * set primal/dual initial guesses for continuous case
  **/
  ALLOW_STD_FEATURE( WARMSTART, true )
  void AddPrimalDualStart(Solution sol0) override;
  /**
  * Specifically, MIP warm start
  **/
  ALLOW_STD_FEATURE( MIPSTART, true )
  void AddMIPStart(
      ArrayRef<double> x0, ArrayRef<int> s0) override;
  /**
  * Obtain inf/unbounded rays
  **/
  ALLOW_STD_FEATURE( RAYS, true )
  ArrayRef<double> Ray() override;
  ArrayRef<double> DRay() override;
  /**
  * Compute the IIS and obtain relevant values
  **/
  ALLOW_STD_FEATURE( IIS, true )
  /// Compute IIS
  void ComputeIIS() override;
  /// Retrieve IIS. Elements correspond to IISStatus
  IIS GetIIS() override;
  /**
  * Get MIP Gap
  **/
  ALLOW_STD_FEATURE( RETURN_MIP_GAP, true )
  /// Should return +inf if not available
  double MIPGap() override;
  /// Should return +inf if not available
  double MIPGapAbs() override;
  /**
  * Get MIP dual bound
  **/
  ALLOW_STD_FEATURE( RETURN_BEST_DUAL_BOUND, true )
  double BestDualBound() override;
  /**
  * Set branch and bound priorities
  **/
  ALLOW_STD_FEATURE( VAR_PRIORITIES, true )
  void VarPriorities(ArrayRef<int> ) override;
  /**
  * Get basis condition value (kappa)
  **/
  ALLOW_STD_FEATURE( KAPPA, true)
  double Kappa() override;
  /**
  * FeasRelax
  **/
  ALLOW_STD_FEATURE( FEAS_RELAX, true)
  /**
  * Report sensitivity analysis suffixes (postsolved)
  **/
  ALLOW_STD_FEATURE( SENSITIVITY_ANALYSIS, true )
  SensRanges GetSensRanges() override;
  /**
  * FixModel - duals, basis, and sensitivity for MIP.
  * No API to overload,
  * Impl should check need_fixed_MIP()
  **/
  ALLOW_STD_FEATURE( FIX_MODEL, true )


  ///////////////////// Model attributes /////////////////////

  /// IsMIP()
  bool IsMIP() const override;
  /// IsQP()
  bool IsQP() const override;
  /// IsQCP()
  bool IsQCP() const override;


  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////// OPTION ACCESSORS ///////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  /// GurobiBackend driver options
  void InitCustomOptions() override;

  /// Chance for the Backend to init solver environment, etc
  void InitOptionParsing() override;
  /// Chance to consider options immediately (open cloud, etc)
  void FinishOptionParsing() override;


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////// SOLVING ACCESSORS ///////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  /// This can actually modify the model -- e.g., suffixes
  void InputExtras() override;

  /// Set the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

  /// Solve, no model modification any more.
  /// Can report intermediate results via ReportIntermediateSolution() during this,
  /// otherwise in ReportResults()
  void Solve() override;


  ///////////////////////////////////////////////////////////////////////////////
  //////////////////// Implementation's internals ///////////////////////////////
  //////////////////// Gurobi methods should contain 'Gurobi' or similar ////////
  //////////////////// to avoid name clashes with the base classes //////////////
  ///////////////////////////////////////////////////////////////////////////////
protected:
  void OpenGurobi();
  void OpenGurobiModel();
  void CloseGurobi();

  void OpenGurobiComputeServer();
  void OpenGurobiCloud();

  void InputGurobiExtras();
  /// .funcpieces etc.
  /// Override the corresponding Gurobi options.
  void InputGurobiFuncApproxParams();
  /// .iisforce
  void InputGurobiIISForceParams();
  /// Functional nonlinear cons attributes need a filter
  /// @return the index of the first used value (-1 if none)
  template <class T>
  int GurobiSetFuncConAttributes(
      const char* attr, const std::vector<T>& vals);
  void DoGurobiFeasRelax();

  void ExportFile(GRBmodel* lp, const std::string& file);

  void PrepareGurobiSolve();
  void SetPartitionValues();

  void DoGurobiTune();

  void WindupGurobiSolve();

  void ReportResults() override;
  void ReportGurobiResults();

  std::pair<int, std::string> ConvertGurobiStatus() const;
  void AddGurobiMessage();

  void ReportGurobiPool();
  /// Creates and solves, marks model_fixed to be used for duals/basis/sens
  void ConsiderGurobiFixedModel();
  /// Return error message if any
  std::string DoGurobiFixedModel();
  ArrayRef<double> CurrentGrbPoolPrimalSolution();
  double CurrentGrbPoolObjectiveValue() const;

  double ObjectiveValue() const;
  mp::ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;

  std::vector<double> GurobiDualSolution_LP();
  std::vector<double> GurobiDualSolution_QCP();

  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  void VarStatii(ArrayRef<int> );
  void ConStatii(ArrayRef<int> );

  ArrayRef<double> Senslbhi() const ;
  ArrayRef<double> Senslblo() const ;
  ArrayRef<double> Sensobjhi() const ;
  ArrayRef<double> Sensobjlo() const ;
  ArrayRef<double> Sensrhshi() const ;
  ArrayRef<double> Sensrhslo() const ;
  ArrayRef<double> Sensubhi() const ;
  ArrayRef<double> Sensublo() const ;

  ArrayRef<int> VarsIIS();
  pre::ValueMapInt ConsIIS();

  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;


protected:  //  ///////////////// Utilities ///////////////////
  /// Swap 0 with -1.
  /// Usefule to convert option / suffix values, where we interpret
  /// them differently, as for IISForce and FuncNonLinear.
    void Swap_0_vs_minus1(std::vector<int>& arr);

private:
  /// Fixed model to produce duals for MIP if requested
  GRBmodel *model_fixed_ = nullptr;

  /// These options are stored in the class as variables
  /// for direct access
  struct Options {
    std::string paramRead_, paramWrite_, exportPresolvedFile_, logFile_;
    std::list<std::string> inlineParams_;

    int nMIPStart_=1;
    int nPoolMode_=2;

    int nFixedMethod_=-2;

    int fFuncPieceSuf_ = 1;
    int fFuncNonlinear_ = 0;
    int nIISForce_ = 1;

    std::string cloudid_, cloudkey_, cloudpool_;
    int cloudpriority_ = 0;

    std::string servers_, server_password_, server_group_, server_router_;
    int server_priority_=0, server_insecure_=0;
    double server_timeout_=10.0;

    std::string tunebase_;
  } storedOptions_;


protected:  //////////// Option accessors ////////////////
  int Gurobi_mipstart() const { return storedOptions_.nMIPStart_; }

  const std::string& paramfile_read() const { return storedOptions_.paramRead_; }
  const std::string& paramfile_write() const { return storedOptions_.paramWrite_; }

  int funcpiecesuf() const { return storedOptions_.fFuncPieceSuf_; }
  int funcnonlinear() const { return storedOptions_.fFuncNonlinear_; }
  int iisforce() const { return storedOptions_.nIISForce_; }

  const std::string& cloudid() const { return storedOptions_.cloudid_; }
  const std::string& cloudkey() const { return storedOptions_.cloudkey_; }
  const std::string& cloudpool() const { return storedOptions_.cloudpool_; }
  int cloudpriority() const { return storedOptions_.cloudpriority_; }

  const std::string& servers() const { return storedOptions_.servers_; }
  const std::string& server_password() const { return storedOptions_.server_password_; }
  const std::string& server_group() const { return storedOptions_.server_group_; }
  const std::string& server_router() const { return storedOptions_.server_router_; }
  int server_priority() const { return storedOptions_.server_priority_; }
  int server_insecure() const { return storedOptions_.server_insecure_; }
  int server_timeout() const { return storedOptions_.server_timeout_; }

  const std::string& tunebase() const { return storedOptions_.tunebase_; }

  const std::string& logfile() const { return storedOptions_.logFile_; }

  const std::string& exportPresolvedFile() const { return storedOptions_.exportPresolvedFile_; }


private: /////////// Suffixes ///////////
  const SuffixDef<int> sufHintPri = { "hintpri", suf::VAR | suf::INPUT };


  /// For "obj:*:method" etc
public:
  using ObjNParamKey = std::pair< std::string, std::string >;
  template <class T>
  using ObjNParam = std::pair< ObjNParamKey, T >;


private:
  std::vector< ObjNParam<int> > objnparam_int_;
  std::vector< ObjNParam<double> > objnparam_dbl_;


protected:
  void GrbSetObjIntParam(const SolverOption& opt, int val);
  void GrbSetObjDblParam(const SolverOption& opt, double val);
  int GrbGetObjIntParam(const SolverOption& opt) const;
  double GrbGetObjDblParam(const SolverOption& opt) const;

  void GrbPlayObjNParams();
};

} // namespace mp

#endif  // MP_GUROBI_BACKEND_H_

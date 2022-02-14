#ifndef MP_GUROBI_BACKEND_H_
#define MP_GUROBI_BACKEND_H_

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

#include "gurobicommon.h"
#include "mp/backend_mip.h"
#include "mp/backend_with_pre.h"

namespace mp {

class GurobiBackend :
    public MIPBackend<GurobiBackend>,
    public BackendWithPresolver,
    public GurobiCommon
{
  using BaseBackend = MIPBackend<GurobiBackend>;

public:
  GurobiBackend();
  ~GurobiBackend();

  ////////////////////////////////////////////////////////////
  //////////////////// PART 1. Accessor API //////////////////
  /// Standard and optional methods to provide or retrieve ///
  /// information to/from or manipulate the solver. Most  ////
  /// of them override placeholders from base classes.   /////
  ////////////////////////////////////////////////////////////

  /// Name displayed in messages
  static const char* GetSolverName() { return "x-Gurobi"; }
  /// Version
  static std::string GetSolverVersion();
  /// Reuse gurobi_options:
  static const char* GetAMPLSolverName() { return "gurobi"; }
  /// In long messages
  static const char* GetAMPLSolverLongName() { return "AMPL-Gurobi"; }

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  USING_STD_FEATURES;
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
   * MULTISOL support
   * No API, use ReportIntermediateSolution()
  **/
  ALLOW_STD_FEATURE( MULTISOL, true )
  /**
  * Set lazy/user cut attributes
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
  void AddMIPStart(ArrayRef<double> x0) override;
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
  double MIPGap() override;
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
  * Report sensitivity analysis suffixes
  **/
  ALLOW_STD_FEATURE( SENSITIVITY_ANALYSIS, true )
  ArrayRef<double> Senslbhi() const override;
  ArrayRef<double> Senslblo() const override;
  ArrayRef<double> Sensobjhi() const override;
  ArrayRef<double> Sensobjlo() const override;
  ArrayRef<double> Sensrhshi() const override;
  ArrayRef<double> Sensrhslo() const override;
  ArrayRef<double> Sensubhi() const override;
  ArrayRef<double> Sensublo() const override;
  /**
  * FixModel - duals, basis, and sensitivity for MIP
  * No API to overload,
  * Impl should check need_fixed_MIP()
  **/
  ALLOW_STD_FEATURE( FIX_MODEL, true )


  ///////////////////// Model attributes /////////////////////
  bool IsMIP() const override;
  bool IsQP() const override;
  bool IsQCP() const override;


  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////// OPTION ACCESSORS ///////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  /// Gurobi-specific options
  void InitCustomOptions() override;

  /// Chance for the Backend to init solver environment, etc
  void InitOptionParsing() override;
  /// Chance to consider options immediately (open cloud, etc)
  void FinishOptionParsing() override;


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////// SOLVING ACCESSORS ///////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  void SetInterrupter(mp::Interrupter* inter) override;

  void SolveAndReportIntermediateResults() override;

  /// Various solution attribute getters.
  Solution GetSolution() override;


  ///////////////////////////////////////////////////////////////////////////////
  //////////////////// PART 2. Implementation's internals ///////////////////////
  //////////////////// Gurobi methods should include name Gurobi or similar /////
  //////////////////// to avoid name clashes with the base classes //////////////
  ///////////////////////////////////////////////////////////////////////////////
protected:
  void OpenGurobi();
  void OpenGurobiModel();
  void CloseGurobi();

  void OpenGurobiComputeServer();
  void OpenGurobiCloud();


  void ExportModel(const std::string& file);

  void PrepareGurobiSolve();
  void DoGurobiFeasRelax();
  void SetPartitionValues();

  void DoGurobiTune();

  void WindupGurobiSolve();
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
  std::vector<double> PrimalSolution();
  pre::ValueMapDbl DualSolution();

  std::vector<double> GurobiDualSolution_LP();
  std::vector<double> GurobiDualSolution_QCP();

  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  void VarStatii(ArrayRef<int> );
  void ConStatii(ArrayRef<int> );

  ArrayRef<int> VarsIIS();
  pre::ValueMapInt ConsIIS();

  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;


private:
  GRBmodel *model_fixed_ = nullptr;

  /// These options are stored in the class as variables
  /// for direct access
  struct Options {
    std::string exportFile_, paramRead_, paramWrite_;

    int nMIPStart_=1;
    int nPoolMode_=2;

    int nFixedMethod_=-2;

    std::string cloudid_, cloudkey_, cloudpool_;
    int cloudpriority_;

    std::string servers_, server_password_, server_group_, server_router_;
    int server_priority_=0, server_insecure_=0;
    double server_timeout_=-1.0;

    std::string tunebase_;
  } storedOptions_;


protected:  //////////// Option accessors ////////////////
  int Gurobi_mipstart() const { return storedOptions_.nMIPStart_; }

  const std::string& paramfile_read() const { return storedOptions_.paramRead_; }
  const std::string& paramfile_write() const { return storedOptions_.paramWrite_; }

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

private: /////////// Suffixes ///////////
  const SuffixDef<int> sufHintPri = { "hintpri", suf::VAR | suf::INPUT };


  /// For "obj:*:method" etc
  /// Should they be handled in the Converter?
public:
  using ObjNParamKey = std::pair< std::string, std::string >;
  template <class T>
  using ObjNParam = std::pair< ObjNParamKey, T >;
private:
  std::vector< ObjNParam<int> > objnparam_int_;
  std::vector< ObjNParam<double> > objnparam_dbl_;
protected:
  /// Assume opt has the * info
  void GrbSetObjIntParam(const SolverOption& opt, int val);
  void GrbSetObjDblParam(const SolverOption& opt, double val);
  int GrbGetObjIntParam(const SolverOption& opt) const;
  double GrbGetObjDblParam(const SolverOption& opt) const;

  void GrbPlayObjNParams();
};

} // namespace mp

#endif  // MP_GUROBI_BACKEND_H_

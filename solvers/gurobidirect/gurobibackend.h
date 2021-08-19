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

extern "C" {
  #include "gurobi_c.h"
}

#if __clang__
# pragma clang diagnostic pop
#elif _MSC_VER
# pragma warning(pop)
#endif

#include <string>

#include "mp/convert/MIP/backend.h"
#include "mp/convert/std_constr.h"

namespace mp {

class GurobiBackend : public MIPBackend<GurobiBackend>
{
  using BaseBackend = MIPBackend<GurobiBackend>;

public:
  GurobiBackend();
  ~GurobiBackend();

  /// Metadata
  static const char* GetSolverName() { return "Gurobi"; }
  static std::string GetSolverVersion();
  static const char* GetSolverInvocationName();
  static const char* GetAMPLSolverLongName() { return nullptr; }
  static const char* GetBackendName();
  static const char* GetBackendLongName() { return nullptr; }

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  ////////////////////////////////////////////////////////////
  /// Default switch, needed if not all std features mentioned
  DEFAULT_STD_FEATURES_TO( false )
  /**
   * MULTIOBJ
  **/
  ALLOW_STD_FEATURE( MULTIOBJ, true )
  std::vector<double> ObjectiveValues() const;
  void ObjPriorities(ArrayRef<int>);
  void ObjWeights(ArrayRef<double>);
  void ObjAbsTol(ArrayRef<double>);
  void ObjRelTol(ArrayRef<double>);
  /**
   * MULTISOL support
  **/
  ALLOW_STD_FEATURE( MULTISOL, true )
  void StartPoolSolutions();
  bool SelectNextPoolSolution();
  void EndPoolSolutions();
  std::vector<double> CurrentPoolPrimalSolution();
  double CurrentPoolObjectiveValue() const;
  /**
  * Get/Set AMPL var/con statii
  **/
  ALLOW_STD_FEATURE( BASIS, true )
  std::vector<int> VarStatii();
  std::vector<int> ConStatii();
  void VarStatii(ArrayRef<int> );
  void ConStatii(ArrayRef<int> );
  /**
  * General warm start, e.g.,
  * set primal/dual initial guesses for continuous case
  **/
  ALLOW_STD_FEATURE( WARMSTART, true )
  void InputPrimalDualStart(ArrayRef<double> x0,
                         ArrayRef<double> pi0);
  /**
  * Specifically, MIP warm start
  **/
  ALLOW_STD_FEATURE( MIPSTART, true )
  void AddMIPStart(ArrayRef<double> x0);
  /**
  * Obtain inf/unbounded rays
  **/
  ALLOW_STD_FEATURE( RAYS, true )
  std::vector<double> Ray();
  std::vector<double> DRay();
  /**
  * Compute the IIS and obtain relevant values
  **/
  ALLOW_STD_FEATURE( IIS, true )
  void ComputeIIS();
  /// Elements correspond to IISStatus
  std::vector<int> VarsIIS();
  std::vector<int> ConsIIS();
  /**
  * Get MIP Gap
  **/
  ALLOW_STD_FEATURE( RETURN_MIP_GAP, true )
  double MIPGap() const;
  /**
  * Get MIP dual bound
  **/
  ALLOW_STD_FEATURE( RETURN_BEST_DUAL_BOUND, true )
  double BestDualBound() const;
  /**
  * Set branch and bound priorities
  **/
  ALLOW_STD_FEATURE( VAR_PRIORITIES, true )
  void VarPriorities(ArrayRef<int> );
  /**
  * Get basis condition value (kappa)
  **/
  ALLOW_STD_FEATURE(KAPPA, true)
  double Kappa() const;
  /**
  * FeasRelax
  **/
  ALLOW_STD_FEATURE(FEAS_RELAX, true)


  /// [[ Prototype an incremental interface ]]
  void InitProblemModificationPhase();
  /// Chance to update the model
  void FinishProblemModificationPhase();

  void AddVariable(Variable var);
  void SetLinearObjective( int iobj, const LinearObjective& lo );
  void SetQuadraticObjective( int iobj, const QuadraticObjective& qo );

  //////////////////////////// GENERAL CONSTRAINTS ////////////////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BaseBackend)

  /// TODO Attributes (lazy/user cut, etc)
  ACCEPT_CONSTRAINT(LinearConstraint, Recommended)
  void AddConstraint(const LinearConstraint& lc);
  ACCEPT_CONSTRAINT(QuadraticConstraint, Recommended)
  void AddConstraint(const QuadraticConstraint& qc);
  ACCEPT_CONSTRAINT(MaximumConstraint, AcceptedButNotRecommended)
  void AddConstraint(const MaximumConstraint& mc);
  ACCEPT_CONSTRAINT(MinimumConstraint, AcceptedButNotRecommended)
  void AddConstraint(const MinimumConstraint& mc);
  ACCEPT_CONSTRAINT(AbsConstraint, AcceptedButNotRecommended)
  void AddConstraint(const AbsConstraint& absc);
  ACCEPT_CONSTRAINT(ConjunctionConstraint, AcceptedButNotRecommended)
  void AddConstraint(const ConjunctionConstraint& cc);
  ACCEPT_CONSTRAINT(DisjunctionConstraint, AcceptedButNotRecommended)
  void AddConstraint(const DisjunctionConstraint& mc);
  /// Enabling built-in indicator for infinite bounds,
  /// but not recommended otherwise --- may be slow
  ACCEPT_CONSTRAINT(IndicatorConstraintLinLE, AcceptedButNotRecommended)
  void AddConstraint(const IndicatorConstraintLinLE& mc);

  /// Nonlinear
  ACCEPT_CONSTRAINT(ExpConstraint, Recommended)
  void AddConstraint(const ExpConstraint& cc);
  ACCEPT_CONSTRAINT(ExpAConstraint, Recommended)
  void AddConstraint(const ExpAConstraint& cc);
  ACCEPT_CONSTRAINT(LogConstraint, Recommended)
  void AddConstraint(const LogConstraint& cc);
  ACCEPT_CONSTRAINT(LogAConstraint, Recommended)
  void AddConstraint(const LogAConstraint& cc);
  ACCEPT_CONSTRAINT(PowConstraint, Recommended)
  void AddConstraint(const PowConstraint& cc);
  ACCEPT_CONSTRAINT(SinConstraint, Recommended)
  void AddConstraint(const SinConstraint& cc);
  ACCEPT_CONSTRAINT(CosConstraint, Recommended)
  void AddConstraint(const CosConstraint& cc);
  ACCEPT_CONSTRAINT(TanConstraint, Recommended)
  void AddConstraint(const TanConstraint& cc);
  ACCEPT_CONSTRAINT(PLConstraint, Recommended)
  void AddConstraint(const PLConstraint& cc);


  ///////////////////// Model attributes /////////////////////
  bool IsMIP() const;
  bool IsQCP() const;

  int NumberOfConstraints() const;
  int NumberOfVariables() const;
  int NumberOfObjectives() const;
  int ModelSense() const;

  void ExportModel(const std::string& file);


  //////////////////////////// SOLVING ///////////////////////////////
  void SetInterrupter(mp::Interrupter* inter);
  void SolveAndReportIntermediateResults();
  std::string ConvertSolutionStatus(
      const mp::Interrupter &interrupter, int &solve_code);

  /// Various solution attribute getters.
  /// Return empty vectors if not available
  std::vector<double> PrimalSolution();
  std::vector<double> DualSolution();
  double ObjectiveValue() const;

  double NodeCount() const;
  double Niterations() const;


  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
public:
  void OpenSolver();
  void CloseSolver();
  void InitCustomOptions();

  static double Infinity() { return GRB_INFINITY; }
  static double MinusInfinity() { return -GRB_INFINITY; }

  /// REMEMBER Gurobi does not update attributes before calling optimize() etc
  /// Scalar attributes. If (flag), set *flag <-> success,
  /// otherwise fail on error
  int GrbGetIntAttr(const char* attr_id, bool *flag=nullptr) const;
  double GrbGetDblAttr(const char* attr_id, bool *flag=nullptr) const;
  /// Vector attributes. Return empty vector on failure
  std::vector<int> GrbGetIntAttrArray(const char* attr_id,
    std::size_t size, std::size_t offset=0) const;
  std::vector<double> GrbGetDblAttrArray(const char* attr_id,
    std::size_t size, std::size_t offset=0) const;

  /// Set attributes.
  /// Return false on failure
  void GrbSetIntAttr(const char* attr_id, int val);
  void GrbSetDblAttr(const char* attr_id, double val);
  ///  Silently ignore empty vector arguments.
  void GrbSetIntAttrArray(const char* attr_id,
                               ArrayRef<int> values, std::size_t start=0);
  void GrbSetDblAttrArray(const char* attr_id,
                               ArrayRef<double> values, std::size_t start=0);
  ///  Silently ignore empty vector arguments.
  void GrbSetIntAttrList(const char* attr_id,
                         const std::vector<int>& idx, const std::vector<int>& val);
  void GrbSetDblAttrList(const char* attr_id,
                         const std::vector<int>& idx, const std::vector<double>& val);

  /// First objective's sense
  void SetMainObjSense(obj::Type s);
  obj::Type GetMainObjSense() const;

protected:
  void PrepareGurobiSolve();
  void DoFeasRelax();

private:
  GRBenv   *env   = NULL;
  GRBmodel *model = NULL;

  /// The sense of the main objective
  obj::Type main_obj_sense_;

private:
  /// These options are stored in the class as variables
  /// for direct access
  struct Options {
    std::string exportFile_;

    int nMIPStart_=1;
    int nPoolMode_=2;
  } storedOptions_;


protected:  //////////// Option accessors ////////////////
  int mipstart() const { return storedOptions_.nMIPStart_; }


private:  /////////// Some working variables ////////////
  int iPoolSolution = -2;          // for SelectNextPoolSolution()


private: /////////// Suffixes ///////////
  const SuffixDef<int> sufHintPri = { "hintpri", suf::VAR | suf::INPUT };


public:  //////////// Wrappers for Get/SetSolverOption()
  int GrbGetIntParam(const char* key) const;
  double GrbGetDblParam(const char* key) const;
  std::string GrbGetStrParam(const char* key) const;
  void GrbSetIntParam(const char* key, int value);
  void GrbSetDblParam(const char* key, double value);
  void GrbSetStrParam(const char* key, const std::string& value);

  /// These methods access Gurobi options. Used by AddSolverOption()
  void GetSolverOption(const char* key, int& value) const;
  void SetSolverOption(const char* key, int value);
  void GetSolverOption(const char* key, double& value) const;
  void SetSolverOption(const char* key, double value);
  void GetSolverOption(const char* key, std::string& value) const;
  void SetSolverOption(const char* key, const std::string& value);

};

} // namespace mp

#endif  // MP_GUROBI_BACKEND_H_

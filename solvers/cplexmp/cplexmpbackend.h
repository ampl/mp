#ifndef MP_CPLEX_BACKEND_H_
#define MP_CPLEX_BACKEND_H_

#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "cplexmpcommon.h"

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
  void ReportCPLEXResults();

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> ConvertCPLEXStatus();
  void AddCPLEXMessages();
  void ReportCPLEXPool();

private:
 
  /// These options are stored in the class
  struct Options {
    std::string exportFile_;
    std::string logFile_;
    int outlev_ = 0;;
    int nPoolMode_=2;
    int populate_ = -1;
    int poolIntensity_ = -1;

  };
  Options storedOptions_;
  // to store IIS
  std::vector<int> iisColIndices, iisColValues,
                    iisRowIndices, iisRowValues;

};

}  // namespace mp

#endif  // MP_CPLEX_BACKEND_H_

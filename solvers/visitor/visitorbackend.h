#ifndef MP_VISITOR_BACKEND_H_
#define MP_VISITOR_BACKEND_H_

#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "visitorcommon.h"

namespace mp {

class VisitorBackend :
    public FlatBackend< MIPBackend<VisitorBackend> >,
    public VisitorCommon
{
  using BaseBackend = FlatBackend< MIPBackend<VisitorBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  /// Construct
  VisitorBackend();
  /// Destruct
  ~VisitorBackend();

  /// Prefix used for the <prefix>_options environment variable
  static const char* GetAMPLSolverName() { return "visitor"; }

  /// AMPL driver name displayed in messages
  static const char* GetAMPLSolverLongName() { return "AMPL-VISITOR"; }
  /// Solver name displayed in messages
  static const char* GetSolverName() { return "x-VISITOR"; }
  /// Version displayed with -v
  std::string GetSolverVersion();
  /// External libraries displayed with -v
  std::string set_external_libs() override { return ""; };
  
  /// Name for diagnostic messages
  static const char* GetBackendName();
  /// "long name", rarely used
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
  * EXPORT PROBLEM
  **/
  ALLOW_STD_FEATURE(WRITE_PROBLEM, true)
  void DoWriteProblem(const std::string& name) override { }

  /**
  * EXPORT SOLUTION
  **/
  ALLOW_STD_FEATURE(WRITE_SOLUTION, true)
  void DoWriteSolution(const std::string& name) override { }


  /**
   * MULTISOL support
   * No API, see ReportIntermediateSolution()
   **/
  ALLOW_STD_FEATURE(MULTISOL, true)

  /**
   * Get/Set AMPL var/con statii
   **/
  ALLOW_STD_FEATURE(BASIS, true)
  /// TODO If getting/setting a basis is supported, implement the
  /// accessor and the setter below.
  /// Should return empty basis if not available
  /// (e.g., not an LP.)
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;

  /**
  * MIP warm start
  **/
  /// TODO If MIP warm start is supported, implement the function below
  /// to set a non-presolved starting solution
  ALLOW_STD_FEATURE(MIPSTART, true)
  void AddMIPStart(ArrayRef<double> x0,
                   ArrayRef<int> sparsity) override;


  /**
  * Get MIP Gap
  **/
  // TODO Implement to return MIP gap
  // (adds option mip:return_gap)
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, true)
  double MIPGap() override;
  double MIPGapAbs() override;

  /**
  * Get MIP dual bound
  **/
  // TODO Implement to return the best dual bound value
  // (adds option mip:bestbound)
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

  /////////////////////////// Model attributes /////////////////////////

  /// Reimplement if the solver gives more information
  /// than just the number of non-fixed integer variables
  /// (e.g., the solver might consider if it has PL expressions.)
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
  void OpenSolver();
  void CloseSolver();

  double ObjectiveValue() const;

	/// PrimalSolution() for flat backends.
	/// @return empty vector if no primal solution.
	ArrayRef<double> PrimalSolution() override;
	/// DualSolution() for flat backends.
	/// @return empty map if no dual solution.
	pre::ValueMapDbl DualSolution() override;
	/// Dual solution for the LP part only.
	/// @return empty vector if none.
  ArrayRef<double> DualSolution_LP();

  void WindupVISITORSolve();

  void ReportResults() override;
  void ReportVISITORResults();

  void ReportVISITORPool();

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> ConvertVISITORStatus();
  void AddVISITORMessages();

  /// Return basis.
  /// @return empty vector if not available.
  ArrayRef<int> VarStatii();
  /// @return empty vector if not available.
  ArrayRef<int> ConStatii();

  /// Set var basis statuses.
  void VarStatii(ArrayRef<int>);
  /// Set con basis statuses.
  void ConStatii(ArrayRef<int>);

  ArrayRef<int> VarsIIS();
  pre::ValueMapInt ConsIIS();


private:
  /// These options are stored in the class
  struct Options {
    std::string option_example_;
    bool flag_option_ = false;
    std::vector<double> list_option_;
  };
  Options storedOptions_;


protected:
  const std::string& get_example_option() const
  { return storedOptions_.option_example_; }
  const std::vector<double>& get_list_option() const
  { return storedOptions_.list_option_; }

};

}  // namespace mp

#endif  // MP_VISITOR_BACKEND_H_

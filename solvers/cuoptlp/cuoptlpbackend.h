#ifndef MP_CUOPTLP_BACKEND_H_
#define MP_CUOPTLP_BACKEND_H_

#include <vector>
#include <list>
#include <string>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "cuoptlpcommon.h"

namespace mp {

class CuoptlpBackend :
    public FlatBackend< MIPBackend<CuoptlpBackend> >,
    public CuoptlpCommon
{
  using BaseBackend = FlatBackend< MIPBackend<CuoptlpBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
  /// Construct
  CuoptlpBackend();
  /// Destruct
  ~CuoptlpBackend();

  /// Prefix used for the <prefix>_options environment variable
  static const char* GetAMPLSolverName() { return "cuoptlp"; }

  /// AMPL driver name displayed in messages
  static const char* GetAMPLSolverLongName() { return "AMPL-CUOPT"; }
  /// Solver name displayed in messages
  static const char* GetSolverName() { return "CUOPT"; }
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
  // For a full list of features possible,
  // grep "STD_FEATURE".
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
   * MULTISOL support.
   *  If (need_multiple_solutions()),
   *  call ReportIntermediateSolution() during solve or after.
   **/
  //ALLOW_STD_FEATURE(MULTISOL, true)

  /**
   * Get/Set AMPL var/con statii
   **/
  //ALLOW_STD_FEATURE(BASIS, true)
  /// TODO If getting/setting a basis is supported, implement the
  /// accessor and the setter below.
  /// Should return empty basis if not available
  /// (e.g., not an LP.)
  //SolutionBasis GetBasis() override;
  //void SetBasis(SolutionBasis) override;

  /**
  * MIP warm start
  **/
  /// TODO If MIP warm start is supported, implement the function below
  /// to set a non-presolved starting solution
  //ALLOW_STD_FEATURE(MIPSTART, true)
  //void AddMIPStart(ArrayRef<double> x0,
  //                 ArrayRef<int> sparsity) override;


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
  //ALLOW_STD_FEATURE(IIS, true)
  /// Compute IIS.
  /// This method can fail (MP_RAISE)
  /// if it discovers a different problem status.
  //void ComputeIIS() override;
  /// Retrieve IIS elements.
  /// Only called if the status was confirmed Infeasible.
  //  IIS GetIIS() override;

  /////////////////////////// Model attributes /////////////////////////

  /// Reimplement if the solver gives more information
  /// than just the number of non-fixed integer variables
  /// (e.g., the solver might consider if it has PL expressions.)
  bool IsMIP() const override;
  //bool IsQCP() const override;

  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;


public:  // public for static polymorphism
  /// Solve, to be overloaded by the solver.
  /// No model modification any more.
  /// @note If using STD_FEATURE( MULTISOL ),
  /// can report intermediate results
  /// via ReportIntermediateSolution() during this
  /// (check if (need_multiple_solutions())),
  /// otherwise afterwards.
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
  /// Dual solution for the QP part only.
  /// @return empty vector if none.
  ArrayRef<double> DualSolution_QP();

  void WindupCUOPTLPSolve();

  void ReportResults() override;
  void ReportCUOPTLPResults();

  void ReportCUOPTLPPool();

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> GetSolveResult() override;
  void AddCUOPTLPMessages();

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

	  std::string paramread_, paramwrite_;
    std::list<std::string> inlineparams_;

    int verbosity_= 1;
  };
  Options storedOptions_;
  void printModelStats();

protected:
  const std::string& get_example_option() const
  { return storedOptions_.option_example_; }
  const std::vector<double>& get_list_option() const
  { return storedOptions_.list_option_; }

};

}  // namespace mp

#endif  // MP_CUOPTLP_BACKEND_H_

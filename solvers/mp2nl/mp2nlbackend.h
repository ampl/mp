#ifndef MP2NL_BACKEND_H_
#define MP2NL_BACKEND_H_

#include <string>
#include <map>

#include "mp/backend-mip.h"
#include "mp/flat/backend_flat.h"
#include "mp2nlcommon.h"

namespace mp {

class MP2NLBackend :
		public FlatBackend< MIPBackend<MP2NLBackend> >,
		public MP2NLCommon
{
	using BaseBackend = FlatBackend< MIPBackend<MP2NLBackend> >;

  //////////////////// [[ The public interface ]] //////////////////////
public:
	MP2NLBackend();
	~MP2NLBackend();

  /// Prefix used for the <prefix>_options environment variable
  static const char* GetAMPLSolverName() { return "mp2nl"; }

  /// AMPL driver name displayed in messages
  static const char* GetAMPLSolverLongName() { return "AMPL-MP2NL"; }
  /// Solver name displayed in messages
  static const char* GetSolverName() { return "MP2NL"; }
  /// Version displayed with -v
  std::string GetSolverVersion();
  /// External libraries displayed with -v
  std::string set_external_libs() override;
  
  /// Name for diagnostic messages
  static const char* GetBackendName();
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
 * MULTISOL support
 * No API, see ReportIntermediateSolution()
**/
  ALLOW_STD_FEATURE(MULTISOL, false)

  /**
  * Get/Set AMPL var/con statii
  **/
  ALLOW_STD_FEATURE(BASIS, false)
  // TODO If getting/setting a basis is supported, implement the 
  // accessor and the setter below
  SolutionBasis GetBasis() override;
  void SetBasis(SolutionBasis) override;
  /**
  * General warm start, e.g.,
  * set primal/dual initial guesses for continuous case
  **/
  ALLOW_STD_FEATURE( WARMSTART, true )
  void AddPrimalDualStart(Solution sol) override;
  /**
  * MIP warm start
  **/
  // If MIP warm start is supported, implement the function below
  // to set a non-presolved starting solution
  ALLOW_STD_FEATURE(MIPSTART, true)
  void AddMIPStart(ArrayRef<double> x0,
									 ArrayRef<int> sparsity) override;

  /**
  * EXPORT PROBLEM
  **/
  ALLOW_STD_FEATURE(WRITE_PROBLEM, false)
  void DoWriteProblem(const std::string& name) override;


 /**
  * Get MIP Gap
  **/
  // return MIP gap
  // (adds option mip:return_gap)
  ALLOW_STD_FEATURE(RETURN_MIP_GAP, false)
  double MIPGap() override;
  double MIPGapAbs() override;
  /**
  * Get MIP dual bound
  **/
  // return the best dual bound value
  // (adds option mip:bestbound)
  ALLOW_STD_FEATURE(RETURN_BEST_DUAL_BOUND, false)
  double BestDualBound() override;

  /////////////////////////// Model attributes /////////////////////////
  
  //////////////////////////// SOLVING ///////////////////////////////

  /// Note the interrupt notifier
  void SetInterrupter(mp::Interrupter* inter) override;

public:  // public for static polymorphism
  /// Solve, no model modification any more (such as feasrelax).
  /// Can report intermediate results via ReportIntermediateSolution() during this,
  /// otherwise/finally via ReportResults()
  void Solve() override;

  /// SOL does not return an obj value
  ArrayRef<double> GetObjectiveValues() override
  { return {}; }

  double Infinity() const { return AMPLInf(); }

  /// Get ini guess
  ArrayRef< std::pair<int, double> > GetInitialGuess()
  { return x0_; }

  /// Get ini dual guess
  ArrayRef<double> GetInitialDualGuess()
  { return y0_; }

	/// Typedef map of default solver configs
	using ConfigMap = std::map<std::string, std::string>;

  //////////////////// [[ Implementation details ]] //////////////////////
  ///////////////////////////////////////////////////////////////////////////////
protected:

  double ObjectiveValue() const;

  /// Solution values. The vectors are emptied if not available
  ArrayRef<double> PrimalSolution() override;
  pre::ValueMapDbl DualSolution() override;
  ArrayRef<double> DualSolution_LP();

  void WindupMP2NLSolve();

  void ReportResults() override;
  void ReportMP2NLResults();

  void ReportSuffixes() override;
  void ReportModelSuffix(const MP2NLModelSuffix& modelsuf);
  void ReportMP2NLPool();

  std::vector<double> getPoolSolution(int i);
  double getPoolObjective(int i);

  /// Solution attributes
  double NodeCount() const;
  double SimplexIterations() const;
  int BarrierIterations() const;

  std::pair<int, std::string> GetSolveResult() override;
  void AddMP2NLMessages();

  ArrayRef<int> VarStatii();
  ArrayRef<int> ConStatii();
  void VarStatii(ArrayRef<int>);
  void ConStatii(ArrayRef<int>);

	int GetOutlev(const SolverOption& ) const
	{ return storedOptions_.outlev_; }
	void SetOutlev(const SolverOption& , int ol);

	/// Only the last one if several solver=... options
	std::string GetSolver(const SolverOption& ) const
	{ return storedOptions_.solver_; }
	void SetSolver(const SolverOption&, fmt::StringRef val);
	std::string ExtractSolverConfigName(fmt::StringRef solver);

	bool GetDummyFlagOption(const SolverOption& ) const
	{ return false; }
	/// Just prints possible configs
	void SetConfigPrintFlag(const SolverOption& , bool );
	/// Only the last one if several config=... options
	std::string GetConfig(const SolverOption& ) const
	{ return storedOptions_.config_; }
	void SetConfig(const SolverOption&, fmt::StringRef val);

	bool DoSetConfig(fmt::StringRef cfg);

private:
  /// These options are stored in the class
  struct Options {
    std::string solver_,
        solver_options_;
		std::string config_, config_attempted_;
    std::string logFile_,
        nlstub_;
    int outlev_ = 0;

    double tilim_ = 1e20;
  };
  Options storedOptions_;
	static ConfigMap config_map_;

  std::unique_ptr<MP2NLSolverQueryCallbacks> p_qc_;

  std::vector< std::pair<int, double> > x0_;   // presolved ini guess
  std::vector<double> y0_;  // presolved dual ini guess, dense
};

}  // namespace mp

#endif  // MP2NL_BACKEND_H_

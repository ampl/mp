#ifndef MODEL_MANAGER_BASE_H
#define MODEL_MANAGER_BASE_H

/**
  * This is an abstract interface to a Model Manager
  * which provides model IO, solution handling and suffixes, etc.
  */

#include <functional>

#include "mp/ampls-ccallbacks.h"

#include "mp/arrayref.h"
#include "mp/solver-base.h"
#include "mp/suffix.h"
#include "mp/obj-option-base.h"

namespace mp {

/// Abstract Model Manager.
/// Standardizes the following tasks:
/// - Input the model;
/// - Access user-provided solution/suffixes;
/// - Report solutions/suffixes.
class BasicModelManager {
public:
  virtual ~BasicModelManager() { }

  /// Setup Model Manager's solver options
  virtual void InitOptions() = 0;

  /// Read NL model
  virtual void ReadNLModel(const std::string& nl_filename,
                           const std::string& filename_no_ext,
                           Checker_AMPLS_ModeltTraits ,
                            std::function<void(char*)> after_header) = 0;

  /// User-provided primal solution
  virtual ArrayRef<double> InitialValues() = 0;
  /// User-provided primal solution: sparsity
  virtual ArrayRef<int> InitialValuesSparsity() = 0;
  /// User-provided dual solution
  virtual ArrayRef<double> InitialDualValues() = 0;
  /// User-provided dual solution: sparsity
  virtual ArrayRef<int> InitialDualValuesSparsity() = 0;

  /// Get suffix names
  virtual std::set<std::string> GetSuffixNames() = 0;

  virtual   SuffixSet& GetSuffixes(suf::Kind kind) = 0;

  /// Read integer suffix
  virtual ArrayRef<int> ReadSuffix(const SuffixDef<int>& suf) = 0;
  /// Read double suffix
  /// @param fint: if not NULL,
  ///   is set to 1 iff the suffix was integer.
  virtual ArrayRef<double> ReadSuffix(
      const SuffixDef<double>& suf, int *fint=nullptr) = 0;

  /// Report integer suffix
  virtual void ReportSuffix(const SuffixDef<int>& suf,
                            ArrayRef<int> values) = 0;
  /// Report double suffix
  virtual void ReportSuffix(const SuffixDef<double>& suf,
                            ArrayRef<double> values) = 0;

  /// Length of a suffix vector of given kind
  virtual size_t GetSuffixSize(int kind) = 0;

  /// Get obj option setter
  virtual BasicObjOptionSetter* GetObjOptionSetter() = 0;

  /// Get solution file name
  virtual void SetSolutionFileName(const std::string& fileName) = 0;

  /// Report final solution
  virtual void HandleSolution(int, fmt::CStringRef,
                              const double *, const double *,
                              double) = 0;
  /// Report intermediate solution
  virtual void HandleFeasibleSolution(
                              int solve_code, fmt::CStringRef,
                              const double *, const double *,
                              double) = 0;

  /// Need and successfully prepared the next solve iteration?
  /// @param get_stt: solution status getter.
  ///   If called, then before get_sol.
  /// @param get_sol: solution getter (for postsolved solution.)
  virtual bool PrepareSolveIteration(
      std::function<sol::Status(void)> get_stt, std::function<Solution(void)> get_sol) = 0;

  /// Objective weights in the 'legacy' format of the obj:multi:weight option
  virtual ArrayRef<double> GetObjWeightsAdapted() = 0;

  /// Integrality flags of the variables in the original instance.
  /// Used for solution rounding
  virtual const std::vector<bool>& IsVarInt() const = 0;

  /// Has unfixed int vars?
  /// This is about the solver-facing instance.
  virtual bool HasUnfixedIntVars() const = 0;
};

} // namespace mp

#endif // MODEL_MANAGER_BASE_H

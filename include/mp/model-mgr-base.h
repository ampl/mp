#ifndef MODEL_MANAGER_BASE_H
#define MODEL_MANAGER_BASE_H

/**
  * This is an abstract interface to a Model Manager
  * which provides model IO, solution handling and suffixes, etc.
  */

#include <functional>

#include "mp/arrayref.h"
#include "mp/solver-base.h"
#include "mp/suffix.h"

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

  /// IO file basename, ideal to know this before option parsing
  virtual void SetBasename(const std::string& filename_no_ext) = 0;

  /// Read NL model
  virtual void ReadNLModel(const std::string& nl_filename,
                           const std::string& filename_no_ext,
                           void (*cb_checkmodel)(size_t, size_t, size_t)) = 0;

  /// User-provided primal solution
  virtual ArrayRef<double> InitialValues() = 0;
  /// User-provided dual solution
  virtual ArrayRef<double> InitialDualValues() = 0;

  /// Read integer suffix
  virtual ArrayRef<int> ReadSuffix(const SuffixDef<int>& suf) = 0;
  /// Read double suffix
  virtual ArrayRef<double> ReadSuffix(const SuffixDef<double>& suf) = 0;

  /// Report integer suffix
  virtual void ReportSuffix(const SuffixDef<int>& suf,
                            ArrayRef<int> values) = 0;
  /// Report double suffix
  virtual void ReportSuffix(const SuffixDef<double>& suf,
                            ArrayRef<double> values) = 0;

  /// Length of a suffix vector of given kind
  virtual size_t GetSuffixSize(int kind) = 0;

  /// Report final solution
  virtual void HandleSolution(int, fmt::CStringRef,
                              const double *, const double *,
                              double) = 0;
  /// Report intermediate solution
  virtual void HandleFeasibleSolution(fmt::CStringRef,
                              const double *, const double *,
                              double) = 0;

  /// Integrality flags of the variables in the original instance.
  /// Used for solution rounding
  virtual const std::vector<bool>& IsVarInt() const = 0;
};

} // namespace mp

#endif // MODEL_MANAGER_BASE_H

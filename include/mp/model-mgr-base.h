#ifndef MODEL_MANAGER_BASE_H
#define MODEL_MANAGER_BASE_H

/**
  * This is an abstract interface to a Model Manager
  * which provides model IO, solution handling and suffixes, etc.
  */

#include "mp/arrayref.h"
#include "mp/solver-base.h"
#include "mp/suffix.h"

namespace mp {

/// Input the model
/// Access/provide the original model solution/suffixes
class BasicModelManager {
public:
  virtual ~BasicModelManager() { }

  /// Setup Model Manager's solver options
  virtual void InitOptions() = 0;

  /// IO file basename, ideal to know this before option parsing
  virtual void SetBasename(const std::string& filename_no_ext) = 0;

  /// Read NL model
  virtual void ReadNLFileAndUpdate(const std::string& nl_filename,
                                   const std::string& filename_no_ext) = 0;

  virtual ArrayRef<double> InitialValues() = 0;
  virtual ArrayRef<double> InitialDualValues() = 0;

  virtual ArrayRef<int> ReadSuffix(const SuffixDef<int>& suf) = 0;
  virtual ArrayRef<double> ReadSuffix(const SuffixDef<double>& suf) = 0;

  virtual void ReportSuffix(const SuffixDef<int>& suf,
                            ArrayRef<int> values) = 0;
  virtual void ReportSuffix(const SuffixDef<double>& suf,
                            ArrayRef<double> values) = 0;

  /// Length of a suffix vector
  virtual size_t GetSuffixSize(int kind) = 0;

  virtual void HandleSolution(int, fmt::CStringRef,
                              const double *, const double *,
                              double) = 0;
  virtual void HandleFeasibleSolution(fmt::CStringRef,
                              const double *, const double *,
                              double) = 0;

  /// Integrality flags of the variables in the original instance
  /// We use this for rounding
  virtual const std::vector<bool>& IsVarInt() const = 0;
};

} // namespace mp

#endif // MODEL_MANAGER_BASE_H

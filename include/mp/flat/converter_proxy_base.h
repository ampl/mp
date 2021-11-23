#ifndef CONVERTER_PROXY_BASE_H
#define CONVERTER_PROXY_BASE_H

/**
  * This is an abstract interface to a 'query class'
  * which provides solution handling and suffixes, etc.
  * To be used by a backend
  */

#include "mp/solver.h"
#include "mp/flat/model.h"
#include "mp/flat/suffix.h"

namespace mp {

class ConverterQuery {
public:
  virtual ~ConverterQuery() { }

  using Model = BasicModel<>;

  virtual ArrayRef<double> InitialValues() = 0;
  virtual ArrayRef<double> InitialDualValues() = 0;

  virtual ArrayRef<int> ReadSuffix(const SuffixDef<int>& suf) = 0;
  virtual ArrayRef<double> ReadSuffix(const SuffixDef<double>& suf) = 0;

  virtual void ReportSuffix(const SuffixDef<int>& suf,
                            ArrayRef<int> values) = 0;
  virtual void ReportSuffix(const SuffixDef<double>& suf,
                            ArrayRef<double> values) = 0;

  virtual void HandleSolution(int, fmt::CStringRef,
                              const double *, const double *,
                              double) = 0;
  virtual void HandleFeasibleSolution(fmt::CStringRef,
                              const double *, const double *,
                              double) = 0;

  virtual const std::vector<bool>& IsVarInt() const = 0;

};

} // namespace mp

#endif // CONVERTER_PROXY_BASE_H

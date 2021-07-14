#ifndef CONVERTER_QUERY_H
#define CONVERTER_QUERY_H

/**
  * This is an abstract interface to a 'query class'
  * which provides solution handling and suffixes, etc.
  * To be used by a backend
  */

#include "mp/solver.h"
#include "mp/convert/model.h"
#include "mp/convert/suffix.h"

namespace mp {

class ConverterQuery {
public:
  virtual ~ConverterQuery() { }

  using Model = BasicModel<>;


  virtual ArrayRef<int> ReadSuffix(const SuffixDef<int>& suf) = 0;

  virtual ArrayRef<double> ReadSuffix(const SuffixDef<double>& suf) = 0;

  virtual void ReportSuffix(const SuffixDef<int>& suf,
                            const std::vector<int>& values) = 0;

  virtual void ReportSuffix(const SuffixDef<double>& suf,
                      const std::vector<double>& values) = 0;

  virtual void HandleSolution(int, fmt::CStringRef,
                              const double *, const double *,
                              double) = 0;

  virtual void HandleFeasibleSolution(fmt::CStringRef,
                              const double *, const double *,
                              double) = 0;

};

} // namespace mp

#endif // CONVERTER_QUERY_H

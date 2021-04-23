#ifndef CONVERTER_QUERY_H
#define CONVERTER_QUERY_H

/**
  * This is an abstract interface to a 'query class'
  * which provides solution handling and suffixes, etc.
  * To be used by a backend
  */

#include "mp/solver.h"
#include "mp/convert/model.h"

namespace mp {

class ConverterQuery {
public:
  virtual ~ConverterQuery() { }

  using Model = BasicModel<>;

  /// Outputs an integer suffix.
  /// name: Suffix name that may not be null-terminated.
  virtual void DeclareAndReportIntSuffix(fmt::StringRef name, int kind,
                      const std::vector<int>& values) = 0;

  virtual void HandleSolution(int, fmt::CStringRef,
      const double *, const double *, double) = 0;

};

} // namespace mp

#endif // CONVERTER_QUERY_H

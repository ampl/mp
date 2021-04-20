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

  using IntSuffixHandler = Model::SuffixHandler<int>;

  /// Adds an integer suffix.
  /// name: Suffix name that may not be null-terminated.
  /// TODO put values right there (1. vector, 2. sparse vector)
  virtual IntSuffixHandler AddIntSuffix(fmt::StringRef name, int kind, int=0) = 0;

  virtual void HandleSolution(int, fmt::CStringRef,
      const double *, const double *, double) = 0;

};

} // namespace mp

#endif // CONVERTER_QUERY_H

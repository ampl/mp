#ifndef CONVERTER_PROXY_BASE_H
#define CONVERTER_PROXY_BASE_H

#include <vector>

#include "mp/arrayref.h"

namespace mp {

/**
  * A basic 'query class' for FlatConverters,
  * providing pre- / postsolve
  * To be used by a backend
  */
class FlatConverterProxy {
public:
  virtual ~FlatConverterProxy() = default;

  /// Convenience methods
  /// E.g., Gurobi reports duals etc separately for linear and QCP constraints

  /// For reporting constraint suffixes
  virtual size_t NumValuedOrigConstr() const = 0;

  /////////////////////////////////////////////////////////////////////////
  /// PRESOLVE ///
  /////////////////////////////////////////////////////////////////////////

  /// From original NL model's suffix or duals
  virtual std::vector<double>
  ExtractLinConValues(ArrayRef<double> allval) = 0;
  virtual std::vector<int>
  ExtractLinConValues(ArrayRef<int> allval) = 0;

  virtual std::vector<double>
  ExtractQCValues(ArrayRef<double> allval) = 0;

  /////////////////////////////////////////////////////////////////////////
  /// POSTSOLVE ///
  /////////////////////////////////////////////////////////////////////////

  /// To original NL model's indexing
  virtual std::vector<double> MakeConstrValuesFromLPAndQCP(
        ArrayRef<double> pi, ArrayRef<double> qcpi) = 0;
  virtual std::vector<int> MakeConstrValuesFromLPAndQCP(
        ArrayRef<int> pi, ArrayRef<int> qcpi) = 0;
};

} // namespace mp

#endif // CONVERTER_PROXY_BASE_H

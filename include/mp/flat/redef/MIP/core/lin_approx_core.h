#ifndef LIN_APPROX_CORE_H
#define LIN_APPROX_CORE_H

/**
 * Converter-independent linear approximation routines
 */

#include "mp/flat/constr_std.h"

namespace mp {

/// Function graph domain
struct FuncGraphDomain {
  /// Bounds for x, y
  double lbx, ubx, lby, uby;
  /// Intersect the domain with another rectangle
  void intersect(const FuncGraphDomain& grDom);
};

/// External parameters
/// for piecewise-linear approximation
struct PLApproxParams {
  /// The domain for approximation
  FuncGraphDomain grDom;
  /// Error upper bound (relative or absolute)
  double ubErr = 1e-5;
};


/// Do the approximation by calling
/// the function-specific approximator
/// @param laPrm: in-out parameter, e.g., bounds can be tightened
/// @return the PL function
template <class FuncCon>
PLPoints PLApproximate(const FuncCon& con,
                             PLApproxParams& laPrm);

}  // namespace mp

#endif // LIN_APPROX_CORE_H

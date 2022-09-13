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

/// Range
struct Range {
  /// Construct
  Range(double l=-1e100, double u=1e100) : lb(l), ub(u) { }
  /// The bounds
  double lb, ub;
};

/// Input and output
/// for piecewise-linear approximation
struct PLApproxParams {
  /////////// INPUT DATA //////////////
  /// Graph domain for approximation
  FuncGraphDomain grDom;
  /// Error upper bound (relative or absolute)
  double ubErr = 1e-5;

  /////////// OUTPUT: RESULT OF APPROXIMATION /////////////
  FuncGraphDomain grDomOut;    // can be tighter
  PLPoints plPoints; // the pl
  bool fUsePeriod;  // whether the approximation is periodic
  Range periodic_factor_range;
  Range period_remainder_range;
};


/// Do the approximation by calling
/// the function-specific approximator
/// @param laPrm: in-out parameter
template <class FuncCon>
void PLApproximate(const FuncCon& con,
                   PLApproxParams& laPrm);

}  // namespace mp

#endif // LIN_APPROX_CORE_H

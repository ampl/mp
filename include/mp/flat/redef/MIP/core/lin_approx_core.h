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
  /// Error upper bound
  /// (relative outside of +-1, absolute inside)
  double ubErr = 1e-5;

  /////////// OUTPUT: RESULT OF APPROXIMATION /////////////
  FuncGraphDomain grDomOut;    // can be tighter than grDom
  PLPoints plPoints; // the pl function
  bool fUsePeriod{false};  // whether the approximation is periodic
  double periodLength;        // the length of the period interval
  Range periodicFactorRange; // range for the n of x = n*periodLength + rmd
  Range periodRemainderRange; // range for the rmd
};


/// Do the approximation by calling
/// the function-specific approximator
/// @param laPrm: in-out parameter
template <class FuncCon>
void PLApproximate(const FuncCon& con,
                   PLApproxParams& laPrm);

}  // namespace mp

#endif // LIN_APPROX_CORE_H

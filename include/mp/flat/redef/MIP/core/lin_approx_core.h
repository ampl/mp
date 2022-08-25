#ifndef LIN_APPROX_CORE_H
#define LIN_APPROX_CORE_H

/**
 * Converter-independent linear approximation routines
 */

#include "mp/flat/constr_std.h"

namespace mp {

/// Function graph domain
struct FuncGraphDomain {
  double lbx, ubx, lby, uby;
  void intersect(const FuncGraphDomain& grDom);
};

/// External parameters
/// for piecewise-linear approximation
struct LinApproxParams {
  FuncGraphDomain grDom;
  double ubErrAbs = 1e-2;
};


/// Do the approximation by calling
/// the function-specific approximator
/// @param laPrm: in-out parameter, e.g., bounds can be tightened
/// @return the PL function
template <class FuncCon>
PLPoints PLApproximate(const FuncCon& con,
                             LinApproxParams& laPrm);

}  // namespace mp

#endif // LIN_APPROX_CORE_H

#ifndef LIN_APPROX_CORE_H
#define LIN_APPROX_CORE_H

/**
 * Converter-independent linear approximation routines
 */

#include "mp/flat/constr_std.h"

namespace mp {

/// Parameters for a generic piecewise-linear approximation
struct LinApproxParams {
  double lbx, ubx, lby, uby;
};

/// Do the approximation
/// @return the PL function
template <class FuncCon>
PLPoints PLApproximate(const FuncCon& con,
                             LinApproxParams laParams) {
  return {};
}

}  // namespace mp

#endif // LIN_APPROX_CORE_H

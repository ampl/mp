#ifndef CONSTR_NONLINEAR_H
#define CONSTR_NONLINEAR_H

/**
 * Nonlinearity traits of flat constraint types
 * (option cvt:pre:boundsbest, values 2-4).
 */

#include <type_traits>

#include "mp/flat/constr_std.h"
#include "mp/flat/nl_expr/constr_nl.h"

namespace mp {

/// Is \a Con a smooth nonlinear functional constraint,
/// i.e., its result is exp, log, pow, ... of its arguments?
template <class Con>
constexpr bool IsSmoothNLFuncCon() {
  return std::is_same_v<Con, DivConstraint>
      || std::is_same_v<Con, MulConstraint>
      || std::is_same_v<Con, ExpConstraint>
      || std::is_same_v<Con, ExpAConstraint>
      || std::is_same_v<Con, LogConstraint>
      || std::is_same_v<Con, LogAConstraint>
      || std::is_same_v<Con, PowConstExpConstraint>
      || std::is_same_v<Con, SignpowConstExpConstraint>
      || std::is_same_v<Con, PowConstraint>
      || std::is_same_v<Con, LogisticConstraint>
      || std::is_same_v<Con, SinConstraint>
      || std::is_same_v<Con, CosConstraint>
      || std::is_same_v<Con, TanConstraint>
      || std::is_same_v<Con, AsinConstraint>
      || std::is_same_v<Con, AcosConstraint>
      || std::is_same_v<Con, AtanConstraint>
      || std::is_same_v<Con, SinhConstraint>
      || std::is_same_v<Con, CoshConstraint>
      || std::is_same_v<Con, TanhConstraint>
      || std::is_same_v<Con, AsinhConstraint>
      || std::is_same_v<Con, AcoshConstraint>
      || std::is_same_v<Con, AtanhConstraint>
      || std::is_same_v<Con, CallConstraint>;
}

/// Is \a Con nonlinear for the solver: a smooth nonlinear or
/// quadratic function, a quadratic row, a cone, or an expression-
/// interface constraint? Logical constraints are not counted.
template <class Con>
constexpr bool IsNonlinearConType() {
  return IsSmoothNLFuncCon<Con>()
      || std::is_same_v<Con, QuadraticFunctionalConstraint>
      || std::is_same_v<Con, QuadConRange>
      || std::is_same_v<Con, QuadConLT>
      || std::is_same_v<Con, QuadConLE>
      || std::is_same_v<Con, QuadConEQ>
      || std::is_same_v<Con, QuadConGE>
      || std::is_same_v<Con, QuadConGT>
      || std::is_same_v<Con, QuadraticConeConstraint>
      || std::is_same_v<Con, RotatedQuadraticConeConstraint>
      || std::is_same_v<Con, ExponentialConeConstraint>
      || std::is_same_v<Con, PowerConeConstraint>
      || std::is_same_v<Con, GeometricConeConstraint>
      || std::is_same_v<Con, NLConstraint>
      || std::is_same_v<Con, NLAssignEQ>
      || std::is_same_v<Con, NLAssignLE>
      || std::is_same_v<Con, NLAssignGE>
      || std::is_same_v<Con, NLComplementarity>;
}

} // namespace mp

#endif // CONSTR_NONLINEAR_H

#ifndef COND_LITERAL_H
#define COND_LITERAL_H

/**
 * Literals of indicator constraints
 * for conversions of conditional comparisons.
 */

#include "mp/flat/constr_std.h"

namespace mp {

/// The binary literal which triggers an implication:
/// (var == val) ==> constraint.
///
/// A conditional comparison's own conversion uses
/// its result variable: (res == 1) ==> body (cmp) rhs
/// in positive context, (res == 0) ==> negation in negative context.
/// Absorbing a disjunction (redef/MIP/disj2ind.h) uses
/// the other disjunct instead.
struct CondLiteral {
  int var_;    ///< the binary variable
  int val_;    ///< the triggering value, 0 or 1
};

/// Can the literal still trigger, i.e., is its variable
/// not fixed to the opposite value?
template <class MC>
bool CanTrigger(const MC& mc, CondLiteral lit) {
  return !mc.is_fixed(lit.var_)
      || lit.val_ == mc.fixed_value(lit.var_);
}

/// Does the literal surely trigger (variable fixed to the value)?
template <class MC>
bool SurelyTriggers(const MC& mc, CondLiteral lit) {
  return mc.is_fixed(lit.var_)
      && lit.val_ == mc.fixed_value(lit.var_);
}

/// Forbid the literal to trigger: fix its variable
/// to the opposite value.
template <class MC>
void ForbidTrigger(MC& mc, CondLiteral lit) {
  mc.NarrowVarBounds(lit.var_, 1 - lit.val_, 1 - lit.val_);
}

/// Linear expression which is 1 iff the literal does NOT trigger:
/// var (if val==0) or 1-var (if val==1).
inline AffineExpr NonTriggerExpr(CondLiteral lit) {
  if (lit.val_)
    return { {{-1.0}, {lit.var_}}, 1.0 };
  return { {{1.0}, {lit.var_}}, 0.0 };
}

} // namespace mp

#endif // COND_LITERAL_H

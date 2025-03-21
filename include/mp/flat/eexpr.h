#ifndef EEXPR_H
#define EEXPR_H

#include <utility>

#include "mp/flat/expr_quadratic.h"

namespace mp {

/// Result expression type for expression conversions
class EExpr : public QuadraticExpr {
public:
  /// Default constructor
  EExpr() = default;
  /// Construct from LinTerms, QuadTerms, const_term
  EExpr(LinTerms lt, QuadTerms qt, double ct)
      : QuadraticExpr( {std::move(lt), std::move(qt) }, ct ) { }
  /// Construct from LinTerms
  EExpr(LinTerms lt): QuadraticExpr( {std::move(lt), {} }, 0.0 ) { }
  /// Constructor from the Constant helper
  EExpr(Constant c) : QuadraticExpr(c) {}
  /// Construct from the Variable helper
  EExpr(Variable v) : QuadraticExpr(v) {}
  /// Construct from c * var[i]
  EExpr(double c, int i) { add_term(c, i); }
};

/// Multiply out two EEXprs.
/// Sort both factors (@todo only in Debug),
/// then produce sorted result
EExpr MultiplyOut(const EExpr& el, const EExpr& er);

} // namespace mp

#endif // EEXPR_H

#ifndef EEXPR_H
#define EEXPR_H

#include <utility>

#include "mp/flat/expr_quadratic.h"

namespace mp {

/// Result expression type for expression conversions
class EExpr : public QuadraticExpr {
public:
  EExpr() = default;
  EExpr(LinTerms lt): QuadraticExpr( {std::move(lt), {} }, 0.0 ) { }
  EExpr(Constant c) : QuadraticExpr(c) {}
  EExpr(Variable v) : QuadraticExpr(v) {}
  EExpr(int i, double c) { add_term(i, c); }
};

} // namespace mp

#endif // EEXPR_H

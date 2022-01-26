#ifndef EEXPR_H
#define EEXPR_H

#include "mp/flat/quad_expr.h"

namespace mp {

/// Result expression type for expression conversions
class EExpr : public QuadExp {
public:
  EExpr() = default;
  EExpr(AffExp ae) noexcept : QuadExp(std::move(ae)) { }
  EExpr(Constant c) : QuadExp(c) {}
  EExpr(Variable v) : QuadExp(v) {}
  EExpr(int i, double c) { add_linear_term(i, c); }
};

} // namespace mp

#endif // EEXPR_H

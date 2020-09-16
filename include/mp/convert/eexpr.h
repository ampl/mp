#ifndef EEXPR_H
#define EEXPR_H

#include "mp/convert/quad_expr.h"

namespace mp {

/// Result expression type for expression conversions
class EExpr : public QuadExpr {
public:
  EExpr() = default;
  EExpr(AffineExpr&& ae) : QuadExpr(std::move(ae)) { }
  EExpr(Constant c) : QuadExpr(c) {}
  EExpr(Variable v) : QuadExpr(v) {}
  EExpr(int i, double c) { AddLinearTerm(i, c); }
};

} // namespace mp

#endif // EEXPR_H

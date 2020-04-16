#ifndef STD_CONSTR_H
#define STD_CONSTR_H

#include "mp/convert/basic_constr.h"
#include "mp/convert/affine_expr.h"

namespace mp {

/// Linear Defining Constraint: r = affine_expr
class LinearDefiningConstraint :
    public BasicConstraint, public DefiningConstraint {
  AffineExpr affine_expr_;
public:
  LinearDefiningConstraint(AffineExpr&& ae, int r) :
    DefiningConstraint(r), affine_expr_(std::move(ae)) { }
  const AffineExpr& GetAffineExpr() const { return affine_expr_; }
};

/// Maximum: r = max(v1, v2, ..., vn)
/// TODO template: using MaximumConstraint = VarArrayDefiningConstraint<...>
class MaximumConstraint :
    public VarArrayArgConstraint, public DefiningConstraint {
public:
  MaximumConstraint(ArgArray&& aa, int r) :
    VarArrayArgConstraint(std::move(aa)), DefiningConstraint(r) { }
};

/// Minimum: r = min(v1, v2, ..., vn)
class MinimumConstraint :
    public VarArrayArgConstraint, public DefiningConstraint {
public:
  MinimumConstraint(ArgArray&& aa, int r) :
    VarArrayArgConstraint(std::move(aa)), DefiningConstraint(r) { }
};

} // namespace mp

#endif // STD_CONSTR_H

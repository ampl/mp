#ifndef STD_CONSTR_H
#define STD_CONSTR_H

#include <vector>

#include "mp/convert/basic_constr.h"
#include "mp/convert/affine_expr.h"

namespace mp {

/// Standard linear constraint
class LinearConstraint : public BasicConstraint {
  const std::vector<double> coefs_;
  const std::vector<int> vars_;
  const double lb_, ub_;
public:
  LinearConstraint(std::vector<double>&& c, std::vector<int>&& v, double l, double u)
    : coefs_(c), vars_(v), lb_(l), ub_(u) { assert(coefs_.size()==vars_.size()); }
  int nnz() const { return (int)coefs_.size(); }
  const double* coefs() const { return coefs_.data(); }
  const int* vars() const { return vars_.data(); }
  double lb() const { return lb_; }
  double ub() const { return ub_; }
};

/// Converting linear expr to 2 vectors.
struct LinearExprUnzipper {
  std::vector<double> c_;
  std::vector<int> v_;
  LinearExprUnzipper(const LinearExpr& e) {
    Reserve(e.num_terms());
    for (LinearExpr::const_iterator it=e.begin(); it!=e.end(); ++it) {
      AddTerm(it->coef(), it->var_index());
    }
  }
  void Reserve(size_t s) { c_.reserve(s); v_.reserve(s); }
  void AddTerm(double c, int v) { c_.push_back(c); v_.push_back(v); }
};

/// Linear Defining Constraint: r = affine_expr
class LinearDefiningConstraint :
    public BasicConstraint, public DefiningConstraint {
  AffineExpr affine_expr_;
public:
  LinearDefiningConstraint(AffineExpr&& ae, int r) :
    DefiningConstraint(r), affine_expr_(std::move(ae)) {
    /// TODO sort elements
  }
  const AffineExpr& GetAffineExpr() const { return affine_expr_; }
  LinearConstraint to_linear_constraint() const {
    const auto& ae = GetAffineExpr();
    LinearExprUnzipper aeu(ae);
    aeu.AddTerm(-1.0, GetResultVar());
    return LinearConstraint(std::move(aeu.c_), std::move(aeu.v_),
                            -ae.constant_term(), -ae.constant_term());
  }
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

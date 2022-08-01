#ifndef EXPR_BOUNDS_H
#define EXPR_BOUNDS_H

/**
 * Bounds and type computation
 * for various expressions
 */

#include "mp/flat/preprocess.h"
#include "mp/flat/expr_quadratic.h"

namespace mp {

/// A mix-in base class
/// providing expression bounds computation
template <class Impl>
class BoundComputations {
public:

  /// ComputeBoundsAndType(LinTerms)
  PreprocessInfoStd ComputeBoundsAndType(const LinTerms& lt) {
    PreprocessInfoStd result;
    result.lb_ = result.ub_ = 0.0;
    result.type_ = var::INTEGER;
    auto& model = MP_DISPATCH( GetModel() );
    for (auto i=lt.size(); i--; ) {
      auto v = lt.var(i);
      auto c = lt.coef(i);
      if (c >= 0.0) {
        result.lb_ += c * model.lb(v);
        result.ub_ += c * model.ub(v);
      } else {
        result.lb_ += c * model.ub(v);
        result.ub_ += c * model.lb(v);
      }
      if (var::INTEGER!=model.var_type(v) || !is_integer(c)) {
        result.type_=var::CONTINUOUS;
      }
    }
    return result;
  }

  /// ComputeBoundsAndType(AlgebraicExpr<>)
  template <class Body>
  PreprocessInfoStd ComputeBoundsAndType(
      const AlgebraicExpression<Body>& ae) {
    PreprocessInfoStd result = ComputeBoundsAndType(ae.GetBody());
    result.lb_ += ae.constant_term();
    result.ub_ += ae.constant_term();
    if (!is_integer(ae.constant_term()))
      result.type_ = var::CONTINUOUS;
    return result;
  }

  /// ComputeBoundsAndType(QuadTerms)
  PreprocessInfoStd ComputeBoundsAndType(const QuadTerms& qt) {
    PreprocessInfoStd result;
    result.lb_ = result.ub_ = 0.0;
    result.type_ = var::INTEGER;
    auto& model = MP_DISPATCH( GetModel() );
    for (auto i=qt.size(); i--; ) {
      auto coef = qt.coef(i);
      auto v1 = qt.var1(i);
      auto v2 = qt.var2(i);
      auto prodBnd = ProductBounds(v1, v2);
      if (coef >= 0.0) {
        result.lb_ += coef * prodBnd.first;
        result.ub_ += coef * prodBnd.second;
      } else {
        result.lb_ += coef * prodBnd.second;
        result.ub_ += coef * prodBnd.first;
      }
      if (var::INTEGER!=model.var_type(v1) ||
          var::INTEGER!=model.var_type(v2) ||
          !is_integer(coef)) {
        result.type_=var::CONTINUOUS;
      }
    }
    return result;
  }

  /// ComputeBoundsAndType(QuadAndLinearTerms)
  PreprocessInfoStd ComputeBoundsAndType(const QuadAndLinTerms& qlt) {
    auto bntLT = ComputeBoundsAndType(qlt.GetLinTerms());
    auto bntQT = ComputeBoundsAndType(qlt.GetQPTerms());
    return AddBoundsAndType(bntLT, bntQT);
  }

  /// Product bounds
  template <class Var>
  std::pair<double, double> ProductBounds(Var x, Var y) const {
    const auto& m = MPCD(GetModel());
    auto lx=m.lb(x), ly=m.lb(y), ux=m.ub(x), uy=m.ub(y);
    std::array<double, 4> pb{lx*ly, lx*uy, ux*ly, ux*uy};
    return { *std::min_element(pb.begin(), pb.end()),
          *std::max_element(pb.begin(), pb.end()) };
  }

  /// Add / merge bounds and type
  static
  PreprocessInfoStd AddBoundsAndType(const PreprocessInfoStd& bnt1,
                                     const PreprocessInfoStd& bnt2) {
    return {
      bnt1.lb()+bnt2.lb(), bnt1.ub()+bnt2.ub(),
      var::INTEGER==bnt1.type() && var::INTEGER==bnt2.type() ?
        var::INTEGER : var::CONTINUOUS
    };
  }
};

} // namespace mp

#endif // EXPR_BOUNDS_H

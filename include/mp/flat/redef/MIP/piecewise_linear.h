#ifndef PIECEWISE_LINEAR_H
#define PIECEWISE_LINEAR_H

/**
 *  Redefinition of PL flat constraint into SOS2 + linear
 */

#include <algorithm>
#include <numeric>

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Convert PLConstraint into SOS2 + linear constraints
template <class ModelConverter>
class PLConverter_MIP :
    public BasicFuncConstrCvt<
      PLConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    PLConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  PLConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = PLConstraint;

  /// Convert in any context
  void Convert(const ItemType& cc, int ) {
    points_ = cc.GetParameters();      // extract PLPoints
    i0=0;                              // first breakpoint
    i1=points_.x_.size()-1;            // last breakpoint
    y = cc.GetResultVar();
    x = cc.GetArguments()[0];
    MP_ASSERT_ALWAYS(i1>=i0, "PL->SOS2: no breakpoints");
    if (i1>i0)                         // Gurobi 9 does this.
      ConsiderExtendingEndSegments();  // Bad when approximating
    ConsiderShorteningPL();
    if (ConsiderDegenerateCases())
      return;
    if (ConsiderConvexity())
      return;
    RedefineInSOS2();
  }


protected:
  /// See if we need to extend the end segments
  /// to cover larger domain
  void ConsiderExtendingEndSegments() {
    if (GetMC().lb(x) < points_.x_.front() &&
        points_.x_.front() > -PLMaxVal())
      ExtendSegTo(0, 1, std::max(GetMC().lb(x), -PLMaxVal()));
    if (GetMC().ub(x) > points_.x_.back() &&
        points_.x_.back() < PLMaxVal())
      ExtendSegTo(i1, i1-1, std::min(GetMC().ub(x), PLMaxVal()));
  }
  /// See if we need to remove/shorten side segments
  void ConsiderShorteningPL() {
    while (i0<i1 && GetMC().lb(x)>=points_.x_[i0+1])
      ++i0;
    while (i0<i1 && GetMC().ub(x)<=points_.x_[i1-1])
      --i1;
    if (i1>i0) {             // don't need otherwise
      if (GetMC().lb(x) > points_.x_[i0])
        ExtendSegTo(i0, i0+1, GetMC().lb(x));
      if (GetMC().ub(x) < points_.x_[i1])
        ExtendSegTo(i1, i1-1, GetMC().ub(x));
    }
  }
  /// See if we obtain a single point or segment
  /// @return true iff that
  /// (a corresponding simpler constraint is added)
  bool ConsiderDegenerateCases() {
    return false;
  }
  /// See if have a convex case
  bool ConsiderConvexity() {
    return false;
  }
  /// Non-convex redefinition into SOS2 + linear
  void RedefineInSOS2() {
    auto lambda = GetMC().AddVars_returnIds(i1-i0+1, 0.0, 1.0);
    std::vector<double> weights(i1-i0+1);
    std::iota(weights.begin(), weights.end(), 1.0);
    GetMC().AddConstraint(      // indicate range of sum(lambda)
          SOS2Constraint(lambda, weights,
                         SOSExtraInfo{ {1.0, 1.0} }) );
    std::fill(weights.begin(), weights.end(), 1.0);
    GetMC().AddConstraint(
          LinConEQ{ {weights, lambda}, {1.0} });
    weights.assign(points_.y_.begin()+i0, points_.y_.begin()+i1+1);
    GetMC().RedefineVariable(y,        // Could just add constraint
                             LinearFunctionalConstraint{
                               {{weights, lambda}, 0.0} });
    weights.assign(points_.x_.begin()+i0, points_.x_.begin()+i1+1);
    weights.push_back(-1.0);
    lambda.push_back(x);
    GetMC().AddConstraint(
          LinConEQ{ {weights, lambda}, {0.0} });
  }

  /// PLMaxVal, currently constant.
  /// Default max abs value of the argument and of the result
  /// of a PL. Applied when the PL is defined on a smaller
  /// domain but the argument allows more.
  double PLMaxVal() const { return 1e6; }

  /// Extend segment (i0, i1) (if i0<i1, then to the left,
  /// otherwise to the right) in PLPoints to start/end
  /// in new x0
  void ExtendSegTo(size_t i0, size_t i1, double x0) {
    assert(i0!=i1);
    assert(2<=points_.x_.size());
    auto& x0_old = points_.x_[i0];
    auto x1 = points_.x_[i1];
    int sign = i0<i1 ? 1 : -1;
    assert(sign*(x1-x0) > 0.0);
    assert(sign*(x1-x0_old)>0.0);
    auto& y0_old = points_.y_[i0];
    auto y1 = points_.y_[i1];
    if (x0 != x0_old) {
      auto slope = (y1-y0_old) / (x1-x0_old);
      y0_old = y1 - slope*(x1-x0);
      x0_old = x0;
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;


private:
  PLPoints points_;
  size_t i0{0};
  size_t i1{0};
  int x{-1};
  int y{-1};
};

} // namespace mp

#endif // PIECEWISE_LINEAR_H

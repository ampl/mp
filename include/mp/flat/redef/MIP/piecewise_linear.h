#ifndef PIECEWISE_LINEAR_H
#define PIECEWISE_LINEAR_H

/*
 *  Redefinition of PL flat constraint into SOS2 + linear
 */

#include <algorithm>
#include <numeric>
#include <cassert>

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
  Context Convert(const ItemType& cc, int ) {
    points_ = cc.GetParameters().GetPLPoints();
    i0=0;                              // first breakpoint
    i1=points_.x_.size()-1;            // last breakpoint
		MP_ASSERT_ALWAYS(i1>=i0, "PL->SOS2: no breakpoints");
		y = cc.GetResultVar();
    x = cc.GetArguments()[0];
		if (ConsiderDegenerateCases(cc))
			return Context::CTX_MIX;
		if (ConsiderConvexity(cc))
			return cc.GetContext();          // convex case
    if (i1>i0)                         // Gurobi 9 does this.
      ConsiderExtendingEndSegments();  // Bad when approximating
    ConsiderShorteningPL();
    RedefineInSOS2();
		return Context::CTX_MIX;           // general case
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
	bool ConsiderDegenerateCases(const ItemType& cc) {
		assert(points_ == cc.GetParameters().GetPLPoints());
		if (2 >= points_.size()) {
			y = cc.GetResultVar();
			x = cc.GetArguments()[0];
			MP_ASSERT_ALWAYS(points_.size(), "PL expression: no breakpoints");
			if (1==points_.size()) {
				GetMC().NarrowVarBounds(x, points_.x_[0], points_.x_[0]);
				GetMC().NarrowVarBounds(y, points_.y_[0], points_.y_[0]);
			} else if (2==points_.size()) {
				MP_ASSERT_ALWAYS(points_.x_[0] < points_.x_[1],
						"Empty 1-segment PL expression");
				long double slope
						= ((long double)(points_.y_[1] - points_.y_[0]))
						/ (points_.x_[1] - points_.x_[0]);
				GetMC().AddConstraint(
							LinConEQ{ {{1.0, double(-slope)}, {y, x}},
							points_.y_[0] - slope*points_.x_[0]});
			}
			return true;
		}
    return false;
  }
  /// See if have a convex case
	bool ConsiderConvexity(const ItemType& cc) {
		if (cc.GetContext().IsPositive()) {
			if (IsConcave(cc.GetParameters())) {
				RedefineConcave(cc);
				return true;
			}
		} else if (cc.GetContext().IsNegative()) {
			if (IsConvex(cc.GetParameters())) {
				RedefineConvex(cc);
				return true;
			}
		}
    return false;
  }
	/// Redefine concave PL
	void RedefineConcave(const ItemType& cc) {
		const auto& slopes = cc.GetParameters().GetPLSlopes().GetSlopes();
		assert(slopes.size() == points_.size()-1);
		for (auto i=slopes.size(); i--; ) {
			auto Xi = points_.x_[i+1];
			auto Yi = points_.y_[i+1];
			GetMC().AddConstraint(
						LinConLE{ {{1.0, -slopes[i]}, {y, x}},
						Yi - slopes[i]*Xi});
		}
	}
	/// Redefine convex PL
	void RedefineConvex(const ItemType& cc) {
		const auto& slopes = cc.GetParameters().GetPLSlopes().GetSlopes();
		assert(slopes.size() == points_.size()-1);
		for (auto i=slopes.size(); i--; ) {
			auto Xi = points_.x_[i+1];
			auto Yi = points_.y_[i+1];
			GetMC().AddConstraint(
						LinConGE{ {{1.0, -slopes[i]}, {y, x}},
						Yi - slopes[i]*Xi});
		}
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
    LinearFunctionalConstraint funccon{ {{weights, lambda}, 0.0} };
    funccon.SetContext( GetMC().GetInitExprContext(y) );
    GetMC().RedefineVariable(y, std::move(funccon));
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

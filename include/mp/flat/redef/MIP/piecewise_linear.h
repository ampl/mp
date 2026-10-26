#ifndef PIECEWISE_LINEAR_H
#define PIECEWISE_LINEAR_H

/*
 *  Redefinition of PL flat constraint into SOS2 + linear
 */

#include <algorithm>
#include <numeric>
#include <cassert>

#include "mp/format.h"
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
    rayL_ = rayR_ = false;
    i0=0;                              // first breakpoint
    i1=points_.x_.size()-1;            // last breakpoint
		MP_ASSERT_ALWAYS(i1>=i0, "PL->SOS2: no breakpoints");
		y = cc.GetResultVar();
    x = cc.GetArguments()[0];
    NarrowBoundsAtConversion(cc);
		if (ConsiderDegenerateCases(cc))
			return Context::CTX_MIX;
		if (ConsiderConvexity(cc))
			return cc.GetContext();          // convex case
    if (i1>i0)                         // Gurobi 9 does this.
      ConsiderExtendingEndSegments();  // Bad when approximating
    ConsiderShorteningPL();
    if (ConsiderShortenedDegenerateCases())
      return Context::CTX_MIX;
    RedefineInSOS2();
		return Context::CTX_MIX;           // general case
  }


protected:
  /// Narrow the bounds of the argument, if it is a linear expression,
  /// and of the result, from the current bounds.
  /// They can be tighter than when this PL was created:
  /// linear constraints can be flattened later (AMPL writes them
  /// after the nonlinear ones), and the PLs this one depends on
  /// are converted earlier (e.g., fixed, see FixAtPoint()).
  void NarrowBoundsAtConversion(const ItemType& cc) {
    if (auto plin = GetMC().template
        GetInitExpressionOfType<LinearFunctionalConstraint>(x)) {
      auto bnt = GetMC().ComputeBoundsAndType(plin->GetAffineExpr());
      GetMC().NarrowVarBestBounds(x, bnt.lb(), bnt.ub());
    }
    const auto range = GetMC().PLResultRange(cc);
    GetMC().NarrowVarBestBounds(y, range.first, range.second);
  }
  /// See if we need to extend the end segments
  /// to cover larger domain.
  /// Unbounded sides with native SOS2 get rays instead.
  void ConsiderExtendingEndSegments() {
    const auto dom = ExtensionDomain();
    if (!rayL_ && GetMC().lb(x) < points_.x_.front() &&
        points_.x_.front() > -dom)
      ExtendSegTo(0, 1, std::max(GetMC().lb(x), -dom));
    if (!rayR_ && GetMC().ub(x) > points_.x_.back() &&
        points_.x_.back() < dom)
      ExtendSegTo(i1, i1-1, std::min(GetMC().ub(x), dom));
  }
  /// Max abs value to extend the end segments to (cvt:pl:domain).
  /// The MIP reformulation of SOS2 amplifies the solver's tolerances
  /// by this value: unless the option is set, it uses the argument's
  /// bounds and fails if one is infinite (as cvt:bigM for indicators).
  /// Native SOS2, unless the option is set, gets exact rays
  /// for infinite sides (see RedefineInSOS2()).
  double ExtensionDomain() {
    const auto dom = GetMC().PLDomain();
    const bool native = !GetMC().GetConstraintKeeper(
          (SOS2Constraint*)nullptr).IfConvertsItems();
    const auto lbx = GetMC().lb(x), ubx = GetMC().ub(x);
    const auto lbExt = lbx < points_.x_.front() ? lbx : 0.0;
    const auto ubExt = ubx > points_.x_.back() ? ubx : 0.0;
    if (dom < 0.0) {
      const bool infL = lbExt <= GetMC().PracticallyMinusInf();
      const bool infR = ubExt >= GetMC().PracticallyInf();
      if (native) {
        rayL_ = infL;
        rayR_ = infR;
      } else if (infL || infR)
        FailUnbounded();
      return INFINITY;
    }
    if (!native &&
        ((lbExt < -dom && points_.x_.front() > -dom)
         || (ubExt > dom && points_.x_.back() < dom)))
      WarnDomainLimited();
    return dom;
  }
  /// Fail: no MIP reformulation for an unbounded argument
  [[noreturn]] void FailUnbounded() {
    throw ConstraintConversionFailure(
          "PLInfBound",
          "The argument of a piecewise-linear expression\n"
          "is unbounded, which prevents its MIP reformulation.\n"
          "Set bounds on the argument,\n"
          "or use option cvt:pl:domain (with caution).\n"
          "See more: mp.ampl.com/modeling-numerics.html");
  }
  /// Warn that the argument's domain is limited by cvt:pl:domain
  void WarnDomainLimited() {
    GetMC().AddWarning("PLDomain", fmt::format(
        "The domain of a piecewise-linear expression's argument\n"
        "has been limited to +-{} (cvt:pl:domain).\n"
        "The MIP reformulation can violate the expression\n"
        "by up to this value times the feasibility tolerance.",
        GetMC().PLDomain()));
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
			MP_ASSERT_ALWAYS(points_.size(), "PL expression: no breakpoints");
			if (1==points_.size()) {
				FixAtPoint(0);
			} else {
				MP_ASSERT_ALWAYS(points_.x_[0] < points_.x_[1],
						"Empty 1-segment PL expression");
				RedefineAsSegment(0);
			}
			return true;
		}
    return false;
  }
  /// After shortening to the argument's domain,
  /// see if a single point or segment remains.
  /// A single breakpoint would otherwise give a 1-variable SOS2.
  /// @return true iff that
  /// (a corresponding simpler constraint is added)
  bool ConsiderShortenedDegenerateCases() {
    if (points_.x_[i0] >= points_.x_[i1])   // fixed argument
      FixAtPoint(i0);
    else if (i0+1 == i1)
      RedefineAsSegment(i0);
    else
      return false;
    return true;
  }
  /// Fix argument and result at breakpoint \a i
  void FixAtPoint(size_t i) {
    GetMC().NarrowVarBounds(x, points_.x_[i], points_.x_[i]);
    GetMC().NarrowVarBounds(y, points_.y_[i], points_.y_[i]);
  }
  /// Redefine as the linear function of segment (i, i+1).
  /// Requires x_[i] < x_[i+1].
  void RedefineAsSegment(size_t i) {
    long double slope
        = ((long double)(points_.y_[i+1] - points_.y_[i]))
        / (points_.x_[i+1] - points_.x_[i]);
    GetMC().AddConstraint(
          LinConEQ{ {{1.0, double(-slope)}, {y, x}},
          double(points_.y_[i] - slope*points_.x_[i])});
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
	/// Non-convex redefinition into SOS2 + linear:
  /// x = sum(lambda_i * x_i), y = sum(lambda_i * y_i), sum(lambda)=1.
  /// An unbounded side gets a ray variable r>=0 as an extra
  /// SOS2 member next to the end lambda: then that lambda is 1,
  /// and (x, y) moves from the end point along the end segment.
  void RedefineInSOS2() {
    auto lambda = GetMC().AddVars_returnIds(i1-i0+1, 0.0, 1.0);
    std::vector<double> weights(lambda.size(), 1.0);
    GetMC().AddConstraint(
          LinConEQ{ {weights, lambda}, {1.0} });
    std::vector<double> cx(points_.x_.begin()+i0, points_.x_.begin()+i1+1);
    std::vector<double> cy(points_.y_.begin()+i0, points_.y_.begin()+i1+1);
    auto sosvars = lambda;
    if (rayL_)
      AddRay(true, sosvars, lambda, cx, cy);
    if (rayR_)
      AddRay(false, sosvars, lambda, cx, cy);
    weights.resize(sosvars.size());
    std::iota(weights.begin(), weights.end(), 1.0);
    GetMC().AddConstraint(      // indicate range of sum(lambda)
          SOS2Constraint(sosvars, weights,
                         (rayL_ || rayR_) ? SOSExtraInfo{}
                                          : SOSExtraInfo{ {1.0, 1.0} }) );
    LinearFunctionalConstraint funccon{ {{cy, lambda}, 0.0} };
    funccon.SetContext( GetMC().GetInitExprContext(y) );
    GetMC().RedefineVariable(y, std::move(funccon));
    cx.push_back(-1.0);
    lambda.push_back(x);
    GetMC().AddConstraint(
          LinConEQ{ {cx, lambda}, {0.0} });
  }
  /// Add a ray variable r>=0 along the left/right end segment:
  /// to \a sosvars at that end, to \a lambda with coefficients
  /// -1/+1 in \a cx and -/+slope in \a cy.
  void AddRay(bool left, std::vector<int>& sosvars,
              std::vector<int>& lambda,
              std::vector<double>& cx, std::vector<double>& cy) {
    const auto r = int( GetMC().AddVar(0.0, INFINITY) );
    const size_t j0 = left ? i0 : i1-1;           // end segment
    const auto slope = (points_.y_[j0+1]-points_.y_[j0])
        / (points_.x_[j0+1]-points_.x_[j0]);
    const double dir = left ? -1.0 : 1.0;
    sosvars.insert(left ? sosvars.begin() : sosvars.end(), r);
    lambda.push_back(r);
    cx.push_back(dir);
    cy.push_back(dir*slope);
  }

  /// Extend segment (i0, i1) (if i0<i1, then to the left,
  /// otherwise to the right) in PLPoints to start/end
  /// in new x0.
  /// x0 can coincide with x1 when shortening both ends of a segment
  /// to a fixed argument; ConsiderShortenedDegenerateCases()
  /// then fixes the result at that point.
  void ExtendSegTo(size_t i0, size_t i1, double x0) {
    assert(i0!=i1);
    assert(2<=points_.x_.size());
    auto& x0_old = points_.x_[i0];
    auto x1 = points_.x_[i1];
    int sign = i0<i1 ? 1 : -1;
    assert(sign*(x1-x0) >= 0.0);
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
  bool rayL_{false};                   // rays for unbounded sides
  bool rayR_{false};
  int x{-1};
  int y{-1};
};

} // namespace mp

#endif // PIECEWISE_LINEAR_H

#ifndef INT_DIV_H
#define INT_DIV_H

#include <algorithm>
#include <cassert>
#include <cmath>

#include "mp/error.h"
#include "mp/flat/redef/redef_base.h"
#include "mp/flat/redef/MIP/mul.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// p = y*q as a MulConstraint, arguments ordered
/// (so that IntDiv and Mod share it)
template <class MC>
int IntDivProduct(MC& mc, int y, int q) {
  assert(y != q);
  return AssignResultVarMixedCtx(mc,
        MulConstraint{ {std::min(y, q), std::max(y, q)} });
}

/// Converts q = x div y = trunc(x/y) for MIP, as MiniZinc's
/// linearization library: with the remainder rem = x - y*q,
///   |rem| <= |y| - eps,  rem has the sign of x (or is 0),
/// eps = 1 for integer rem, cvt:mip:eps otherwise.
/// For a variable y, y*q is a MulConstraint (exact for integer q
/// with finite bounds, see MulConverter_MIP), and y != 0 is enforced.
template <class ModelConverter>
class IntDivConverter_MIP :
    public BasicFuncConstrCvt<
      IntDivConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    IntDivConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  IntDivConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = IntDivConstraint;

  /// Convert in any context
  Context Convert(const ItemType& dc, int ) {
    const int x = dc.GetArguments()[0], y = dc.GetArguments()[1];
    const int q = dc.GetResultVar();
    assert(var::INTEGER == GetMC().var_type(q));  // see preprocessing
    assert(!GetMC().is_fixed(y) || 0.0 != GetMC().fixed_value(y)); // prepro
    GetMC().MarkVarTypeRequired(q);
    // The remainder enters the constraints as linear terms:
    // a result variable used only in bounds would be dropped as unused
    const LinTerms rem = Remainder(x, y, q);
    BoundRemainderByDivisor(y, rem);
    SignOfRemainder(x, rem);
    return Context::CTX_MIX;
  }

  /// rem = x - y*q, as linear terms
  LinTerms Remainder(int x, int y, int q) {
    if (GetMC().is_fixed(y))
      return { {1.0, -GetMC().fixed_value(y)}, {x, q} };
    return { {1.0, -1.0}, {x, IntDivProduct(GetMC(), y, q)} };
  }

  /// |rem| <= |y| - eps, y != 0
  void BoundRemainderByDivisor(int y, const LinTerms& rem) {
    const double eps = GetMC().ComparisonEps(
          GetMC().ComputeBoundsAndType(rem));
    if (GetMC().is_fixed(y)) {
      const double b = std::fabs(GetMC().fixed_value(y)) - eps;
      GetMC().AddConstraint(LinConLE(LinTerms(rem), {b}));
      GetMC().AddConstraint(LinConGE(LinTerms(rem), {-b}));
      return;
    }
    const int w = AssignResultVarMixedCtx(GetMC(),  // w = |y|
          AbsConstraint{ {y} } );
    LinTerms rem_minus_w = rem, rem_plus_w = rem;
    rem_minus_w.add_term(-1.0, w);
    rem_plus_w.add_term(1.0, w);
    GetMC().AddConstraint(LinConLE(std::move(rem_minus_w), {-eps}));
    GetMC().AddConstraint(LinConGE(std::move(rem_plus_w), {eps}));
    GetMC().AddConstraint(LinConGE({{1.0}, {w}},       // y != 0
                                   {GetMC().ComparisonEps(y)}));
  }

  /// rem has the sign of x: x >= 0 ==> rem >= 0, x <= 0 ==> rem <= 0
  void SignOfRemainder(int x, const LinTerms& rem) {
    if (GetMC().lb(x) >= 0.0)
      GetMC().AddConstraint(LinConGE(LinTerms(rem), {0.0}));
    else if (GetMC().ub(x) <= 0.0)
      GetMC().AddConstraint(LinConLE(LinTerms(rem), {0.0}));
    else {
      const int b = (int)GetMC().AddVar(0.0, 1.0, var::INTEGER);
      GetMC().AddConstraint(IndicatorConstraintLinGE(b, 1, {{{1.0}, {x}}, 0.0}));
      GetMC().AddConstraint(IndicatorConstraintLinGE(b, 1, {LinTerms(rem), 0.0}));
      GetMC().AddConstraint(IndicatorConstraintLinLE(b, 0, {{{1.0}, {x}}, 0.0}));
      GetMC().AddConstraint(IndicatorConstraintLinLE(b, 0, {LinTerms(rem), 0.0}));
    }
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// Converts r = x mod y = x - y*(x div y) for MIP,
/// via q = x div y (IntDivConstraint).
template <class ModelConverter>
class ModConverter_MIP :
    public BasicFuncConstrCvt<
      ModConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    ModConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  ModConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = ModConstraint;

  /// Convert in any context
  Context Convert(const ItemType& mc, int ) {
    const int x = mc.GetArguments()[0], y = mc.GetArguments()[1];
    const int r = mc.GetResultVar();
    const int q = AssignResultVarMixedCtx(GetMC(),
                                          IntDivConstraint{ {x, y} });
    if (GetMC().is_fixed(y))                       // r - x + c*q = 0
      GetMC().AddConstraint(LinConEQ(
          { {1.0, -1.0, GetMC().fixed_value(y)}, {r, x, q} }, {0.0}));
    else                                           // r - x + y*q = 0
      GetMC().AddConstraint(LinConEQ(
          { {1.0, -1.0, 1.0}, {r, x, IntDivProduct(GetMC(), y, q)} }, {0.0}));
    return Context::CTX_MIX;
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // INT_DIV_H

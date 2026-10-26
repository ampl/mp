#ifndef ROUNDING_H
#define ROUNDING_H

#include <cassert>
#include <cmath>

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Converts r = floor(x) for MIP:
/// r integer, r <= x <= r + 1 - eps,
/// with eps = cvt:mip:eps (x is not integer: then preprocessing
/// uses x itself).
template <class ModelConverter>
class FloorConverter_MIP :
    public BasicFuncConstrCvt<
      FloorConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    FloorConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  FloorConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = FloorConstraint;

  /// Convert in any context
  Context Convert(const ItemType& fc, int ) {
    const int x = fc.GetArguments()[0];
    const int r = fc.GetResultVar();
    assert(!GetMC().is_var_integer(x));   // preprocessing: r is x then
    GetMC().MarkVarTypeRequired(r);       // integrality defines r
    GetMC().AddConstraint(LinConLE({{1.0, -1.0}, {r, x}}, {0.0}));
    GetMC().AddConstraint(LinConLE({{1.0, -1.0}, {x, r}},
                                   {1.0 - GetMC().ComparisonEps(x)}));
    return Context::CTX_MIX;
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// Converts r = ceil(x) for MIP:
/// r integer, r - 1 + eps <= x <= r.
template <class ModelConverter>
class CeilConverter_MIP :
    public BasicFuncConstrCvt<
      CeilConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    CeilConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  CeilConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = CeilConstraint;

  /// Convert in any context
  Context Convert(const ItemType& cc, int ) {
    const int x = cc.GetArguments()[0];
    const int r = cc.GetResultVar();
    assert(!GetMC().is_var_integer(x));
    GetMC().MarkVarTypeRequired(r);
    GetMC().AddConstraint(LinConLE({{1.0, -1.0}, {x, r}}, {0.0}));
    GetMC().AddConstraint(LinConLE({{1.0, -1.0}, {r, x}},
                                   {1.0 - GetMC().ComparisonEps(x)}));
    return Context::CTX_MIX;
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};


/// Converts r = round(x) (ties to even, as AMPL) for MIP:
/// r = 2m + p (m integer, p binary),
/// |x - r| <= 1/2 - eps*p,
/// so that an odd r is only possible strictly inside
/// the half-unit band, and x = k + 1/2 rounds to the even neighbor.
template <class ModelConverter>
class RoundConverter_MIP :
    public BasicFuncConstrCvt<
      RoundConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    RoundConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  RoundConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = RoundConstraint;

  /// Convert in any context
  Context Convert(const ItemType& rc, int ) {
    const int x = rc.GetArguments()[0];
    const int r = rc.GetResultVar();
    assert(!GetMC().is_var_integer(x));
    GetMC().MarkVarTypeRequired(r);
    const int m = (int)GetMC().AddVar(             // r = 2m + p
          std::floor((GetMC().lb(r) - 1.0) / 2.0),
          std::ceil(GetMC().ub(r) / 2.0), var::INTEGER);
    const int p = (int)GetMC().AddVar(0.0, 1.0, var::INTEGER);
    GetMC().AddConstraint(LinConEQ({{1.0, -2.0, -1.0}, {r, m, p}}, {0.0}));
    const double eps = GetMC().ComparisonEps(x);
    GetMC().AddConstraint(LinConLE({{1.0, -1.0, eps}, {x, r, p}}, {0.5}));
    GetMC().AddConstraint(LinConLE({{1.0, -1.0, eps}, {r, x, p}}, {0.5}));
    return Context::CTX_MIX;
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // ROUNDING_H

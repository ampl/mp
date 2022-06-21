#ifndef POWER_CONST_H
#define POWER_CONST_H

#include <cmath>
#include <cassert>

#include "mp/flat/redef/redef_base.h"
#include "mp/flat/constr_std.h"

namespace mp {

/// Converts PowConstraint (const exponent) for MIP
template <class ModelConverter>
class PowConstExponentConverter_MIP :
    public BasicFuncConstrCvt<
      PowConstExponentConverter_MIP<ModelConverter>, ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
    PowConstExponentConverter_MIP<ModelConverter>, ModelConverter>;
  /// Constructor
  PowConstExponentConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = PowConstraint;

  /// Check whether the constraint
  /// needs to be converted despite being accepted by ModelAPI.
  bool IfNeedsConversion(const ItemType& con, int ) {
    auto pwr = con.GetParameters()[0];
    return GetMC().lb(con.GetArguments()[0]) < 0.0 &&
        GetMC().is_integer_value(pwr) && pwr>=0.0;  // TODO also < 0
  }

  /// Convert in any context
  void Convert(const ItemType& con, int ) {
    assert(!con.GetContext().IsNone());
    auto pwr = con.GetParameters()[0];
    if (!GetMC().is_integer_value(pwr) || pwr < 0.0)
      throw ConstraintConversionFailure( "PowConNegOrFracExp",
          "Not converting PowConstraint with negative "
          "or fractional exponent");
    auto arg = con.GetArguments()[0];
    auto arg1 = arg, arg2 = arg;         // new variables
    if (2.0<pwr) {
      auto pwr1 =   // for even values, split into even subpowers
            GetMC().is_integer_value(pwr / 2.0) ?
              std::floor(pwr/4.0) * 2 :
              std::floor(pwr/2.0);
      auto pwr2 = pwr-pwr1;
      assert(pwr2>0.0);
      arg1 = GetMC().AssignResultVar2Args(
            PowConstraint{ {arg}, DblParamArray1{pwr1} });
      arg2 = GetMC().AssignResultVar2Args(
            PowConstraint{ {arg}, DblParamArray1{pwr2} });
    }
    GetMC().RedefineVariable(con.GetResultVar(),
          QuadraticFunctionalConstraint(
            QuadraticExpr(
              QuadAndLinTerms( { }, { {1.0}, {arg1}, {arg2} } ),
              0.0) ));
    GetMC(). PropagateResultOfInitExpr(     // propagate ctx into new constr
          con.GetResultVar(), con.GetContext());
  }

protected:

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // POWER_CONST_H

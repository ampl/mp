#ifndef POWER_CONST_H
#define POWER_CONST_H

#include <cmath>
#include <cassert>

#include "mp/flat/redef/MIP/lin_approx.h"

namespace mp {

/// Converts PowConstraint (const exponent) for MIP
template <class ModelConverter>
class PowConstExponentConverter_MIP :
    public FuncConConverter_MIP_CRTP<      // Derive from PL Approximator
      PowConstExponentConverter_MIP<ModelConverter>,
      ModelConverter,
      PowConstraint> {
public:
  /// Base class
  using Base = FuncConConverter_MIP_CRTP<
    PowConstExponentConverter_MIP<ModelConverter>,
    ModelConverter, PowConstraint>;
  /// Constructor
  PowConstExponentConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = PowConstraint;

  /// Check whether the constraint
  /// needs to be converted despite being accepted by ModelAPI.
  /// This covers cases not accepted by GenConstrPow in Gurobi 10:
  /// negative lower bound for x while positive integer exponent.
  /// But we do even more.
  /// Note that ^2 has been quadratized in ProblemFlattener;
  /// ^3, ^4, ^6, etc, will be quadratized only with lb(x) < 0.
  bool IfNeedsConversion(const ItemType& con, int ) {
    auto pwr = con.GetParameters()[0];
    return
        GetMC().lb(con.GetArguments()[0]) < 0.0
        && pwr>0.0 && std::floor(pwr) == pwr;
  }

  /// Convert in any context
  void Convert(const ItemType& con, int i) {
    assert(!con.GetContext().IsNone());
    auto pwr = con.GetParameters()[0];
    if (GetMC().IfQuadratizePowConstPosIntExp() &&
        GetMC().is_integer_value(pwr) && pwr > 0.0)
      Convert2Quadratics(con, i);
    else
      Convert2PL(con, i);
  }


protected:
  /// Convert into quadratics
  void Convert2Quadratics(const ItemType& con, int ) {
    auto pwr = con.GetParameters()[0];
    assert(GetMC().is_integer_value(pwr) && pwr > 0.0);
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
    /// propagate ctx into new constr,
    /// particularly into the arguments which are new constraints
    GetMC().PropagateResultOfInitExpr(
          con.GetResultVar(), con.GetContext());
  }

  /// PL approximate
  void Convert2PL(const ItemType& con, int i) {
    Base::Convert(con, i);
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // POWER_CONST_H

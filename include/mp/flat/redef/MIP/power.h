#ifndef POWER_MIP_H
#define POWER_MIP_H

#include <cmath>
#include <cassert>

#include "mp/flat/redef/MIP/lin_approx.h"
#include "mp/flat/redef/MIP/mul.h"

namespace mp {

/// Converts PowConstExpConstraint for MIP;
/// for PowConstraint see below.
template <class ModelConverter>
class PowConstExponentConverter_MIP
    : public FuncConConverter_MIP_CRTP<      // Derive from PL Approximator
          PowConstExponentConverter_MIP<ModelConverter>,
          ModelConverter,
          PowConstExpConstraint> {
public:
  /// Base class
  using Base = FuncConConverter_MIP_CRTP<
    PowConstExponentConverter_MIP<ModelConverter>,
    ModelConverter, PowConstExpConstraint>;
  /// Constructor
  PowConstExponentConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = PowConstExpConstraint;

  /// Check whether the constraint
  /// needs to be converted despite being accepted by ModelAPI.
  /// This used to cover cases not accepted by GenConstrPow in Gurobi <=10:
  /// negative lower bound for x while positive integer exponent.
  /// But we did this even for all integer powers >=2
  /// to avoid PL approximation in Gurobi.
  /// Cancelled in Gurobi 11.
  /// Note that ^2 might be quadratized in ProblemFlattener.
  bool IfNeedsConversion(const ItemType& , int ) {
    return false;
  }

  /// Convert in any context
  Context Convert(const ItemType& con, int i) {
    assert(!con.GetContext().IsNone());
    auto pwr = con.GetParameters()[0];
    const bool f_int_pwr = GetMC().is_integer_value(pwr) && pwr > 0.0;
    if (GetMC().IfQuadratizePowConstPosIntExp() && f_int_pwr)
      Convert2Quadratics(con, i);
    else if (f_int_pwr && pwr >= 2.0 && FactorToExpand(
               GetMC(), con.GetArguments()[0], con.GetArguments()[0]) >= 0)
      Convert2Products(con);             // exact for integer arguments
    else
      Convert2PL(con, i);
    return Context::CTX_MIX;
  }


protected:
  /// x^pwr = x^pwr1 * x^pwr2, pwr integer: the two factors
  std::pair<int, int> SplitPower(const ItemType& con) {
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
            PowConstExpConstraint{ {arg}, DblParamArray1{pwr1} });
      arg2 = GetMC().AssignResultVar2Args(
            PowConstExpConstraint{ {arg}, DblParamArray1{pwr2} });
    }
    return { arg1, arg2 };
  }

  /// Convert into products r = x^pwr1 * x^pwr2 (MulConstraint),
  /// for an integer argument: MulConverter_MIP expands them exactly
  void Convert2Products(const ItemType& con) {
    const auto args = SplitPower(con);
    GetMC().RedefineVariable(con.GetResultVar(),
          MulConstraint{ {std::min(args.first, args.second),
                          std::max(args.first, args.second)} });
    GetMC().PropagateResultOfInitExpr(
          con.GetResultVar(), con.GetContext());
  }

  /// Convert into quadratics
  void Convert2Quadratics(const ItemType& con, int ) {
    const auto [arg1, arg2] = SplitPower(con);
    GetMC().RedefineVariable(con.GetResultVar(),
          QuadraticFunctionalConstraint(
            QuadraticExpr(
              QuadAndLinTerms( { }, { {1.0}, {arg1}, {arg2} } ),
              0.0) ));
    /// propagate ctx into new constr,
    /// particularly into the arguments
    /// which are new expressions when pwr>2
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


/// Converts PowConstraint for MIP
template <class ModelConverter>
class PowConverter_MIP
    : public BasicFuncConstrCvt<
          PowConverter_MIP<ModelConverter>,
          ModelConverter> {
public:
  /// Base class
  using Base = BasicFuncConstrCvt<
      PowConverter_MIP<ModelConverter>,
      ModelConverter>;
  /// Constructor
  PowConverter_MIP(ModelConverter& mc) : Base(mc) { }
  /// Converted item type
  using ItemType = PowConstraint;

  /// Convert in any context
  Context Convert(const ItemType& con, int i) {
    assert(!con.GetContext().IsNone());
    Convert2ExpLog(con, i);
    return Context::CTX_MIX;
  }


protected:
  void Convert2ExpLog(const ItemType& con, int ) {
    auto x = con.GetArguments()[0];
    auto y = con.GetArguments()[1];
    auto logx = GetMC().AssignResultVar2Args(
        LogConstraint{ {x} });
    auto y_logx = GetMC().AssignResultVar2Args(
        QuadraticFunctionalConstraint(
            QuadraticExpr(
                QuadAndLinTerms( { }, { {1.0}, {y}, {logx} } ),
                0.0) ));
    GetMC().RedefineVariable(con.GetResultVar(),
                             ExpConstraint({y_logx}));

    /// propagate ctx into new constr,
    /// particularly into the arguments which are new constraints
    GetMC().PropagateResultOfInitExpr(
        con.GetResultVar(), con.GetContext());
  }

  /// Reuse the stored ModelConverter
  using Base::GetMC;
};

} // namespace mp

#endif // POWER_MIP_H
